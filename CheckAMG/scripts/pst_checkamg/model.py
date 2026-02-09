from collections import defaultdict

import pst
import torch
from pst.typing import OptTensor
from torch_geometric.nn import MLP
from torchmetrics.functional.classification import multiclass_accuracy

from CheckAMG.scripts.pst_checkamg.config import ModelConfig
from CheckAMG.scripts.pst_checkamg.constants import NUM_CLASSES, TRUE_AVG_IDX
from CheckAMG.scripts.pst_checkamg.loss import AVGViewsTripletLoss, TripletLoss


class Attention(torch.nn.Module):
    def __init__(self, embed_dim: int, num_heads: int, **kwargs) -> None:
        super().__init__()

        self.attn = torch.nn.MultiheadAttention(
            embed_dim=embed_dim, num_heads=num_heads, **kwargs
        )
        self.attn_norm = torch.nn.LayerNorm(embed_dim)
        self.lin = MLP(
            in_channels=embed_dim,
            out_channels=embed_dim,
            hidden_channels=embed_dim,
            num_layers=2,
        )
        self.lin_norm = torch.nn.LayerNorm(embed_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x_orig = x

        x = self.attn_norm(x)
        # attn
        x = self.attn(x, x, x)[0]
        x_post_attn = x_orig + x

        x = self.lin_norm(x_post_attn)
        x = self.lin(x)

        return x + x_post_attn


class CheckAMGPST(pst.BaseProteinSetTransformerEncoder[ModelConfig]):
    def __init__(self, config: ModelConfig):
        super().__init__(config)

        embedding_out_dim = self.config.out_dim

        self.global_attn = torch.nn.Sequential(
            *[
                Attention(
                    embed_dim=embedding_out_dim,
                    num_heads=self.config.n_attn_heads,
                    batch_first=True,
                )
                for _ in range(self.config.n_attn_layers)
            ]
        )

        self.final_norm = torch.nn.LayerNorm(embedding_out_dim, bias=False)

        # number of batches skipped for each epoch
        self._skipped_batches: dict[str, defaultdict[int, int]] = {
            "train": defaultdict(int),
            "val": defaultdict(int),
        }

    def setup_objective(
        self,
        margin: float,
        max_ap_pairs: int | None,
        positive_mining_strategy: str,
        negative_mining_strategy: str,
        context_size: int,
        **kwargs,
    ):
        tl_loss_module = TripletLoss(
            margin=margin,
            max_ap_pairs=max_ap_pairs,
            positive_mining_strategy=positive_mining_strategy,  # type: ignore
            negative_mining_strategy=negative_mining_strategy,  # type: ignore
        )

        avg_views_tl_loss_model = AVGViewsTripletLoss(
            context_size=context_size,
            margin=margin,
            pst_encoder=self,
        )

        return tl_loss_module, avg_views_tl_loss_model

    def concatenated_embeddings(self, batch: pst.GenomeGraphBatch) -> torch.Tensor:
        return self.internal_embeddings(batch)[0]

    def databatch_forward(
        self, batch: pst.GenomeGraphBatch, x: OptTensor = None
    ) -> torch.Tensor:
        x_embed = super().databatch_forward(batch=batch, x=x).out

        for attn in self.global_attn:
            x_embed = attn(x_embed)

        return self.final_norm(x_embed)

    def forward(self, batch: pst.GenomeGraphBatch, **kwargs) -> torch.Tensor:
        x_cat = self.concatenated_embeddings(batch)
        x_embed = self.databatch_forward(batch=batch, x=x_cat)
        return x_embed

    def _knn_classify(
        self, embed_dist: torch.Tensor, labels: torch.Tensor
    ) -> torch.Tensor:
        # add 1 and subset first hit since that will be self
        k = min(self.config.knn + 1, embed_dist.shape[0])
        closest_idx = embed_dist.topk(k=k, largest=False).indices[:, 1:]

        closest_labels = labels[closest_idx]

        # majority voting
        return torch.mode(closest_labels, -1).values

    def _train_val_step(
        self, batch: pst.GenomeGraphBatch, stage: str, **kwargs
    ) -> dict[str, torch.Tensor] | None:
        y_target: torch.Tensor = batch.y  # type: ignore

        if len(y_target.unique()) == 1:
            # cannot do triplet loss since we need more than 1 class
            self._skipped_batches[stage][self.current_epoch] += 1
            return None

        x_embed: torch.Tensor = self(batch)

        embed_dist = torch.cdist(x_embed, x_embed)

        tl_loss_module: TripletLoss
        avg_views_tl_loss_module: AVGViewsTripletLoss
        tl_loss_module, avg_views_tl_loss_module = self.criterion  # type: ignore

        tl_loss: torch.Tensor = tl_loss_module(
            class_labels=y_target,
            dist=embed_dist,
            class_weights=batch.weight[0],  # type: ignore
        )

        avg_tl_loss: OptTensor
        if torch.any(y_target == TRUE_AVG_IDX):
            avg_tl_loss = avg_views_tl_loss_module(
                batch=batch,
                x_embed=x_embed,
                embed_dist=embed_dist,
            )
        else:
            avg_tl_loss = None

        y_pred = self._knn_classify(embed_dist.detach(), y_target)

        # this is a balanced accuracy that averages recall across each class
        accuracy = multiclass_accuracy(
            preds=y_pred,
            target=y_target,
            num_classes=NUM_CLASSES,
            average="macro",
        )

        results = {}

        if avg_tl_loss is not None:
            loss: torch.Tensor = tl_loss + avg_tl_loss
            results["loss"] = loss
            results["tl_loss"] = tl_loss
            results["AVG_loss"] = avg_tl_loss
        else:
            loss = results["loss"] = results["tl_loss"] = tl_loss

        results["acc"] = accuracy

        batch_size = batch.num_proteins.numel()
        self.log_dict(
            {f"{stage}_{k}": v for k, v in results.items()},
            prog_bar=True,
            logger=True,
            on_step=self.config.verbose,
            on_epoch=True,
            batch_size=batch_size,
            sync_dist=True,
        )

        return results

    def training_step(self, batch: pst.GenomeGraphBatch, **kwargs) -> OptTensor:
        results = self._train_val_step(batch, stage="train", **kwargs)

        if results is not None:
            return results["loss"]

        return results

    def validation_step(self, batch: pst.GenomeGraphBatch, **kwargs) -> OptTensor:
        results = self._train_val_step(batch, stage="val", **kwargs)

        if results is not None:
            return results["loss"]

        return results

    @torch.no_grad()
    def predict_step(self, batch: pst.GenomeGraphBatch, **kwargs) -> torch.Tensor:
        self.eval()

        x: torch.Tensor = self(batch)
        embed_dist = torch.cdist(x, x)
        return self._knn_classify(embed_dist, batch.y)  # type: ignore
