from typing import TYPE_CHECKING, Literal, cast

import pst
import torch
from pst.typing import OptTensor, PairTensor
from torch_geometric.data import Batch

from CheckAMG.scripts.pst_checkamg.anchored_view import AnchoredView
from CheckAMG.scripts.pst_checkamg.constants import TRUE_AVG_IDX

if TYPE_CHECKING:
    from CheckAMG.scripts.pst_checkamg.model import CheckAMGPST


class TripletLoss(torch.nn.Module):
    def __init__(
        self,
        margin: float,
        max_ap_pairs: int | None = None,
        positive_mining_strategy: Literal["all", "hard"] = "all",
        negative_mining_strategy: Literal["hard", "semihard", "both"] = "semihard",
    ):
        super().__init__()

        self.margin = margin
        self.max_ap_pairs = max_ap_pairs
        self.negative_mining_strat = negative_mining_strategy
        self.positive_mining_strat = positive_mining_strategy

    def _get_same_label(self, class_labels: torch.Tensor) -> torch.Tensor:
        return class_labels.unsqueeze(-2) == class_labels.unsqueeze(-1)

    def _get_all_positive_pairs(self, same_label: torch.Tensor) -> PairTensor:
        all_anchor_idx, all_pos_idx = (
            # get the upper triangle only, minus diagonal (ie ignore self-self triplets)
            torch.triu(same_label, diagonal=1).nonzero(as_tuple=True)
        )

        return all_anchor_idx, all_pos_idx

    def _negative_mining_condition(
        self, triplet_margin: torch.Tensor, same_label: torch.Tensor
    ) -> torch.Tensor:
        # triplet_margin = ap_dist - an_dist
        # thus, hard negatives will be closer, so an_dist < ap_dist
        # so hard negatives are: an_dist < ap_dist
        # semihard negatives are ap_dist < an_dist < ap_dist + margin

        diff_label = ~same_label
        lte_margin = triplet_margin > -self.margin
        if self.negative_mining_strat == "both":
            return lte_margin & diff_label

        # NOTE: i think the conditions were backwards....
        if self.negative_mining_strat == "hard":
            return lte_margin & (triplet_margin > 0.0) & diff_label

        # semihard
        return lte_margin & (triplet_margin <= 0.0) & diff_label

    def _mine_negative(
        self, triplet_margin: torch.Tensor, same_label: torch.Tensor
    ) -> torch.Tensor:
        # shape: [AP-pairs, num_proteins]
        meets_condition = self._negative_mining_condition(triplet_margin, same_label)

        masked_triplet_margin = triplet_margin.masked_fill(
            ~meets_condition, -torch.inf
        )

        # just choose any negative from a different class
        no_neg = meets_condition.any(-1).logical_not()
        if no_neg.any():
            mask = no_neg.unsqueeze(-1) & (~same_label)
            masked_triplet_margin.masked_fill_(mask, 1.0)

        scaled_masked_triplet_margin = (
            masked_triplet_margin - masked_triplet_margin.amax(-1, keepdim=True)
        )

        # the highest weights will be the largest margins
        # for hard negatives, triplet_margin > 0
        # so prioritize the hardest negative (ie greatest margin)
        # for semihard negatives, -self.margin > triplet_margin > 0
        # but we still want the value closest to 0 (ie the max), so this should still be fine
        sampling_probs = scaled_masked_triplet_margin.softmax(-1)

        # shape: [AP-pairs]
        # if there's only 1 pair, then make sure this is atleast_1d
        return torch.atleast_1d(sampling_probs.multinomial(1).squeeze())

    @torch.no_grad()
    def _hard_positive_mining(
        self, dist: torch.Tensor, same_label: torch.Tensor
    ) -> PairTensor:
        masked_dist = (
            # mask other classes
            dist.masked_fill(~same_label, -torch.inf)
            # mask self
            .fill_diagonal_(-torch.inf)
        )

        anchor_idx = torch.arange(dist.shape[0], device=dist.device)
        pos_idx = masked_dist.argmax(-1)

        has_pos_pair = torch.all(masked_dist == -torch.inf, -1).logical_not()

        return anchor_idx[has_pos_pair], pos_idx[has_pos_pair]

    @torch.no_grad()
    def _mine_triplets(
        self,
        dist: torch.Tensor,
        all_anchor_idx: torch.Tensor,
        all_pos_idx: torch.Tensor,
        same_label: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        num_pairs = all_anchor_idx.shape[0]
        if self.max_ap_pairs is None:
            subset_size = num_pairs
        else:
            subset_size = min(self.max_ap_pairs, num_pairs)

        tl_indices: dict[str, list[torch.Tensor]] = {
            "anchor": [],
            "pos": [],
            "neg": [],
        }
        for start_idx in range(0, num_pairs, subset_size):
            end_idx = min(start_idx + subset_size, num_pairs)

            subset_anchor_idx = all_anchor_idx[start_idx:end_idx]
            subset_pos_idx = all_pos_idx[start_idx:end_idx]

            ap_dist = dist[subset_anchor_idx, subset_pos_idx]
            an_dist = dist[subset_anchor_idx]

            triplet_margin = ap_dist.unsqueeze(-1) - an_dist
            subset_neg_idx = self._mine_negative(
                triplet_margin=triplet_margin,
                same_label=same_label[subset_anchor_idx],
            )
            tl_indices["anchor"].append(subset_anchor_idx)
            tl_indices["pos"].append(subset_pos_idx)
            tl_indices["neg"].append(subset_neg_idx)

        cat_tl_indices: dict[str, torch.Tensor] = {
            k: torch.cat(v) for k, v in tl_indices.items()
        }

        return cat_tl_indices

    def _validate_indices(
        self, triplet_indices: dict[str, torch.Tensor], class_labels: torch.Tensor
    ):
        anchor_idx = triplet_indices["anchor"]
        pos_idx = triplet_indices["pos"]
        neg_idx = triplet_indices["neg"]

        if torch.any(anchor_idx == pos_idx):
            raise RuntimeError(
                f"Anchor is positive at: {torch.where(anchor_idx == pos_idx)[0]}"
            )

        if torch.any(anchor_idx == neg_idx):
            where = torch.where(anchor_idx == neg_idx)[0]

            raise RuntimeError(
                f"Anchor is negative at: {where} {anchor_idx[where]} {neg_idx[where]}"
            )

        y_anchor = class_labels[anchor_idx]
        y_pos = class_labels[pos_idx]
        y_neg = class_labels[neg_idx]

        if not torch.all(y_anchor == y_pos):
            raise RuntimeError(
                "Not all anchor and positive pairs belong to the same class"
            )

        if torch.any(y_anchor == y_neg):
            raise RuntimeError("Some negatives are in the same class as the anchor")

    def forward(
        self,
        class_labels: torch.Tensor,
        x: OptTensor = None,
        dist: OptTensor = None,
        class_weights: OptTensor = None,
    ) -> torch.Tensor:
        # should skip loss computation if only 1 class
        if class_labels.unique().numel() == 1:
            raise ValueError(
                "Cannot compute triplet loss with members of only 1 class"
            )

        if dist is None:
            if x is not None:
                dist = torch.cdist(x, x)
            else:
                raise ValueError(
                    "Either a embedding distance matrix OR the embeddings must be passed."
                )

        same_label = self._get_same_label(class_labels)
        if self.positive_mining_strat == "hard":
            all_anchor_idx, all_pos_idx = self._hard_positive_mining(
                dist, same_label
            )
        else:
            all_anchor_idx, all_pos_idx = self._get_all_positive_pairs(same_label)

        triplet_indices = self._mine_triplets(
            dist=dist,
            all_anchor_idx=all_anchor_idx,
            all_pos_idx=all_pos_idx,
            same_label=same_label,
        )

        self._validate_indices(triplet_indices, class_labels)

        anchor_pos_dist = dist[triplet_indices["anchor"], triplet_indices["pos"]]
        anchor_neg_dist = dist[triplet_indices["anchor"], triplet_indices["neg"]]

        loss = anchor_pos_dist + self.margin - anchor_neg_dist

        if class_weights is not None:
            # class_weights shape: [num_classes]
            # class weights are just a global inverse frequency
            # theoretically we could also scale this based on the number of pairs (ie quadratically instead of linearly)
            anchor_labels = class_labels[triplet_indices["anchor"]]
            per_elem_weight = class_weights[anchor_labels]
            loss = loss * per_elem_weight

        return loss.relu().mean()


class AVGViewsTripletLoss(torch.nn.Module):
    def __init__(
        self,
        context_size: int,
        margin: float,
        pst_encoder: "CheckAMGPST",
    ):
        super().__init__()
        self.context_size = context_size
        self.margin = margin
        self.encoder = pst_encoder

    @torch.no_grad()
    def _mine_negatives(
        self,
        protein_is_avg: torch.Tensor,
        avg_idx: torch.Tensor,
        embed_dist: torch.Tensor,
    ) -> torch.Tensor:
        # shape: [num_AVGs, num_proteins]
        # need to choose negatives without autograd so that there is only one autograd
        # fn at the end for indexing
        avg_dist = embed_dist[avg_idx]

        neg_sampling_weights = (
            avg_dist
            # mask dist to all other AVGs
            .masked_fill(
                mask=protein_is_avg.repeat((avg_idx.numel(), 1)),
                value=torch.inf,
                # value=0.0,
            )
            # mul -1 so that we choose the closest ones more frequently
            .mul(-1)
        )

        neg_sampling_weights = neg_sampling_weights - neg_sampling_weights.amax(
            -1, keepdim=True
        )

        neg_sampling_weights = neg_sampling_weights.softmax(-1)

        # shape: [num_avgs]
        neg_idx = neg_sampling_weights.multinomial(num_samples=1).squeeze()
        return torch.atleast_1d(neg_idx)

    def _create_view_graphs(
        self, batch: pst.GenomeGraphBatch, avg_idx: torch.Tensor
    ) -> pst.GenomeGraphBatch:
        view_graphs: list[pst.GenomeGraph] = []
        device = batch.x.device
        for curr_avg_idx in avg_idx:
            curr_avg_idx = int(curr_avg_idx)

            scaffold_idx = int(batch.batch[curr_avg_idx])
            anchor_idx = int(curr_avg_idx - batch.ptr[scaffold_idx])
            nptns = int(batch.num_proteins[scaffold_idx])
            view = AnchoredView(
                anchor_idx=anchor_idx,
                num_proteins=nptns,
                context_size=self.context_size,
            )

            graphslice = slice(
                int(batch.ptr[scaffold_idx]), int(batch.ptr[scaffold_idx + 1])
            )

            if view.total_views > 0:
                for edge_index in view.all_anchored_views():
                    unique: torch.Tensor
                    rel_edge_index: torch.Tensor
                    unique, rel_edge_index = edge_index.to(device).unique(
                        sorted=True, return_inverse=True
                    )

                    num_ptns = unique.numel()

                    x = batch.x[graphslice][unique]
                    strand = batch.strand[graphslice][unique]
                    pos = batch.pos[graphslice][unique]
                    y = batch.y[graphslice][unique]  # type: ignore
                    scaffold_label = int(batch.scaffold_label[scaffold_idx])
                    genome_label = int(batch.genome_label[scaffold_idx])
                    is_anchor = unique == view.anchor_idx

                    # this makes a mini "scaffold" that only contains the proteins in the view
                    minigraph = pst.GenomeGraph(
                        x=x,
                        strand=strand,
                        num_proteins=num_ptns,
                        pos=pos,
                        scaffold_label=scaffold_label,
                        genome_label=genome_label,
                        edge_index=rel_edge_index,
                        avg_idx=curr_avg_idx,
                        y=y,
                        is_anchor=is_anchor,
                    )

                    view_graphs.append(minigraph)

        avgs_view_batch: Batch = Batch.from_data_list(view_graphs).to(device)  # type: ignore

        return cast(pst.GenomeGraphBatch, avgs_view_batch)

    def _filter_scaffolds_with_fewer_proteins_than_context(
        self, batch: pst.GenomeGraphBatch, avg_idx: torch.Tensor
    ) -> torch.Tensor:
        # these are the scaffold indices for all AVGs
        scaffold_has_avg_indices = batch.batch[avg_idx]

        num_proteins_per_scaffolds = batch.num_proteins[scaffold_has_avg_indices]

        mask = num_proteins_per_scaffolds > self.context_size
        return avg_idx[mask]

    def _embed_views(self, avg_batch: pst.GenomeGraphBatch) -> torch.Tensor:
        all_embed: torch.Tensor = self.encoder(batch=avg_batch)

        return all_embed[avg_batch.is_anchor]  # type: ignore

    def _embed(
        self, batch: pst.GenomeGraphBatch, avg_idx: torch.Tensor
    ) -> PairTensor:
        avg_batch = self._create_view_graphs(batch=batch, avg_idx=avg_idx)
        num_views: torch.Tensor = avg_batch.avg_idx.unique(return_counts=True)[  # type: ignore
            1
        ]

        avg_view_embed = self._embed_views(avg_batch=avg_batch)

        return avg_view_embed, num_views

    def forward(
        self,
        batch: pst.GenomeGraphBatch,
        x_embed: torch.Tensor,
        embed_dist: OptTensor = None,
    ) -> OptTensor:
        if embed_dist is None:
            embed_dist = torch.cdist(x_embed, x_embed)

        y: torch.Tensor = batch.y  # type: ignore
        protein_is_avg = y == TRUE_AVG_IDX
        avg_idx: torch.Tensor = torch.atleast_1d(protein_is_avg.nonzero().squeeze())

        # we need to ignore AVGs whose scaffold encodes a number of ptns <= self.context_size
        # since in those cases, the short ctx view will be identical long ctx view and there will
        # be nothing to optimize
        avg_idx = self._filter_scaffolds_with_fewer_proteins_than_context(
            batch=batch, avg_idx=avg_idx
        )

        if len(avg_idx) == 0:
            return None

        avg_neg_idx = self._mine_negatives(
            protein_is_avg=protein_is_avg, avg_idx=avg_idx, embed_dist=embed_dist
        )

        avg_views_embed, num_views = self._embed(batch=batch, avg_idx=avg_idx)

        # repeat to match the number of views per AVG
        avg_idx = avg_idx.repeat_interleave(num_views)
        avg_neg_idx = avg_neg_idx.repeat_interleave(num_views)

        avg_embed = x_embed[avg_idx]

        ap_dist = torch.square(avg_embed - avg_views_embed).sum(-1).sqrt()
        an_dist = embed_dist[avg_idx, avg_neg_idx]
        loss = ap_dist - an_dist + self.margin

        return loss.relu().mean()
