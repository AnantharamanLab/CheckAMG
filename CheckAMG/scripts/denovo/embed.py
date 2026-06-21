import pst
import tables as tb
import torch
from pst.typing import FilePath
from torch.utils.data import Subset
from tqdm import tqdm

from CheckAMG.scripts.denovo.cli import EmbedArgs
from CheckAMG.scripts.denovo.constants import EMBED_FIELD, LABEL_FIELD
from CheckAMG.scripts.denovo.data import GenomeDataModule
from CheckAMG.scripts.denovo.model import CheckAMGPST
from CheckAMG.scripts.denovo.utils import load_model

FILTERS = tb.Filters(complevel=4, complib="blosc2:lz4hc")  # type: ignore

_DatasetT = Subset[pst.GenomeGraphBatch] | pst.LazyGenomeDataset | pst.GenomeDataset


def _peek_num_proteins(dataset: _DatasetT) -> int:
    if isinstance(dataset, Subset):
        return _peek_num_proteins(dataset.dataset)  # type: ignore

    return dataset.num_proteins


@torch.no_grad()
def embed(
    model: CheckAMGPST,
    dataset: _DatasetT,
    output: FilePath,
    batch_size: int = 1,
    verbose: bool = True,
):
    with tb.open_file(
        output,  # type: ignore
        "w",
        filters=FILTERS,
    ) as fp:
        num_proteins = _peek_num_proteins(dataset)
        out_dim = model.config.out_dim

        storage = fp.create_earray(
            where=fp.root,
            name=EMBED_FIELD,
            shape=(0, out_dim),
            expectedrows=num_proteins,
            atom=tb.Float32Atom(),
        )

        # label storage is created on first batch if labels are present
        label_storage = None

        dataloader = pst.ScaffoldDataLoader(
            dataset=dataset,  # type: ignore
            batch_size=batch_size,
            shuffle=False,
        )

        batch: pst.GenomeGraphBatch
        for batch in tqdm(
            dataloader, disable=not verbose, desc="Embedding proteins"
        ):
            batch = batch.to(model.device)  # type: ignore
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=model.device.type == "cuda"):
                x: torch.Tensor = model(batch=batch)

            if x.ndim == 1:
                x = x.unsqueeze(0)
            storage.append(x.float().cpu().numpy()) # cast back to fp32 before writing

            y = batch.y # type: ignore
            if y is None:
                continue

            if label_storage is None:
                label_storage = fp.create_earray(
                    where=fp.root,
                    name=LABEL_FIELD,
                    shape=(0,),
                    expectedrows=num_proteins,
                    atom=tb.Int64Atom(),
                )

            if y.ndim > 1:
                y = y.squeeze()

            label_storage.append(y.cpu().numpy())


def _main():
    args = EmbedArgs.parse_from_cli()

    model = load_model(args.model_ckpt)
    datamodule = GenomeDataModule(test_file=args.query_file)

    datamodule.setup()

    embed(
        model=model,
        dataset=datamodule.test_dataset,  # type: ignore
        output=args.output_file,
        batch_size=args.batch_size,
        verbose=args.verbose,
    )

    datamodule.close()


if __name__ == "__main__":
    _main()