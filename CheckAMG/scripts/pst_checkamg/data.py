from typing import Literal

import numpy as np
import pst
import torch
from pst.typing import EdgeIndexStrategy, FilePath, KwargType
from torch.utils.data import Subset
from torch.utils.data.sampler import WeightedRandomSampler

from CheckAMG.scripts.pst_checkamg.constants import (
    CHUNK_SIZE,
    EDGE_STRATEGY,
    FRAGMENT_SIZE,
    THRESHOLD,
)


def register_class_label_and_weights(
    dataset: pst.LazyGenomeDataset, class_weighting: str = "class_freq"
):
    protein_is_AVG = dataset.get_registered_feature("protein_is_AVG")
    protein_is_viral = dataset.get_registered_feature("protein_is_viral")

    condlist = [
        (protein_is_AVG & protein_is_viral),  # 3 True True
        (protein_is_AVG & ~protein_is_viral),  # 2 True False
        (~protein_is_AVG & protein_is_viral),  # 1 False True
        # (~protein_is_AVG & ~protein_is_viral) # 0 False False
    ]

    choicelist = [3, 2, 1]
    default = 0

    # 0 = nonAMG/nonviral
    # 1 = nonAMG/viral
    # 2 = AMG/nonviral
    # 3 = AMG/viral

    label = torch.from_numpy(np.select(condlist, choicelist, default=default))

    dataset.register_feature(
        name="y",
        data=label,
        feature_level="protein",
    )

    counts = label.bincount()

    if class_weighting == "class_freq":
        weights = counts.sum() / counts
        weights /= weights.amin()
    elif class_weighting == "pair_freq":
        pairs = (counts * (counts) - 1) // 2
        weights = pairs.sum() / pairs
        weights /= weights.amin()
    else:
        raise ValueError(
            f"Invalid class_weighting option: {class_weighting}. Should be one of the following: class_freq, pair_freq"
        )

    # now this is a hack to register weights with the dataset
    # just expand it be the shape of the dataset, and it will always be part
    # of the minibatches
    weights = weights.unsqueeze(0).repeat((len(dataset), 1)).unsqueeze(-2)
    dataset.register_feature(name="weight", data=weights, feature_level="scaffold")


def register_scaffold_target_category(dataset: pst.LazyGenomeDataset):
    scaffold_has_AVG = dataset.get_registered_feature("scaffold_has_AVG")
    scaffold_has_nonviral = dataset.get_registered_feature("scaffold_has_nonviral")

    condlist = [
        (scaffold_has_AVG & scaffold_has_nonviral),  # True True
        (scaffold_has_AVG & ~scaffold_has_nonviral),  # True False
        (~scaffold_has_AVG & scaffold_has_nonviral),  # False True
        # (~scaffold_has_AVG & ~scaffold_has_nonviral) # False False
    ]

    choicelist = [3, 2, 1]
    default = 0

    # 0 = nonAMG/nonviral
    # 1 = nonAMG/viral
    # 2 = AMG/nonviral
    # 3 = AMG/viral

    scaffold_category = torch.from_numpy(
        np.select(condlist, choicelist, default=default)
    )

    dataset.register_feature(
        "scaffold_category",
        scaffold_category,
        feature_level="scaffold",
    )


def get_dataset_hparams(ckpt: FilePath) -> KwargType:
    loaded_ckpt = torch.load(ckpt, map_location="cpu")
    data_hparams = loaded_ckpt["datamodule_hyper_parameters"]
    keys = ["edge_strategy", "threshold", "chunk_size", "fragment_size"]

    return {k: data_hparams[k] for k in keys}


_SubsetT = Subset[pst.GenomeGraphBatch]
_TrainValSubsets = tuple[_SubsetT, _SubsetT]


def train_val_split(dataset: pst.LazyGenomeDataset) -> _TrainValSubsets:
    scaffold_val_mask = dataset.get_registered_feature("scaffold_val_mask")

    val_idx = scaffold_val_mask.nonzero(as_tuple=True)[0]
    train_idx = scaffold_val_mask.logical_not().nonzero(as_tuple=True)[0]

    train_dataset = Subset(dataset, train_idx.tolist())
    val_dataset = Subset(dataset, val_idx.tolist())

    return train_dataset, val_dataset


def _get_sampling_weights(train_val_dataset: pst.LazyGenomeDataset) -> torch.Tensor:
    scaffold_category = train_val_dataset.get_registered_feature("scaffold_category")
    category_counts = scaffold_category.bincount()
    category_weights = 1.0 / category_counts
    scaffold_category_weights = category_weights[scaffold_category]

    scaffold_val_mask = train_val_dataset.get_registered_feature("scaffold_val_mask")

    return scaffold_category_weights[~scaffold_val_mask]


def get_training_data_sampler(
    train_dataset: _SubsetT, weights: torch.Tensor
) -> WeightedRandomSampler:
    return WeightedRandomSampler(
        weights=weights,  # type: ignore
        num_samples=len(train_dataset),
        replacement=True,
    )


# this is meant to mimic a lightning DataModule
class GenomeDataModule:
    def __init__(
        self,
        train_file: FilePath | None = None,
        test_file: FilePath | None = None,
        class_weighting: Literal["class_freq", "pair_freq"] = "class_freq",
        # these are the defaults from PST-TL-P__large.ckpt
        edge_strategy: EdgeIndexStrategy = EDGE_STRATEGY,
        chunk_size: int = CHUNK_SIZE,
        threshold: int = THRESHOLD,
        fragment_size: int = FRAGMENT_SIZE,
    ):
        if train_file is None and test_file is None:
            raise ValueError("At least 1 training or test file must be provided.")

        self.dataset_init_kwargs: KwargType = dict(
            edge_strategy=edge_strategy,
            chunk_size=chunk_size,
            threshold=threshold,
            fragment_size=fragment_size,
        )
        self.train_file = train_file
        self.test_file = test_file
        self.class_weighting = class_weighting

        self._test_file_is_inference = False

    def setup(self):
        if self.train_file is not None:
            self.train_val_dataset = pst.LazyGenomeDataset(
                file=self.train_file, **self.dataset_init_kwargs
            )

            print(f"weighting classes by: {self.class_weighting}")
            register_class_label_and_weights(
                self.train_val_dataset, class_weighting=self.class_weighting
            )
            print(
                f"class weights: {self.train_val_dataset.get_registered_feature('weight')[0, 0]}"
            )

            register_scaffold_target_category(self.train_val_dataset)

            self.train_dataset, self.val_dataset = train_val_split(
                self.train_val_dataset
            )

            train_sampling_weights = _get_sampling_weights(self.train_val_dataset)
            self.train_sampler = get_training_data_sampler(
                self.train_dataset, train_sampling_weights
            )

        else:
            self.train_val_dataset = None
            self.train_dataset = None
            self.val_dataset = None

        if self.test_file is not None:
            self.test_dataset = pst.LazyGenomeDataset(
                file=self.test_file, **self.dataset_init_kwargs
            )
            try:
                register_class_label_and_weights(
                    self.test_dataset, class_weighting=self.class_weighting
                )
            except KeyError:
                # for inference datasets, no labels are present by definition, so can't register these
                self._test_file_is_inference = True
        else:
            self.test_dataset = None

    def _check_dataset_exists_before_loading(self, stage: str):
        dataset: pst.LazyGenomeDataset | None = getattr(self, f"{stage}_dataset")

        if dataset is None:
            raise ValueError(
                f"Dataset at stage {stage} does not exist since the file to load this data was not provided upon initialization."
            )

    def _dataloader(
        self,
        stage: str,
        dataset: pst.LazyGenomeDataset | None,
        batch_size: int = 1,
        **kwargs,
    ) -> pst.ScaffoldDataLoader:
        self._check_dataset_exists_before_loading(stage)

        return pst.ScaffoldDataLoader(
            dataset=dataset,  # type: ignore
            batch_size=batch_size,
            **kwargs,
        )

    def train_dataloader(
        self, batch_size: int = 1, **kwargs
    ) -> pst.ScaffoldDataLoader:
        return self._dataloader(
            stage="train",
            dataset=self.train_dataset,  # type: ignore
            batch_size=batch_size,
            shuffle=None,
            sampler=self.train_sampler,
            **kwargs,
        )

    def val_dataloader(
        self, batch_size: int = 1, **kwargs
    ) -> pst.ScaffoldDataLoader:
        return self._dataloader(
            stage="val",
            dataset=self.val_dataset,  # type: ignore
            batch_size=batch_size,
            shuffle=True,
            **kwargs,
        )

    def test_dataloader(
        self, batch_size: int = 1, **kwargs
    ) -> pst.ScaffoldDataLoader:
        return self._dataloader(
            stage="test",
            dataset=self.test_dataset,  # type: ignore
            batch_size=batch_size,
            shuffle=False,
            **kwargs,
        )

    def close(self, stage: str | list[str] | None = None):
        if stage is None:  # close all
            stage = ["train", "test"]
        elif isinstance(stage, str):
            stage = [stage]

        for key in stage:
            if key not in {"train", "test"}:
                raise ValueError(
                    f"Invalid stage to close: {key}. Expecting only values from: ['train', 'test']"
                )

            file: str | None = getattr(self, f"{key}_file")
            if file is not None:
                dataset: pst.LazyGenomeDataset = getattr(
                    self,
                    f"{key}_dataset" if key == "test" else f"{key}_val_dataset",
                )
                dataset._file.close()
