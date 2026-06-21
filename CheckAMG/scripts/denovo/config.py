import attrs
import pst
from attrs import validators
from pst.utils.attrs.validators import positive_float, positive_int

_MARGIN_FIELD = attrs.field(default=1.0, validator=positive_float)

_KNN_FIELD = attrs.field(default=5, validator=positive_int)

_POSITIVE_MINING_FIELD = attrs.field(default="esm_easy")
_NEGATIVE_MINING_FIELD = attrs.field(default="semihard")
_OPPOSITE_CLASS_NEGATIVES_FIELD = attrs.field(default=True)


def divisor_of_n_validator(n: int):
    def _inner(instance, attribute: attrs.Attribute, value: int):
        if n % value != 0:
            raise ValueError(f"Value ({value}) must be a divisor of {n}")

    return _inner


_N_ATTN_HEADS_FIELD = attrs.field(
    default=8,
    validator=validators.and_(positive_int, divisor_of_n_validator(n=800)),
)
_CONTEXT_SIZE_FIELD = attrs.field(
    default=4, validator=validators.and_(positive_int, validators.le(15))
)


def optional_positive_int(instance, attribute: attrs.Attribute, value: int | None):
    if value is not None and value <= 0:
        raise ValueError(
            f"If value is not None, then it must be a positive integer. Got: {value}"
        )


_MAX_AP_PAIRS_FIELD = attrs.field(
    default=None,
    validator=optional_positive_int,
)


@attrs.define
class LossConfig(pst.BaseLossConfig):
    margin: float = _MARGIN_FIELD
    max_ap_pairs: int | None = _MAX_AP_PAIRS_FIELD
    positive_mining_strategy: str = _POSITIVE_MINING_FIELD
    negative_mining_strategy: str = _NEGATIVE_MINING_FIELD
    opposite_class_negatives: bool = _OPPOSITE_CLASS_NEGATIVES_FIELD
    context_size: int = _CONTEXT_SIZE_FIELD


@attrs.define
class ModelConfig(pst.BaseModelConfig):
    knn: int = _KNN_FIELD
    verbose: bool = attrs.field(default=False)
    loss: LossConfig = attrs.field(factory=LossConfig)
