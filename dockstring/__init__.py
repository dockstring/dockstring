try:
    # Python >=3.8
    from importlib.metadata import PackageNotFoundError, version  # pyright: ignore[reportMissingImports]
except ModuleNotFoundError:
    # Python 3.7
    from importlib_metadata import PackageNotFoundError, version  # type: ignore[no-redef]

from .errors import (
    CanonicalizationError,
    DockingError,
    DockstringError,
    EmbeddingError,
    FormatConversionError,
    OutputError,
    ParsingError,
    PoseProcessingError,
    ProtonationError,
    SanityError,
    StructureOptimizationError,
    VinaError,
)
from .target import list_all_target_names, load_target
from .utils import setup_logger

__all__ = [
    'load_target', 'list_all_target_names', 'setup_logger', 'DockstringError', 'CanonicalizationError', 'ParsingError',
    'OutputError', 'SanityError', 'EmbeddingError', 'StructureOptimizationError', 'FormatConversionError',
    'ProtonationError', 'PoseProcessingError', 'VinaError', 'DockingError'
]

try:
    __version__ = version("dockstring")
except PackageNotFoundError:
    # package is not installed
    pass
