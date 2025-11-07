# This file is intentionally minimal.
# The actual module is built from Rust using PyO3.
# When installed, the compiled module will be available as 'orphos'.

from .orphos import (
    OrphosOptions,
    OrphosResult,
    analyze_sequence,
    analyze_file,
    __version__,
)

__all__ = [
    "OrphosOptions",
    "OrphosResult",
    "analyze_sequence",
    "analyze_file",
    "__version__",
]
