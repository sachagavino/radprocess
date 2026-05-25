"""
Shared fixtures and import setup for radprocess tests.

pymses is an optional dependency not available in CI.
We mock it here before any test module imports radprocess.
"""
import sys
import types


def _make_mock(name):
    mod = types.ModuleType(name)
    mod.__path__ = []  # make it a package
    return mod


# Create mock pymses hierarchy
_pymses = _make_mock("pymses")
_pymses.RamsesOutput = None  # attribute used in Convert.py

_utils = _make_mock("pymses.utils")
_utils.constants = types.ModuleType("pymses.utils.constants")
_pymses.utils = _utils

_filters = _make_mock("pymses.filters")
_filters.CellsToPoints = None
_pymses.filters = _filters

_sources = _make_mock("pymses.sources")
_sources_ramses = _make_mock("pymses.sources.ramses")
_sources_ramses.output = _make_mock("pymses.sources.ramses.output")
_sources_ramses.filename_utils = _make_mock("pymses.sources.ramses.filename_utils")
_sources.ramses = _sources_ramses
_pymses.sources = _sources

_analysis = _make_mock("pymses.analysis")
_analysis.visualization = _make_mock("pymses.analysis.visualization")
_pymses.analysis = _analysis

# Register all in sys.modules
for name, mod in [
    ("pymses", _pymses),
    ("pymses.utils", _utils),
    ("pymses.utils.constants", _utils.constants),
    ("pymses.filters", _filters),
    ("pymses.sources", _sources),
    ("pymses.sources.ramses", _sources_ramses),
    ("pymses.sources.ramses.output", _sources_ramses.output),
    ("pymses.sources.ramses.filename_utils", _sources_ramses.filename_utils),
    ("pymses.analysis", _analysis),
    ("pymses.analysis.visualization", _analysis.visualization),
]:
    sys.modules[name] = mod