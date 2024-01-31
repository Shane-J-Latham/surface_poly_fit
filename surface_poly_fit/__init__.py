"""
Fit polynomial surface to mesh vertices.
"""
from . import _version
from ._spf_cgal import PolyhedralSurface, MongeJetFitter  # noqa: F401

__version__ = _version.get_versions()['version']
