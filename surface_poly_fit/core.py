"""
Polynomial fitting classes and functions.


.. autosummary::
   :toctree: generated/

"""
from ._spf_cgal import PolyhedralSurface as _PolyhedralSurface   # noqa: F401
from ._spf_cgal import MongeJetFitter as _MongeJetFitter


class PolyhedralSurface(_PolyhedralSurface):
    pass


class MongeJetFitter(_MongeJetFitter):
    pass


__all__ = [s for s in dir() if not s.startswith('_')]
