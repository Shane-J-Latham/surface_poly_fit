"""
Fit polynomial surface to mesh vertices.


.. autosummary::
   :toctree: generated/

   cli - Command line interface.
   core - Polynomial fitting classes and functions.
   results_to_mesh - Command line interface for converting :samp:`".npz"` to :samp:`".vtu"`.
   tests - :mod:`surface_poly_fit` unit-tests.

"""
from . import _version
from .core import PolyhedralSurface, MongeJetFitter  # noqa: F401

__version__ = _version.get_versions()['version']

__all__ = [s for s in dir() if not s.startswith('_')]
