"""
Tests for :mod:`surface_poly_fit`.
"""
import numpy as _np
import unittest as _unittest


have_trimesh = False
try:
    import trimesh as _trimesh
    have_trimesh = True
except Exception:
    pass


class SurfacePolyFitTest(_unittest.TestCase):

    """
    Base class for polynomial surface fitting.
    """

    def setUp(self):
        """
        """
        _np.random.seed(54317953)


class SurfacePolyFitImportTest(SurfacePolyFitTest):

    def test_import_surface_poly_fit(self):
        import surface_poly_fit

        self.assertIsNotNone(surface_poly_fit)

    def test_import_surface_poly_fit_spf_cgal(self):
        from surface_poly_fit import _spf_cgal

        self.assertIsNotNone(_spf_cgal)


__all__ = [s for s in dir() if not s.startswith('_')]

_unittest.main(__name__)
