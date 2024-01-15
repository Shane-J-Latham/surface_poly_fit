"""
Tests for :mod:`surface_poly_fit`.
"""
import numpy as _np
import unittest as _unittest


# have_trimesh = False
# try:
#     import trimesh as _trimesh
#     have_trimesh = True
# except Exception:
#     pass


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

    def test_surface_poly_fit_spf_cgal_polyhedral_surface(self):
        from surface_poly_fit import _spf_cgal

        self.assertTrue(hasattr(_spf_cgal, "PolyhedralSurface"))
        self.assertIsNotNone(_spf_cgal.PolyhedralSurface())


class PolyhedralSurfaceTest(SurfacePolyFitTest):

    def test_construct(self):
        from trimesh.primitives import Capsule
        from surface_poly_fit._spf_cgal import PolyhedralSurface

        self.assertIsNotNone(PolyhedralSurface())
        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)
        self.assertEqual(len(trimesh_mesh.vertices), poly_surf.num_vertices)
        self.assertEqual(len(trimesh_mesh.faces), poly_surf.num_faces)


__all__ = [s for s in dir() if not s.startswith('_')]

_unittest.main(__name__)
