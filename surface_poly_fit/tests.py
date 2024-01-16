"""
Tests for :mod:`surface_poly_fit`.
"""
import numpy as _np
import unittest as _unittest
import logging as _logging


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
        # self.logger = logging.getLogger(__name__ + "." + self.__class__.__name__)
        self.logger = _logging.getLogger()


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

        self.logger.info("")
        self.logger.info("len(poly surf faces) = %s", len(poly_surf.get_faces()))
        self.logger.info("poly surf faces = %s", poly_surf.get_faces())
        self.assertSequenceEqual(
            trimesh_mesh.vertices.tolist(),
            poly_surf.get_vertices().tolist()
        )
        poly_surf_faces = poly_surf.get_faces()
        self.assertSequenceEqual(
            trimesh_mesh.faces.tolist(),
            poly_surf_faces.tolist() if hasattr(poly_surf_faces, "tolist") else poly_surf_faces
        )


__all__ = [s for s in dir() if not s.startswith('_')]

_logging.basicConfig()
_unittest.main(__name__)
