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
        self.logger.info("")
        self.logger.info("")
        self.logger.info("trimesh_mesh.face_normals:")
        for i in range(0, 10):
            self.logger.info(trimesh_mesh.face_normals.tolist()[i])
        self.logger.info("\npoly_surf.get_face_normals():")
        for i in range(0, 10):
            self.logger.info(poly_surf.get_face_normals().tolist()[i])
        self.assertTrue(
            _np.allclose(
                trimesh_mesh.face_normals,
                poly_surf.get_face_normals()
            )
        )

        self.logger.info("")
        self.logger.info("")
        self.logger.info("trimesh_mesh.vertex_normals:")
        for i in range(0, 10):
            self.logger.info(trimesh_mesh.vertex_normals.tolist()[i])
        self.logger.info("poly_surf.get_vertex_normals():")
        for i in range(0, 10):
            self.logger.info(poly_surf.get_vertex_normals().tolist()[i])

        nrml_angles = \
            _np.rad2deg(
                _np.arccos(
                    _np.sum(
                        trimesh_mesh.vertex_normals
                        *
                        poly_surf.get_vertex_normals(),
                        axis=1
                    )
                )
            )
        self.logger.info("")
        self.logger.info("nrml_angles:\n%s", nrml_angles.tolist()[0:10])
        self.logger.info("min nrml_angles:\n%s", nrml_angles.min())
        self.logger.info("max nrml_angles:\n%s", nrml_angles.max())
        self.logger.info("nrml_angles.shape = %s", nrml_angles.shape)

        self.assertTrue(
            _np.all(
                nrml_angles
                <
                5.0
            )
        )


__all__ = [s for s in dir() if not s.startswith('_')]

if __name__ == "__main__":
    _logging.basicConfig(level=_logging.DEBUG)
    _unittest.main(__name__)
