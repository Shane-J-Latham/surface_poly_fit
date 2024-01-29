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

    def export_mesh(self, file_name, poly_surface):
        """
        """
        import trimesh
        from trimesh.exchange.export import export_mesh

        self.logger.debug("vertices=%s", poly_surface.get_vertices())
        self.logger.debug("faces=%s", poly_surface.get_faces())
        mesh = \
            trimesh.Trimesh(
                vertices=poly_surface.get_vertices(),
                faces=poly_surface.get_faces(),
                process=False,
                validate=False
            )
        export_mesh(mesh, file_name)


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

    def test_set_vertex_normals(self):
        from trimesh.primitives import Capsule
        from surface_poly_fit._spf_cgal import PolyhedralSurface

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)
        self.assertTrue(_np.any(poly_surf.get_vertex_normals() != trimesh_mesh.vertex_normals))
        poly_surf.set_vertex_normals(trimesh_mesh.vertex_normals)
        self.assertTrue(_np.all(poly_surf.get_vertex_normals() == trimesh_mesh.vertex_normals))

        self.assertRaises(Exception, poly_surf.set_vertex_normals, trimesh_mesh.vertex_normals[1:])

    def test_create_ring_patch(self):
        poly_surface = create_monge_surface()
        origin_vertex_index = poly_surface.num_vertices // 2
        self.assertTrue(
            _np.allclose(
                (0.0, 0.0, 0.0),
                poly_surface.get_vertices()[origin_vertex_index].tolist()
            )
        )
        # self.export_mesh("monge_surface.ply", poly_surface)

        patch_surface = poly_surface.create_ring_patch(origin_vertex_index, 1)
        # self.export_mesh("monge_surface_patch_ring1.ply", patch_surface)
        self.assertEqual(5, patch_surface.num_vertices)
        self.assertEqual(4, patch_surface.num_faces)

        patch_surface = poly_surface.create_ring_patch(origin_vertex_index, 2)
        # self.export_mesh("monge_surface_patch_ring2.ply", patch_surface)
        self.assertEqual(21, patch_surface.num_vertices)
        self.assertEqual(28, patch_surface.num_faces)

        patch_surface = poly_surface.create_ring_patch(origin_vertex_index, 3)
        # self.export_mesh("monge_surface_patch_ring3.ply", patch_surface)
        self.assertEqual(45, patch_surface.num_vertices)
        self.assertEqual(68, patch_surface.num_faces)

        patch_surface = poly_surface.create_ring_patch(origin_vertex_index, 4)
        # self.export_mesh("monge_surface_patch_ring4.ply", patch_surface)
        self.assertEqual(77, patch_surface.num_vertices)
        self.assertEqual(124, patch_surface.num_faces)


class MongePolynomial:
    def __init__(self, k, b, c):
        self._k = _np.asanyarray(k).copy()
        self._b = _np.asanyarray(b).copy()
        self._c = _np.asanyarray(c).copy()

    @property
    def k(self):
        return self._k

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    def evaluate(self, xy):
        xy = _np.asanyarray(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        k = self.k
        b = self.b
        c = self.c

        z = \
            (
                0.5 * (k[0] * (x**2) + k[1] * (y**2))
                +
                (1.0 / 6.0) * (
                    b[0] * (x**3) + 3 * b[1] * (x**2) * y + 3 * b[2] * x * (y**2) + b[3] * y**3
                )
                +
                (1.0 / 24.0) * (
                    c[0] * (x**4) + 4 * c[1] * (x**3) * y + 6 * c[2] * (x**2) * (y**2)
                    +
                    4 * c[3] * x * (y**3) + c[4] * (y**4)
                )
            )
        return z

    def __call__(self, xy):
        return self.evaluate(xy)


def create_monge_surface(monge_polynomial=None, xy=None):
    """
    """
    from scipy.spatial import Delaunay
    from surface_poly_fit._spf_cgal import PolyhedralSurface

    if monge_polynomial is None:
        monge_polynomial = \
            MongePolynomial(
                k=(1.0, 0.5),
                b=(0.45, -0.2, 0.4, -0.50),
                c=(-0.125, 0.100, 0.05, -0.08, 0.075)
            )
    if xy is None:
        xy = \
            _np.meshgrid(
                _np.linspace(-4.0, 4.0, 33),
                _np.linspace(-4.0, 4.0, 33),
                indexing="ij"
            )
    xy = _np.asarray([xy[0].flatten(), xy[1].flatten()]).T.copy()
    z = monge_polynomial(xy)
    z = z.reshape((len(z), 1))
    dlny = Delaunay(xy)
    faces = dlny.simplices
    vertices = _np.hstack((xy, z))

    return PolyhedralSurface(vertices=vertices, faces=faces)


class MongeJetFitterTest(SurfacePolyFitTest):

    def test_construct(self):
        from trimesh.primitives import Capsule
        from surface_poly_fit._spf_cgal import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        fitter = MongeJetFitter(poly_surf)
        self.assertIsNotNone(fitter)

        fitter = MongeJetFitter(poly_surf, 4)
        self.assertIsNotNone(fitter)

        fitter = MongeJetFitter(poly_surf, 4, 4)
        self.assertIsNotNone(fitter)

    def test_properties(self):
        from trimesh.primitives import Capsule
        from surface_poly_fit._spf_cgal import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        fitter = MongeJetFitter(poly_surf, 4, 2)
        self.assertEqual(4, fitter.degree_poly_fit)
        self.assertEqual(2, fitter.degree_monge)
        self.assertEqual(15, fitter.min_num_fit_points)

    def test_fit(self):
        # import trimesh
        from surface_poly_fit._spf_cgal import MongeJetFitter

        monge_polynomial = \
            MongePolynomial(
                k=(0.50, 0.25),
                b=(0.0, 0.0, 0.0, 0.0),
                c=(0.0, 0.0, 0.0, 0.0, 0.0)
            )
        poly_surface = create_monge_surface(monge_polynomial)
        self.logger.info("")
        self.logger.info(
            "poly_surface (min_z, max_z) = (%s, %s).",
            poly_surface.get_vertices()[:, 2].min(),
            poly_surface.get_vertices()[:, 2].max(),
        )
        # tmesh = \
        #    trimesh.Trimesh(vertices=poly_surface.get_vertices(), faces=poly_surface.get_faces())
        # trimesh.exchange.export.export_mesh(tmesh, "monge_poly_surface.ply")

        monge_origin_vertex_index = poly_surface.num_vertices // 2
        self.assertTrue(
            _np.allclose(
                (0.0, 0.0, 0.0),
                poly_surface.get_vertices()[monge_origin_vertex_index].tolist()
            )
        )
        degree_monge = 2
        degree_poly_fit = 2
        fitter = MongeJetFitter(poly_surface, degree_poly_fit, degree_monge)
        result = fitter.fit_at_vertex(monge_origin_vertex_index, num_rings=16)
        self.logger.info("Result = %s", result)
        self.assertEqual(monge_origin_vertex_index, result["vertex_index"][0])
        self.assertEqual(degree_monge, result["degree_monge"][0])
        self.assertEqual(degree_poly_fit, result["degree_poly_fit"][0])
        self.assertAlmostEqual(monge_polynomial.k[0], result["k"][0][0])
        self.assertAlmostEqual(monge_polynomial.k[1], result["k"][0][1])
        self.assertTrue(
            _np.allclose((0.0, 0.0, 0.0), result["origin"][0])
        )
        self.assertTrue(
            _np.allclose(_np.identity(3, dtype=_np.float64), _np.absolute(result["direction"][0]))
        )
        if degree_monge < 3:
            self.assertTrue(
                _np.all(result["b"][0] == 0.0)
            )
        if degree_monge < 4:
            self.assertTrue(
                _np.all(result["c"][0] == 0.0)
            )


__all__ = [s for s in dir() if not s.startswith('_')]

if __name__ == "__main__":
    _logging.basicConfig(level=_logging.DEBUG)
    _unittest.main(__name__)
