"""
Tests for :mod:`surface_poly_fit`.

.. autosummary::
   :toctree: generated/

"""
import numpy as _np
import unittest as _unittest
import logging as _logging


class SurfacePolyFitTest(_unittest.TestCase):

    """
    Base class for polynomial surface fitting.
    """

    def setUp(self):
        """
        Set-up, initialise logger.
        """
        self.logger = _logging.getLogger()

    def export_mesh(self, file_name, poly_surface):
        """
        Export a :obj:`surface_poly_fit.core.PolyhedralSurface` to file.

        :type file_name: :obj:`str`
        :param file_name: File path for mesh export.
        :type poly_surface: :obj:`surface_poly_fit.core.PolyhedralSurface`
        :param poly_surface: Mesh to export.
        """
        import meshio

        mio_mesh = poly_surface.as_meshio_mesh()
        meshio.write(file_name, mio_mesh)


class SurfacePolyFitImportTest(SurfacePolyFitTest):
    """
    Basic import tests, make sure the pybind11 extension can be imported.
    """

    def test_import_surface_poly_fit(self):
        """
        Test import of :mod:`surface_poly_fit` package.
        """
        import surface_poly_fit

        self.assertIsNotNone(surface_poly_fit)

    def test_import_surface_poly_fit_spf_cgal(self):
        """
        Test import of pybind11 extension.
        """
        from surface_poly_fit import _spf_cgal

        self.assertIsNotNone(_spf_cgal)

    def test_surface_poly_fit_spf_cgal_polyhedral_surface(self):
        """
        Test pybind11 class export (can construct an instance).
        """
        from surface_poly_fit import _spf_cgal

        self.assertTrue(hasattr(_spf_cgal, "PolyhedralSurface"))
        self.assertIsNotNone(_spf_cgal.PolyhedralSurface())


class PolyhedralSurfaceTest(SurfacePolyFitTest):
    """
    Tests for :obj:`surface_poly_fit.core.PolyhedralSurface`.
    """

    def test_construct(self):
        """
        Test default construction, construction with vertices and faces.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface

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
        """
        Test :meth:`surface_poly_fit.core.PolyhedralSurface.set_vertex_normals`
        and :meth:`surface_poly_fit.core.PolyhedralSurface.get_vertex_normals`.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)
        self.assertTrue(_np.any(poly_surf.get_vertex_normals() != trimesh_mesh.vertex_normals))
        poly_surf.set_vertex_normals(trimesh_mesh.vertex_normals)
        self.assertTrue(_np.all(poly_surf.get_vertex_normals() == trimesh_mesh.vertex_normals))

        self.assertRaises(Exception, poly_surf.set_vertex_normals, trimesh_mesh.vertex_normals[1:])

    def test_create_ring_patch(self):
        """
        Test :meth:`surface_poly_fit.core.PolyhedralSurface.create_ring_patch`.
        """
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

    def test_vertex_normals(self):
        """
        Test that :meth:`surface_poly_fit.core.PolyhedralSurface.get_vertex_normals`
        normals are valid.
        """
        poly_surface = create_monge_surface()
        poly_surface_nrmls = poly_surface.get_vertex_normals()
        for nrml in poly_surface_nrmls:
            self.logger.debug("nrml=%s", nrml.tolist())
            self.assertTrue(_np.all(~_np.isnan(nrml)))
            self.assertAlmostEqual(1.0, _np.linalg.norm(nrml))

    def test_face_normals(self):
        """
        Test that :meth:`surface_poly_fit.core.PolyhedralSurface.get_face_normals`
        normals are valid.
        """
        poly_surface = create_monge_surface()
        poly_surface_nrmls = poly_surface.get_face_normals()
        for nrml in poly_surface_nrmls:
            self.logger.debug("nrml=%s", nrml.tolist())
            self.assertTrue(_np.all(~_np.isnan(nrml)))
            self.assertAlmostEqual(1.0, _np.linalg.norm(nrml))

    def test_to_meshio_mesh(self):
        """
        Test :meth:`surface_poly_fit.core.PolyhedralSurface.to_meshio_mesh`.
        """
        poly_surface = create_monge_surface()
        mio_mesh = poly_surface.to_meshio_mesh()

        self.assertSequenceEqual(
            poly_surface.get_vertices().tolist(),
            mio_mesh.points.tolist()
        )
        self.assertSequenceEqual(
            poly_surface.get_faces().tolist(),
            mio_mesh.cells[0].data.tolist()
        )
        self.assertSequenceEqual(
            poly_surface.get_vertex_normals().tolist(),
            mio_mesh.point_data["Normals"].tolist(),
        )

    def test_from_meshio_mesh(self):
        """
        Test :meth:`surface_poly_fit.core.PolyhedralSurface.from_meshio_mesh`.
        """
        from surface_poly_fit.core import PolyhedralSurface
        poly_surface_orig = create_monge_surface()
        mio_mesh = poly_surface_orig.to_meshio_mesh()
        poly_surface_from_mio = PolyhedralSurface.from_meshio_mesh(mio_mesh)

        self.assertSequenceEqual(
            poly_surface_orig.get_vertices().tolist(),
            poly_surface_from_mio.get_vertices().tolist(),
        )
        # Note: face order is only guaranteed to be preserved for single-degree face meshes
        # i.e. all-triangle-face meshes, all-quad-face meshes.
        self.assertSequenceEqual(
            poly_surface_orig.get_faces().tolist(),
            poly_surface_from_mio.get_faces().tolist()
        )
        self.assertSequenceEqual(
            poly_surface_orig.get_vertex_normals().tolist(),
            poly_surface_from_mio.get_vertex_normals().tolist(),
        )
        self.assertSequenceEqual(
            poly_surface_orig.get_face_normals().tolist(),
            poly_surface_from_mio.get_face_normals().tolist(),
        )


class MongePolynomial:
    """
    Evaluates Monge polynomial.
    """
    def __init__(self, k, b, c):
        self._k = _np.asanyarray(k).copy()
        self._b = _np.asanyarray(b).copy()
        self._c = _np.asanyarray(c).copy()

    @property
    def k(self):
        """
        Principal curvature coefficients, a :samp:`(2,)` shaped :obj:`numpy.ndarray`.
        """
        return self._k

    @property
    def b(self):
        """
        Curvature derivative coefficients, a :samp:`(4,)` shaped :obj:`numpy.ndarray`.
        """
        return self._b

    @property
    def c(self):
        """
        Second order curvature derivative coefficients, a :samp:`(5,)` shaped :obj:`numpy.ndarray`.
        """
        return self._c

    def evaluate(self, xy):
        """
        Evaluates polynomial at 2D coordinates :samp:`{xy}`.

        :type xy: :obj:`numpy.ndarray`
        :param xy: A :samp:`(N, 2)` shaped array of 2D coordinates.
        :rtype: :obj:`numpy.ndarray`
        :return: A :samp:`(N,)` shaped array of polynomial surface *heights*.
        """
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
        """
        See :meth:`evaluate`.
        """
        return self.evaluate(xy)


def create_monge_surface(monge_polynomial=None, xy=None):
    """
    Create a Monge polynomial shaped :obj:`surface_poly_fit.core.PolyhedralSurface`.
    Performs *Delaunay* triangulation of the :samp:`{xy}` coordinates to form
    the triangular mesh faces.

    :type monge_polynomial: :obj:`MongePolynomial`
    :param monge_polynomial: Defines the polynomial surface.
    :type xy: :obj:`numpy.ndarray`
    :param xy: A :samp:`(N, 2)` shaped array of 2D coordinates at which polynomial is evaluated.
    :rtype: :obj:`surface_poly_fit.core.PolyhedralSurface`
    :return: A Monge polynomial :obj:`surface_poly_fit.core.PolyhedralSurface`.
    """
    from scipy.spatial import Delaunay
    from surface_poly_fit.core import PolyhedralSurface

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
    """
    Tests for :obj:`surface_poly_fit.core.MongeJetFitter`.
    """
    def test_construct(self):
        """
        Test :obj:`surface_poly_fit.core.MongeJetFitter` instance creation.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        fitter = MongeJetFitter(poly_surf)
        self.assertIsNotNone(fitter)

        fitter = MongeJetFitter(poly_surf, 4)
        self.assertIsNotNone(fitter)

        fitter = MongeJetFitter(poly_surf, 4, 4)
        self.assertIsNotNone(fitter)

    def test_properties(self):
        """
        Test property get and set.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surf = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        fitter = MongeJetFitter(poly_surf, 4, 2)
        self.assertEqual(4, fitter.degree_poly_fit)
        self.assertEqual(2, fitter.degree_monge)
        self.assertEqual(15, fitter.min_num_fit_points)
        fitter.ring_normal_gaussian_sigma = 4.0
        self.assertEqual(4.0, fitter.ring_normal_gaussian_sigma)

    def test_fit(self):
        """
        Test :meth:`surface_poly_fit.core.MongeJetFitter.fit_at_vertex`.
        """
        from surface_poly_fit.core import MongeJetFitter

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
        # self.export_mesh("monge_poly_surface.ply", poly_surface)

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
        self.assertTrue(
            _np.allclose(
                _np.identity(3, dtype=_np.float64),
                _np.absolute(result["poly_fit_basis"][0])
            )
        )
        self.assertTrue(
            _np.allclose(
                0.0,
                _np.absolute(result["poly_fit_residual_stats"].tolist()[0])
            )
        )

        if degree_monge < 3:
            self.assertTrue(
                _np.all(result["b"][0] == 0.0)
            )
        if degree_monge < 4:
            self.assertTrue(
                _np.all(result["c"][0] == 0.0)
            )

    def test_oriented_fit(self):
        """
        Test :meth:`surface_poly_fit.core.MongeJetFitter.fit_at_vertex` for
        polynomial which is *not aligned to coordinate axes*.
        """
        import trimesh
        from scipy.spatial.transform import Rotation
        from surface_poly_fit.core import MongeJetFitter, PolyhedralSurface

        monge_polynomial = \
            MongePolynomial(
                k=(0.50, 0.25),
                b=(0.0, 0.0, 0.0, 0.0),
                c=(0.0, 0.0, 0.0, 0.0, 0.0)
            )
        poly_surface = create_monge_surface(monge_polynomial)
        ax = _np.asarray((1.0, 2.0, 0.25))
        ax /= _np.linalg.norm(ax)
        R = Rotation.from_rotvec(ax * 55, degrees=True).as_matrix()
        t = (2.0, 4.0, 8.0)
        M = _np.identity(4, dtype=_np.float64)
        M[:3, :3] = R
        M[:3, 3] = t
        tmesh = \
            trimesh.Trimesh(vertices=poly_surface.get_vertices(), faces=poly_surface.get_faces())
        tmesh.apply_transform(M)
        poly_surface = PolyhedralSurface(vertices=tmesh.vertices, faces=tmesh.faces)

        self.logger.info("")
        self.logger.info(
            "poly_surface (min_z, max_z) = (%s, %s).",
            poly_surface.get_vertices()[:, 2].min(),
            poly_surface.get_vertices()[:, 2].max(),
        )
        # self.export_mesh("monge_surface_oriented.ply", poly_surface)

        monge_origin_vertex_index = poly_surface.num_vertices // 2
        degree_monge = 2
        degree_poly_fit = 2
        fitter = MongeJetFitter(poly_surface, degree_poly_fit, degree_monge)
        result = fitter.fit_at_vertex(monge_origin_vertex_index, num_rings=16)

        self.logger.info("R=%s", R.tolist())
        self.logger.info("Result = %s", result)
        self.assertEqual(monge_origin_vertex_index, result["vertex_index"][0])
        self.assertEqual(degree_monge, result["degree_monge"][0])
        self.assertEqual(degree_poly_fit, result["degree_poly_fit"][0])
        self.assertAlmostEqual(monge_polynomial.k[0], result["k"][0][0])
        self.assertAlmostEqual(monge_polynomial.k[1], result["k"][0][1])
        self.assertTrue(
            _np.allclose(t, result["origin"][0])
        )
        self.assertTrue(
            _np.allclose(_np.absolute(R), _np.absolute(result["direction"][0]))
        )
        self.assertTrue(
            _np.allclose(
                _np.dot(R, _np.asarray((0.0, 0.0, 1.0)).reshape((3, 1))).T,
                _np.dot(result["poly_fit_basis"][0], _np.asarray((0.0, 0.0, 1.0)).reshape((3, 1))).T
            )
        )
        self.assertTrue(
            _np.allclose(
                0.0,
                _np.absolute(result["poly_fit_residual_stats"].tolist()[0])
            )
        )

        if degree_monge < 3:
            self.assertTrue(
                _np.all(result["b"][0] == 0.0)
            )
        if degree_monge < 4:
            self.assertTrue(
                _np.all(result["c"][0] == 0.0)
            )

    def test_fit_all(self):
        """
        Test :meth:`surface_poly_fit.core.MongeJetFitter.fit_all`.
        """
        from surface_poly_fit.core import MongeJetFitter

        poly_surface = create_monge_surface()

        degree_monge = 4
        degree_poly_fit = 4
        fitter = MongeJetFitter(poly_surface, degree_poly_fit, degree_monge)
        result = fitter.fit_all(num_rings=8)
        self.logger.info("Result = %s", result)
        self.assertEqual(poly_surface.num_vertices, len(result))
        self.assertTrue(_np.all(result["degree_monge"] == degree_monge))
        self.assertTrue(_np.all(result["degree_poly_fit"] == degree_poly_fit))

        if degree_monge < 3:
            self.assertTrue(
                _np.all(result["b"] == 0.0)
            )
        if degree_monge < 4:
            self.assertTrue(
                _np.all(result["c"] == 0.0)
            )

    def test_fit_all_bounding_area(self):
        """
        Checks on the bounding-area (:samp:`"poly_fit_bounding_area"` field) results
        returned by :meth:`surface_poly_fit.core.MongeJetFitter.fit_all`.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surface = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        degree_monge = 2
        degree_poly_fit = 2
        fitter = MongeJetFitter(poly_surface, degree_poly_fit, degree_monge)
        results = []
        for nr in [2, 4, 8, 16]:
            results.append(fitter.fit_all(num_rings=nr))
        for i in 0, 1, 2, 3:
            self.assertTrue(
                _np.all(
                    results[i]["poly_fit_bounding_area"]["ellipse_min_radius"]
                    <=
                    results[i]["poly_fit_bounding_area"]["ellipse_max_radius"]
                )
            )
            self.assertTrue(
                _np.all(
                    results[i]["poly_fit_bounding_area"]["rectangle_min_side_length"]
                    <=
                    results[i]["poly_fit_bounding_area"]["rectangle_max_side_length"]
                )
            )
            self.assertTrue(
                _np.all(
                    results[i]["poly_fit_bounding_area"]["ellipse_min_radius"]
                    <=
                    results[i]["poly_fit_bounding_area"]["circle_radius"]
                )
            )
            self.assertTrue(
                _np.all(
                    results[i]["poly_fit_bounding_area"]["ellipse_max_radius"]
                    <=
                    results[i]["poly_fit_bounding_area"]["circle_radius"]
                )
            )
            if i > 0:
                self.assertTrue(
                    _np.all(
                        results[i]["poly_fit_bounding_area"]["circle_radius"]
                        >
                        results[i - 1]["poly_fit_bounding_area"]["circle_radius"]
                    )
                )
                self.assertTrue(
                    _np.all(
                        results[i]["poly_fit_bounding_area"]["ellipse_min_radius"]
                        >
                        results[i - 1]["poly_fit_bounding_area"]["ellipse_min_radius"]
                    )
                )
                self.assertTrue(
                    _np.all(
                        results[i]["poly_fit_bounding_area"]["ellipse_max_radius"]
                        >
                        results[i - 1]["poly_fit_bounding_area"]["ellipse_max_radius"]
                    )
                )
                self.assertTrue(
                    _np.all(
                        results[i]["poly_fit_bounding_area"]["rectangle_min_side_length"]
                        >
                        results[i - 1]["poly_fit_bounding_area"]["rectangle_min_side_length"]
                    )
                )
                self.assertTrue(
                    _np.all(
                        results[i]["poly_fit_bounding_area"]["rectangle_max_side_length"]
                        >
                        results[i - 1]["poly_fit_bounding_area"]["rectangle_max_side_length"]
                    )
                )

    def test_to_meshio_mesh(self):
        """
        Test :meth:`surface_poly_fit.core.MongeJetFitter.to_meshio_mesh`.
        """
        from trimesh.primitives import Capsule
        from surface_poly_fit.core import PolyhedralSurface, MongeJetFitter

        trimesh_mesh = Capsule()
        poly_surface = PolyhedralSurface(vertices=trimesh_mesh.vertices, faces=trimesh_mesh.faces)

        degree_monge = 2
        degree_poly_fit = 2
        fitter = MongeJetFitter(poly_surface, degree_poly_fit, degree_monge)
        result_array = fitter.fit_all(num_rings=8)

        mio_mesh = fitter.to_meshio_mesh(result_array)
        self.assertSequenceEqual(
            fitter.poly_surface.get_vertices().tolist(),
            mio_mesh.points.tolist()
        )
        self.assertSequenceEqual(
            result_array["vertex_index"].tolist(),
            mio_mesh.point_data["vertex_index"].tolist()
        )
        self.assertSequenceEqual(
            result_array["num_fitting_points"].tolist(),
            mio_mesh.point_data["num_fitting_points"].tolist()
        )
        self.assertSequenceEqual(
            result_array["k"][:, 0].tolist(),
            mio_mesh.point_data["k0"].tolist()
        )
        self.assertSequenceEqual(
            result_array["k"][:, 1].tolist(),
            mio_mesh.point_data["k1"].tolist()
        )
        self.assertSequenceEqual(
            result_array["direction"][:, :, 0].tolist(),
            mio_mesh.point_data["k0_dir"].tolist()
        )
        self.assertSequenceEqual(
            result_array["direction"][:, :, 1].tolist(),
            mio_mesh.point_data["k1_dir"].tolist()
        )
        self.assertSequenceEqual(
            result_array["direction"][:, :, 2].tolist(),
            mio_mesh.point_data["k_normal"].tolist()
        )


__all__ = [s for s in dir() if not s.startswith('_')]

if __name__ == "__main__":
    logging_format = "%(asctime)s|%(process)-8s|%(name)-8s|%(levelname)-8s|%(message)s"
    _logging.basicConfig(format=logging_format, level=_logging.INFO)
    _unittest.main(__name__)
