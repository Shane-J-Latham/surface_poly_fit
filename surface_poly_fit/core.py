"""
Polynomial fitting classes and functions.


.. autosummary::
   :toctree: generated/

"""
import numpy as _np
from ._spf_cgal import PolyhedralSurface as _PolyhedralSurface   # noqa: F401
from ._spf_cgal import MongeJetFitter as _MongeJetFitter


def _fitting_basis_type_doc_string():
    doc_str = (
        "\n"
        +
        ".. _fitting-basis-type-description:\n\n"
        +
        "`FittingBasisType` Descriptions\n"
        +
        "+++++++++++++++++++++++++++++++\n\n"
        +
        "The :meth:`MongeJetFitter.fit_at_vertex` and :meth:`MongeJetFitter.fit_all`"
        +
        " methods take a :samp:`{fit_basis_type}` argument which specifies how the"
        +
        " polynomial fitting coordinate system (fitting-basis) is calculated"
        +
        " for a polyhedral-surface-patch:\n\n"
        +
        "".join(
            list(
                ":attr:`MongeJetFitter.%s`\n   %s\n" % (k, v[1])
                for (k, v) in _MongeJetFitter.FittingBasisType.__entries.items()
            )
        )
        +
        "\n\n"
    )
    return doc_str


class PolyhedralSurface(_PolyhedralSurface):
    """
    A mesh consisting of polygonal faces.
    """

    def create_ring_patch(self, vertex_index, num_rings):
        """
        Creates a polyhedral surface patch about a specified vertex.

        :type vertex_index: :obj:`int`
        :param vertex_index: The index of the vertex about which the patch is created.
        :type num_rings: :obj:`int`
        :param num_rings: The number of edge-hops from vertex :samp:`{vertex_index}`
           used to form the patch.
        :rtype: :obj:`PolyhedralSurface`
        :return: A new instance :obj:`PolyhedralSurface` patch.
        """

        return \
            self._create_child_ring_patch(
                vertex_index=vertex_index,
                num_rings=num_rings,
                child_class=self.__class__
            )

    def to_meshio_mesh(self, float_dtype=_np.float64):
        """
        Return a :obj:`meshio.Mesh` version of this surface.

        :type float_dtype: :obj:`numpy.dtype`
        :param float_dtype: The floating point type for vertices and normals of the returned mesh.
        :rtype: :obj:`meshio.Mesh`
        :return: This mesh converted to a :obj:`meshio.Mesh` instance.
        """
        import meshio
        from collections import defaultdict

        points = self.get_vertices().astype(float_dtype)
        faces = self.get_faces()
        face_normals = self.get_face_normals().astype(float_dtype)

        # Convert face list to list-of-tuples (group faces of same degree/number-of-vertices).
        cell_type_dict = defaultdict(lambda: "polygon")
        cell_type_dict.update({3: "triangle", 4: "quad"})
        cell_dict = defaultdict(list)
        cell_nrmls_dict = defaultdict(list)
        for face_idx in range(len(faces)):
            face = faces[face_idx]
            cell_dict[len(face)].append(face)
            cell_nrmls_dict[len(face)].append(face_normals[face_idx])
        cells = \
            list(
                (cell_type_dict[num_verts], _np.asarray(cell_dict[num_verts]))
                for num_verts in sorted(list(cell_dict.keys()))
            )
        del cell_dict

        point_data = {"Normals": self.get_vertex_normals().astype(float_dtype)}
        cell_data = {
            "Normals": list(
                cell_nrmls_dict[num_verts] for num_verts in sorted(list(cell_nrmls_dict.keys()))
            )
        }

        return meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data)

    @classmethod
    def from_meshio_mesh(cls, meshio_mesh):
        """
        Return a :obj:`PolyhedralSurface` version of the :samp:`{meshio_mesh}` mesh.

        :type meshio_mesh: :obj:`meshio.Mesh`
        :param meshio_mesh: Convert this mesh to a :obj:`PolyhedralSurface`.
        :rtype: :obj:`PolyhedralSurface`
        :return: The :samp:`{meshio_mesh}` mesh converted to a :obj:`PolyhedralSurface` instance.
        """
        faces = \
            sum(
                list(_np.asanyarray(cells.data).tolist() for cells in meshio_mesh.cells),
                list()
            )
        ps = PolyhedralSurface(vertices=meshio_mesh.points, faces=faces)
        if "Normals" in meshio_mesh.point_data:
            ps.set_vertex_normals(meshio_mesh.point_data["Normals"])

        return ps


class MongeJetFitter(_MongeJetFitter):
    """
    Fits Monge polynomial to :obj:`PolyhedralSurface` patch vertices.
    """

    def get_field_data(self, result_array):
        """
        Returns a :obj:`dict` which can be used as :attr:`meshio.Mesh.field_data`.

        :type result_array: :obj:`numpy.ndarray`
        :param result_array: Polynomial fitting result array, e.g. as returned by :meth:`fit_all`.
        :rtype: :obj:`dict`
        :return: A dictionary of :samp:`(str, obj)` *field-data* entries.
        """
        degree_monge = _np.unique(result_array["degree_monge"])
        if len(degree_monge) == 1:
            degree_monge = degree_monge[0]
        else:
            degree_monge = str(tuple(degree_monge))

        degree_poly_fit = _np.unique(result_array["degree_poly_fit"])
        if len(degree_poly_fit) == 1:
            degree_poly_fit = degree_poly_fit[0]
        else:
            degree_poly_fit = str(tuple(degree_poly_fit))

        num_rings = _np.unique(result_array["num_rings"])
        if len(num_rings) == 1:
            num_rings = num_rings[0]
        else:
            num_rings = str(tuple(num_rings))

        return \
            {
                "degree_monge": degree_monge,
                "degree_poly_fit": degree_poly_fit,
                "num_rings": num_rings,
            }

    def convert_result_to_point_data(self, result_array, float_dtype=_np.float64):
        """
        Convert the :samp:`{result_array}` record array to a *point-data* :obj:`dict`
        suitable for :attr:`meshio.Mesh.point_data`.

        :type float_dtype: :obj:`numpy.dtype`
        :param float_dtype: The floating point type used for float fields in the returned
           point-data :obj:`dict`.
        :type result_array: :obj:`numpy.ndarray`
        :param result_array: Polynomial fitting result array, e.g. as returned by :meth:`fit_all`.
        :rtype: :obj:`dict`
        :return: A dictionary of :samp:`(str, numpy.ndarray)` entries.
        """
        import numpy as np

        point_data_dict = dict()

        max_vertex_index = _np.max(result_array["vertex_index"])
        vi_dtype = _np.uint64
        if max_vertex_index < np.iinfo(_np.uint8).max:
            vi_dtype = _np.uint8
        elif max_vertex_index < np.iinfo(_np.uint16).max:
            vi_dtype = _np.uint16
        elif max_vertex_index < np.iinfo(_np.uint32).max:
            vi_dtype = _np.uint32

        point_data_dict["vertex_index"] = result_array["vertex_index"].astype(vi_dtype)

        max_num_fit_coordinates = _np.max(result_array["num_fitting_points"])
        nfp_dtype = _np.uint32
        if max_num_fit_coordinates < np.iinfo(_np.uint8).max:
            nfp_dtype = _np.uint8
        if max_num_fit_coordinates < np.iinfo(_np.uint16).max:
            nfp_dtype = _np.uint16
        point_data_dict["num_fitting_points"] = result_array["num_fitting_points"].astype(nfp_dtype)

        point_data_dict["pf_condition_number"] = \
            result_array["poly_fit_condition_number"].astype(float_dtype)

        ra_pfrs = result_array["poly_fit_residual_stats"]
        point_data_dict["pfrs_min_abs"] = ra_pfrs["min_abs"].astype(float_dtype)
        point_data_dict["pfrs_max_abs"] = ra_pfrs["max_abs"].astype(float_dtype)
        point_data_dict["pfrs_mean_abs"] = ra_pfrs["mean_abs"].astype(float_dtype)
        point_data_dict["pfrs_median_abs"] = ra_pfrs["median_abs"].astype(float_dtype)
        point_data_dict["pfrs_stdd"] = ra_pfrs["stdd"].astype(float_dtype)

        point_data_dict["k0_dir"] = result_array["direction"][:, :, 0].astype(float_dtype)
        point_data_dict["k1_dir"] = result_array["direction"][:, :, 1].astype(float_dtype)
        point_data_dict["k_normal"] = result_array["direction"][:, :, 2].astype(float_dtype)
        point_data_dict["k0"] = result_array["k"][:, 0].astype(float_dtype)
        point_data_dict["k1"] = result_array["k"][:, 1].astype(float_dtype)

        k_abs = np.absolute(result_array["k"]).astype(float_dtype)
        k_abs_min_idx = tuple(np.argmin(k_abs, axis=1))
        k_abs_max_idx = tuple(np.argmax(k_abs, axis=1))
        k_idx = tuple(np.arange(k_abs.shape[0]))
        point_data_dict["k_abs_min"] = k_abs[k_idx, k_abs_min_idx].astype(float_dtype)
        point_data_dict["k_abs_max"] = k_abs[k_idx, k_abs_max_idx].astype(float_dtype)
        point_data_dict["k_abs_min_dir"] = \
            np.asarray(
                list(
                    result_array["direction"][k_idx, i, k_abs_min_idx] for i in (0, 1, 2)
                )
            ).T.astype(float_dtype)
        point_data_dict["k_abs_max_dir"] = \
            np.asarray(
                list(
                    result_array["direction"][k_idx, i, k_abs_max_idx] for i in (0, 1, 2)
                )
            ).T.astype(float_dtype)

        ra_pfba = result_array["poly_fit_bounding_area"]
        point_data_dict["ba_rect_minor"] = ra_pfba["rectangle_min_side_length"].astype(float_dtype)
        point_data_dict["ba_rect_major"] = ra_pfba["rectangle_max_side_length"].astype(float_dtype)
        point_data_dict["ba_ellip_minor"] = ra_pfba["ellipse_min_radius"].astype(float_dtype)
        point_data_dict["ba_ellip_major"] = ra_pfba["ellipse_max_radius"].astype(float_dtype)
        point_data_dict["ba_circle"] = ra_pfba["circle_radius"].astype(float_dtype)

        return point_data_dict

    def to_meshio_mesh(self, result_array, float_dtype=_np.float64):
        """
        Return a :obj:`meshio.Mesh` version of the surface with *point-data*
        fields assigned from the :samp:`{result_array}`.

        :type float_dtype: :obj:`numpy.dtype`
        :param float_dtype: The floating point type used for vertices, normals and
           float fields in the returned :obj:`meshio.Mesh`.
        :rtype: :obj:`meshio.Mesh`
        :return: A :obj:`meshio.Mesh` instance with point-data assigned
           from the :samp:`{result_array}`.
        """
        if len(result_array) != self.poly_surface.num_vertices:
            raise ValueError(
                f"Got len(result_array)={len(result_array)}"
                +
                ", expected len(result_array)=self.poly_surface.num_vertices="
                +
                f"{self.poly_surface.num_vertices}"
            )
        mio_mesh = self.poly_surface.to_meshio_mesh(float_dtype=float_dtype)
        mio_mesh.point_data.update(
            self.convert_result_to_point_data(result_array, float_dtype=float_dtype)
        )
        mio_mesh.field_data.update(self.get_field_data(result_array))

        return mio_mesh


MongeJetFitter.__doc__ += _fitting_basis_type_doc_string()

MongeJetFitter.__doc__ += \
"""

.. _fitting-result-array-description:

Fitting Result Array Field Definitions
++++++++++++++++++++++++++++++++++++++

The :meth:`MongeJetFitter.fit_at_vertex` and :meth:`MongeJetFitter.fit_all` return
a :mod:`numpy` `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_
with fields:


`"vertex_index"` (:obj:`numpy.int64`)
   The index of the vertex at which polynomial fitting was performed.
`"degree_monge"` (:obj:`numpy.uint8`)
   The *degree* of the Monge polynomial (converted from the fitting polynomial).
`"degree_poly_fit"` (:obj:`numpy.uint8`)
   The *degree* of the fitting polynomial.
`"num_rings"` (:obj:`numpy.int64`)
   The number of ring neighbours defining the polynomial fitting coordinates.
`"num_fitting_points"` (:obj:`numpy.int64`)
   The number of neighbour coordinates (neighbours) in the `"num_rings"` neighbourhood.
   This is the number of coordinates used to fit the polynomial.
`"poly_fit_condition_number"` (:obj:`numpy.float64`)
   The polynomial fitting matrix condition number.
`"pca_eigenvalues"` (:obj:`numpy.float64`)
   A :samp:`(3,)` shaped array. If :attr:`MongeJetFitter.PCA` is used for the fitting
   basis determination, then these are the eigenvalues of the principal component
   analysis (otherwise all values are :obj:`numpy.nan`).
`"poly_fit_basis"` (:obj:`numpy.float64`)
   A :samp:`(3, 3)` shaped rotation matrix, indicating the orientation of
   the polynomial fitting coordinate-system (basis).
`"poly_fit_residual_stats"` (:obj:`numpy.dtype`)
   A *struct* containing statistics on the polynomial fitting residual values.
   See :ref:`residual stats description <fitting-result-array-poly-fit-residual-stats-description>`.
`"poly_fit_bounding_area"` (:obj:`numpy.dtype`)
   A *struct* containing info about 2D bounding areas of the fitting coordinates.
   See :ref:`bounding area description <fitting-result-array-poly-fit-bounding-area-description>`.
`"origin"` (:obj:`numpy.float64`)
   A :samp:`(3,)` shaped array indicating the Monge polynomial origin coordinate.
   The Monge polynomial intersects this coordinate. Note that the `"vertex_index"`
   coordinate is the *origin* of the fitting-coordinate-system.
`"direction"` (:obj:`numpy.float64`)
   A :samp:`(3, 3)` shaped rotation matrix :samp:`R`, where columns :samp:`R[:, 0]`
   and :samp:`R[:, 1]` are the principal curvature directions, and where :samp:`R[:, 2]`
   is *outward facing* and orthogonal to :samp:`R[:, 0]` and :samp:`R[:, 1]`.
`"k"` (:obj:`numpy.float64`)
   A :samp:`(2,)` shaped array of the principal curvature magnitudes.
`"b"` (:obj:`numpy.float64`)
   A :samp:`(4,)` shaped array, :samp:`b[0]` and :samp:`b[3]`
   are the directional derivatives of :samp:`k[0]` and :samp:`k[1]`
   along their respective curvature line. :samp:`b[1]` and :samp:`b[2]` are the
   directional derivatives of :samp:`k[0]` and :samp:`k[1]` along the other curvature lines.
`"c"` (:obj:`numpy.float64`)
   A :samp:`(5,)` shaped array of higher order derivatives of curvature.

.. _fitting-result-array-poly-fit-residual-stats-description:

The `"poly_fit_residual_stats"` structure contains residual statistics and has fields:

`"min"` (:obj:`numpy.float64`)
   Minimum residual value.
`"max"` (:obj:`numpy.float64`)
   Maximum residual value.
`"mean"` (:obj:`numpy.float64`)
   Mean residual value.
`"median"` (:obj:`numpy.float64`)
   Median residual value.
`"min_abs"` (:obj:`numpy.float64`)
   Minimum absolute residual value.
`"max_abs"` (:obj:`numpy.float64`)
   Maximum absolute residual value.
`"mean_abs"` (:obj:`numpy.float64`)
   Mean absolute residual value.
`"median_abs"` (:obj:`numpy.float64`)
   Median absolute residual value.
`"stdd"` (:obj:`numpy.float64`)
   Standard deviation residual value.


.. _fitting-result-array-poly-fit-bounding-area-description:

The `"poly_fit_bounding_area"` structure has 2D bounding area info in fields:

`"rectangle_min_side_length"` (:obj:`numpy.float64`)
   Minimum oriented bounding rectangle smallest side length.
`"rectangle_max_side_length"` (:obj:`numpy.float64`)
   Minimum oriented bounding rectangle largest side length.
`"circle_radius"` (:obj:`numpy.float64`)
   Minimum bounding circle radius.
`"ellipse_min_radius"` (:obj:`numpy.float64`)
   Minimum oriented bounding ellipse smallest radius.
`"ellipse_max_radius"` (:obj:`numpy.float64`)
   Minimum oriented bounding ellipse largest radius.

"""  # noqa: E122

__all__ = [s for s in dir() if not s.startswith('_')]
