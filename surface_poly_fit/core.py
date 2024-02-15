"""
Polynomial fitting classes and functions.


.. autosummary::
   :toctree: generated/

"""
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

    def to_meshio_mesh(self):
        """
        Return a :obj:`meshio.Mesh` version of this surface.

        :rtype: :obj:`meshio.Mesh`
        :return: This mesh converted to a :obj:`meshio.Mesh` instance.
        """
        import numpy as np
        import meshio
        from collections import defaultdict

        points = self.get_vertices()
        faces = self.get_faces()
        face_normals = self.get_face_normals()

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
                (cell_type_dict[num_verts], np.asarray(cell_dict[num_verts]))
                for num_verts in sorted(list(cell_dict.keys()))
            )
        del cell_dict

        point_data = {"Normals": self.get_vertex_normals()}
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
        import numpy as np

        faces = sum(list(np.asanyarray(cells.data).tolist() for cells in meshio_mesh.cells), list())
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
        import numpy as np

        degree_monge = np.unique(result_array["degree_monge"])
        if len(degree_monge) == 1:
            degree_monge = degree_monge[0]
        else:
            degree_monge = str(tuple(degree_monge))

        degree_poly_fit = np.unique(result_array["degree_poly_fit"])
        if len(degree_poly_fit) == 1:
            degree_poly_fit = degree_poly_fit[0]
        else:
            degree_poly_fit = str(tuple(degree_poly_fit))

        num_rings = np.unique(result_array["num_rings"])
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

    def convert_result_to_point_data(self, result_array):
        """
        Convert the :samp:`{result_array}` record array to a *point-data* :obj:`dict`
        suitable for :attr:`meshio.Mesh.point_data`.

        :type result_array: :obj:`numpy.ndarray`
        :param result_array: Polynomial fitting result array, e.g. as returned by :meth:`fit_all`.
        :rtype: :obj:`dict`
        :return: A dictionary of :samp:`(str, numpy.ndarray)` entries.
        """
        point_data_dict = dict()
        point_data_dict["vertex_index"] = result_array["vertex_index"]
        point_data_dict["num_fitting_points"] = result_array["num_fitting_points"]
        point_data_dict["pf_condition_number"] = result_array["poly_fit_condition_number"]

        ra_pfrs = result_array["poly_fit_residual_stats"]
        point_data_dict["pfrs_min_abs"] = ra_pfrs["min_abs"]
        point_data_dict["pfrs_max_abs"] = ra_pfrs["max_abs"]
        point_data_dict["pfrs_mean_abs"] = ra_pfrs["mean_abs"]
        point_data_dict["pfrs_median_abs"] = ra_pfrs["median_abs"]
        point_data_dict["pfrs_stdd"] = ra_pfrs["stdd"]

        point_data_dict["k0_dir"] = result_array["direction"][:, :, 0]
        point_data_dict["k1_dir"] = result_array["direction"][:, :, 1]
        point_data_dict["k_normal"] = result_array["direction"][:, :, 2]
        point_data_dict["k0"] = result_array["k"][:, 0]
        point_data_dict["k1"] = result_array["k"][:, 1]

        ra_pfba = result_array["poly_fit_bounding_area"]
        point_data_dict["ba_rect_minor"] = ra_pfba["rectangle_min_side_length"]
        point_data_dict["ba_rect_major"] = ra_pfba["rectangle_max_side_length"]
        point_data_dict["ba_ellip_minor"] = ra_pfba["ellipse_min_radius"]
        point_data_dict["ba_ellip_major"] = ra_pfba["ellipse_max_radius"]
        point_data_dict["ba_circle"] = ra_pfba["circle_radius"]

        return point_data_dict

    def to_meshio_mesh(self, result_array):
        """
        Return a :obj:`meshio.Mesh` version of the surface with *point-data*
        assigned from the :samp:`{result_array}`.

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
        mio_mesh = self.poly_surface.to_meshio_mesh()
        mio_mesh.point_data.update(self.convert_result_to_point_data(result_array))
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
