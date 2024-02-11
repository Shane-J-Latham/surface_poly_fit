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


class MongeJetFitter(_MongeJetFitter):
    """
    Fits Monge polynomial to :obj:`PolyhedralSurface` patch vertices.
    """
    pass


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
