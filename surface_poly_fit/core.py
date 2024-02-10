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
        "+++++++++++++++++++++++++++++++\n"
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
    pass


class MongeJetFitter(_MongeJetFitter):
    """
    Fits Monge polynomial to :obj:`PolyhedralSurface` patch vertices.
    """
    pass


MongeJetFitter.__doc__ += _fitting_basis_type_doc_string()


__all__ = [s for s in dir() if not s.startswith('_')]
