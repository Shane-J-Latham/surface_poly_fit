"""
Command line interface for surface polynomial fitting.
"""
from surface_poly_fit._spf_cgal import MongeJetFitter


def read_polyhedral_surface(file_name):
    """
    Read mesh from file. Uses :mod:`meshio` to parse file.

    :type file_name: :obj:`str`
    :param file_name: Mesh file path.
    :rtype: :obj:`surface_poly_fit.PolyhedralSurface`
    """
    import numpy as np
    import meshio
    from . import PolyhedralSurface

    mesh = meshio.read(file_name)
    faces = sum(list(np.asanyarray(cells.data).tolist() for cells in mesh.cells), list())
    ps = PolyhedralSurface(vertices=mesh.points, faces=faces)
    return ps


def write_result_array(output_file_name, result_ary, polyhedral_surface=None):
    """
    Writes polynomial fit results and :samp:`{polyhedral_surface}` vertices, face-lists
    and vertex-normals to :samp:`.npz` file. The *fields* within the written :samp:`.npz`
    are:

    `"surface_poly_fit_results"`
       The :samp:`{result_ary}`.

    `"vertices"`
       The :samp:`{polyhedral_surface}.get_vertices()`.

    `"faces"`
       The :samp:`{polyhedral_surface}.get_faces()`.

    `"vertex_normals"`
       The :samp:`{polyhedral_surface}.get_vertex_normals()`.


    :type output_file_name: :obj:`str`
    :param output_file_name: Output file path for :func:`numpy.savez_compressed` file.
    :type result_ary: :obj:`numpy.ndarray`
    :param result_ary: Polynomial fitting result
       array (e.g. as returned by :meth:`surface_poly_fit.MongeJetFitter.fit_all`).
    :type polyhedral_surface: :obj:`surface_poly_fit.PolyhedralSurface`
    :param polyhedral_surface: If not :obj:`None`, write vertices, faces and vertex-normals
        the :samp:`.npz` file.
    """
    import numpy as np

    vertices = None
    faces = None
    vertex_normals = None
    if polyhedral_surface is not None:
        vertices = polyhedral_surface.get_vertices()
        faces = polyhedral_surface.get_faces()
        vertex_normals = polyhedral_surface.get_vertex_normals()

    np.savez_compressed(
        output_file_name,
        surface_poly_fit_results=result_ary,
        vertices=vertices,
        faces=faces,
        vertex_normals=vertex_normals
    )


def surface_poly_fit_cli(args):
    """
    Command line interface for reading mesh file, fitting polynomials and writing results to file.

    :type args: :obj:`types.SimpleNamespace`
    :param args: Parsed command line arguments (e.g. as returned
       by :samp:`surface_poly_fit.get_argument_parser().parse_args()`.
    """
    import logging
    import os
    import numpy as np
    from . import MongeJetFitter as Fitter

    logging.basicConfig(level=getattr(logging, args.log_level))
    logger = logging.getLogger()

    output_file_name = args.output_file
    if output_file_name is None:
        output_file_name = os.path.splitext(os.path.split(args.mesh_file)[1])[0]
        output_file_name += "_surface_poly_fit" + os.path.extsep + "npz"
    num_rings_list = eval(args.num_rings)
    # Convert the argument string to enum type.
    args.poly_fit_basis_type = MongeJetFitter.FittingBasisType.__members__[args.poly_fit_basis_type]
    logger.info("Reading mesh from file %s...", args.mesh_file)
    polyhedral_surface = read_polyhedral_surface(args.mesh_file)
    fitter = \
        Fitter(
            polyhedral_surface,
            degree_poly_fit=args.degree_poly_fit,
            degree_monge=args.degree_monge
        )
    fitter.ring_normal_gaussian_sigma = args.poly_fit_basis_gaussian_sigma
    logger.info("poly_fit_basis_type=%s...", args.poly_fit_basis_type)
    logger.info("Fitting for num_rings=%s...", num_rings_list[0])
    result_ary = \
        fitter.fit_all(
            num_rings=num_rings_list[0],
            fit_basis_type=args.poly_fit_basis_type
        )
    for num_rings in num_rings_list[1:]:
        logger.info("Fitting for num_rings=%s...", num_rings)
        result_ary = np.hstack(fitter.fit_all(num_rings=num_rings))

    logger.info("Writing fitting result to file %s...", output_file_name)
    write_result_array(output_file_name, result_ary, polyhedral_surface)


def get_argument_parser():
    """
    Returns a :obj:`argparse.ArgumentParser` to handle command line option processing.

    :rtype :obj:`argparse.ArgumentParser`
    :return: Argument parser.
    """
    import argparse
    from . import MongeJetFitter
    ap = \
        argparse.ArgumentParser(
            "surface_poly_fit",
            description=(
                "Fit a Monge polynomial to vertex neighbourhoods (for all vertices of a mesh)."
            ),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    ap.add_argument(
        "-l", "--log_level",
        action='store',
        help="Amount of logging.",
        choices=("NONE", "CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO"
    )
    ap.add_argument(
        "-o", "--output_file",
        action='store',
        help="Name of the .npz file containing per-vertex poly-fit info.",
        default=None
    )
    ap.add_argument(
        "-p", "--degree_poly_fit",
        action='store',
        help="Degree of polynomial fit to the surface.",
        type=int,
        default=2
    )
    ap.add_argument(
        "-m", "--degree_monge",
        action='store',
        help="Degree of monge-polynomial form.",
        type=int,
        default=2
    )
    ap.add_argument(
        "-r", "--num_rings",
        action='store',
        help="List of neighbourhood-ring sizes.",
        default="[2, 4, 8]"
    )
    ap.add_argument(
        "--poly_fit_basis_type",
        choices=tuple(MongeJetFitter.FittingBasisType.__members__.keys()),
        help=(
            "How the polynomial fitting basis is determined from the neighbourhood-ring"
            +
            " of vertex-coorindates and/or the vertex-normals."
        ),
        default="RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN"
    )
    ap.add_argument(
        "--poly_fit_basis_gaussian_sigma",
        action="store",
        help=(
            "The sigma value used when calculating Gaussian weights for "
            +
            " RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA polynomial fitting basis."
        ),
        type=float,
        default=2.0
    )
    ap.add_argument(
        "mesh_file",
        help="Mesh file path."
    )

    return ap


def surface_poly_fit_main_cli():
    """
    Entry-point function for command line interface.
    """
    surface_poly_fit_cli(get_argument_parser().parse_args())


if __name__ == "__main__":
    surface_poly_fit_main_cli()
