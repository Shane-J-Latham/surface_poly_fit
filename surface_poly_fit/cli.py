"""
Command line interface for surface polynomial fitting.
"""


def read_polyhedral_surface(file_name):
    """
    Read mesh from file. Uses :mod:`meshio` to parse file.

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
    Write polynomial fit results and :samp:`{polyhedral_surface}` vertices and face-lists
    to file.
    """
    import numpy as np

    np.savez(
        output_file_name,
        surface_poly_fit_results=result_ary,
        vertices=polyhedral_surface.get_vertices() if polyhedral_surface is not None else None,
        faces=polyhedral_surface.get_faces() if polyhedral_surface is not None else None
    )


def surface_poly_fit_cli(args):
    """
    Command line interface for reading mesh file, fitting polynomials and writing results to file.
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
    logger.info("Reading mesh from file %s...", args.mesh_file)
    ps = read_polyhedral_surface(args.mesh_file)
    fitter = Fitter(ps, degree_poly_fit=args.degree_poly_fit, degree_monge=args.degree_monge)
    logger.info("Fitting for num_rings=%s...", num_rings_list[0])
    result_ary = fitter.fit_all(num_rings_list[0])
    for num_rings in num_rings_list[1:]:
        logger.info("Fitting for num_rings=%s...", num_rings)
        result_ary = np.hstack(fitter.fit_all(num_rings=num_rings))

    logger.info("Writing fitting result to file %s...", output_file_name)
    write_result_array(output_file_name, result_ary)


def get_argument_parser():
    """
    Returns a :obj:`argparse.ArgumentParser` to handle command line option processing.
    """
    import argparse
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
        "mesh_file",
        help="Mesh file path."
    )

    return ap


def surface_poly_fit_main_cli():
    """
    Entry-point for command line interface.
    """
    surface_poly_fit_cli(get_argument_parser().parse_args())


if __name__ == "__main__":
    surface_poly_fit_main_cli()
