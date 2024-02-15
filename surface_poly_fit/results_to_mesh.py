"""
Command line interface for converting polynomial-fitting results :samp:`".npz"`
(:mod:`surface_poly_fit.cli`) to a VTK :samp:`".vtu"` file for visualization.

.. autosummary::
   :toctree: generated/

"""


def read_results(file_name):
    """
    Read polynomial fitting results and corresponding polyhedral surface from :samp:`".npz"`
    file.

    :type file_name: :obj:`str`
    :param file_name: Path to input polynomial fitting results
       file (i.e. as written by :mod:`surface_poly_fit.cli`).
    :rtype: :obj:`tuple`
    :return: A :samp:`(poly_surface, result_array)` pair.
    """
    import numpy
    from .core import PolyhedralSurface

    npz_dict = numpy.load(file_name)
    results_array = npz_dict["surface_poly_fit_results"]
    ps = PolyhedralSurface(vertices=npz_dict["vertices"], faces=npz_dict["faces"])
    ps.set_vertex_normals(npz_dict["vertex_normals"])

    return (ps, results_array)


def write_result_mesh(output_file_name, polyhedral_surface, result_array):
    """
    Writes mesh file with fit results point-data.

    :type output_file_name: :obj:`str`
    :param output_file_name: Output file path for :func:`meshio.write` file.
    :type polyhedral_surface: :obj:`surface_poly_fit.core.PolyhedralSurface`
    :param polyhedral_surface: Defines vertices, faces and vertex normals written
      to mesh file.
    :type result_array: :obj:`numpy.ndarray`
    :param result_array: Polynomial fitting result
       array (e.g. as returned by :meth:`surface_poly_fit.core.MongeJetFitter.fit_all`).
       These results are written as *point-data* in the output mesh file (if supported).
    """
    from pathlib import Path
    import meshio
    from .core import MongeJetFitter

    output_file_path = Path(output_file_name)
    if not output_file_path.parent.exists():
        output_file_path.parent.mkdir(parents=True, exists_ok=True)

    fitter = MongeJetFitter(polyhedral_surface)
    mio_mesh = fitter.to_meshio_mesh(result_array)
    meshio.write(output_file_name, mio_mesh)


def results_to_mesh_cli(args):
    """
    Command line interface for reading polynomial fitting results file
    and generating a mesh file containing the original mesh along with
    the fitting-results *point data*.

    :type args: :obj:`types.SimpleNamespace`
    :param args: Parsed command line arguments (e.g. as returned
       by :samp:`surface_poly_fit.results_to_mesh.get_argument_parser().parse_args()`.
    """
    import logging
    from pathlib import Path
    import numpy as np

    logging_format = '%(asctime)s|%(process)-8s|%(name)-8s|%(levelname)-8s|%(message)s'
    logging.basicConfig(format=logging_format, level=getattr(logging, args.log_level))
    logger = logging.getLogger()

    output_file_path = args.output_file
    if output_file_path is None:
        output_file_path = \
            Path(Path(args.surface_poly_fit_results_file).stem + "_nr%(num_rings)03d")
        output_file_path = output_file_path.with_suffix(".vtu")
    else:
        output_file_path = Path(output_file_path)

    logger.info("Loading result data from file %s...", args.surface_poly_fit_results_file)
    poly_surface, result_array = read_results(args.surface_poly_fit_results_file)

    num_rings_list = sorted(np.unique(result_array["num_rings"]).tolist())
    for num_rings in num_rings_list:
        nr_output_file_path = Path(str(output_file_path) % {"num_rings": num_rings})
        logger.info("Writing mesh and result data to file %s...", nr_output_file_path)
        write_result_mesh(
            nr_output_file_path,
            poly_surface,
            result_array[result_array["num_rings"] == num_rings]
        )


def get_argument_parser():
    """
    Returns a :obj:`argparse.ArgumentParser` to handle command line option processing.

    :rtype: :obj:`argparse.ArgumentParser`
    :return: Argument parser.
    """
    import argparse
    ap = \
        argparse.ArgumentParser(
            "surface_poly_fit_r2m",
            description=(
                "Write mesh and polynomial-fit results point-data to mesh-file."
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
        help=(
            "Substitution string used for the name of the output 'mesh' files"
            +
            " to be generated with polynomial fitting results *point data*."
            +
            " Default is to replace the extension and add suffixes to input '.npz' file name."
            +
            " One mesh-file per 'num_rings', e.g. 'output_mesh_file_name_nr%(num_rings)%03ds.vtu'."
        ),
        default=None
    )
    ap.add_argument(
        "surface_poly_fit_results_file",
        help="A polynomial fit results '.npz' file."
    )

    return ap


def results_to_mesh_main_cli():
    """
    Entry-point function for command line interface.
    """
    results_to_mesh_cli(get_argument_parser().parse_args())


if __name__ == "__main__":
    results_to_mesh_main_cli()
