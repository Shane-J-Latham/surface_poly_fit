<!-- Badges-Start-->
![ReadTheDocs](https://readthedocs.org/projects/surface-poly-fit/badge/?version=latest&style=plastic)
![Python Test](https://github.com/Shane-J-Latham/surface_poly_fit/actions/workflows/python-test.yml/badge.svg)
![Python Test (vcpkg build)](https://github.com/Shane-J-Latham/surface_poly_fit/actions/workflows/python-test-vcpkg.yml/badge.svg)
![Release](https://img.shields.io/github/v/release/Shane-J-Latham/surface_poly_fit.svg?style=flat)
[![GitHub version](https://badge.fury.io/gh/Shane-J-Latham%2Fsurface_poly_fit.svg)](https://badge.fury.io/gh/Shane-J-Latham%2Fsurface_poly_fit)
<!-- Badges-Finish-->


# `surface_poly_fit`

Python package for fitting polynomials to mesh surface patches.
Under the hood the implementation uses [CGAL's](https://cgal.org) [Local Differential Properties Estimation](https://doc.cgal.org/latest/Jet_fitting_3/index.html#Chapter_Estimation_of_Local_Differential_Properties_of_Point-Sampled_Surfaces)
via [Monge Jet Fitting](https://doc.cgal.org/latest/Jet_fitting_3/classCGAL_1_1Monge__via__jet__fitting.html).

## Quick Start

### Python Package

```python
import meshio
from surface_poly_fit import PolyhedralSurface, MongeJetFitter

# Read in the mesh
mesh = meshio.read("my_mesh.ply")

# Create the polyhedral surface mesh data structure.
faces = sum(list(np.asanyarray(cells.data).tolist() for cells in mesh.cells), list())
polyhedral_surface = PolyhedralSurface(vertices=mesh.points, faces=faces)

# Create the object which performs the polynomial fitting to surface patches.
degree_poly_fit = 2  # Degree of fitting polynomial
degree_monge = 2  # Degree of Monge polynomial
fitter = MongeJetFitter(polyhedral_surface, degree_poly_fit, degree_monge)

# Fit a polynomial to the vertices in the 4-ring neighbourhood of a single vertex.
# Number of rings (num_rings) is the number of mesh-edge-hops,
# the more rings, the larger the surface patch.
result_array_vtx = fitter.fit_at_vertex(polyhedral_surface.num_vertices // 2, num_rings=4)

# Fit a polynomial in the 4-ring neighbourhood of each vertex.
result_array = fitter.fit_all(num_rings=4)

# Principle curvatures for all vertices of the surface-mesh:
result_array["k"]
```

### Command Line

Perform polynomial fitting at all vertices of mesh for different
neighbourhood sizes (``num_rings=2,4,8,16``):

```console
$ surface_poly_fit                           \
    --degree_poly_fit=2                      \
    --degree_monge=2                         \
    --num_rings=2,4,8,16                     \
    --output_file=my_mesh_surf_poly_fit.npz  \
    my_mesh.ply
```

The output file ``my_mesh_surf_poly_fit.npz`` is a [numpy](https://numpy.org)
[.npz](https://numpy.org/doc/stable/reference/generated/numpy.savez_compressed.html)
file containing entries: ``"surface_poly_fit_result"``, ``"vertices"``, ``"faces"`` and ``"vertex_normals"``. 

### Parallelism

[OpenMP](https://www.openmp.org/) threads are used for coarse parallelism
(each thread runs a *chunk* of vertex fitting). The number of threads is
controlled with the `OMP_NUM_THREADS` environment variable, a *maximum* number
of threads (depends on execution environment) is used by default.

## Installation

Install from latest github source:

```console
$ python -m pip install --user git+https://github.com/Shane-J-Latham/surface_poly_fit.git#egg=surface_poly_fit
```

or from source directory:

```console
$ git clone git@github.com/Shane-J-Latham/surface_poly_fit.git
$ cd surface_poly_fit
$ python -m pip install --user .
```

If you're on windows, you can use [vcpkg](https://github.com/microsoft/vcpkg) to
manage non-python dependencies (can also be used on Linux and MacOS):

```powershell
PS > git clone https://github.com/microsoft/vcpkg
PS > .\vcpkg\bootstrap-vcpkg.bat
PS > $Env:VCPKG_ROOT=$(Resolve-Path ./vcpkg)
PS > git clone git@github.com/Shane-J-Latham/surface_poly_fit.git
PS > cd surface_poly_fit
PS > python -m pip install --user .
```

You also still need to have build tools installed (some kind of C/C++ compiler).
One way to achieve this is to install Visual Studio Build tools. Visual studio
build tools likely require the installation of visual studio community edition first.
This link should (hopefully) get you started:

 https://visualstudio.microsoft.com/downloads/


## Requirements

Build requires:

- [pybind11](https://github.com/pybind)
- [eigen3](https://eigen.tuxfamily.org/)
- [boost](https://boost.org)
- [CGAL](https://cgal.org/)

Runtime requires:

- [python-3](https://www.python.org/doc/)
- [numpy](http://www.numpy.org/)
- [meshio](https://pypi.org/project/meshio/)

Unit-tests additionally require:

- [trimesh](https://trimesh.org)
- [scipy](https://scipy.org)


## Testing

Run unit-tests:

```console
$ python -m surface_poly_fit.tests
```

## Latest source code

Source at github:

https://github.com/Shane-J-Latham/surface_poly_fit

## Documentation

Documentation at [surface_poly_fit read-the-docs](https://surface-poly-fit.readthedocs.io/en/latest/).

Generate sphinx docs from the source tree using `sphinx_build` directly::

```console
$ sphinx_build -j 4 docs/source docs/_build/html
```

and browse at ``docs/_build/html/index.html``.

Documentation build requires:

- [mistune<1](https://pypi.org/project/mistune/0.8.4/) (version < 1 required for [m2r](https://pypi.org/project/m2r/))
- [m2r](https://pypi.org/project/m2r/)
- [sphinx_rtd_theme](https://pypi.org/project/sphinx_rtd_theme/)
- [sphinx-argparse](https://pypi.org/project/sphinx-argparse/)
- [myst_parser](https://pypi.org/project/myst-parser/)


## License information

See the file [LICENSE](https://github.com/Shane-J-Latham/surface_poly_fit/blob/main/LICENSE)
for terms & conditions, for usage and a DISCLAIMER OF ALL WARRANTIES.

