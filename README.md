<!-- Badges-Start-->
![Python Test](https://github.com/Shane-J-Latham/surface_poly_fit/actions/workflows/python-test.yml/badge.svg)
![Python Test (vcpkg build)](https://github.com/Shane-J-Latham/surface_poly_fit/actions/workflows/python-test-vcpkg.yml/badge.svg)
<!-- Badges-Finish-->





# `surface_poly_fit`

Python package for fitting polynomials to mesh surface patches.
Under the hood the implementation uses [CGAL's](https://cgal.org) [Local Differential Properties Estimation](https://doc.cgal.org/latest/Jet_fitting_3/index.html#Chapter_Estimation_of_Local_Differential_Properties_of_Point-Sampled_Surfaces)
via [Monge Jet Fitting](https://doc.cgal.org/latest/Jet_fitting_3/classCGAL_1_1Monge__via__jet__fitting.html).

## Quick Start

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

