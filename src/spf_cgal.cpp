
#include <surface_poly_fit/spf_cgal.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/property_map.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>

#include <vector>
#include <stdexcept>
#include <fstream>
#include <iostream>
// types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace spf
{

}
