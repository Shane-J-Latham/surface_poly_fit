
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <surface_poly_fit/spf_cgal.h>
#include <surface_poly_fit/polyhedral_surface_py.h>
#include <surface_poly_fit/monge_jet_fitting_py.h>

namespace spf
{
}

PYBIND11_MODULE(_spf_cgal, m)
{
    spf::export_polyhedral_surface(m);
    spf::export_monge_jet_fitter(m);
}
