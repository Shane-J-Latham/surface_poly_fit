
#include "surface_poly_fit/monge_jet_fitting_py.h"
#include "surface_poly_fit/monge_jet_fitting.h"
#include "surface_poly_fit/polyhedral_surface.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace spf
{

namespace py = pybind11;

/// Python wrapper for MongeJetFitter.
class MongeJetFitterPy
{
public:
  typedef MongeJetFitter<PolyhedralSurface> MongeFitter;

  MongeJetFitterPy(py::object poly_surface, std::size_t degree_poly_fit, std::size_t degree_monge)
    : monge_fitter_(degree_poly_fit, degree_monge),
      poly_surface_obj_()
  {
  }

protected:
  MongeFitter monge_fitter_;
  py::object poly_surface_obj_;
};


/// Export MongeJetFitterPy class to specified python module.
void export_monge_jet_fitter(pybind11::module_ m)
{
  py::class_<MongeJetFitterPy>(m, "MongeJetFitter")
    .def(
        py::init<py::object, std::size_t, std::size_t>(),
        py::arg("poly_surface"),
        py::arg("degree_poly_fit")=2,
        py::arg("degree_monge")=2
    )
  ;
}

}
