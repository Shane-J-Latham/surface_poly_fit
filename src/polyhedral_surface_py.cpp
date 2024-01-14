
#include "surface_poly_fit/polyhedral_surface_py.h"
#include "surface_poly_fit/polyhedral_surface.h"

namespace spf
{

namespace py = pybind11;

class PolyhedralSurfacePy
{
public:
    PolyhedralSurfacePy()
      : surface_()
    {
    }

    PolyhedralSurfacePy(py::object vertices, py::object faces)
      : surface_()
    {
    }

    PolyhedralSurface surface_;
};


void export_polyhedral_surface(pybind11::module_ m)
{
    py::class_<PolyhedralSurfacePy>(m, "PolyhedralSurface")
        .def(py::init<>())
        .def(py::init<py::object, py::object>(), py::arg("vertices"), py::arg("faces"))
    ;
}

}
