
#include "surface_poly_fit/polyhedral_surface_py.h"
#include "surface_poly_fit/polyhedral_surface.h"

namespace spf
{

namespace py = pybind11;

class PolyhedraSurfacePy
{
public:
    PolyhedraSurfacePy()
      : surface_()
    {
    }

    PolyhedralSurface surface_;
};


void export_polyhedral_surface(pybind11::module_ m)
{
    py::class_<PolyhedraSurfacePy>(m, "PolyhedralSurface")
        .def(py::init<>())
    ;
}

}
