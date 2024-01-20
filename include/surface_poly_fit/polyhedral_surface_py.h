#ifndef SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_PY_H
#define SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_PY_H

#include <cstdint>
#include <pybind11/pybind11.h>
#include "surface_poly_fit/polyhedral_surface.h"

namespace spf
{
namespace py = pybind11;

/// Python wrapper for PolyhedralSurface.
class PolyhedralSurfacePy
{
public:

  PolyhedralSurfacePy();

  PolyhedralSurfacePy(py::object vertices, py::object faces);

  std::int64_t get_num_vertices() const;

  py::object get_vertices();

  py::object get_vertex_normals();

  std::int64_t get_num_faces() const;

  py::object get_faces();

  py::object get_face_normals();

  PolyhedralSurface surface_;
};

void export_polyhedral_surface(py::module_ m);

}

#endif
