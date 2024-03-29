#ifndef SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_PY_H
#define SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_PY_H

#include <cstdint>
#include <memory>
#include <pybind11/pybind11.h>
#include "surface_poly_fit/polyhedral_surface.h"

namespace spf
{
namespace py = pybind11;

/// Python wrapper for PolyhedralSurface.
class PolyhedralSurfacePy
{
public:
  typedef std::unique_ptr<PolyhedralSurface> PolyhedralSurfacePtr;
  typedef std::unique_ptr<PolyhedralSurfacePy> PolyhedralSurfacePyPtr;

  PolyhedralSurfacePy();

  PolyhedralSurfacePy(PolyhedralSurfacePy const & other);

  PolyhedralSurfacePy(py::object vertices, py::object faces);

  std::int64_t get_num_vertices() const;

  py::object get_vertices();

  py::object get_vertex_normals();

  void set_vertex_normals(py::object normals);

  std::int64_t get_num_faces() const;

  py::object get_faces();

  py::object get_face_normals();

  PolyhedralSurfacePyPtr create_ring_patch(const std::int64_t vertex_index, const std::int64_t num_rings);

  py::object create_child_ring_patch(const std::int64_t vertex_index, const std::int64_t num_rings, py::object child_class);

  PolyhedralSurfacePtr surface_;
};

void export_polyhedral_surface(py::module_ m);

}

#endif
