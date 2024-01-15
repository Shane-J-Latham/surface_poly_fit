
#include "surface_poly_fit/polyhedral_surface_py.h"
#include "surface_poly_fit/polyhedral_surface.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace spf
{

namespace py = pybind11;

/// Populates a half-edge data stucture from vertices and faces
/// representation of a mesh.
template <typename THalfEdgeDS>
class PolyhedralSurfaceBuilderPy: public CGAL::Modifier_base<THalfEdgeDS>
{
public:
  typedef CGAL::Modifier_base<THalfEdgeDS> Inherited;

  PolyhedralSurfaceBuilderPy(py::object vertices, py::object faces)
    : Inherited(),
      vertices_(vertices),
      faces_(faces)
  {

  }

  void operator()(THalfEdgeDS & hds)
  {
    typedef CGAL::Polyhedron_incremental_builder_3<THalfEdgeDS> Builder;

    const std::int64_t num_vertices = py::len(this->vertices_);
    const std::int64_t num_faces = py::len(this->faces_);

    Builder bldr(hds, true);

    bldr.begin_surface(num_vertices, num_faces, 0, Builder::ABSOLUTE_INDEXING);

    {
      // Add vertices.
      py::object j_objs[] = {py::cast(long(0)), py::cast(long(1)), py::cast(long(2))};
      for (std::int64_t i = 0; i < num_vertices; ++i)
      {
          typename THalfEdgeDS::Traits::Kernel::FT p[3];
          py::object i_obj(py::cast(i));
          for (std::int64_t j = 0; j < 3; ++j)
          {
              p[j] = py::cast<typename THalfEdgeDS::Traits::Kernel::FT>(this->vertices_[i_obj][j_objs[j]]);
          }
          bldr.add_vertex(typename THalfEdgeDS::Vertex::Point(p[0], p[1], p[2]));
      }
    }

    // Add facets.
    for (std::int64_t i = 0; i < num_faces; ++i)
    {
      py::object i_obj(py::cast(i));
      py::object face_obj = this->faces_[i_obj];
      const std::int64_t num_face_verts = py::len(face_obj);

      bldr.begin_facet();
      for (std::size_t j = 0; j < num_face_verts; ++j)
      {
        const std::int64_t vertex_index = py::cast<std::int64_t>(face_obj[py::cast(j)]);
        bldr.add_vertex_to_facet(vertex_index);
      }
      bldr.end_facet();
    }

    bldr.end_surface();
  }

  py::object vertices_;
  py::object faces_;
};


/// Python wrapper for PolyhedralSurface.
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
    PolyhedralSurfaceBuilderPy<PolyhedralSurface::HalfedgeDS> bldr(vertices, faces);
    this->surface_.delegate(bldr);
  }

  std::int64_t get_num_vertices() const
  {
    return std::int64_t(std::distance(this->surface_.vertices_begin(), this->surface_.vertices_end()));
  }

  std::int64_t get_num_faces() const
  {
    return std::int64_t(std::distance(this->surface_.facets_begin(), this->surface_.facets_end()));
  }

  PolyhedralSurface surface_;
};


/// Export PolyhedralSurfacePy class to specified python module.
void export_polyhedral_surface(pybind11::module_ m)
{
  py::class_<PolyhedralSurfacePy>(m, "PolyhedralSurface")
    .def(py::init<>())
    .def(py::init<py::object, py::object>(), py::arg("vertices"), py::arg("faces"))
    .def_property_readonly("num_vertices", &PolyhedralSurfacePy::get_num_vertices)
    .def_property_readonly("num_faces", &PolyhedralSurfacePy::get_num_faces)
  ;
}

}
