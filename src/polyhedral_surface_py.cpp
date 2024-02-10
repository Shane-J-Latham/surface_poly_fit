
#include "surface_poly_fit/polyhedral_surface_py.h"
#include "surface_poly_fit/polyhedral_surface.h"

#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

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
    typedef typename THalfEdgeDS::Vertex Vertex;
    typedef typename THalfEdgeDS::Vertex::Point Point;


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
          auto vertex_handle = bldr.add_vertex(Point(p[0], p[1], p[2]));
          vertex_handle->index = i;
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


PolyhedralSurfacePy::PolyhedralSurfacePy()
  : surface_(std::make_unique<PolyhedralSurface>())
{
}

PolyhedralSurfacePy::PolyhedralSurfacePy(PolyhedralSurfacePy const & other)
  : surface_(std::make_unique<PolyhedralSurface>(*(other.surface_)))
{
}

PolyhedralSurfacePy::PolyhedralSurfacePy(py::object vertices, py::object faces)
  : surface_(std::make_unique<PolyhedralSurface>())
{
  PolyhedralSurfaceBuilderPy<PolyhedralSurface::HalfedgeDS> bldr(vertices, faces);
  this->surface_->delegate(bldr);
  this->surface_->update();
}

std::int64_t PolyhedralSurfacePy::get_num_vertices() const
{
  return std::int64_t(std::distance(this->surface_->vertices_begin(), this->surface_->vertices_end()));
}

py::object PolyhedralSurfacePy::get_vertices()
{
  typedef PolyhedralSurface::HalfedgeDS::Vertex Vertex;
  typedef Vertex::Point Point;

  std::size_t const shape[2]{std::size_t(this->get_num_vertices()), 3};
  auto ary = py::array_t<PolyhedralSurface::HalfedgeDS::Traits::Kernel::FT>(shape);

  {
    auto vtx_it=this->surface_->vertices_begin();
    std::size_t i = 0;
    for ( ; i < shape[0]; ++vtx_it, ++i)
    {
        Point const & pt = vtx_it->point();
        for (std::size_t j = 0; j < 3; j++)
        {
            ary.mutable_at(i, j) = pt[j];
        }
    }
  }
  return py::object(ary);
}

py::object PolyhedralSurfacePy::get_vertex_normals()
{
  typedef PolyhedralSurface::HalfedgeDS::Vertex Vertex;
  typedef PolyhedralSurface::Traits::Vector_3 Vector;
  typedef Vertex::Point Point;

  std::size_t const shape[2]{std::size_t(this->get_num_vertices()), 3};
  auto ary = py::array_t<PolyhedralSurface::HalfedgeDS::Traits::Kernel::FT>(shape);

  {
    auto vtx_it = this->surface_->vertices_begin();
    std::size_t i = 0;
    for ( ; i < shape[0]; ++vtx_it, ++i)
    {
        Vector const & nrml = vtx_it->normal;
        for (std::size_t j = 0; j < 3; j++)
        {
            ary.mutable_at(i, j) = nrml[j];
        }
    }
  }

  return py::object(ary);
}

void PolyhedralSurfacePy::set_vertex_normals(py::object normalsObj)
{
  typedef PolyhedralSurface::Traits::Vector_3 Vector;

  const std::size_t num_normals = py::len(normalsObj);
  if (num_normals != this->get_num_vertices())
  {
    std::stringstream msg;
    msg
      << "Got num_normals="  << num_normals
      << " != num_vertices=" << this->get_num_vertices();
    throw std::runtime_error(msg.str());
  }
  py::object j_objs[] = {py::cast(long(0)), py::cast(long(1)), py::cast(long(2))};
  auto vtx_it = this->surface_->vertices_begin();
  std::size_t i = 0;
  for ( ; (vtx_it != this->surface_->vertices_end()) && (i < num_normals); ++vtx_it, ++i)
  {
    PolyhedralSurface::Traits::Kernel::FT nrml[3];
    py::object i_obj(py::cast(i));
    for (std::int64_t j = 0; j < 3; ++j)
    {
      nrml[j] = py::cast<typename PolyhedralSurface::Traits::Kernel::FT>(normalsObj[i_obj][j_objs[j]]);
    }
    vtx_it->normal = Vector(nrml[0], nrml[1], nrml[2]);
  }
}

std::int64_t PolyhedralSurfacePy::get_num_faces() const
{
  return std::int64_t(std::distance(this->surface_->facets_begin(), this->surface_->facets_end()));
}

py::object PolyhedralSurfacePy::get_faces()
{
  typedef PolyhedralSurface::HalfedgeDS::Face Face;
  typedef PolyhedralSurface::HalfedgeDS::Vertex Vertex;

  bool faces_all_same_degree = true;
  const std::size_t face_degree = this->surface_->facets_begin()->facet_degree();
  {
    auto face_it = this->surface_->facets_begin();
    for (++face_it; face_it != this->surface_->facets_end(); ++face_it)
    {
      if (face_it->facet_degree() != face_degree)
      {
        faces_all_same_degree = false;
        break;
      }
    }
  }

  py::object return_ary;
  if (faces_all_same_degree)
  {
    std::size_t const shape[2]{std::size_t(this->get_num_faces()), face_degree};
    auto ary = py::array_t<std::int64_t>(shape);
    auto face_it = this->surface_->facets_begin();
    for (std::size_t i = 0; i < shape[0]; ++i, ++face_it)
    {
      auto halfedge_it = face_it->facet_begin();
      for (std::size_t j = 0; j < shape[1]; ++j, ++halfedge_it)
      {
        ary.mutable_at(i, j) = halfedge_it->vertex()->index;
      }
    }
    return_ary = py::object(ary);
  }
  else
  {
    py::list face_list;
    auto face_it = this->surface_->facets_begin();
    for (std::size_t i = 0; i < this->get_num_faces(); ++i, ++face_it)
    {
      py::list vertex_list;
      auto halfedge_it = face_it->facet_begin();
      for (std::size_t j = 0; j < face_it->facet_degree(); ++j, ++halfedge_it)
      {
        vertex_list.append(halfedge_it->vertex()->index);
      }
      face_list.append(vertex_list);
    }
    return_ary = py::object(face_list);
  }

  return return_ary;
}

py::object PolyhedralSurfacePy::get_face_normals()
{
  typedef PolyhedralSurface::HalfedgeDS::Face Face;
  typedef PolyhedralSurface::HalfedgeDS::Vertex Vertex;
  typedef PolyhedralSurface::HalfedgeDS::Traits::Vector_3 Vector;

  std::size_t const shape[2]{std::size_t(this->get_num_faces()), 3};
  auto ary = py::array_t<PolyhedralSurface::HalfedgeDS::Traits::Kernel::FT>(shape);

  {
    auto face_it = this->surface_->facets_begin();
    std::size_t i = 0;
    for (; face_it != this->surface_->facets_end(); ++face_it, ++i)
    {
      Vector const & nrml = this->surface_->facet_prop_map_[face_it];
      for (std::size_t j = 0; j < 3; j++)
      {
          ary.mutable_at(i, j) = nrml[j];
      }
    }
  }

  return py::object(ary);
}

PolyhedralSurfacePy::PolyhedralSurfacePyPtr PolyhedralSurfacePy::create_ring_patch(
    const std::int64_t vertex_index,
    const std::int64_t num_rings
)
{
  PolyhedralSurface::PolyhedralSurfacePtr ps_ptr = this->surface_->create_ring_patch(vertex_index, num_rings);

  PolyhedralSurfacePyPtr ret_psp_ptr = std::make_unique<PolyhedralSurfacePy>();
  ret_psp_ptr->surface_.swap(ps_ptr);

  return ret_psp_ptr;
}


/// Export PolyhedralSurfacePy class to specified python module.
void export_polyhedral_surface(pybind11::module_ m)
{
  py::class_<PolyhedralSurfacePy> ps_cls(m, "PolyhedralSurface", "A mesh consisting of polygonal faces.");

  ps_cls
    .def(py::init<>(), "Default construction, no vertices, no faces.")
    .def(
        py::init<py::object, py::object>(),
        py::arg("vertices"),
        py::arg("faces"),
        "Construct, specify vertices and faces.\n\n"
        ":type vertices: :obj:`numpy.ndarray`\n"
        ":param vertices: A :samp:`(N, 3)` shaped array of vertex coordinates.\n"
        ":type faces: :obj:`list` of :obj:`list`\n"
        ":param faces: Each element of :samp:`{faces}` is an *ordered* sequence"
        " of vertex indices defining a face.\n"
    )
    .def_property_readonly(
        "num_vertices",
        &PolyhedralSurfacePy::get_num_vertices,
        "An :obj:`int` indicating the number of vertices in this polyhedral surface.\n"
    )
    .def_property_readonly(
        "num_faces",
        &PolyhedralSurfacePy::get_num_faces,
        "An :obj:`int` indicating the number of faces in this polyhedral surface.\n"
    )
    .def(
        "get_vertices",
        &PolyhedralSurfacePy::get_vertices,
        "Returns the vertex coordinates of this polyhedral surface.\n\n"
        ":rtype: :obj:`numpy.ndarray`\n"
        ":return: A :samp:`(N, 3)` shaped :obj:`numpy.ndarray` of coordinates.\n"
    )
    .def(
        "get_vertex_normals",
        &PolyhedralSurfacePy::get_vertex_normals,
        "Returns the vertex normals of this polyhedral surface.\n\n"
        ":rtype: :obj:`numpy.ndarray`\n"
        ":return: A :samp:`(N, 3)` shaped :obj:`numpy.ndarray` of normals.\n"
    )
    .def(
        "set_vertex_normals",
        &PolyhedralSurfacePy::set_vertex_normals,
        py::arg("normals"),
        "Sets vertex normals.\n\n"
        ":type normals: :obj:`numpy.ndarray`\n"
        ":param normals: A :samp:`(N, 3)` shaped :obj:`numpy.ndarray`"
        " of vertex-normals. The :samp:`{normals}[i]` direction is assigned"
        " to vertex :samp:`i`.\n"
    )
    .def(
        "get_faces",
        &PolyhedralSurfacePy::get_faces,
        "Returns the faces of this polyhedral surface.\n\n"
        ":rtype: :obj:`list`\n"
        ":return: Each element of the returned list is an *ordered* sequence of"
        " vertex integers indicating the vertices (vertex indexes) which form"
        " the face."
    )
    .def(
        "get_face_normals",
        &PolyhedralSurfacePy::get_face_normals,
        "Returns the face normals of this polyhedral surface.\n\n"
        ":rtype: :obj:`numpy.ndarray`\n"
        ":return: A :samp:`(N, 3)` shaped :obj:`numpy.ndarray` of normals.\n"
    )
    .def(
        "create_ring_patch",
        &PolyhedralSurfacePy::create_ring_patch,
        py::arg("vertex_index"), py::arg("num_rings"),
        "Creates a polyhedral surface patch about a specified vertex.\n\n"
        ":type vertex_index: :obj:`int`\n"
        ":param vertex_index: The index of the vertex about which the patch is created.\n"
        ":type num_rings: :obj:`int`\n"
        ":param num_rings: The number of edge-hops from vertex :samp:`{vertex_index}`"
        " used to form the patch.\n"
        ":rtype: :obj:`PolyhedralSurface`\n"
        ":return: A new instance :obj:`PolyhedralSurface` patch.\n"
    )
  ;
}

}
