#ifndef SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_H
#define SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_H

#include "surface_poly_fit/polyhedral_surface_ops.h"
#include "surface_poly_fit/polyhedral_surface_rings.h"

#include <CGAL/Handle_hash_function.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/property_map.h>
#include <boost/graph/properties.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>


#include <algorithm>
#include <vector>
#include <list>
#include <memory>

#include <cstdlib>
#include <cstdio>
#include <functional>

namespace spf
{

//---------------------------------------------------------------- A
//redefined items class for the Polyhedron_3 with a refined vertex
//class that contains nothing more! (the _ring_tag is given by an
//externa std::map; a refined facet with a normal vector instead of
//the plane equation (this is an alternative solution instead of using
//Polyhedron_traits_with_normals_3). edges with the length
//----------------------------------------------------------------

template < class Refs, class Tag, class Pt, class FGeomTraits >
class PsVertex:public CGAL::HalfedgeDS_vertex_base < Refs, Tag, Pt >
{
typedef typename FGeomTraits::Point_3 Point_3;
typedef typename FGeomTraits::Vector_3 Vector_3;

public:
 PsVertex(const std::int64_t idx, const Point_3 & pt):
   CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt),
   index(idx),
   normal(0.0, 0.0, 0.0)
  {
  }

 PsVertex(const Point_3 & pt):
   CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt),
   index(-1),
   normal(0.0, 0.0, 0.0)
  {
  }

  PsVertex() :
    CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (),
    index(-1),
    normal(0.0, 0.0, 0.0)
  {
  }

  // Vertex index.
  std::int64_t index;
  // Vertex normal
  Vector_3 normal;
};

//----------------------------------------------------------------
// Facet with normal and possibly more types. types are recovered
//from the FGeomTraits template arg
//----------------------------------------------------------------
template < class Refs, class Tag, class FGeomTraits >
class PsFacet:public CGAL::HalfedgeDS_face_base < Refs, Tag >
{
public:
 typedef typename FGeomTraits::Vector_3 Vector_3;

protected:
  Vector_3 normal;
  //int ring_index;

public:
  const Vector_3& get_unit_normal() const { return normal; }
  Vector_3& get_unit_normal() { return normal; }

  //PsFacet(): ring_index(-1) {}
  //void setNormal(Vector_3  n) { normal = n; }
//   //this is for collecting i-th ring neighbors
//   void setRingIndex(int i) { ring_index = i; }
//   int getRingIndex() { return ring_index; }
//   void resetRingIndex() { ring_index = -1; }
};


template <class TPoly>
class Facet_PM :
  public boost::put_get_helper<typename TPoly::Traits::Vector_3&, Facet_PM<TPoly> >
{
public:

  //read_write
  typedef boost::lvalue_property_map_tag category;
  typedef typename TPoly::Facet key_type;
  typedef typename TPoly::Traits::Vector_3 value_type;
  typedef typename TPoly::Traits::Vector_3& reference;

  Facet_PM() {}
  reference operator[](key_type f) const {return f.get_unit_normal();}
};

} // namespace spf

//XFC: we should have Facet instead of PsVertex!
namespace boost{
  enum vertex_attribute_t        { vertex_attribute        = 1111 };
  //BOOST_INSTALL_PROPERTY(facet, attribute);
  BOOST_INSTALL_PROPERTY(vertex, attribute);

}

namespace spf
{

template <class TPoly>
Facet_PM<TPoly> get_fpm(boost::vertex_attribute_t, TPoly& ) {return Facet_PM<TPoly>();}



//----------------------------------------------------------------
// Halfedge
//----------------------------------------------------------------
 //  int ring_index;
//   PsHalfEdge(): ring_index(-1) {}
//   void setRingIndex(int i) {        ring_index = i;    }
//   int getRingIndex() {return ring_index;    }
//   void resetRingIndex() {ring_index = -1;    }

template < class Refs, class Tprev, class Tvertex, class Tface>
class PsHalfEdge:public CGAL::HalfedgeDS_halfedge_base < Refs, Tprev, Tvertex, Tface >
{
public:
  double len;
public:
  PsHalfEdge(): len(-1) {}
  double& get_length()  { return len; }
};

/*XFC: tentative ... failed so far...*/
//property map associated to the half edge
template <class TPoly>
class HEdge_PM :
  public boost::put_get_helper<typename TPoly::Traits::FT&, HEdge_PM<TPoly> >//double
{
public:
  typedef boost::lvalue_property_map_tag category;
  typedef typename TPoly::Halfedge key_type;
  typedef typename TPoly::Traits::FT value_type;
  typedef typename TPoly::Traits::FT& reference;

  HEdge_PM() {}
  reference operator[](key_type h) const {return h.len;}
};

//use the std edge_weight_t tag...
template <class TPoly>
HEdge_PM<TPoly> get_hepm(boost::edge_weight_t, TPoly& )
{return HEdge_PM<TPoly>();}


//------------------------------------------------
// Wrappers [Vertex, Face, Halfedge]
//------------------------------------------------
struct Wrappers_VFH:public CGAL::Polyhedron_items_3 {
  // wrap vertex
  template < class Refs, class Traits > struct Vertex_wrapper {
    typedef struct FGeomTraits {
    public:
      typedef typename Traits::Point_3 Point_3;
      typedef typename Traits::Vector_3 Vector_3;
    } FGeomTraits;
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::Vector_3 Vector_3;
    typedef PsVertex < Refs, CGAL::Tag_true, Point_3, FGeomTraits > Vertex;
  };

  // wrap face
  //NOTE: [HDS, Face] renamed [Polyhedron, Facet]
  template < class Refs, class Traits > struct Face_wrapper {
    //typedef typename Traits::Vector_3 Vector_3;
    //all types needed by the facet...
    typedef struct FGeomTraits {
    public:
      typedef typename Traits::Vector_3 Vector_3;
    } FGeomTraits;
    //custom type instantiated...
    typedef PsFacet < Refs, CGAL::Tag_true, FGeomTraits > Face;
  };

  // wrap halfedge
  template < class Refs, class Traits > struct Halfedge_wrapper {
   typedef PsHalfEdge < Refs,
      CGAL::Tag_true,
      CGAL::Tag_true, CGAL::Tag_true>  Halfedge;
  };
};

template <typename THalfedgeDS>
class PolyhedralPatchBuilder: public CGAL::Modifier_base<THalfedgeDS>
{
public:
  typedef CGAL::Modifier_base<THalfedgeDS> Inherited;
  typedef THalfedgeDS HalfedgeDS;
  typedef typename HalfedgeDS::Vertex Vertex;
  typedef typename HalfedgeDS::Vertex_handle Vertex_handle;
  typedef typename HalfedgeDS::Vertex::Point Point;
  typedef std::vector < Vertex const * > VertexStlVec;
  typedef typename HalfedgeDS::Halfedge Halfedge;
  typedef typename Halfedge::Facet Facet;
  typedef typename Halfedge::Facet_handle Facet_handle;
  typedef boost::unordered_map<std::int64_t, std::int64_t> VertexIndexToNewIndexMap;
  typedef std::vector<Facet_handle> FacetHandleStlVec;
  typedef CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> Builder;
  typedef boost::unordered_set<Facet const *, CGAL::Handle_hash_function> FacetHandleSet;

  PolyhedralPatchBuilder(
      VertexStlVec const & vertices
  )
    : Inherited(),
      vertices_(vertices)
  {
  }

  void operator()(HalfedgeDS & hds)
  {
    const std::int64_t num_vertices = this->vertices_.size();

    std::vector<std::vector<std::int64_t> > faces;
    // Calculate the number of faces.

    {
      VertexIndexToNewIndexMap vtxIdxToNewIdxMap;

      // initialise the vtxIdxToNewIdxMap.
      for (std::int64_t i = 0; i < num_vertices; ++i)
      {
        auto orig_vtx_handle = this->vertices_[i];
        vtxIdxToNewIdxMap[orig_vtx_handle->index] = i;
      }

      // Work out the faces.
      FacetHandleSet facetHdlSet;
      std::vector<std::int64_t> vertex_indices;

      for (std::int64_t i = 0; i < num_vertices; ++i)
      {
        auto orig_vtx_handle = this->vertices_[i];
        auto he_it = orig_vtx_handle->vertex_begin();
        do
        {
          auto facetHdl = he_it->facet();
          if (facetHdlSet.find(&(*facetHdl)) == facetHdlSet.end())
          {
            facetHdlSet.insert(&(*facetHdl));
            auto fct_he_it = facetHdl->facet_begin();
            vertex_indices.clear();
            do
            {
              auto new_vtx_idx_it = vtxIdxToNewIdxMap.find(fct_he_it->vertex()->index);
              if (new_vtx_idx_it != vtxIdxToNewIdxMap.end())
              {
                vertex_indices.push_back(new_vtx_idx_it->second);
              }
              else
              {
                break;
              }
            }
            while (++fct_he_it != facetHdl->facet_begin());

            if (facetHdl->facet_degree() == vertex_indices.size())
            {
              faces.push_back(vertex_indices);
            }
          }
        } while (++he_it != orig_vtx_handle->vertex_begin());
      }
    }

    Builder bldr(hds, true);

    bldr.begin_surface(num_vertices, faces.size(), 0, Builder::ABSOLUTE_INDEXING);
    {
      // Add vertices.
      for (std::int64_t i = 0; i < num_vertices; ++i)
      {
        auto orig_vtx_handle = this->vertices_[i];
        auto new_vertex_handle = bldr.add_vertex(orig_vtx_handle->point());
        new_vertex_handle->index = i;
        new_vertex_handle->normal = orig_vtx_handle->normal;
      }
    }

    // Add facets.
    for (auto face_it = faces.begin(); face_it != faces.end(); ++face_it)
    {
      bldr.add_facet(face_it->begin(), face_it->end());
    }
    bldr.end_surface();
  }

  VertexStlVec const & vertices_;
};



//------------------------------------------------
//PolyhedralSurface
//------------------------------------------------
typedef double                DFT;
typedef CGAL::Simple_cartesian<DFT>  Data_Kernel;
typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH,  CGAL::HalfedgeDS_vector > Polyhedron;
typedef Data_Kernel::Vector_3 Vector_3;

inline
std::size_t hash_value(boost::graph_traits<Polyhedron>::face_descriptor const & fd)
{
  return CGAL::Handle_hash_function()(fd);
}

} // namespace spf


namespace boost
{

inline
std::size_t hash_value(boost::graph_traits<spf::Polyhedron>::face_descriptor const & fd)
{
  return CGAL::Handle_hash_function()(fd);
}

}

namespace std
{
// Custom specialization of std::hash can be injected in namespace std.
template<>
struct hash<boost::graph_traits<spf::Polyhedron>::face_descriptor>
{
  std::size_t operator()(const boost::graph_traits<spf::Polyhedron>::face_descriptor& fd) const noexcept
  {
    return CGAL::Handle_hash_function()(fd);
  }
};

}

namespace spf
{

class PolyhedralSurface: public Polyhedron
{
public:
  typedef std::unique_ptr<PolyhedralSurface> PolyhedralSurfacePtr;

  typedef Traits::Point_3 Point_3;

  struct Hedge_cmp{
    bool operator()(Halfedge_handle a,  Halfedge_handle b) const{
      return &*a < &*b;
    }
  };

  struct Facet_cmp{
    bool operator()(Facet_handle a, Facet_handle b) const{
      return &*a < &*b;
    }
  };

  //Vertex property map, with boost::unordered_map
  typedef std::unordered_map<Vertex const *, int, CGAL::Handle_hash_function> Vertex2int_map_type;
  typedef boost::associative_property_map< Vertex2int_map_type > Vertex_PM_type;
  typedef PolyhedralSurfaceRings<PolyhedralSurface, Vertex_PM_type > Poly_rings;

  //Hedge property map, with boost::unordered_map
  typedef boost::unordered_map<Halfedge_handle, double, CGAL::Handle_hash_function> Hedge2double_map_type;
  typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
  typedef PolyhedralSurfaceHedgeOps<PolyhedralSurface, Hedge_PM_type> Poly_hedge_ops;

  //Facet property map, with boost::unordered_map
  typedef boost::unordered_map<Facet_handle, Vector_3, CGAL::Handle_hash_function> Facet2normal_map_type;
  typedef boost::associative_property_map<Facet2normal_map_type> Facet_PM_type;
  typedef PolyhedralSurfaceFacetOps<PolyhedralSurface, Facet_PM_type> Poly_facet_ops;

  PolyhedralSurface() :
    facet_map_(),
    facet_prop_map_(facet_map_),
    hedge_map_(),
    hedge_prop_map_(hedge_map_),
    vertex_map_(),
    vertex_prop_map_(vertex_map_)
  {
  }

  void reset_vertex_prop_map()
  {
    //initialize the tag of all vertices to -1
    auto vitb = this->vertices_begin();
    auto vite = this->vertices_end();
    this->vertex_map_.clear();
    CGAL_For_all(vitb, vite) put(this->vertex_prop_map_, &(*vitb), -1);
  }

  void update_edge_lengths()
  {
    Poly_hedge_ops::compute_edges_length(*this, this->hedge_prop_map_);
  }

  void update_vertex_normals()
  {
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

    boost::unordered_map<vertex_descriptor, Traits::Vector_3, CGAL::Handle_hash_function> vnormals;
    CGAL::Polygon_mesh_processing::compute_vertex_normals(
      *(dynamic_cast<Polyhedron *>(this)),
      boost::make_assoc_property_map(vnormals)
    );

    for (auto vtx_it = this->vertices_begin(); vtx_it != this->vertices_end(); ++vtx_it)
    {
      vtx_it->normal = vnormals[vtx_it];
    }
  }

  void update_face_normals()
  {
    CGAL::Polygon_mesh_processing::compute_face_normals(
      *(dynamic_cast<Polyhedron *>(this)),
      this->facet_prop_map_
    );
  }

  void update_vertex_and_face_normals(const bool updateVertexNormals=true)
  {
    this->update_face_normals();
    if (updateVertexNormals)
    {
      this->update_vertex_normals();
    }
  }

  void update(const bool updateVertexNormals=true)
  {
    this->update_edge_lengths();
    this->update_vertex_and_face_normals(updateVertexNormals);
    this->reset_vertex_prop_map();
  }

  static
  void gather_fitting_points(
      Vertex const * v,
      std::size_t num_rings,
      std::vector<Point_3> & in_points,
      std::vector<Vector_3> & in_normals,
      std::vector<std::int32_t> & in_ring,
      Vertex_PM_type & vpm
  )
  {
    //container to collect vertices of v on the PolyhedralSurface
    std::vector<Vertex const *> gathered;
    //initialize
    in_points.clear();
    in_normals.clear();
    in_ring.clear();

    Poly_rings::collect_i_rings(v, num_rings, gathered, vpm, false);

    //store the gathered points
    in_points.reserve(gathered.size());
    in_normals.reserve(gathered.size());
    in_ring.reserve(gathered.size());
    std::vector<Vertex const *>::iterator
      itb = gathered.begin(), ite = gathered.end();
    CGAL_For_all(itb, ite)
    {
      in_points.push_back((*itb)->point());
      in_normals.push_back((*itb)->normal);
      in_ring.push_back(vpm[*itb]);
    }
    Poly_rings::reset_ring_indices(gathered, vpm);
  }

  PolyhedralSurfacePtr create_ring_patch(
      std::size_t vertex_index,
      std::size_t num_rings
  )
  {
    auto vtx_it = this->vertices_begin() + vertex_index;
    std::vector<Vertex const *> gathered;
    Poly_rings::collect_i_rings(&(*vtx_it), num_rings, gathered, this->vertex_prop_map_);

    PolyhedralPatchBuilder<HalfedgeDS> bldr(gathered);
    PolyhedralSurfacePtr polyPatchPtr = std::make_unique<PolyhedralSurface>();
    polyPatchPtr->delegate(bldr);
    const bool updateVertexNormals = false;
    polyPatchPtr->update(updateVertexNormals);

    return polyPatchPtr;
  }

  Facet2normal_map_type facet_map_;
  Facet_PM_type facet_prop_map_;
  Hedge2double_map_type hedge_map_;
  Hedge_PM_type hedge_prop_map_;
  Vertex2int_map_type vertex_map_;
  Vertex_PM_type vertex_prop_map_;
};

} // namespace spf

#endif
