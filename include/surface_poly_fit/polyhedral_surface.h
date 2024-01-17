#ifndef SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_H
#define SURFACE_POLY_FIT_POLYHEDRAL_SURFACE_H

#include "surface_poly_fit/polyhedral_surface_ops.h"
#include "surface_poly_fit/polyhedral_surface_rings.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/property_map.h>
#include <boost/graph/properties.hpp>


#include <algorithm>
#include <vector>
#include <list>

#include <cstdlib>
#include <cstdio>

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

//------------------------------------------------
//PolyhedralSurface
//------------------------------------------------


typedef double                DFT;
typedef CGAL::Simple_cartesian<DFT>  Data_Kernel;
typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH,  CGAL::HalfedgeDS_vector > Polyhedron;
typedef Data_Kernel::Vector_3 Vector_3;

class PolyhedralSurface: public Polyhedron
{
public:

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

  //Vertex property map, with std::map
  typedef std::map<Vertex*, int> Vertex2int_map_type;
  typedef boost::associative_property_map< Vertex2int_map_type > Vertex_PM_type;
  typedef PolyhedralSurfaceRings<PolyhedralSurface, Vertex_PM_type > Poly_rings;

  //Hedge property map, with std::map
  typedef std::map<Halfedge_handle, double, Hedge_cmp> Hedge2double_map_type;
  typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
  typedef PolyhedralSurfaceHedgeOps<PolyhedralSurface, Hedge_PM_type> Poly_hedge_ops;

  //Facet property map, with std::map
  typedef std::map<Facet_handle, Vector_3, Facet_cmp> Facet2normal_map_type;
  typedef boost::associative_property_map<Facet2normal_map_type> Facet_PM_type;
  typedef PolyhedralSurfaceFacetOps<PolyhedralSurface, Facet_PM_type> Poly_facet_ops;

  PolyhedralSurface() :
    facet_map_(),
    facet_prop_map_(facet_map_),
    hedge_map_(),
    hedge_prop_map_(hedge_map_)
  {
  }

  void update_edge_lengths()
  {
    Poly_hedge_ops::compute_edges_length(*this, this->hedge_prop_map_);
  }

  void update_vertex_and_face_normals()
  {
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

    std::map<face_descriptor,Traits::Vector_3> fnormals;
    std::map<vertex_descriptor,Traits::Vector_3> vnormals;
    CGAL::Polygon_mesh_processing::compute_normals(
      *(dynamic_cast<Polyhedron *>(this)),
      boost::make_assoc_property_map(vnormals),
      this->facet_prop_map_
    );

    for (auto vtx_it = this->vertices_begin(); vtx_it != this->vertices_end(); ++vtx_it)
    {
      vtx_it->normal = vnormals[vtx_it];
    }
  }

  void update()
  {
    this->update_edge_lengths();
    this->update_vertex_and_face_normals();
  }

  Facet2normal_map_type facet_map_;
  Facet_PM_type facet_prop_map_;
  Hedge2double_map_type hedge_map_;
  Hedge_PM_type hedge_prop_map_;
};

} // namespace spf

#endif
