#ifndef SURFACE_POLY_FIT_MONGE_JET_FITTING_H
#define SURFACE_POLY_FIT_MONGE_JET_FITTING_H

#include "surface_poly_fit/polyhedral_surface.h"
#include <Eigen/Dense>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Eigen_svd.h>
#include <cstdint>

namespace spf
{

template <typename TPoly, typename TLocalFloatType=double>
class MongeJetFitter
{
public:
  typedef TPoly PolySurf;
  typedef typename PolySurf::Point_3 Point_3;
  typedef TLocalFloatType LocalFloatType;
  typedef CGAL::Simple_cartesian<TLocalFloatType> LocalKernel;
  typedef typename TPoly::Traits::Kernel DataKernel;
  typedef CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel> MongeViaJetFitting;
  typedef typename LocalKernel::Vector_3 Vector_3;
  typedef Eigen::Matrix<typename LocalKernel::FT, 3, 3> Matrix_3x3;

  struct MongeForm: public MongeViaJetFitting::Monge_form
  {
  public:
    typedef typename MongeViaJetFitting::Monge_form Inherited;

    MongeForm() :
      Inherited(),
      vertex_index_(-1),
      num_rings_(-1),
      poly_fit_condition_number_(0.0),
      pca_eigenvalues_(0.0, 0.0, 0.0),
      pca_eigenvectors_(Matrix_3x3::Zero())
    {
      this->m_origin_pt = typename Inherited::Point_3(0.0, 0.0, 0.0);
      this->m_d1 = typename Inherited::Vector_3(0.0, 0.0, 0.0);
      this->m_d2 = typename Inherited::Vector_3(0.0, 0.0, 0.0);
      this->m_n = typename Inherited::Vector_3(0.0, 0.0, 0.0);
      this->m_coefficients = std::vector<typename Inherited::FT>(11, 0.0);
    }

    MongeForm(Inherited const & mf) :
      Inherited(mf),
      vertex_index_(-1),
      num_rings_(-1),
      poly_fit_condition_number_(0.0),
      pca_eigenvalues_(0.0, 0.0, 0.0),
      pca_eigenvectors_(Matrix_3x3::Zero())
    {
    }

    std::int64_t vertex_index_;
    std::uint8_t degree_monge_;
    std::uint8_t degree_poly_fit_;
    std::int64_t num_rings_;
    std::int64_t num_fitting_points_;
    LocalFloatType poly_fit_condition_number_;
    Vector_3 pca_eigenvalues_;
    Matrix_3x3 pca_eigenvectors_;
  };

  MongeJetFitter(const std::size_t degree_poly_fit=2, const std::size_t degree_monge=2) :
    in_points_(),
    degree_poly_fit_(degree_poly_fit),
    degree_monge_(degree_monge)
  {
  }

  std::size_t get_min_num_fit_points() const
  {
    return ((this->degree_poly_fit_ + 1) * (this->degree_poly_fit_ + 2)) / 2;
  }

  MongeForm fit_at_vertex(
      PolySurf & poly_surface,
      const std::int64_t vertex_index,
      const std::int64_t num_rings
  )
  {
    auto vtx_it = poly_surface.vertices_begin() + vertex_index;
    poly_surface.gather_fitting_points(&(*vtx_it), num_rings,  this->in_points_, poly_surface.vertex_prop_map_);

    MongeForm monge_form;
    MongeViaJetFitting monge_fit;
    if (this->in_points_.size() >= this->get_min_num_fit_points())
    {
      monge_form =
          MongeForm(
              monge_fit(
                  this->in_points_.begin(),
                  this->in_points_.end(),
                  this->degree_poly_fit_,
                  this->degree_monge_
              )
          );
      monge_form.comply_wrt_given_normal(vtx_it->normal);
      monge_form.poly_fit_condition_number_ = monge_fit.condition_number();
      monge_form.pca_eigenvalues_ =
          Vector_3(
              monge_fit.pca_basis(0).first,
              monge_fit.pca_basis(1).first,
              monge_fit.pca_basis(2).first
          );
    }
    monge_form.vertex_index_ = vertex_index;
    monge_form.degree_monge_ = this->degree_monge_;
    monge_form.degree_poly_fit_ = this->degree_poly_fit_;
    monge_form.num_rings_ = num_rings;
    monge_form.num_fitting_points_ = this->in_points_.size();

    return monge_form;
  }

  std::vector<Point_3> in_points_;  //container for data points
  std::size_t degree_poly_fit_;
  std::size_t degree_monge_;
};

}

#endif
