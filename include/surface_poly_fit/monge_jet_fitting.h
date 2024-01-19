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
      pca_eigenvalues_(0.0),
      pca_eigenvectors_(0.0)
    {
    }

    LocalFloatType poly_fit_condition_number_;
    Vector_3 pca_eigenvalues_;
    Matrix_3x3 pca_eigenvectors_;
  };

  MongeJetFitter(const std::size_t degree_poly_fit=2, const std::size_t degree_monge=2) :
    degree_poly_fit_(degree_poly_fit),
    degree_monge_(degree_monge)
  {
  }

  std::size_t min_num_fit_points() const
  {
    return ((this->degree_poly_fit_ + 1) * (this->degree_poly_fit_ + 2)) / 2;
  }

  std::size_t degree_poly_fit_;
  std::size_t degree_monge_;
};

}

#endif
