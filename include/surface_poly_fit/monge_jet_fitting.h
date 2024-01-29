#ifndef SURFACE_POLY_FIT_MONGE_JET_FITTING_H
#define SURFACE_POLY_FIT_MONGE_JET_FITTING_H

#include "surface_poly_fit/polyhedral_surface.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
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
  typedef typename LocalKernel::Vector_3 Vector_3;
  typedef Eigen::Matrix<typename LocalKernel::FT, 3, 3> Matrix_3x3;
  typedef Eigen::Matrix<typename LocalKernel::FT, 3, 1> Matrix_3x1;
  typedef Eigen::Quaternion<typename LocalKernel::FT> Quaternion;

  enum FittingBasisType {
    PCA = 0,
    VERTEX_NORMAL = 1,
    RING_NORMAL_MEAN = 2,
    RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN = 3,
    RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA = 4
  };

  struct MongeForm: public CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel>::Monge_form
  {
  public:
    typedef typename CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel>::Monge_form Inherited;

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

  class MongeViaJetFitting: public CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel>
  {
  public:
    typedef CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel> Inherited;
    typedef typename Inherited::Monge_form Monge_form;
    typedef typename Inherited::LAMatrix LAMatrix;
    typedef typename Inherited::LAVector LAVector;
    typedef typename Inherited::Point_3 Point_3;
    typedef typename Inherited::Vector_3 Vector_3;
    typedef typename Inherited::Aff_transformation Aff_transformation;

    template <class InputIterator>
    std::pair<Matrix_3x3, Matrix_3x1> calculate_PCA_basis(InputIterator begin, InputIterator end) const
    {
      int n = std::distance(begin, end);
      typename LocalKernel::FT x, y, z,
        sumX = 0., sumY = 0., sumZ = 0.,
        sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
        sumXY = 0., sumXZ = 0., sumYZ = 0.,
        xx, yy, zz, xy, xz, yz;

      for (; begin != end; begin++)
        {
          Point_3 lp = D2L_converter(*begin);
          x = lp.x();
          y = lp.y();
          z = lp.z();
          sumX += x / n;
          sumY += y / n;
          sumZ += z / n;
          sumX2 += x * x / n;
          sumY2 += y * y / n;
          sumZ2 += z * z / n;
          sumXY += x * y / n;
          sumXZ += x * z / n;
          sumYZ += y * z / n;
        }
      xx = sumX2 - sumX * sumX;
      yy = sumY2 - sumY * sumY;
      zz = sumZ2 - sumZ * sumZ;
      xy = sumXY - sumX * sumY;
      xz = sumXZ - sumX * sumZ;
      yz = sumYZ - sumY * sumZ;

      // assemble covariance matrix as a
      // semi-definite matrix.
      // Matrix numbering:
      // 0 1 2
      //   3 4
      //     5
      std::array<typename LocalKernel::FT, 6> covariance = {{ xx,xy,xz,yy,yz,zz }};
      std::array<typename LocalKernel::FT, 3> eigen_values = {{ 0., 0., 0. }};
      std::array<typename LocalKernel::FT, 9> eigen_vectors = {{ 0., 0., 0., 0., 0., 0., 0., 0., 0. }};

      // solve for eigenvalues and eigenvectors.
      // eigen values are sorted in ascending order,
      // eigen vectors are sorted in accordance.
      CGAL::Default_diagonalize_traits<typename LocalKernel::FT, 3>::diagonalize_selfadjoint_covariance_matrix
        (covariance, eigen_values, eigen_vectors);

      //store in evals
      Matrix_3x1 evals;
      for (int i=0; i<3; i++)
      {
        evals(i) =  eigen_values[2-i];
      }

      Vector_3 v1(eigen_vectors[6],eigen_vectors[7],eigen_vectors[8]);
      Vector_3 v2(eigen_vectors[3],eigen_vectors[4],eigen_vectors[5]);
      Vector_3 v3(eigen_vectors[0],eigen_vectors[1],eigen_vectors[2]);
      this->switch_to_direct_orientation(v1, v2, v3);

      return
          std::make_pair(
            Matrix_3x3(
                v1[0], v2[0], v3[0],
                v1[1], v2[1], v3[1],
                v1[2], v2[2], v3[2]
            ),
            evals
          );
    }

    Matrix_3x3 calc_rotation_matrix_from_normal(Vector_3 const & normal) const
    {
      const Matrix_3x1 z_dir(0.0, 0.0, 0.0);
      const Matrix_3x1 nrml_dir(normal[0], normal[1], normal[2]);
      Matrix_3x3 R;
      R = Quaternion().setFromTwoVectors(z_dir, nrml_dir).toRotationMatrix();

      return R;
    }

    void set_world_to_fitting_from_normal(Vector_3 const & fit_normal)
    {
      this->set_world_to_fitting_from_rotation(
          this->calc_rotation_matrix_from_normal(fit_normal)
      );
    }

    void set_world_to_fitting_from_rotation(Matrix_3x3 const & R)
    {
      Matrix_3x3 const Rt = R.transpose();

      Aff_transformation
        change_basis (Rt(0, 0), Rt(0, 1), Rt(0, 2),
                      Rt(1, 0), Rt(1, 1), Rt(1, 2),
                      Rt(2, 0), Rt(2, 1), Rt(2, 2));

      this->change_world2fitting = change_basis;
    }

    template <class InputIterator>
    MongeForm operator()(
        InputIterator begin,
        InputIterator end,
        std::size_t d,
        std::size_t dprime,
        const bool doPCA=true
    )
    {
      // precondition: on the degrees, jet and monge
      CGAL_precondition( (d >=1) && (dprime >= 1)
                         && (dprime <= 4) && (dprime <= d) );
      this->deg = static_cast<int>(d);
      this->deg_monge = static_cast<int>(dprime);
      this->nb_d_jet_coeff = static_cast<int>((d+1)*(d+2)/2);
      this->nb_input_pts = static_cast<int>(end - begin);
      // precondition: solvable ls system
      CGAL_precondition( this->nb_input_pts >= this->nb_d_jet_coeff );

      //Initialize
      Monge_form monge_form;
      monge_form.set_up(dprime);
      //for the system MA=Z
      LAMatrix M(this->nb_input_pts, this->nb_d_jet_coeff);
      LAVector Z(this->nb_input_pts);

      if (doPCA)
      {
        this->compute_PCA(begin, end);
      }
      this->fill_matrix(begin, end, d, M, Z);//with precond
      this->solve_linear_system(M, Z);  //correct with precond
      this->compute_Monge_basis(Z.vector(), monge_form);
      if ( dprime >= 3) this->compute_Monge_coefficients(Z.vector(), dprime, monge_form);

      MongeForm ret_monge_form(monge_form);
      ret_monge_form.poly_fit_condition_number_ = this->condition_number();
      if (doPCA)
      {
        ret_monge_form.pca_eigenvalues_ =
            Vector_3(
                this->pca_basis(0).first,
                this->pca_basis(1).first,
                this->pca_basis(2).first
            );
      }

      return ret_monge_form;
    }

    template <class InputIterator>
    MongeForm operator()(
        InputIterator begin,
        InputIterator end,
        std::size_t d,
        std::size_t dprime,
        Vector_3 const & fit_normal
    )
    {
      this->set_world_to_fitting_from_normal(fit_normal);
      return (*this)(begin, end, d, dprime, false);
    }

  };


  MongeJetFitter(const std::size_t degree_poly_fit=2, const std::size_t degree_monge=2) :
    in_points_(),
    in_normals_(),
    in_ring_(),
    degree_poly_fit_(degree_poly_fit),
    degree_monge_(degree_monge)
  {
  }

  std::size_t get_min_num_fit_points() const
  {
    return ((this->degree_poly_fit_ + 1) * (this->degree_poly_fit_ + 2)) / 2;
  }

  void set_monge_form_data(
      MongeForm & mf,
      const std::int64_t vertex_index,
      const std::int64_t num_rings
  )
  {
    mf.vertex_index_ = vertex_index;
    mf.degree_monge_ = this->degree_monge_;
    mf.degree_poly_fit_ = this->degree_poly_fit_;
    mf.num_rings_ = num_rings;
    mf.num_fitting_points_ = this->in_points_.size();
  }

  void gather_fitting_points(
      PolySurf & poly_surface,
      const std::int64_t vertex_index,
      const std::int64_t num_rings
  )
  {
    auto vtx_it = poly_surface.vertices_begin() + vertex_index;
    poly_surface.gather_fitting_points(
        &(*vtx_it),
        num_rings,
        this->in_points_,
        this->in_normals_,
        this->in_ring_,
        poly_surface.vertex_prop_map_
    );
  }

  MongeForm fit_at_vertex(
      PolySurf & poly_surface,
      const std::int64_t vertex_index,
      const std::int64_t num_rings
  )
  {
    this->gather_fitting_points(poly_surface, vertex_index, num_rings);

    MongeForm monge_form;
    MongeViaJetFitting monge_fit;
    if (this->in_points_.size() >= this->get_min_num_fit_points())
    {
      monge_form =
        monge_fit(
          this->in_points_.begin(),
          this->in_points_.end(),
          this->degree_poly_fit_,
          this->degree_monge_,
          true
        );
      monge_form.comply_wrt_given_normal(this->in_normals_[0]);
    }
    this->set_monge_form_data(monge_form, vertex_index, num_rings);

    return monge_form;
  }

  std::vector<Point_3> in_points_;  //container for data points
  std::vector<Vector_3> in_normals_;  //container for data point normals
  std::vector<std::int32_t> in_ring_;  //container for point ring-number
  std::size_t degree_poly_fit_;
  std::size_t degree_monge_;
};

}

#endif
