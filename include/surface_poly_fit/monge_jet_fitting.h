#ifndef SURFACE_POLY_FIT_MONGE_JET_FITTING_H
#define SURFACE_POLY_FIT_MONGE_JET_FITTING_H

#include "surface_poly_fit/polyhedral_surface.h"
#if SPF_HAVE_OPENMP
#include <omp.h>
#endif
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Eigen_svd.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/convex_hull_constructive_traits_2.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_2.h>

#include <cstdint>
#include <memory>
#include <vector>

namespace spf
{

template <typename TPoly, typename TLocalFloatType=double>
class MongeJetFitter
{
public:
  typedef TPoly PolySurf;
  typedef typename PolySurf::Vertex Vertex;
  typedef typename PolySurf::Vertex_handle Vertex_handle;
  typedef typename PolySurf::Point_3 Point_3;
  typedef TLocalFloatType LocalFloatType;
  typedef CGAL::Simple_cartesian<TLocalFloatType> LocalKernel;
  typedef typename TPoly::Traits::Kernel DataKernel;
  typedef typename LocalKernel::Vector_3 Vector_3;
  typedef Eigen::Matrix<typename LocalKernel::FT, 3, 3> Matrix_3x3;
  typedef Eigen::Matrix<typename LocalKernel::FT, 3, 1> Matrix_3x1;
  typedef Eigen::Quaternion<typename LocalKernel::FT> Quaternion;
  typedef CGAL::Eigen_svd CGALSvdTraits;

  enum FittingBasisType {
    PCA = 0,
    VERTEX_NORMAL = 1,
    RING_NORMAL_MEAN = 2,
    RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN = 3,
    RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA = 4
  };

  struct BoundingArea
  {
    BoundingArea() :
      rectangle_min_side_length(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      rectangle_max_side_length(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      circle_radius(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      ellipse_min_radius(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      ellipse_max_radius(std::numeric_limits<LocalFloatType>::quiet_NaN())
    {
    }

    template <typename TPoly2d>
    static
    BoundingArea from_rectangle(TPoly2d const & rect)
    {
      BoundingArea return_ba;
      return_ba.rectangle_min_side_length = std::numeric_limits<LocalFloatType>::max();
      return_ba.rectangle_max_side_length = std::numeric_limits<LocalFloatType>::min();
      auto prev_it = rect.vertices_begin();
      for (auto it = prev_it + 1; it != rect.vertices_end(); ++it, ++prev_it)
      {
        auto dist_sqrd = (*it - *prev_it).squared_length();
        if (dist_sqrd < return_ba.rectangle_min_side_length)
        {
          return_ba.rectangle_min_side_length = dist_sqrd;
        }
        if (dist_sqrd > return_ba.rectangle_max_side_length)
        {
          return_ba.rectangle_max_side_length = dist_sqrd;
        }

      }
      return_ba.rectangle_min_side_length = std::sqrt(return_ba.rectangle_min_side_length);
      return_ba.rectangle_max_side_length = std::sqrt(return_ba.rectangle_max_side_length);
      return return_ba;
    }

    LocalFloatType rectangle_min_side_length;
    LocalFloatType rectangle_max_side_length;
    LocalFloatType circle_radius;
    LocalFloatType ellipse_min_radius;
    LocalFloatType ellipse_max_radius;
  };

  struct ResidualStats
  {
    ResidualStats() :
      min(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      max(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      max_abs(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      mean(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      median(std::numeric_limits<LocalFloatType>::quiet_NaN()),
      stdd(std::numeric_limits<LocalFloatType>::quiet_NaN())
    {
    }

    LocalFloatType min;
    LocalFloatType max;
    LocalFloatType max_abs;
    LocalFloatType mean;
    LocalFloatType median;
    LocalFloatType stdd;
  };

  struct MongeForm: public CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel, CGALSvdTraits>::Monge_form
  {
  public:
    typedef typename CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel, CGALSvdTraits>::Monge_form Inherited;

    MongeForm() :
      Inherited(),
      vertex_index_(-1),
      num_rings_(-1),
      poly_fit_condition_number_(0.0),
      pca_eigenvalues_(0.0, 0.0, 0.0),
      fitting_basis_(Matrix_3x3::Zero()),
      fitting_residual_stats_()
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
      fitting_basis_(Matrix_3x3::Zero())
    {
    }

    std::int64_t vertex_index_;
    std::uint8_t degree_monge_;
    std::uint8_t degree_poly_fit_;
    std::int64_t num_rings_;
    std::int64_t num_fitting_points_;
    LocalFloatType poly_fit_condition_number_;
    Vector_3 pca_eigenvalues_;
    Matrix_3x3 fitting_basis_;
    ResidualStats fitting_residual_stats_;
    BoundingArea fitting_bounding_area_;
  };
  typedef std::vector<MongeForm> MongeFormStlVec;
  typedef std::unique_ptr<MongeFormStlVec> MongeFormStlVecPtr;

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

    MongeViaJetFitting() :
      Inherited(),
      residual_stats_(),
      bounding_area_()
    {
    }

    static
    void switch_to_direct_orientation(
      Vector_3 & v1,
      Vector_3 const & v2,
      Vector_3 const & v3)
    {
      if (CGAL::orientation (v1, v2, v3) == CGAL::NEGATIVE)
        v1 = -v1;
    }

    static
    void switch_to_direct_orientation(
      Matrix_3x1 & v1,
      Matrix_3x1 const & v2,
      Matrix_3x1 const & v3)
    {
      if (CGAL::orientation (v1, v2, v3) == CGAL::NEGATIVE)
        v1 = -v1;
    }

    template <class InputIterator>
    static
    std::pair<Matrix_3x3, Matrix_3x1> calculate_PCA_basis(InputIterator begin, InputIterator end)
    {
      CGAL::Cartesian_converter<DataKernel, LocalKernel> D2L_converter;
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
      Matrix_3x1 evals(eigen_values[2], eigen_values[1], eigen_values[0]);

      Vector_3 v1(eigen_vectors[6],eigen_vectors[7],eigen_vectors[8]);
      Vector_3 v2(eigen_vectors[3],eigen_vectors[4],eigen_vectors[5]);
      Vector_3 v3(eigen_vectors[0],eigen_vectors[1],eigen_vectors[2]);
      switch_to_direct_orientation(v1, v2, v3);

      Matrix_3x3 R;
      R <<
        v1[0], v2[0], v3[0],
        v1[1], v2[1], v3[1],
        v1[2], v2[2], v3[2]
        ;

      return std::make_pair(R, evals);
    }

    static
    Matrix_3x3 calc_rotation_matrix_from_normal(Vector_3 const & normal)
    {
      const Matrix_3x1 z_dir(0.0, 0.0, 1.0);
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
    void
    fill_matrix(InputIterator begin, InputIterator end,
                std::size_t d, LAMatrix &M, LAVector& Z)
    {
      typedef typename LocalKernel::FT FT;

      typedef CGAL::Simple_cartesian<FT> K;
      typedef typename K::Point_2        Point2d;
      typedef CGAL::Polygon_2<K>         Polygon2d;
      typedef CGAL::Min_sphere_of_points_d_traits_2<K,FT>    MinSphereTraits;
      typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
      typedef CGAL::Approximate_min_ellipsoid_d_traits_2<K, CGAL::MP_Float> MinEllipseTraits;
      typedef CGAL::Approximate_min_ellipsoid_d<MinEllipseTraits>           ApproxMinEllipse;
      using CGAL::fact;

      CGAL::Cartesian_converter<DataKernel, LocalKernel> D2L_converter;
      //origin of fitting coord system = first input data point
      Point_3 point0 = D2L_converter(*begin);
      //transform coordinates of sample points with a
      //translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
      Point_3 orig(0.,0.,0.);
      Vector_3 v_point0_orig(orig - point0);
      Aff_transformation transl(CGAL::TRANSLATION, v_point0_orig);
      this->translate_p0 = transl;
      Aff_transformation transf_points = this->change_world2fitting *
        this->translate_p0;

      //compute and store transformed points
      std::vector<Point_3> pts_in_fitting_basis;
      pts_in_fitting_basis.reserve(this->nb_input_pts);
      CGAL_For_all(begin,end){
        Point_3 cur_pt = transf_points(D2L_converter(*begin));
        pts_in_fitting_basis.push_back(cur_pt);
      }

      //Compute preconditionning
      FT precond = 0.;
      typename std::vector<Point_3>::iterator itb = pts_in_fitting_basis.begin(),
        ite = pts_in_fitting_basis.end();
      CGAL_For_all(itb,ite) precond += CGAL::abs(itb->x()) + CGAL::abs(itb->y());
      precond /= 2*this->nb_input_pts;
      this->preconditionning = precond;

      {
        // Calculate minimum bounding rectangle of 2D points.
        std::vector<Point2d> pts2d;
        pts2d.reserve(this->nb_input_pts);
        itb = pts_in_fitting_basis.begin();
        CGAL_For_all(itb,ite) {
          pts2d.push_back(Point2d(itb->x(), itb->y()));
        }
        std::vector<Point2d> ch2d;
        ch2d.reserve(this->nb_input_pts / 2);
        CGAL::ch_graham_andrew( pts2d.rbegin(), pts2d.rend(), std::back_inserter(ch2d), CGAL::Convex_hull_constructive_traits_2<K>());

        Polygon2d bounding_rect;
        CGAL::min_rectangle_2(ch2d.begin(), ch2d.end(), std::back_inserter(bounding_rect));
        this->bounding_area_ = BoundingArea::from_rectangle(bounding_rect);
        Min_sphere ms(ch2d.begin(), ch2d.end());
        this->bounding_area_.circle_radius = ms.radius();
        ApproxMinEllipse me(0.002, ch2d.begin(), ch2d.end());
        this->bounding_area_.ellipse_max_radius = 0.5 * *(me.axes_lengths_begin());
        this->bounding_area_.ellipse_min_radius = 0.5 * *(me.axes_lengths_begin() + 1);
        if (this->bounding_area_.ellipse_max_radius < this->bounding_area_.ellipse_min_radius)
        {
          std::swap(this->bounding_area_.ellipse_max_radius, this->bounding_area_.ellipse_min_radius);
        }
      }

      //fill matrices M and Z
      itb = pts_in_fitting_basis.begin();
      int line_count = 0;
      FT x, y;

      itb = pts_in_fitting_basis.begin();
      CGAL_For_all(itb,ite) {
        x = itb->x();
        y = itb->y();

        //  Z[line_count] = itb->z();
        Z.set(line_count,itb->z());
        for (std::size_t k=0; k <= d; k++) {
          for (std::size_t i=0; i<=k; i++) {
            M.set(line_count, k*(k+1)/2+i,
                  std::pow(x,static_cast<int>(k-i))
                  * std::pow(y,static_cast<int>(i))
                  /( fact(static_cast<unsigned int>(i)) *
                     fact(static_cast<unsigned int>(k-i))
                     *std::pow(this->preconditionning,static_cast<int>(k))));
          }
        }
        line_count++;
      }
    }

    void update_residual_stats(LAVector const & residuals)
    {
      boost::accumulators::accumulator_set<
          LocalFloatType,
          boost::accumulators::stats<
            boost::accumulators::tag::min,
            boost::accumulators::tag::max,
            boost::accumulators::tag::mean,
            boost::accumulators::tag::median,
            boost::accumulators::tag::variance
          >
        > acc;

      for (auto it = residuals.begin(); it != residuals.end(); ++it)
      {
        acc(*it);
      }
      this->residual_stats_.min = boost::accumulators::min(acc);
      this->residual_stats_.max = boost::accumulators::max(acc);
      this->residual_stats_.max_abs =
        std::max(std::fabs(this->residual_stats_.min), std::fabs(this->residual_stats_.max));
      this->residual_stats_.mean = boost::accumulators::mean(acc);
      this->residual_stats_.median = boost::accumulators::median(acc);
      this->residual_stats_.stdd = std::sqrt(boost::accumulators::variance(acc));
    }

    void solve_linear_system(LAMatrix &M, LAVector& Z)
    {
      {
        LAVector residuals(Z);
        this->condition_nb = CGALSvdTraits::solve(M, Z);

        LAVector Y(M.cols());
        for (std::size_t i = 0; i < M.cols(); ++i) Y.set(i, Z(i));
        residuals = (M * Y - residuals);
        this->update_residual_stats(residuals);
      }
      for (int k=0; k <= this->deg; k++) for (int i=0; i<=k; i++)
        // Z[k*(k+1)/2+i] /= std::pow(this->preconditionning,k);
        Z.set( k*(k+1)/2+i, Z(k*(k+1)/2+i) / std::pow(this->preconditionning,k) );
    }

    template <class InputIterator>
    MongeForm operator()(
        InputIterator begin,
        InputIterator end,
        std::size_t d,
        std::size_t dprime,
        Matrix_3x3 const & fitting_basis
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

      this->set_world_to_fitting_from_rotation(fitting_basis);

      this->fill_matrix(begin, end, d, M, Z);//with precond
      this->solve_linear_system(M, Z);  //correct with precond
      this->compute_Monge_basis(Z.vector(), monge_form);
      if ( dprime >= 3) this->compute_Monge_coefficients(Z.vector(), dprime, monge_form);

      MongeForm ret_monge_form(monge_form);
      ret_monge_form.poly_fit_condition_number_ = this->condition_number();
      ret_monge_form.fitting_residual_stats_ = this->residual_stats_;
      ret_monge_form.fitting_bounding_area_ = this->bounding_area_;

      return ret_monge_form;
    }

    ResidualStats residual_stats_;
    BoundingArea bounding_area_;
  };


  MongeJetFitter(const std::size_t degree_poly_fit=2, const std::size_t degree_monge=2) :
    num_rings_(0),
    in_points_(),
    in_normals_(),
    in_ring_(),
    degree_poly_fit_(degree_poly_fit),
    degree_monge_(degree_monge),
    pca_eigenvals_(
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN()
    ),
    ring_normal_gaussian_sigma_(1.6)
  {
  }

  std::size_t get_min_num_fit_points() const
  {
    return ((this->degree_poly_fit_ + 1) * (this->degree_poly_fit_ + 2)) / 2;
  }

  void set_monge_form_data(
      MongeForm & mf,
      const std::int64_t vertex_index,
      const std::int64_t num_rings,
      Matrix_3x3 const & fitting_basis
  )
  {
    mf.vertex_index_ = vertex_index;
    mf.degree_monge_ = this->degree_monge_;
    mf.degree_poly_fit_ = this->degree_poly_fit_;
    mf.num_rings_ = num_rings;
    mf.num_fitting_points_ = this->in_points_.size();
    mf.pca_eigenvalues_ = this->pca_eigenvals_;
    mf.fitting_basis_ = fitting_basis;
  }

  void gather_fitting_points(
      Vertex const * vertex,
      std::int64_t const num_rings,
      typename PolySurf::Vertex_PM_type & vertex_prop_map
  )
  {
    PolySurf::gather_fitting_points(
        vertex,
        num_rings,
        this->in_points_,
        this->in_normals_,
        this->in_ring_,
        vertex_prop_map
    );
  }

  void gather_fitting_points(
      PolySurf & poly_surface,
      std::int64_t const vertex_index,
      std::int64_t const num_rings
  )
  {
    auto vtx_it = poly_surface.vertices_begin() + vertex_index;
    return this->gather_fitting_points(&(*vtx_it), num_rings, poly_surface.vertex_prop_map_);
  }

  Matrix_3x3 calculate_fit_basis(FittingBasisType const & fit_basis_type)
  {
    typedef typename LocalKernel::FT FT;
    Matrix_3x3 R;
    switch (fit_basis_type)
    {
    case PCA:
    {
      auto R_evals_pair = MongeViaJetFitting::calculate_PCA_basis(this->in_points_.begin(), this->in_points_.end());
      R = R_evals_pair.first;
      this->pca_eigenvals_ =
        Vector_3(
          R_evals_pair.second[0],
          R_evals_pair.second[1],
          R_evals_pair.second[2]
        );
      break;
    }
    case VERTEX_NORMAL:
    {
      R = MongeViaJetFitting::calc_rotation_matrix_from_normal(*(this->in_normals_.begin()));
      break;
    }
    case RING_NORMAL_MEAN:
    {
      Vector_3 dir(0.0, 0.0, 0.0);
      for (auto it = this->in_normals_.begin(); it != this->in_normals_.end(); ++it)
      {
        dir += *it;
      }
      dir /= FT(this->in_normals_.size());
      R = MongeViaJetFitting::calc_rotation_matrix_from_normal(dir);
      break;
    }
    case RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN:
    {
      Vector_3 dir(0.0, 0.0, 0.0);
      FT weight_sum = 0.0;
      FT const sigma = FT(this->num_rings_) / 3.0;
      FT const sigma_sqrd = sigma * sigma;
      auto rit = this->in_ring_.begin();
      for (auto nit = this->in_normals_.begin(); nit != this->in_normals_.end(); ++nit, ++rit)
      {
        FT const r = FT(*rit);
        FT const weight = std::exp(-r*r / sigma_sqrd);
        weight_sum += weight;
        dir += weight * (*nit);
      }
      dir /= weight_sum;
      R = MongeViaJetFitting::calc_rotation_matrix_from_normal(dir);
      break;
    }
    case RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA:
    {
      Vector_3 dir(0.0, 0.0, 0.0);
      FT weight_sum = 0.0;
      FT const sigma = this->ring_normal_gaussian_sigma_;
      FT const sigma_sqrd = sigma * sigma;
      auto rit = this->in_ring_.begin();
      for (auto nit = this->in_normals_.begin(); nit != this->in_normals_.end(); ++nit, ++rit)
      {
        FT const r = FT(*rit);
        FT const weight = std::exp(-r*r / sigma_sqrd);
        weight_sum += weight;
        dir += weight * (*nit);
      }
      dir /= weight_sum;
      R = MongeViaJetFitting::calc_rotation_matrix_from_normal(dir);
      break;
    }

    default:
    {
      std::stringstream msg;
      msg << "Unhandled fit_basis_type=" << int(fit_basis_type);
      throw std::runtime_error(msg.str());
      break;
    }
    }

    return R;
  }

  MongeForm fit_at_vertex(
      PolySurf & poly_surface,
      const std::int64_t vertex_index,
      const std::int64_t num_rings,
      FittingBasisType const & fit_basis_type
  )
  {
    this->num_rings_ = num_rings;
    this->pca_eigenvals_ =
      Vector_3(
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN()
      );
    this->gather_fitting_points(poly_surface, vertex_index, num_rings);
    Matrix_3x3 const fit_basis = this->calculate_fit_basis(fit_basis_type);
    MongeForm monge_form = this->fit_at_vertex(fit_basis);
    this->set_monge_form_data(monge_form, vertex_index, num_rings, fit_basis);
    return monge_form;
  }

  MongeFormStlVecPtr fit_all(
      PolySurf const & poly_surface,
      const std::int64_t num_rings,
      FittingBasisType const & fit_basis_type
  )
  {
    MongeFormStlVecPtr monge_forms_ptr =
      std::make_unique<MongeFormStlVec>(
        std::distance(poly_surface.vertices_begin(), poly_surface.vertices_end())
      );

    this->num_rings_ = num_rings;
    this->pca_eigenvals_ =
      Vector_3(
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN(),
        std::numeric_limits<typename LocalKernel::FT>::quiet_NaN()
      );

    const std::int64_t num_vertices = std::int64_t(monge_forms_ptr->size());

#pragma omp parallel shared(poly_surface, num_rings, fit_basis_type, monge_forms_ptr, num_vertices)
    {
      std::int64_t vtx_idx = 0;
      typename PolySurf::Vertex2int_map_type vertex_map;
      typename PolySurf::Vertex_PM_type vertex_prop_map(vertex_map);
      auto vitb = poly_surface.vertices_begin();
      auto vite = poly_surface.vertices_end();
      vertex_map.clear();
      vertex_map.reserve(num_vertices);
      CGAL_For_all(vitb, vite) put(vertex_prop_map, &(*vitb), -1);
      MongeJetFitter fitter(*this);
#pragma omp for schedule(dynamic, 128)
      for (vtx_idx=0; vtx_idx < num_vertices; ++vtx_idx)
      {
        auto vtx_it = poly_surface.vertices_begin() + vtx_idx;
        fitter.gather_fitting_points(&(*vtx_it), num_rings, vertex_prop_map);
        Matrix_3x3 const fit_basis = fitter.calculate_fit_basis(fit_basis_type);
        MongeForm monge_form = fitter.fit_at_vertex(fit_basis);
        fitter.set_monge_form_data(monge_form, vtx_idx, num_rings, fit_basis);
        (*monge_forms_ptr)[vtx_idx] = monge_form;
      }
    }
    return monge_forms_ptr;
  }

  MongeForm fit_at_vertex(
      Matrix_3x3 const & fit_basis
  ) const
  {
    MongeViaJetFitting monge_fit;
    MongeForm monge_form;
    if (this->in_points_.size() >= this->get_min_num_fit_points())
    {
      monge_form =
        monge_fit(
          this->in_points_.begin(),
          this->in_points_.end(),
          this->degree_poly_fit_,
          this->degree_monge_,
          fit_basis
        );
      monge_form.comply_wrt_given_normal(this->in_normals_[0]);
    }

    return monge_form;
  }

  std::int64_t num_rings_; // number of vertex neighbourhood rings
  std::vector<Point_3> in_points_;  //container for data points
  std::vector<Vector_3> in_normals_;  //container for data point normals
  std::vector<std::int32_t> in_ring_;  //container for point ring-number
  std::size_t degree_poly_fit_; // Degree of the fitting polynomial
  std::size_t degree_monge_; // Degree of the Monge polynomial
  Vector_3 pca_eigenvals_; // Degree of the Monge polynomial
  typename LocalKernel::FT ring_normal_gaussian_sigma_;
};

}

#endif
