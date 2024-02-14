
#include "surface_poly_fit/monge_jet_fitting_py.h"
#include "surface_poly_fit/monge_jet_fitting.h"
#include "surface_poly_fit/polyhedral_surface_py.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace spf
{

namespace py = pybind11;

using Array2 = std::array<double, 2>;
using Array3 = std::array<double, 3>;
using Array4 = std::array<double, 4>;
using Array5 = std::array<double, 5>;
using Array3x3 = std::array<std::array<double, 3>, 3>;

typedef MongeJetFitter<PolyhedralSurface> MongeFitter;
typedef MongeFitter::MongeForm MongeForm;
typedef MongeFitter::MongeFormStlVec MongeFormStlVec;
typedef MongeFitter::MongeFormStlVecPtr MongeFormStlVecPtr;


struct BoundingAreaNumpy: public MongeFitter::BoundingArea
{
public:
  BoundingAreaNumpy() :
    MongeFitter::BoundingArea()
  {
  }

  BoundingAreaNumpy(MongeFitter::BoundingArea const & ba) :
    MongeFitter::BoundingArea(ba)
  {
  }

};

struct ResidualStatsNumpy: public MongeFitter::ResidualStats
{
public:
  ResidualStatsNumpy() :
    MongeFitter::ResidualStats()
  {
  }

  ResidualStatsNumpy(MongeFitter::ResidualStats const & rs) :
    MongeFitter::ResidualStats(rs)
  {
  }

};

#pragma pack(1)
struct MongeFormNumpy
{
  MongeFormNumpy()
  {
  }

  explicit
  MongeFormNumpy(MongeForm const & mf)
  {
    this->vertex_index = mf.vertex_index_;
    this->degree_monge = mf.degree_monge_;
    this->degree_poly_fit = mf.degree_poly_fit_;
    this->num_rings = mf.num_rings_;
    this->num_fitting_points = mf.num_fitting_points_;
    this->poly_fit_condition_number = mf.poly_fit_condition_number_;

    for (std::size_t i = 0; i < 3; ++i)
    {
      this->pca_eigenvalues[i] = mf.pca_eigenvalues_[i];
    }

    for (std::size_t i = 0; i < 3; ++i)
    {
      for (std::size_t j = 0; j < 3; ++j)
      {
        this->poly_fit_basis[i][j] = mf.fitting_basis_(i, j);
      }
    }

    this->poly_fit_residual_stats = ResidualStatsNumpy(mf.fitting_residual_stats_);
    this->poly_fit_bounding_area = BoundingAreaNumpy(mf.fitting_bounding_area_);

    for (std::size_t i = 0; i < 3; ++i)
    {
      this->origin[i] = mf.origin()[i];
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
      this->direction[i][0] = mf.maximal_principal_direction()[i];
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
      this->direction[i][1] = mf.minimal_principal_direction()[i];
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
      this->direction[i][2] = mf.normal_direction()[i];
    }
    if (this->degree_monge > 1)
    {
      this->k[0] = mf.principal_curvatures(0);
      this->k[1] = mf.principal_curvatures(1);
    }
    else
    {
      this->k[0] = 0.0;
      this->k[1] = 0.0;
    }
    if (this->degree_monge > 2)
    {
      this->b[0] = mf.third_order_coefficients(0);
      this->b[1] = mf.third_order_coefficients(1);
      this->b[2] = mf.third_order_coefficients(2);
      this->b[3] = mf.third_order_coefficients(3);
    }
    else
    {
      this->b[0] = 0.0;
      this->b[1] = 0.0;
      this->b[2] = 0.0;
      this->b[3] = 0.0;
    }
    if (this->degree_monge > 3)
    {
      this->c[0] = mf.fourth_order_coefficients(0);
      this->c[1] = mf.fourth_order_coefficients(1);
      this->c[2] = mf.fourth_order_coefficients(2);
      this->c[3] = mf.fourth_order_coefficients(3);
      this->c[4] = mf.fourth_order_coefficients(4);
    }
    else
    {
      this->c[0] = 0.0;
      this->c[1] = 0.0;
      this->c[2] = 0.0;
      this->c[3] = 0.0;
      this->c[4] = 0.0;
    }
  }

  std::int64_t vertex_index;
  std::uint8_t degree_monge;
  std::uint8_t degree_poly_fit;
  std::int64_t num_rings;
  std::int64_t num_fitting_points;
  double poly_fit_condition_number;
  Array3 pca_eigenvalues;
  Array3x3 poly_fit_basis;
  ResidualStatsNumpy poly_fit_residual_stats;
  BoundingAreaNumpy poly_fit_bounding_area;
  Array3 origin;
  Array3x3 direction;
  Array2 k;
  Array4 b;
  Array5 c;
};
#pragma pack(0)

/// Python wrapper for MongeJetFitter.
class MongeJetFitterPy
{
public:

  MongeJetFitterPy(py::object poly_surface, std::size_t degree_poly_fit, std::size_t degree_monge)
    : monge_fitter_(degree_poly_fit, degree_monge),
      poly_surface_obj_(poly_surface)
  {
  }

  std::int64_t get_degree_poly_fit() const
  {
    return std::int64_t(this->monge_fitter_.degree_poly_fit_);
  }

  std::int64_t get_degree_monge() const
  {
    return std::int64_t(this->monge_fitter_.degree_monge_);
  }

  std::int64_t get_min_num_fit_points() const
  {
    return std::int64_t(this->monge_fitter_.get_min_num_fit_points());
  }

  double get_ring_normal_gaussian_sigma() const
  {
    return this->monge_fitter_.ring_normal_gaussian_sigma_;
  }

  void set_ring_normal_gaussian_sigma(const double sigma)
  {
    this->monge_fitter_.ring_normal_gaussian_sigma_ = sigma;
  }

  py::object monge_form_to_array(MongeForm const & monge_form) const
  {
    std::size_t const shape[1] = {1};
    auto ary = py::array_t<MongeFormNumpy>(shape);
    auto req = ary.request();
    auto ptr = static_cast<MongeFormNumpy *>(req.ptr);

    MongeFormNumpy const ary_elem(monge_form);
    ptr[0] = ary_elem;

    return py::object(ary);
  }

  py::object monge_form_to_array(MongeFormStlVec const & monge_forms) const
  {
    std::size_t const shape[1] = {monge_forms.size()};
    auto ary = py::array_t<MongeFormNumpy>(shape);
    auto req = ary.request();
    auto ptr = static_cast<MongeFormNumpy *>(req.ptr);

    for (std::size_t i = 0; i < shape[0]; ++i)
    {
      MongeFormNumpy const ary_elem(monge_forms[i]);
      ptr[i] = ary_elem;
    }
    return py::object(ary);
  }

  py::object fit_at_vertex(std::int64_t vertex_index, std::int64_t num_rings, MongeFitter::FittingBasisType fit_basis_type)
  {
    MongeForm mf_result =
        this->monge_fitter_.fit_at_vertex(
            *(py::cast<PolyhedralSurfacePy &>(this->poly_surface_obj_).surface_),
            vertex_index,
            num_rings,
            fit_basis_type
        );

    return this->monge_form_to_array(mf_result);
  }

  py::object fit_all(std::int64_t num_rings, MongeFitter::FittingBasisType fit_basis_type)
  {
    MongeFormStlVecPtr mf_result =
        this->monge_fitter_.fit_all(
            *(py::cast<PolyhedralSurfacePy &>(this->poly_surface_obj_).surface_),
            num_rings,
            fit_basis_type
        );

    return this->monge_form_to_array(*mf_result);
  }

  py::object get_poly_surface()
  {
    return this->poly_surface_obj_;
  }

protected:
  MongeFitter monge_fitter_;
  py::object poly_surface_obj_;
};


/// Export MongeJetFitterPy class to specified python module.
void export_monge_jet_fitter(pybind11::module_ m)
{
  PYBIND11_NUMPY_DTYPE(
      BoundingAreaNumpy,
      rectangle_min_side_length,
      rectangle_max_side_length,
      circle_radius,
      ellipse_min_radius,
      ellipse_max_radius
  );

  PYBIND11_NUMPY_DTYPE(
      ResidualStatsNumpy,
      min,
      max,
      mean,
      median,
      min_abs,
      max_abs,
      mean_abs,
      median_abs,
      stdd
  );

  PYBIND11_NUMPY_DTYPE(
      MongeFormNumpy,
      vertex_index,
      degree_monge,
      degree_poly_fit,
      num_rings,
      num_fitting_points,
      poly_fit_condition_number,
      pca_eigenvalues,
      poly_fit_basis,
      poly_fit_residual_stats,
      poly_fit_bounding_area,
      origin,
      direction,
      k,
      b,
      c
  );

  py::class_<MongeJetFitterPy> mjf_class(
      m,
      "MongeJetFitter"
  );

  py::enum_<MongeFitter::FittingBasisType>(
    mjf_class,
    "FittingBasisType",
    "How the polynomial fitting coordinate system (fitting-basis) is calculated"
    " for a polyhedral-surface-patch."
  )
    .value(
        "PCA", MongeFitter::FittingBasisType::PCA,
        "Uses Principal Component Analysis of the patch coordinates"
        " with the smallest component assigned to the :math:`z` direction"
        " and the largest component the :math:`x` direction."
    )
    .value(
        "VERTEX_NORMAL", MongeFitter::FittingBasisType::VERTEX_NORMAL,
        "Use the vertex-normal as the polynomial fitting coordinate system :math:`z` direction."
    )
    .value(
        "RING_NORMAL_MEAN", MongeFitter::FittingBasisType::RING_NORMAL_MEAN,
        "Use mean of all surface-patch vertex-normals as the polynomial fitting coordinate system :math:`z` direction."
    )
    .value(
        "RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN", MongeFitter::FittingBasisType::RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN,
        "Use a Gaussian-weighted-mean of all surface-patch vertex-normals as the polynomial fitting coordinate system :math:`z` direction."
        " The sigma used for calculating the Gaussian weights is :samp:`num_rings / 3.0`."

    )
    .value(
        "RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA", MongeFitter::FittingBasisType::RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA,
        "Use a Gaussian-weighted-mean of all surface-patch vertex-normals as the polynomial fitting coordinate system :math:`z` direction."
        " The sigma used for calculating the Gaussian weights is :attr:`MongeJetFitter.ring_normal_gaussian_sigma`."
    )
    .export_values()
    ;

  mjf_class
    .def(
        py::init<py::object, std::size_t, std::size_t>(),
        py::arg("poly_surface"),
        py::arg("degree_poly_fit")=2,
        py::arg("degree_monge")=2,
        "Construct.\n\n"
        ":type poly_surface: :obj:`PolyhedralSurface`\n"
        ":param poly_surface: Fit polynomial surfaces to sub-patches of this surface mesh.\n"
        ":type degree_poly_fit: :obj:`int`\n"
        ":param degree_poly_fit: The degree (highest exponent) of the fitting polynomial."
        " :samp:`{degree_poly_fit} >= 1`.\n"
        ":type degree_monge: :obj:`int`\n"
        ":param degree_monge: The degree (highest exponent) of the Monge polynomial."
        " :samp:`1 <= {degree_monge} <= {degree_poly_fit} <= 4`.\n"
    )
    .def_property_readonly(
        "poly_surface",
        &MongeJetFitterPy::get_poly_surface,
        "Polynomial patch vertex fitting performed on this :obj:`PolyhedralSurface` .\n"
    )
    .def_property_readonly(
        "degree_poly_fit",
        &MongeJetFitterPy::get_degree_poly_fit,
        "An :obj:`int` indicating the *degree* of the fitting polynomial.\n"
    )
    .def_property_readonly(
        "degree_monge",
        &MongeJetFitterPy::get_degree_monge,
        "An :obj:`int` indicating the *degree* of the Monge polynomial.\n"
    )
    .def_property_readonly(
        "min_num_fit_points",
        &MongeJetFitterPy::get_min_num_fit_points,
        "An :obj:`int` indicating the minimum number of points (vertex coordinates)"
        " required to fit a polynomial of degree :attr:`degree_poly_fit`.\n"
    )
    .def_property(
      "ring_normal_gaussian_sigma",
      &MongeJetFitterPy::get_ring_normal_gaussian_sigma,
      &MongeJetFitterPy::set_ring_normal_gaussian_sigma,
      "The *sigma* :obj:`float` value used to calculate Gaussian weights"
      " when determining the polynomial fitting basis :math:`z` direction.\n"
    )
    .def(
        "fit_at_vertex",
        &MongeJetFitterPy::fit_at_vertex,
        py::arg("vertex_index"),
        py::arg("num_rings")=1,
        py::arg("fit_basis_type")=MongeFitter::FittingBasisType::VERTEX_NORMAL,
        "Fits a polynomial surface to the neigbourhood of vertices of vertex"
        " with index :samp:`{vertex_index}`.\n\n"
        ":type vertex_index: :obj:`int`\n"
        ":param vertex_index: Fit polynomial to neighbourhood vertices of"
        " this vertex.\n"
        ":type num_rings: :obj:`int`\n"
        ":param num_rings: Size the neighbourhood. Vertices within :samp:`{num_rings}`"
        " edge-hops of the :samp:`vertex_index` vertex are used in the polynomial-surface"
        " fitting.\n"
        ":type fit_basis_type: :obj:`MongeJetFitter.FittingBasisType`\n"
        ":param fit_basis_type: This specifies the method used to calculate"
        " the :math:`z` direction of the *fitting coordinate system*."
        " See :ref:`fitting-basis-type-description`.\n"
        ":rtype: :obj:`numpy.ndarray`\n"
        ":return: A :samp:`(1,)` shaped `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_"
        " containing fitting results."
        " See :ref:`fitting-result-array-description`.\n"
     )
     .def(
         "fit_all",
         &MongeJetFitterPy::fit_all,
         py::arg("num_rings")=1,
         py::arg("fit_basis_type")=MongeFitter::FittingBasisType::VERTEX_NORMAL,
         "Fits a polynomial surface to the neigbourhoods of all vertices.\n\n"
         ":type num_rings: :obj:`int`\n"
         ":param num_rings: Size the neighbourhood. Vertices within :samp:`{num_rings}`"
         " edge-hops of the :samp:`vertex_index` vertex are used in the polynomial-surface"
         " fitting.\n"
         ":type fit_basis_type: :obj:`MongeJetFitter.FittingBasisType`\n"
         ":param fit_basis_type: This specifies the method used to calculate"
         " the :math:`z` direction of the *fitting coordinate system*."
         " See :ref:`fitting-basis-type-description`.\n"
         ":rtype: :obj:`numpy.ndarray`\n"
         ":return: A :samp:`(N,)` shaped `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_"
         " containing fitting results."
         " See :ref:`fitting-result-array-description`.\n"
      )
  ;

}

}
