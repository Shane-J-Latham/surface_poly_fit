
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

protected:
  MongeFitter monge_fitter_;
  py::object poly_surface_obj_;
};


/// Export MongeJetFitterPy class to specified python module.
void export_monge_jet_fitter(pybind11::module_ m)
{
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
      origin,
      direction,
      k,
      b,
      c
  );

  py::class_<MongeJetFitterPy> mjf_class(m, "MongeJetFitter");

  py::enum_<MongeFitter::FittingBasisType>(mjf_class, "FittingBasisType")
    .value("PCA", MongeFitter::FittingBasisType::PCA)
    .value("VERTEX_NORMAL", MongeFitter::FittingBasisType::VERTEX_NORMAL)
    .value("RING_NORMAL_MEAN", MongeFitter::FittingBasisType::RING_NORMAL_MEAN)
    .value("RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN", MongeFitter::FittingBasisType::RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN)
    .value("RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA", MongeFitter::FittingBasisType::RING_NORMAL_GAUSSIAN_WEIGHTED_MEAN_SIGMA)
    .export_values();

  mjf_class
    .def(
        py::init<py::object, std::size_t, std::size_t>(),
        py::arg("poly_surface"),
        py::arg("degree_poly_fit")=2,
        py::arg("degree_monge")=2
    )
    .def_property_readonly("degree_poly_fit", &MongeJetFitterPy::get_degree_poly_fit)
    .def_property_readonly("degree_monge", &MongeJetFitterPy::get_degree_monge)
    .def_property_readonly("min_num_fit_points", &MongeJetFitterPy::get_min_num_fit_points)
    .def_property(
      "ring_normal_gaussian_sigma",
      &MongeJetFitterPy::get_ring_normal_gaussian_sigma,
      &MongeJetFitterPy::set_ring_normal_gaussian_sigma
    )
    .def(
        "fit_at_vertex",
        &MongeJetFitterPy::fit_at_vertex,
        py::arg("vertex_index"),
        py::arg("num_rings")=1,
        py::arg("fit_basis_type")=MongeFitter::FittingBasisType::VERTEX_NORMAL
     )
     .def(
         "fit_all",
         &MongeJetFitterPy::fit_all,
         py::arg("num_rings")=1,
         py::arg("fit_basis_type")=MongeFitter::FittingBasisType::VERTEX_NORMAL
      )
  ;

}

}
