#ifndef SURFACE_POLY_FIT_MONGE_JET_FITTING_PY_H
#define SURFACE_POLY_FIT_MONGE_JET_FITTING_PY_H

#include <pybind11/pybind11.h>


namespace spf
{

void export_monge_jet_fitter(pybind11::module_ m);

}

#endif
