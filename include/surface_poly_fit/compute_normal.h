
#ifndef SURFACE_POLY_FIT_COMPUTE_NORMAL_H
#define SURFACE_POLY_FIT_COMPUTE_NORMAL_H

#include <CGAL/version.h>
#if CGAL_VERSION_NR >= 1050600000
#include "surface_poly_fit/compute_normal_CGAL_5pt6.h"
#elif CGAL_VERSION_NR >= 1050500000
#include "surface_poly_fit/compute_normal_CGAL_5pt5.h"
#else
#include "surface_poly_fit/compute_normal_CGAL_5pt4.h"
#endif
#endif // SURFACE_POLY_FIT_COMPUTE_NORMAL_H
