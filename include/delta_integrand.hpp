#ifndef __DELTA_INTEGRAND_HPP__
#define __DELTA_INTEGRAND_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "projection.hpp"
//
#include "shared_definitions.hpp"
#include "R_mn.hpp"
#include "integrand_projection_1d.hpp"
#include "integrand_projection_2d.hpp"
#include "integrand_projection_3d.hpp"

// Definitions

// Functions
complex_t delta_2d_2d(const basis_2d_t b_m, const basis_2d_t b_n, const complex_t k, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag);
    
#endif
