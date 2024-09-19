#ifndef __TESTBENCH_2D_HPP__
#define __TESTBENCH_2D_HPP__

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
#include "engine.hpp"

// Definitions

// Functions
void test_engine_2d_2d();
void test_engine_2d_sphere_RCS();
void test_engine_2d_sphere_near_field();
void test_engine_2d_debug();
void test_engine_2d_sheet_near_field();
void test_engine_2d_box_RCS();
void test_engine_2d_box_near_field();
void test_engine_2d_sphere_near_field_2d();

#endif