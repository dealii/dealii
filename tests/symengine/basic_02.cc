// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check that SymEngine can perform symbolic differentiation correctly

// References:
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_basic.cpp
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_number.cpp
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_subs.cpp
// https://github.com/symengine/symengine/blob/master/symengine/symengine_casts.h
// https://github.com/symengine/symengine/blob/master/symengine/symengine_rcp.h
// https://github.com/symengine/symengine/blob/master/symengine/derivative.h
// https://github.com/symengine/symengine/blob/master/symengine/subs.h

#include "../tests.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"
#include <symengine/complex_double.h>
#include <symengine/derivative.h>
#include <symengine/eval_double.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>
#include <symengine/symengine_casts.h>
#include <symengine/visitor.h>
#pragma clang diagnostic pop

#include <fstream>
#include <iomanip>

namespace SE = SymEngine;

SE::RCP<const SE::Number>
make_symengine_rcp(const int &val)
{
  return SE::integer(val);
}

SE::RCP<const SE::Number>
make_symengine_rcp(const float &val)
{
  return SE::real_double(val);
}

SE::RCP<const SE::Number>
make_symengine_rcp(const double &val)
{
  return SE::real_double(val);
}

template <typename NumberType>
SE::RCP<const SE::Number>
make_symengine_rcp(const std::complex<NumberType> &val)
{
  // Build complex from two SymEngine numbers
  return SE::Complex::from_two_nums(*make_symengine_rcp(val.real()),
                                    *make_symengine_rcp(val.imag()));
}

SE::RCP<const SE::Symbol>
make_symengine_rcp(const std::string &name)
{
  return SE::symbol(name);
}

template <typename NumberType>
void
test_number()
{
  SE::RCP<const SE::Number> a = make_symengine_rcp(NumberType(3.1));
  SE::RCP<const SE::Basic>  b = make_symengine_rcp(NumberType(7.5));
  deallog << "a: " << *a << std::endl;
  deallog << "b: " << *b << std::endl;

  SE::RCP<const SE::Symbol> x(make_symengine_rcp("x"));
  SE::RCP<const SE::Basic>  y(make_symengine_rcp("y"));
  deallog << "x: " << *x << std::endl;
  deallog << "y: " << *y << std::endl;

  SE::RCP<const SE::Basic> c;
  // deallog << "c (no constructor): " << *c << std::endl; // Fails

  // Construction of symbolic function
  c = SE::mul(y, SE::mul(SE::sub(y, b), SE::add(a, x)));
  deallog << "c = y*(y-b)*(a+x): " << *c << std::endl;

  // Perform symbolic differentiation
  SE::RCP<const SE::Basic> dc_dx = c->diff(x);
  // SE::RCP<const SE::Basic> dc_dy = diff(c,SE::implicit_cast<const
  // SE::RCP<const SE::Symbol> &>(y));
  SE::RCP<const SE::Basic> dc_dy =
    diff(c, SE::rcp_static_cast<const SE::Symbol>(y));
  SE::RCP<const SE::Basic> dc_dy_2 = sdiff(c, y);

  deallog << "dc_dx = a*y*(y-b): " << *dc_dx << std::endl;
  deallog << "dc_dy = 2*y*(a+x): " << *dc_dy << std::endl;
  deallog << "dc_dy = 2*y*(a+x): " << *dc_dy_2 << std::endl;

  // Substitute values
  SE::map_basic_basic sub_vals;
  sub_vals[x] = make_symengine_rcp(NumberType(1));
  sub_vals[y] = make_symengine_rcp(NumberType(2.2));
  deallog << "dc_dx(x=1,y=2.2): " << *(dc_dx->subs(sub_vals)) << std::endl;

  sub_vals[x] = make_symengine_rcp(NumberType(0.75));
  sub_vals[y] = make_symengine_rcp(NumberType(-2.08));
  deallog << "dc_dy(x=0.75,y=-1.08): " << *(dc_dy->subs(sub_vals)) << std::endl;

  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog.push("Integer");
  test_number<int>();
  deallog.pop();

  deallog.push("Float");
  test_number<float>();
  deallog.pop();

  deallog.push("Double");
  test_number<double>();
  deallog.pop();

  // Not available yet
  // SymEngine::SymEngineException: Invalid Format: Expected Integer or Rational
  // deallog << "Complex double" << std::endl;
  // test_number<std::complex<double>>();
  // deallog.pop();

  deallog << "OK" << std::endl;
}
