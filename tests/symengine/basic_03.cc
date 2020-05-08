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


// Check that SymEngine can perform symbolic differentiation using function
// symbols correctly

// References:
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_functions.cpp
// https://github.com/symengine/symengine/blob/master/symengine/function.h
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

SE::RCP<const SE::Basic>
make_symengine_rcp(const std::string &                          name,
                   const std::vector<SE::RCP<const SE::Basic>> &args)
{
  return SE::function_symbol(name, args);
}

SE::RCP<const SE::Basic>
make_symengine_rcp(const std::string &name, const SE::RCP<const SE::Basic> &arg)
{
  return SE::function_symbol(name, arg);
}

template <typename NumberType>
void
test_number()
{
  // Modified snippet from SymEngine test
  // {
  //   SE::RCP<const SE::Symbol> x = SE::symbol("x");
  //   SE::RCP<const SE::Symbol> _xi_1 = SE::symbol("_xi_1");
  //   SE::RCP<const SE::Basic> f = SE::function_symbol("f", x);
  //   SE::RCP<const SE::Basic> r1, r2;
  //
  //   f = function_symbol("f", SE::pow(x, SE::integer(2)));
  //   deallog << "f: " << *f << std::endl;
  //   r1 = f->diff(x);
  //   deallog << "r1: " << *r1 << std::endl;
  //   r2 = SE::Derivative::create(SE::function_symbol("f", _xi_1), {_xi_1});
  //   deallog << "r2: " << *r2 << std::endl;
  //   r2 = SE::Subs::create(r2, {{_xi_1, SE::pow(x, SE::integer(2))}});
  //   deallog << "r2: " << *r2 << std::endl;
  // }

  // Normal symbols
  SE::RCP<const SE::Symbol> x(make_symengine_rcp("x"));
  SE::RCP<const SE::Basic>  y(make_symengine_rcp("y"));
  deallog << "x: " << *x << std::endl;
  deallog << "y: " << *y << std::endl;

  // Function symbols
  SE::RCP<const SE::Basic> f(make_symengine_rcp("f", x));
  SE::RCP<const SE::Basic> g(make_symengine_rcp("g", {x, y}));
  SE::RCP<const SE::Basic> h(
    make_symengine_rcp("h", {SE::add(x, y), SE::mul(x, y)}));
  deallog << "f: " << *f << std::endl;
  deallog << "g: " << *g << std::endl;
  deallog << "h: " << *h << std::endl;

  // Perform symbolic differentiation
  SE::RCP<const SE::Basic> df_dx = f->diff(x);
  SE::RCP<const SE::Basic> dg_dx = g->diff(x);
  SE::RCP<const SE::Basic> dh_dx = h->diff(x);
  deallog << "df_dx: " << *df_dx << std::endl;
  deallog << "dg_dx: " << *dg_dx << std::endl;
  deallog << "dh_dx: " << *dh_dx << std::endl;
  SE::RCP<const SE::Basic> df_dy =
    f->diff(SE::rcp_static_cast<const SE::Symbol>(y));
  SE::RCP<const SE::Basic> dg_dy =
    g->diff(SE::rcp_static_cast<const SE::Symbol>(y));
  SE::RCP<const SE::Basic> dh_dy =
    h->diff(SE::rcp_static_cast<const SE::Symbol>(y));
  deallog << "df_dy: " << *df_dy << std::endl;
  deallog << "dg_dy: " << *dg_dy << std::endl;
  deallog << "dh_dy: " << *dh_dy << std::endl;

  // Substitute values
  SE::map_basic_basic sub_vals;
  sub_vals[x] = make_symengine_rcp(NumberType(1));
  sub_vals[y] = make_symengine_rcp(NumberType(2.2));
  sub_vals[f] = SE::pow(x, SE::integer(2));
  sub_vals[g] = SE::add(x, y);
  sub_vals[h] = SE::sub(x, y);

  deallog << "f(x=1,y=2.2): " << *(f->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  deallog << "g(x=1,y=2.2): " << *(g->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  deallog << "h(x=1,y=2.2): " << *(h->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  deallog << "Eval: f(x=1,y=2.2): "
          << eval_double(*(f->subs(sub_vals)->subs(sub_vals))) << std::endl;
  deallog << "Eval: g(x=1,y=2.2): "
          << eval_double(*(g->subs(sub_vals)->subs(sub_vals))) << std::endl;
  deallog << "Eval: h(x=1,y=2.2): "
          << eval_double(*(h->subs(sub_vals)->subs(sub_vals))) << std::endl;

  // Not yet implemented
  deallog << "df_dx(x=1,y=2.2): " << *(df_dx->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  deallog << "dg_dx(x=1,y=2.2): " << *(dg_dx->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  deallog << "dh_dx(x=1,y=2.2): " << *(dh_dx->subs(sub_vals))
          << std::endl; // ->subs(sub_vals)
  // deallog << "Eval: df_dx(x=1,y=2.2): " <<
  // eval_double(*(df_dx->subs(sub_vals)->subs(sub_vals))) << std::endl; deallog
  // << "Eval: dg_dx(x=1,y=2.2): " <<
  // eval_double(*(dg_dx->subs(sub_vals)->subs(sub_vals))) << std::endl; deallog
  // << "Eval: dh_dx(x=1,y=2.2): " <<
  // eval_double(*(dh_dx->subs(sub_vals)->subs(sub_vals))) << std::endl;
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
