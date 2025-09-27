// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check consistency of function implementations

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/base/function_bessel.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <boost/core/demangle.hpp>

#include <string>
#include <typeinfo>
#include <vector>

#include "../tests.h"

#include "functions.h"

template <typename FunctionType, typename... Args>
void
check(Args... args)
{
  FunctionType f(args...);
  std::string  function_name = boost::core::demangle(typeid(f).name());
  // remove the preceding dealii:: to match older output
  function_name.erase(0, sizeof("dealii::") - 1 /*remember the NUL*/);
  deallog << function_name << std::endl;
  check_function_value_consistency(f, 5);
  check_function_gradient_consistency(f, 5);
  check_gradient(f, 5);
}



int
main()
{
  initlog();

  check<Functions::SquareFunction<1>>();
  check<Functions::SquareFunction<2>>();
  check<Functions::SquareFunction<3>>();
  check<Functions::CosineFunction<1>>();
  check<Functions::CosineFunction<2>>();
  check<Functions::CosineFunction<3>>();
  check<Functions::CosineGradFunction<1>>();
  check<Functions::CosineGradFunction<2>>();
  check<Functions::CosineGradFunction<3>>();
  check<Functions::ExpFunction<1>>();
  check<Functions::ExpFunction<2>>();
  check<Functions::ExpFunction<3>>();
  // Bessel1 is only implemented in 2D right now
  // check<Functions::Bessel1<1>>(3, 2.0, Point<1>(1.0));
  check<Functions::Bessel1<2>>(3, 2.0, Point<2>(1.0, 2.0));
  // check<Functions::Bessel1<3>>(3, 2.0, Point<2>(1.0, 2.0));
  // check(Functions::Bessel1<3>);
}
