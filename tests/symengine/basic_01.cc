// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that all fundamental operations with primitive SymEngine number
// types work as expected

// References:
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_basic.cpp
// https://github.com/symengine/symengine/blob/master/symengine/tests/basic/test_number.cpp

#include "../tests.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"
#include <symengine/complex_double.h>
#include <symengine/eval_double.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
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

template <typename NumberType>
void
test_number()
{
  SE::RCP<const SE::Number> a = make_symengine_rcp(NumberType(4.2));
  deallog << "a: " << *a << std::endl;

  SE::RCP<const SE::Basic> b(make_symengine_rcp(NumberType(2.1)));
  deallog << "b: " << *b << std::endl;

  SE::RCP<const SE::Basic> c;
  // deallog << "c (no constructor): " << *c << std::endl; // Fails

  c = SE::add(a, b);
  deallog << "c = a+b: " << *c << std::endl;

  c = SE::sub(a, b);
  deallog << "c = a-b: " << *c << std::endl;

  c = SE::mul(a, b);
  deallog << "c = a*b: " << *c << std::endl;

  c = SE::div(a, b);
  deallog << "c = a/b: " << *c << std::endl;

  c = a;
  deallog << "c = a: " << *c << std::endl;

  c = SE::add(c, a);
  deallog << "c += a: " << *c << std::endl;

  c = SE::sub(c, a);
  deallog << "c -= a: " << *c << std::endl;

  c = SE::mul(c, a);
  deallog << "c *= a: " << *c << std::endl;

  c = SE::div(c, a);
  deallog << "c /= a: " << *c << std::endl;
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
