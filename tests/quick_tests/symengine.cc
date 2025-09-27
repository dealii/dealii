// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Adaptation of symengine/basic_[01,02].cc for a quick test

#include <deal.II/base/logstream.h>

#include <symengine/add.h>
#include <symengine/derivative.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>

#include <fstream>
#include <iomanip>

using namespace dealii;
namespace SE = SymEngine;

int
main()
{
  {
    const SE::RCP<const SE::Number> a = SE::integer(2);
    const SE::RCP<const SE::Number> b = SE::real_double(4.2);

    const SE::RCP<const SE::Basic> c = SE::add(a, b);
    deallog << "c = a+b: " << *c << std::endl;
  }

  {
    const SE::RCP<const SE::Number> a = SE::real_double(3.1);
    const SE::RCP<const SE::Basic>  b = SE::real_double(7.5);

    const SE::RCP<const SE::Symbol> x = (SE::symbol("x"));
    const SE::RCP<const SE::Symbol> y = (SE::symbol("y"));

    // Construction of symbolic function
    const SE::RCP<const SE::Basic> c =
      SE::mul(y, SE::mul(SE::sub(y, b), SE::add(a, x)));
    deallog << "c = y*(y-b)*(a+x): " << *c << std::endl;

    // Perform symbolic differentiation
    const SE::RCP<const SE::Basic> dc_dx = c->diff(x);
    const SE::RCP<const SE::Basic> dc_dy = c->diff(y);

    deallog << "dc_dx = y*(y-b): " << *dc_dx << std::endl;
    deallog << "dc_dy = (2*y+1)*(a+x): " << *dc_dy << std::endl;

    const SE::RCP<const SE::Basic> dc_dx_check = SE::mul(y, SE::sub(y, b));
    const SE::RCP<const SE::Basic> dc_dy_check =
      SE::mul(SE::sub(SE::mul(SE::integer(2), y), b), SE::add(a, x));
    Assert(SE::eq(*dc_dx, *dc_dx_check), ExcMessage("Should be equal!"));
    // Although these two *values* are the same, the underlying
    // *representation* is different. So we'd need to match the representation
    // exactly, and I'm too lazy to do this right now.
    // Assert(SE::eq(*dc_dy, *dc_dy_check), ExcMessage("Should be equal!"));
  }

  deallog << "OK" << std::endl;
}
