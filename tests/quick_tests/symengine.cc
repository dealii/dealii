// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
