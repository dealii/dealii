// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_math.h>

#  include <symengine/add.h>
#  include <symengine/functions.h>
#  include <symengine/mul.h>
#  include <symengine/pow.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;

    /* --------------------------- Math functions ------------------------- */

    Expression
    pow(const Expression &base, const Expression &exponent)
    {
      return SE::pow(base.get_RCP(), exponent.get_RCP());
    }


    Expression
    sqrt(const Expression &x)
    {
      return SE::sqrt(x.get_RCP());
    }


    Expression
    cbrt(const Expression &x)
    {
      return SE::cbrt(x.get_RCP());
    }


    Expression
    exp(const Expression &exponent)
    {
      return SE::exp(exponent.get_RCP());
    }


    Expression
    log(const Expression &x)
    {
      return SE::log(x.get_RCP());
    }


    Expression
    log(const Expression &x, const Expression &base)
    {
      return SE::log(x.get_RCP(), base.get_RCP());
    }


    Expression
    log10(const Expression &x)
    {
      return log(x.get_RCP(), Expression(10.0));
    }


    Expression
    sin(const Expression &x)
    {
      return SE::sin(x.get_RCP());
    }


    Expression
    cos(const Expression &x)
    {
      return SE::cos(x.get_RCP());
    }


    Expression
    tan(const Expression &x)
    {
      return SE::tan(x.get_RCP());
    }


    Expression
    csc(const Expression &x)
    {
      return SE::csc(x.get_RCP());
    }


    Expression
    sec(const Expression &x)
    {
      return SE::sec(x.get_RCP());
    }


    Expression
    cot(const Expression &x)
    {
      return SE::cot(x.get_RCP());
    }


    Expression
    asin(const Expression &x)
    {
      return SE::asin(x.get_RCP());
    }


    Expression
    acos(const Expression &x)
    {
      return SE::acos(x.get_RCP());
    }


    Expression
    atan(const Expression &x)
    {
      return SE::atan(x.get_RCP());
    }


    Expression
    atan2(const Expression &y, const Expression &x)
    {
      return SE::atan2(y.get_RCP(), x.get_RCP());
    }


    Expression
    acsc(const Expression &x)
    {
      return SE::acsc(x.get_RCP());
    }


    Expression
    asec(const Expression &x)
    {
      return SE::asec(x.get_RCP());
    }


    Expression
    acot(const Expression &x)
    {
      return SE::acot(x.get_RCP());
    }


    Expression
    sinh(const Expression &x)
    {
      return SE::sinh(x.get_RCP());
    }


    Expression
    cosh(const Expression &x)
    {
      return SE::cosh(x.get_RCP());
    }


    Expression
    tanh(const Expression &x)
    {
      return SE::tanh(x.get_RCP());
    }


    Expression
    csch(const Expression &x)
    {
      return SE::csch(x.get_RCP());
    }


    Expression
    sech(const Expression &x)
    {
      return SE::sech(x.get_RCP());
    }


    Expression
    coth(const Expression &x)
    {
      return SE::coth(x.get_RCP());
    }


    Expression
    asinh(const Expression &x)
    {
      return SE::asinh(x.get_RCP());
    }


    Expression
    acosh(const Expression &x)
    {
      return SE::acosh(x.get_RCP());
    }


    Expression
    atanh(const Expression &x)
    {
      return SE::atanh(x.get_RCP());
    }


    Expression
    acsch(const Expression &x)
    {
      return SE::acsch(x.get_RCP());
    }


    Expression
    asech(const Expression &x)
    {
      return SE::asech(x.get_RCP());
    }


    Expression
    acoth(const Expression &x)
    {
      return SE::acoth(x.get_RCP());
    }


    Expression
    abs(const Expression &x)
    {
      return SE::abs(x.get_RCP());
    }


    Expression
    fabs(const Expression &x)
    {
      return SE::abs(x.get_RCP());
    }


    Expression
    sign(const Expression &x)
    {
      return SE::sign(x.get_RCP());
    }


    Expression
    copysign(const Expression &value, const Expression &sign)
    {
      return value * Expression(SE::sign(sign.get_RCP()));
    }


    Expression
    floor(const Expression &x)
    {
      return SE::floor(x.get_RCP());
    }


    Expression
    ceil(const Expression &x)
    {
      return SE::ceiling(x.get_RCP());
    }


    Expression
    max(const Expression &a, const Expression &b)
    {
      return SE::max({a.get_RCP(), b.get_RCP()});
    }


    Expression
    min(const Expression &a, const Expression &b)
    {
      return SE::min({a.get_RCP(), b.get_RCP()});
    }


    Expression
    erf(const Expression &x)
    {
      return SE::erf(x.get_RCP());
    }


    Expression
    erfc(const Expression &x)
    {
      return SE::erfc(x.get_RCP());
    }

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
