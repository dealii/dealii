// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Check that the wrapper for symengine expressions can be used within
// math functions

#include <deal.II/differentiation/sd.h>

#include <math.h>

#include <fstream>
#include <iomanip>
#include <limits>

#include "../tests.h"

namespace SD = Differentiation::SD;
namespace SE = SymEngine;

template <typename NumberType1, typename NumberType2>
bool
same_number(const NumberType1 &val1,
            const NumberType2 &val2,
            const bool         print   = true,
            const double       tol_eps = 10.0)
{
  if (print == true)
    std::cout << "Diff: " << std::abs(val1 - val2) << std::endl;
  const long double eps1 = std::numeric_limits<NumberType1>::epsilon();
  const long double eps2 = std::numeric_limits<NumberType2>::epsilon();
  return std::abs(val1 - val2) <= tol_eps * std::max(eps1, eps2);
}

template <typename NumberType>
void
test_number_functions()
{
  typedef SD::Expression                SD_number_t;
  typedef SE::RCP<const SE::RealDouble> SE_double_t;
  // typedef SD_number_t result_t;
  typedef NumberType result_t;
  typedef double     float_result_t;

  // Primitive numbers
  const NumberType x_(NumberType(2.2));
  const NumberType a_(NumberType(3.3));
  const NumberType b_(NumberType(-4.4));
  const NumberType inv_a_(NumberType(1.0 / a_));
  // Raw SymEngine number pointers
  const SE_double_t x__(SE::real_double(x_));
  const SE_double_t a__(SE::real_double(a_));
  const SE_double_t b__(SE::real_double(b_));
  const SE_double_t inv_a__(SE::real_double(inv_a_));
  // Our wrapped numbers
  const SD_number_t x(x_);
  const SD_number_t a(a_);
  const SD_number_t b(b_);
  const SD_number_t inv_a(inv_a_);

  // SD_number_t f; // Temporary value; Check default constructor
  // SD_number_t::substitution_map_t sub_vals;

  // --- Values ---
  // NOTE: The returned results from SymEngineWrapper functions are actually
  // symbolic in nature, so we cast them to the number in order to evaluate
  // them.

  deallog.push("Power functions");
  {
    deallog << "pow(a_,x_): " << std::pow(a_, x_) << std::endl;
    deallog << "pow(a__,x__): " << *SE::pow(a__, x__) << std::endl;
    deallog << "pow(a,x):   " << static_cast<result_t>(pow(a, x)) << std::endl;
    deallog << "pow(a,x_):  " << static_cast<result_t>(pow(a, x_)) << std::endl;
    deallog << "pow(a_,x):  " << static_cast<result_t>(pow(a_, x)) << std::endl;
    Assert(same_number(static_cast<result_t>(pow(a, x)), std::pow(a_, x_)),
           ExcMessage("Incorrect result from pow function"));
    // --------------
    deallog << "sqrt(a_): " << std::sqrt(a_) << std::endl;
    deallog << "sqrt(a__): " << *SE::sqrt(a__) << std::endl;
    deallog << "sqrt(a):  " << sqrt(a) << std::endl;
    deallog << "sqrt(a):  " << static_cast<float_result_t>(sqrt(a))
            << std::endl;
    Assert(same_number(static_cast<result_t>(sqrt(a)), std::sqrt(a_)),
           ExcMessage("Incorrect result from sqrt function"));
    // --------------
    deallog << "cbrt(a_): " << std::cbrt(a_) << std::endl;
    deallog << "cbrt(a__): " << *SE::cbrt(a__) << std::endl;
    deallog << "cbrt(a):  " << cbrt(a) << std::endl;
    deallog << "cbrt(a):  " << static_cast<float_result_t>(cbrt(a))
            << std::endl;
    Assert(same_number(static_cast<result_t>(cbrt(a)), std::cbrt(a_)),
           ExcMessage("Incorrect result from cbrt function"));
    // --------------
    deallog << "exp(a_): " << std::exp(a_) << std::endl;
    deallog << "exp(a__): " << *SE::exp(a__) << std::endl;
    deallog << "exp(a):  " << exp(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "exp(a):  " << static_cast<float_result_t>(exp(a)) << std::endl;
    Assert(same_number(static_cast<result_t>(exp(a)), std::exp(a_)),
           ExcMessage("Incorrect result from exp function"));
    // --------------
    deallog << "log(a_): " << std::log(a_) << std::endl;
    deallog << "log(a__): " << *SE::log(a__) << std::endl;
    deallog << "log(a):  " << log(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "log(a):  " << static_cast<float_result_t>(log(a)) << std::endl;
    Assert(same_number(static_cast<result_t>(log(a)), std::log(a_)),
           ExcMessage("Incorrect result from log function"));
  }
  deallog.pop();

  deallog.push("Other functions");
  {
    deallog << "abs(b_): " << std::abs(b_) << std::endl;
    deallog << "abs(b):  " << static_cast<result_t>(abs(b)) << std::endl;
    Assert(same_number(static_cast<result_t>(abs(b)), std::abs(b_)),
           ExcMessage("Incorrect result from abs function"));
    // --------------
    deallog << "copysign(1.0, b_): " << std::copysign(1.0, b_) << std::endl;
    deallog << "sign(b):  " << static_cast<result_t>(sign(b)) << std::endl;
    Assert(same_number(static_cast<result_t>(sign(b)), std::copysign(1.0, b_)),
           ExcMessage("Incorrect result from sign function"));
    // --------------
    deallog << "copysign(a_, b_): " << std::copysign(a_, b_) << std::endl;
    deallog << "copysign(a, b):  " << static_cast<result_t>(copysign(a, b))
            << std::endl;
    Assert(same_number(static_cast<result_t>(copysign(a, b)),
                       std::copysign(a_, b_)),
           ExcMessage("Incorrect result from copysign function"));
    // --------------
    deallog << "floor(b_): " << std::floor(b_) << std::endl;
    deallog << "floor(b):  " << static_cast<result_t>(floor(b)) << std::endl;
    Assert(same_number(static_cast<result_t>(floor(b)), std::floor(b_)),
           ExcMessage("Incorrect result from floor function"));
    // --------------
    deallog << "ceil(b_): " << std::ceil(b_) << std::endl;
    deallog << "ceil(b):  " << static_cast<result_t>(ceil(b)) << std::endl;
    Assert(same_number(static_cast<result_t>(ceil(b)), std::ceil(b_)),
           ExcMessage("Incorrect result from ceil function"));
    // --------------
    deallog << "max(a_,x_): " << std::max(a_, x_) << std::endl;
    deallog << "max(a__,x__): " << *SE::max({a__, x__}) << std::endl;
    deallog << "max(a,x):   " << static_cast<result_t>(max(a, x)) << std::endl;
    deallog << "max(a,x_):  " << static_cast<result_t>(max(a, x_)) << std::endl;
    deallog << "max(a_,x):  " << static_cast<result_t>(max(a_, x)) << std::endl;
    Assert(same_number(static_cast<result_t>(max(a, x)), std::max(a_, x_)),
           ExcMessage("Incorrect result from max function"));
    // --------------
    deallog << "min(a_,x_): " << std::min(a_, x_) << std::endl;
    deallog << "min(a__,x__): " << *SE::min({a__, x__}) << std::endl;
    deallog << "min(a,x):   " << static_cast<result_t>(min(a, x)) << std::endl;
    deallog << "min(a,x_):  " << static_cast<result_t>(min(a, x_)) << std::endl;
    deallog << "min(a_,x):  " << static_cast<result_t>(min(a_, x)) << std::endl;
    Assert(same_number(static_cast<result_t>(min(a, x)), std::min(a_, x_)),
           ExcMessage("Incorrect result from max function"));
    // --------------
    deallog << "erf(a_): " << std::erf(a_) << std::endl;
    deallog << "erf(a__): " << *SE::erf(a__) << std::endl;
    deallog << "erf(a):  " << erf(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "erf(a):  " << static_cast<float_result_t>(erf(a)) << std::endl;
    Assert(same_number(static_cast<result_t>(erf(a)), std::erf(a_)),
           ExcMessage("Incorrect result from erf function"));
    // --------------
    deallog << "erfc(a_): " << std::erfc(a_) << std::endl;
    deallog << "erfc(a__): " << *SE::erfc(a__) << std::endl;
    deallog << "erfc(a):  " << erfc(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "erfc(a):  " << static_cast<float_result_t>(erfc(a))
            << std::endl;
    Assert(same_number(static_cast<result_t>(erfc(a)), std::erfc(a_)),
           ExcMessage("Incorrect result from erfc function"));
  }
  deallog.pop();

  // References:
  // https://en.wikipedia.org/wiki/Trigonometric_functions
  // https://en.wikipedia.org/wiki/Hyperbolic_function
  // http://mathworld.wolfram.com/InverseCosecant.html
  deallog.push("Trig functions");
  {
    deallog << "sin(a_): " << std::sin(a_) << std::endl;
    deallog << "sin(a__): " << *SE::sin(a__) << std::endl;
    deallog << "sin(a):  " << sin(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "sin(a):  " << static_cast<float_result_t>(sin(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(sin(a)), std::sin(a_)),
           ExcMessage("Incorrect result from sin function"));
    // --------------
    deallog << "cos(a_): " << std::cos(a_) << std::endl;
    deallog << "cos(a__): " << *SE::cos(a__) << std::endl;
    deallog << "cos(a):  " << cos(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "cos(a):  " << static_cast<float_result_t>(cos(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(cos(a)), std::cos(a_)),
           ExcMessage("Incorrect result from cos function"));
    // --------------
    deallog << "tan(a_): " << std::tan(a_) << std::endl;
    deallog << "tan(a__): " << *SE::tan(a__) << std::endl;
    deallog << "tan(a):  " << tan(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "tan(a):  " << static_cast<float_result_t>(tan(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(tan(a)), std::tan(a_)),
           ExcMessage("Incorrect result from tan function"));
  }
  deallog.pop();

  deallog.push("Reciprocal trig functions");
  {
    deallog << "csc(a_): " << (NumberType(1.0) / std::sin(a_)) << std::endl;
    deallog << "csc(a__): " << *SE::csc(a__) << std::endl;
    deallog << "csc(a):  " << csc(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "csc(a):  " << static_cast<float_result_t>(csc(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(csc(a)),
                       (NumberType(1.0) / std::sin(a_))),
           ExcMessage("Incorrect result from csc function"));
    // --------------
    deallog << "sec(a_): " << (NumberType(1.0) / std::cos(a_)) << std::endl;
    deallog << "sec(a__): " << *SE::sec(a__) << std::endl;
    deallog << "sec(a):  " << sec(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "sec(a):  " << static_cast<float_result_t>(sec(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(sec(a)),
                       (NumberType(1.0) / std::cos(a_))),
           ExcMessage("Incorrect result from sec function"));
    // --------------
    deallog << "cot(a_): " << (NumberType(1.0) / std::tan(a_)) << std::endl;
    deallog << "cot(a__): " << *SE::cot(a__) << std::endl;
    deallog << "cot(a):  " << cot(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "cot(a):  " << static_cast<float_result_t>(cot(a)) << std::endl;
    Assert(same_number(static_cast<float_result_t>(cot(a)),
                       (NumberType(1.0) / std::tan(a_))),
           ExcMessage("Incorrect result from cot function"));
  }
  deallog.pop();

  deallog.push("Inverse trig functions");
  {
    deallog << "asin(inv_a_): " << std::asin(inv_a_) << std::endl;
    deallog << "asin(inv_a__): " << *SE::asin(inv_a__) << std::endl;
    deallog << "asin(inv_a):  " << asin(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "asin(inv_a):  " << static_cast<float_result_t>(asin(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(asin(inv_a)),
                       std::asin(inv_a_)),
           ExcMessage("Incorrect result from asin function"));
    // --------------
    deallog << "acos(inv_a_): " << std::acos(inv_a_) << std::endl;
    deallog << "acos(inv_a__): " << *SE::acos(inv_a__) << std::endl;
    deallog << "acos(inv_a):  " << acos(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acos(inv_a):  " << static_cast<float_result_t>(acos(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acos(inv_a)),
                       std::acos(inv_a_)),
           ExcMessage("Incorrect result from acos function"));
    // --------------
    deallog << "atan(inv_a_): " << std::atan(inv_a_) << std::endl;
    deallog << "atan(inv_a__): " << *SE::atan(inv_a__) << std::endl;
    deallog << "atan(inv_a):  " << atan(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "atan(inv_a):  " << static_cast<float_result_t>(atan(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(atan(inv_a)),
                       std::atan(inv_a_)),
           ExcMessage("Incorrect result from atan function"));
    // --------------
    deallog << "atan2(x_,a_): " << std::atan2(x_, a_) << std::endl;
    deallog << "atan2(x__,a__): " << *SE::atan2(x__, a__) << std::endl;
    deallog << "atan2(x,a):  " << atan2(x, a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "atan2(x,a):  " << static_cast<float_result_t>(atan2(x, a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(atan2(x, a)),
                       std::atan2(x_, a_)),
           ExcMessage("Incorrect result from atan2 function"));
  }
  deallog.pop();

  deallog.push("Reciprocal inverse trig functions");
  {
    deallog << "acsc(a_): " << std::asin(NumberType(1.0) / a_) << std::endl;
    deallog << "acsc(a__): " << *SE::acsc(a__) << std::endl;
    deallog << "acsc(a):  " << acsc(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acsc(a):  " << static_cast<float_result_t>(acsc(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acsc(a)),
                       std::asin(NumberType(1.0) / a_)),
           ExcMessage("Incorrect result from acsc function"));
    // --------------
    deallog << "asec(a_): " << std::acos(NumberType(1.0) / a_) << std::endl;
    deallog << "asec(a__): " << *SE::asec(a__) << std::endl;
    deallog << "asec(a):  " << asec(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "asec(a):  " << static_cast<float_result_t>(asec(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(asec(a)),
                       std::acos(NumberType(1.0) / a_)),
           ExcMessage("Incorrect result from asec function"));
    // --------------
    deallog << "acot(a_): " << std::atan(NumberType(1.0) / a_) << std::endl;
    deallog << "acot(a__): " << *SE::acot(a__) << std::endl;
    deallog << "acot(a):  " << acot(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acot(a):  " << static_cast<float_result_t>(acot(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acot(a)),
                       std::atan(NumberType(1.0) / a_)),
           ExcMessage("Incorrect result from acot function"));
  }
  deallog.pop();

  deallog.push("Hyperbolic trig functions");
  {
    deallog << "sinh(a_): " << std::sinh(a_) << std::endl;
    deallog << "sinh(a__): " << *SE::sinh(a__) << std::endl;
    deallog << "sinh(a):  " << sinh(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "sinh(a):  " << static_cast<float_result_t>(sinh(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(sinh(inv_a)),
                       std::sinh(inv_a_)),
           ExcMessage("Incorrect result from sinh function"));
    // --------------
    deallog << "cosh(a_): " << std::cosh(a_) << std::endl;
    deallog << "cosh(a__): " << *SE::cosh(a__) << std::endl;
    deallog << "cosh(a):  " << cosh(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "cosh(a):  " << static_cast<float_result_t>(cosh(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(cosh(inv_a)),
                       std::cosh(inv_a_)),
           ExcMessage("Incorrect result from cosh function"));
    // --------------
    deallog << "tanh(a_): " << std::tanh(a_) << std::endl;
    deallog << "tanh(a__): " << *SE::tanh(a__) << std::endl;
    deallog << "tanh(a):  " << tanh(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "tanh(a):  " << static_cast<float_result_t>(tanh(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(tanh(inv_a)),
                       std::tanh(inv_a_)),
           ExcMessage("Incorrect result from tanh function"));
  }
  deallog.pop();

  deallog.push("Reciprocal hyperbolic trig functions");
  {
    deallog << "csch(a_): " << (NumberType(1.0) / std::sinh(a_)) << std::endl;
    deallog << "csch(a__): " << *SE::csch(a__) << std::endl;
    deallog << "csch(a):  " << csch(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "csch(a):  " << static_cast<float_result_t>(csch(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(csch(a)),
                       (NumberType(1.0) / std::sinh(a_))),
           ExcMessage("Incorrect result from csch function"));
    // --------------
    deallog << "sech(a_): " << (NumberType(1.0) / std::cosh(a_)) << std::endl;
    deallog << "sech(a__): " << *SE::sech(a__) << std::endl;
    deallog << "sech(a):  " << sech(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "sech(a):  " << static_cast<float_result_t>(sech(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(sech(a)),
                       (NumberType(1.0) / std::cosh(a_))),
           ExcMessage("Incorrect result from sech function"));
    // --------------
    deallog << "coth(a_): " << (NumberType(1.0) / std::tanh(a_)) << std::endl;
    deallog << "coth(a__): " << *SE::coth(a__) << std::endl;
    deallog << "coth(a):  " << coth(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "coth(a):  " << static_cast<float_result_t>(coth(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(coth(a)),
                       (NumberType(1.0) / std::tanh(a_))),
           ExcMessage("Incorrect result from coth function"));
  }
  deallog.pop();

  deallog.push("Inverse hyperbolic trig functions");
  {
    deallog << "asinh(inv_a_): " << std::asinh(inv_a_) << std::endl;
    deallog << "asinh(inv_a__): " << *SE::asinh(inv_a__) << std::endl;
    deallog << "asinh(inv_a):  " << asinh(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "asinh(inv_a):  " << static_cast<float_result_t>(asinh(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(asinh(inv_a)),
                       std::asinh(inv_a_)),
           ExcMessage("Incorrect result from asinh function"));
    // --------------
    deallog << "acosh(a_): " << std::acosh(a_) << std::endl;
    deallog << "acosh(a__): " << *SE::acosh(a__) << std::endl;
    deallog << "acosh(a):  " << acosh(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acosh(a):  " << static_cast<float_result_t>(acosh(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acosh(a)), std::acosh(a_)),
           ExcMessage("Incorrect result from acosh function"));
    // --------------
    deallog << "atanh(inv_a_): " << std::atanh(inv_a_) << std::endl;
    deallog << "atanh(inv_a__): " << *SE::atanh(inv_a__) << std::endl;
    deallog << "atanh(inv_a):  " << atanh(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "atanh(inv_a):  " << static_cast<float_result_t>(atanh(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(atanh(inv_a)),
                       std::atanh(inv_a_)),
           ExcMessage("Incorrect result from atanh function"));
  }
  deallog.pop();

  deallog.push("Reciprocal inverse hyperbolic trig functions");
  {
    deallog << "acsch(a_): " << std::asinh(NumberType(1.0) / a_) << std::endl;
    deallog << "acsch(a__): " << *SE::acsch(a__) << std::endl;
    deallog << "acsch(a):  " << acsch(a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acsch(a):  " << static_cast<float_result_t>(acsch(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acsch(a)),
                       std::asinh(NumberType(1.0) / a_)),
           ExcMessage("Incorrect result from acsch function"));
    // --------------
    deallog << "asech(inv_a_): " << std::acosh(a_) << std::endl;
    deallog << "asech(inv_a__): " << *SE::asech(inv_a__) << std::endl;
    deallog << "asech(inv_a):  " << asech(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "asech(inv_a):  " << static_cast<float_result_t>(asech(inv_a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(asech(inv_a)),
                       std::acosh(a_)),
           ExcMessage("Incorrect result from asech function"));
    // --------------
    deallog << "acoth(a_): " << std::atanh(NumberType(1.0) / a_) << std::endl;
    deallog << "acoth(a__): " << *SE::acoth(a__) << std::endl;
    deallog << "acoth(a):  " << acoth(inv_a)
            << std::endl; // Returns a SymEngine function, not a number
    deallog << "acoth(a):  " << static_cast<float_result_t>(acoth(a))
            << std::endl;
    Assert(same_number(static_cast<float_result_t>(acoth(a)),
                       std::atanh(NumberType(1.0) / a_)),
           ExcMessage("Incorrect result from acoth function"));
  }
  deallog.pop();
}

int
main()
{
  initlog();

  // Most of the functions tested return irrational values. Making a test
  // work for both integer and floating point values is hard!
  // deallog.push("Integer");
  // test_number_functions<int>();

  deallog.push("Float");
  test_number_functions<float>();
  deallog.pop();

  deallog.push("Double");
  test_number_functions<double>();
  deallog.pop();

  // Not available yet
  // SymEngine::SymEngineException: Invalid Format: Expected Integer or Rational
  // deallog.push("Complex double");
  // test_number_functions<std::complex<double>>();
  // deallog.pop();

  deallog << "OK" << std::endl;
}
