// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mu_parser_internal_h
#define dealii_mu_parser_internal_h

// This file contains functions used internally by the FunctionParser
// and the TensorFunctionParser class.

#include <deal.II/base/config.h>

#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN



#ifdef DEAL_II_WITH_MUPARSER

namespace internal
{
  namespace FunctionParser
  {
    int
    mu_round(double val);

    double
    mu_if(double condition, double thenvalue, double elsevalue);

    double
    mu_or(double left, double right);

    double
    mu_and(double left, double right);

    double
    mu_int(double value);

    double
    mu_ceil(double value);

    double
    mu_floor(double value);

    double
    mu_cot(double value);

    double
    mu_csc(double value);

    double
    mu_sec(double value);

    double
    mu_log(double value);

    double
    mu_pow(double a, double b);

    double
    mu_erfc(double value);

    // returns a random value in the range [0,1] initializing the generator
    // with the given seed
    double
    mu_rand_seed(double seed);

    // returns a random value in the range [0,1]
    double
    mu_rand();

    extern std::vector<std::string> function_names;

  } // namespace FunctionParser

} // namespace internal
#endif



DEAL_II_NAMESPACE_CLOSE

#endif
