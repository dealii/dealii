// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/tensor_function_parser.h>
#include <deal.II/base/utilities.h>

#include <map>

DEAL_II_NAMESPACE_OPEN


template <int rank, int dim, typename Number>
const std::vector<std::string> &
TensorFunctionParser<rank, dim, Number>::get_expressions() const
{
  return this->expressions;
}



template <int rank, int dim, typename Number>
TensorFunctionParser<rank, dim, Number>::TensorFunctionParser(
  const double initial_time)
  : TensorFunction<rank, dim, Number>(initial_time)
  , n_components(Utilities::pow(dim, rank))
{}


template <int rank, int dim, typename Number>
TensorFunctionParser<rank, dim, Number>::TensorFunctionParser(
  const std::string &expression,
  const std::string &constants,
  const std::string &variable_names)
  : TensorFunction<rank, dim, Number>()
  , n_components(Utilities::pow(dim, rank))
{
  auto constants_map = Patterns::Tools::Convert<ConstMap>::to_value(
    constants,
    Patterns::Map(Patterns::Anything(),
                  Patterns::Double(),
                  0,
                  Patterns::Map::max_int_value,
                  ",",
                  "="));
  initialize(variable_names,
             expression,
             constants_map,
             Utilities::split_string_list(variable_names, ",").size() ==
               dim + 1);
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string                   &variables,
  const std::vector<std::string>      &expressions,
  const std::map<std::string, double> &constants,
  const bool                           time_dependent)
{
  AssertThrow(this->n_components == expressions.size(),
              ExcInvalidExpressionSize(this->n_components, expressions.size()));
  internal::FunctionParser::ParserImplementation<dim, Number>::initialize(
    variables, expressions, constants, time_dependent);
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string                   &vars,
  const std::string                   &expression,
  const std::map<std::string, double> &constants,
  const bool                           time_dependent)
{
  initialize(vars,
             Utilities::split_string_list(expression, ';'),
             constants,
             time_dependent);
}



template <int rank, int dim, typename Number>
Tensor<rank, dim, Number>
TensorFunctionParser<rank, dim, Number>::value(const Point<dim> &p) const
{
  std::array<Number, Tensor<rank, dim, Number>::n_independent_components>
       values;
  auto values_view = make_array_view(values.begin(), values.end());
  this->do_all_values(p, this->get_time(), values_view);

  return Tensor<rank, dim, Number>(
    make_array_view(values.begin(), values.end()));
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::value_list(
  const std::vector<Point<dim>>          &p,
  std::vector<Tensor<rank, dim, Number>> &values) const
{
  Assert(p.size() == values.size(),
         ExcDimensionMismatch(p.size(), values.size()));

  for (unsigned int i = 0; i < p.size(); ++i)
    {
      values[i] = value(p[i]);
    }
}

// explicit instantiations
#include "base/tensor_function_parser.inst"

DEAL_II_NAMESPACE_CLOSE
