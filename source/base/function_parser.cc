// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/function_parser.h>
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>

#include <map>

DEAL_II_NAMESPACE_OPEN


template <int dim>
const std::vector<std::string> &
FunctionParser<dim>::get_expressions() const
{
  return this->expressions;
}



template <int dim>
FunctionParser<dim>::FunctionParser(const unsigned int n_components,
                                    const double       initial_time,
                                    const double       h)
  : AutoDerivativeFunction<dim>(h, n_components, initial_time)
{}


template <int dim>
FunctionParser<dim>::FunctionParser(const std::string &expression,
                                    const std::string &constants,
                                    const std::string &variable_names,
                                    const double       h)
  : AutoDerivativeFunction<dim>(
      h,
      Utilities::split_string_list(expression, ';').size())
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



template <int dim>
void
FunctionParser<dim>::initialize(const std::string              &variables,
                                const std::vector<std::string> &expressions,
                                const std::map<std::string, double> &constants,
                                const bool time_dependent)
{
  AssertThrow(this->n_components == expressions.size(),
              ExcInvalidExpressionSize(this->n_components, expressions.size()));
  internal::FunctionParser::ParserImplementation<dim, double>::initialize(
    variables, expressions, constants, time_dependent);
}



template <int dim>
void
FunctionParser<dim>::initialize(const std::string                   &vars,
                                const std::string                   &expression,
                                const std::map<std::string, double> &constants,
                                const bool time_dependent)
{
  initialize(vars,
             Utilities::split_string_list(expression, ';'),
             constants,
             time_dependent);
}



template <int dim>
double
FunctionParser<dim>::value(const Point<dim>  &p,
                           const unsigned int component) const
{
  return this->do_value(p, this->get_time(), component);
}

// Explicit Instantiations.

template class FunctionParser<1>;
template class FunctionParser<2>;
template class FunctionParser<3>;

DEAL_II_NAMESPACE_CLOSE
