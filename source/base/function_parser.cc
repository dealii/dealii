// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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


#include <deal.II/base/function_parser.h>
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>

#include <cmath>
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
FunctionParser<dim>::initialize(const std::string &             variables,
                                const std::vector<std::string> &expressions,
                                const std::map<std::string, double> &constants,
                                const bool time_dependent)
{
  this->parser_data.clear(); // this will reset all thread-local objects

  this->constants   = constants;
  this->var_names   = Utilities::split_string_list(variables, ',');
  this->expressions = expressions;
  AssertThrow(((time_dependent) ? dim + 1 : dim) == this->var_names.size(),
              ExcMessage("Wrong number of variables"));

  // We check that the number of components of this function matches the
  // number of components passed in as a vector of strings.
  AssertThrow(this->n_components == expressions.size(),
              ExcInvalidExpressionSize(this->n_components, expressions.size()));

  // Now we define how many variables we expect to read in.  We distinguish
  // between two cases: Time dependent problems, and not time dependent
  // problems. In the first case the number of variables is given by the
  // dimension plus one. In the other case, the number of variables is equal
  // to the dimension. Once we parsed the variables string, if none of this is
  // the case, then an exception is thrown.
  if (time_dependent)
    this->n_vars = dim + 1;
  else
    this->n_vars = dim;

  // create a parser object for the current thread we can then query in
  // value() and vector_value(). this is not strictly necessary because a user
  // may never call these functions on the current thread, but it gets us
  // error messages about wrong formulas right away
  this->init_muparser();

  // finally set the initialization bit
  this->initialized = true;
}



template <int dim>
void
FunctionParser<dim>::initialize(const std::string &                  vars,
                                const std::string &                  expression,
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
FunctionParser<dim>::value(const Point<dim> & p,
                           const unsigned int component) const
{
  return this->do_value(p, this->get_time(), component);
}

// Explicit Instantiations.

template class FunctionParser<1>;
template class FunctionParser<2>;
template class FunctionParser<3>;

DEAL_II_NAMESPACE_CLOSE
