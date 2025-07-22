// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_symbolic_function_templates_h
#define dealii_symbolic_function_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/symbolic_function.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
#ifdef DEAL_II_WITH_SYMENGINE

  namespace internal
  {
    inline std::vector<Differentiation::SD::Expression>
    to_expressions(const std::vector<std::string> &strings)
    {
      std::vector<Differentiation::SD::Expression> ret(strings.size());
      for (unsigned int i = 0; i < strings.size(); ++i)
        ret[i] = Differentiation::SD::Expression(strings[i], true);
      return ret;
    }
  } // namespace internal



  template <int dim, typename RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>::SymbolicFunction(
    const std::vector<Differentiation::SD::Expression>    &functions,
    const Tensor<1, dim, Differentiation::SD::Expression> &argument,
    const Differentiation::SD::Expression                 &time,
    const Differentiation::SD::types::substitution_map &user_substitution_map)
    : Function<dim, RangeNumberType>(functions.size(), 0.0)
    , user_function(functions)
    , user_substitution_map(user_substitution_map)
    , coordinate_symbols(argument)
    , time_symbol(time)
  {}



  template <int dim, typename RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>::SymbolicFunction(
    const std::string &expressions)
    : Function<dim, RangeNumberType>(
        Utilities::split_string_list(expressions, ';').size(),
        0.0)
    , user_function(internal::to_expressions(
        Utilities::split_string_list(expressions, ';')))
    , coordinate_symbols(get_default_coordinate_symbols())
    , time_symbol(Differentiation::SD::make_symbol("t"))
  {}



  template <int dim, typename RangeNumberType>
  void
  SymbolicFunction<dim, RangeNumberType>::update_values() const
  {
    // This is only necessary if the function vector is empty. A call to
    // update_user_substitution_map() will clear all internal vectors.
    if (function.empty())
      {
        function.resize(user_function.size());
        for (unsigned int i = 0; i < user_function.size(); ++i)
          {
            function[i] = user_function[i].substitute(user_substitution_map);
          }
      }
  }



  template <int dim, typename RangeNumberType>
  void
  SymbolicFunction<dim, RangeNumberType>::update_first_derivatives() const
  {
    // This is only necessary if the gradient vector is empty. A call to
    // update_user_substitution_map() will clear all internal vectors.
    if (function_gradient.empty())
      {
        update_values();
        function_gradient.resize(user_function.size());
        for (unsigned int i = 0; i < user_function.size(); ++i)
          {
            function_gradient[i] =
              Differentiation::SD::differentiate(function[i],
                                                 coordinate_symbols);
          }
      }
  }



  template <int dim, typename RangeNumberType>
  void
  SymbolicFunction<dim, RangeNumberType>::update_second_derivatives() const
  {
    // This is only necessary if the gradient vector is empty. A call to
    // update_user_substitution_map() will clear all internal vectors.
    if (function_hessian.empty())
      {
        update_first_derivatives();
        function_hessian.resize(user_function.size());
        function_laplacian.resize(user_function.size());

        for (unsigned int i = 0; i < user_function.size(); ++i)
          {
            function_hessian[i] =
              Differentiation::SD::differentiate(function_gradient[i],
                                                 coordinate_symbols);
            function_laplacian[i] = trace(function_hessian[i]);
          }
      }
  }


  template <int dim, typename RangeNumberType>
  void
  SymbolicFunction<dim, RangeNumberType>::update_user_substitution_map(
    const Differentiation::SD::types::substitution_map &substitutions)
  {
    user_substitution_map = substitutions;

    // Now reset all internal vectors.
    function.resize(0);
    function_gradient.resize(0);
    function_hessian.resize(0);
    function_laplacian.resize(0);
  }



  template <int dim, typename RangeNumberType>
  void
  SymbolicFunction<dim, RangeNumberType>::set_additional_function_arguments(
    const Differentiation::SD::types::substitution_map &arguments)
  {
    additional_function_arguments = arguments;
  }



  template <int dim, typename RangeNumberType>
  Tensor<1, dim, Differentiation::SD::Expression>
  SymbolicFunction<dim, RangeNumberType>::get_default_coordinate_symbols()
  {
    static_assert(dim <= 3, "Not implemented yet.");
    const std::vector<std::string> names = {"x", "y", "z"};

    Tensor<1, dim, Differentiation::SD::Expression> x;
    for (unsigned int i = 0; i < dim; ++i)
      x[i] = Differentiation::SD::make_symbol(names[i]);
    return x;
  }



  template <int dim, typename RangeNumberType>
  const Tensor<1, dim, Differentiation::SD::Expression> &
  SymbolicFunction<dim, RangeNumberType>::get_coordinate_symbols() const
  {
    return coordinate_symbols;
  }



  template <int dim, typename RangeNumberType>
  const Differentiation::SD::Expression &
  SymbolicFunction<dim, RangeNumberType>::get_time_symbol() const
  {
    return time_symbol;
  }



  template <int dim, typename RangeNumberType>
  const std::vector<Differentiation::SD::Expression> &
  SymbolicFunction<dim, RangeNumberType>::get_symbolic_function_expressions()
    const
  {
    return user_function;
  }



  template <int dim, typename RangeNumberType>
  const Differentiation::SD::types::substitution_map &
  SymbolicFunction<dim, RangeNumberType>::get_user_substitution_map() const
  {
    return user_substitution_map;
  }



  template <int dim, typename RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>::time_derivative() const
  {
    std::vector<Differentiation::SD::Expression> df_dt(user_function.size());
    for (unsigned int i = 0; i < user_function.size(); ++i)
      df_dt[i] = user_function[i].differentiate(time_symbol);
    return SymbolicFunction<dim, RangeNumberType>(df_dt,
                                                  coordinate_symbols,
                                                  time_symbol,
                                                  user_substitution_map);
  }



  template <int dim, typename RangeNumberType>
  RangeNumberType
  SymbolicFunction<dim, RangeNumberType>::value(
    const Point<dim>  &p,
    const unsigned int component) const
  {
    update_values();
    AssertIndexRange(component, function.size());
    return Differentiation::SD::substitute_and_evaluate<RangeNumberType>(
      function[component], create_evaluation_substitution_map(p));
  }



  template <int dim, typename RangeNumberType>
  Tensor<1, dim, RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>::gradient(
    const Point<dim>  &p,
    const unsigned int component) const
  {
    update_first_derivatives();
    AssertIndexRange(component, function.size());
    return Differentiation::SD::substitute_and_evaluate<RangeNumberType>(
      function_gradient[component], create_evaluation_substitution_map(p));
  }



  template <int dim, typename RangeNumberType>
  RangeNumberType
  SymbolicFunction<dim, RangeNumberType>::laplacian(
    const Point<dim>  &p,
    const unsigned int component) const
  {
    update_second_derivatives();
    AssertIndexRange(component, function.size());
    return Differentiation::SD::substitute_and_evaluate<RangeNumberType>(
      function_laplacian[component], create_evaluation_substitution_map(p));
  }



  template <int dim, typename RangeNumberType>
  SymmetricTensor<2, dim, RangeNumberType>
  SymbolicFunction<dim, RangeNumberType>::hessian(
    const Point<dim>  &p,
    const unsigned int component) const
  {
    update_second_derivatives();
    AssertIndexRange(component, function.size());
    return SymmetricTensor<2, dim, RangeNumberType>(
      Differentiation::SD::substitute_and_evaluate<RangeNumberType>(
        function_hessian[component], create_evaluation_substitution_map(p)));
  }



  template <int dim, typename RangeNumberType>
  Differentiation::SD::types::substitution_map
  SymbolicFunction<dim, RangeNumberType>::create_evaluation_substitution_map(
    const Point<dim> &p) const
  {
    auto map = additional_function_arguments;
    Differentiation::SD::add_to_substitution_map(
      map,
      std::make_pair(time_symbol, this->get_time()),
      std::make_pair(coordinate_symbols, p));
    return map;
  }
#endif
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
