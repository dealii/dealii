// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim>
  ParsedFunction<dim>::ParsedFunction(const unsigned int n_components,
                                      const double       h)
    : AutoDerivativeFunction<dim>(h, n_components)
    , function_object(n_components)
  {}



  template <int dim>
  void
  ParsedFunction<dim>::declare_parameters(ParameterHandler  &prm,
                                          const unsigned int n_components,
                                          const std::string &input_expr)
  {
    Assert(n_components > 0, ExcZero());

    std::string vnames;
    switch (dim)
      {
        case 1:
          vnames = "x,t";
          break;
        case 2:
          vnames = "x,y,t";
          break;
        case 3:
          vnames = "x,y,z,t";
          break;
        default:
          AssertThrow(false, ExcNotImplemented());
          break;
      }
    prm.declare_entry(
      "Variable names",
      vnames,
      Patterns::Anything(),
      "Names of the independent variables, separated by commas.");

    // The expression of the function
    // If the string is an empty string, 0 is set for each components.
    std::string expr = input_expr;
    if (expr == "")
      {
        expr = "0";
        for (unsigned int i = 1; i < n_components; ++i)
          expr += "; 0";
      }
    else
      {
        // If the user specified an input expr, the number of component
        // specified need to match n_components.
        AssertDimension((std::count(expr.begin(), expr.end(), ';') + 1),
                        n_components);
      }


    prm.declare_entry(
      "Function expression",
      expr,
      Patterns::Anything(),
      "Semicolon-separated formulas for the function components. Supports standard "
      "operations, functions `sin`, `cos`, etc., and conditionals `if(x>0, 1, -1)`.");
    prm.declare_entry(
      "Function constants",
      "",
      Patterns::Anything(),
      "Symbolic constants for the function expression, in the form "
      "'var1=value1, var2=value2, ...'. Note: 'pi', 'Pi', 'PI', and 'E' are predefined.");
  }



  template <int dim>
  void
  ParsedFunction<dim>::parse_parameters(ParameterHandler &prm)
  {
    std::string vnames         = prm.get("Variable names");
    std::string expression     = prm.get("Function expression");
    std::string constants_list = prm.get("Function constants");

    std::vector<std::string> const_list =
      Utilities::split_string_list(constants_list, ',');
    std::map<std::string, double> constants;

    // set pi, Pi, and PI as synonyms for the corresponding value, and set the
    // Nepero number E.
    constants["pi"] = numbers::PI;
    constants["Pi"] = numbers::PI;
    constants["PI"] = numbers::PI;
    constants["E"]  = numbers::E;

    for (const auto &constant : const_list)
      {
        std::vector<std::string> this_c =
          Utilities::split_string_list(constant, '=');
        AssertThrow(this_c.size() == 2,
                    ExcMessage("The list of constants, <" + constants_list +
                               ">, is not a comma-separated list of "
                               "entries of the form 'name=value'."));
        constants[this_c[0]] = Utilities::string_to_double(this_c[1]);
      }

    const unsigned int nn = (Utilities::split_string_list(vnames)).size();
    switch (nn)
      {
        case dim:
          // Time independent function
          function_object.initialize(vnames, expression, constants);
          break;
        case dim + 1:
          // Time dependent function
          function_object.initialize(vnames, expression, constants, true);
          break;
        default:
          AssertThrow(false,
                      ExcMessage(
                        "The list of variables specified is <" + vnames +
                        "> which is a list of length " +
                        Utilities::int_to_string(nn) +
                        " but it has to be a list of length equal to" +
                        " either dim (for a time-independent function)" +
                        " or dim+1 (for a time-dependent function)."));
      }
  }



  template <int dim>
  void
  ParsedFunction<dim>::vector_value(const Point<dim> &p,
                                    Vector<double>   &values) const
  {
    function_object.vector_value(p, values);
  }



  template <int dim>
  double
  ParsedFunction<dim>::value(const Point<dim> &p, unsigned int comp) const
  {
    return function_object.value(p, comp);
  }



  template <int dim>
  void
  ParsedFunction<dim>::set_time(const double newtime)
  {
    function_object.set_time(newtime);
    AutoDerivativeFunction<dim>::set_time(newtime);
  }


  // Explicit instantiations
  template class ParsedFunction<1>;
  template class ParsedFunction<2>;
  template class ParsedFunction<3>;
} // namespace Functions
DEAL_II_NAMESPACE_CLOSE
