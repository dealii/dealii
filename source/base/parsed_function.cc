// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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

#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <cstdio>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim>
  ParsedFunction<dim>::ParsedFunction (const unsigned int n_components, const double h)
    :
    AutoDerivativeFunction<dim>(h, n_components),
    function_object(n_components)
  {}



  template <int dim>
  void
  ParsedFunction<dim>::declare_parameters(ParameterHandler  &prm,
                                          const unsigned int n_components)
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
    prm.declare_entry("Variable names", vnames, Patterns::Anything(),
                      "The name of the variables as they will be used in the "
                      "function, separated by commas. By default, the names of variables "
                      "at which the function will be evaluated is `x' (in 1d), `x,y' (in 2d) or "
                      "`x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then "
                      "use these variable names in your function expression and they will be "
                      "replaced by the values of these variables at which the function is "
                      "currently evaluated. However, you can also choose a different set "
                      "of names for the independent variables at which to evaluate your function "
                      "expression. For example, if you work in spherical coordinates, you may "
                      "wish to set this input parameter to `r,phi,theta,t' and then use these "
                      "variable names in your function expression.");

    // The expression of the function
    std::string expr = "0";
    for (unsigned int i=1; i<n_components; ++i)
      expr += "; 0";

    prm.declare_entry("Function expression", expr, Patterns::Anything(),
                      "The formula that denotes the function you want to evaluate for "
                      "particular values of the independent variables. This expression "
                      "may contain any of the usual operations such as addition or "
                      "multiplication, as well as all of the common functions such as "
                      "`sin' or `cos'. In addition, it may contain expressions like "
                      "`if(x>0, 1, -1)' where the expression evaluates to the second "
                      "argument if the first argument is true, and to the third argument "
                      "otherwise. For a full overview of possible expressions accepted "
                      "see the documentation of the muparser library at http://muparser.beltoforion.de/."
                      "\n\n"
                      "If the function you are describing represents a vector-valued "
                      "function with multiple components, then separate the expressions "
                      "for individual components by a semicolon.");
    prm.declare_entry("Function constants", "", Patterns::Anything(),
                      "Sometimes it is convenient to use symbolic constants in the "
                      "expression that describes the function, rather than having to "
                      "use its numeric value everywhere the constant appears. These "
                      "values can be defined using this parameter, in the form "
                      "`var1=value1, var2=value2, ...'."
                      "\n\n"
                      "A typical example would be to set this runtime parameter to "
                      "`pi=3.1415926536' and then use `pi' in the expression of the "
                      "actual formula. (That said, for convenience this class actually "
                      "defines both `pi' and `Pi' by default, but you get the idea.)");
  }



  template <int dim>
  void ParsedFunction<dim>::parse_parameters(ParameterHandler &prm)
  {
    std::string vnames = prm.get("Variable names");
    std::string expression = prm.get("Function expression");
    std::string constants_list = prm.get("Function constants");

    std::vector<std::string> const_list =
      Utilities::split_string_list(constants_list, ',');
    std::map<std::string, double> constants;
    for (unsigned int i = 0; i < const_list.size(); ++i)
      {
        std::vector<std::string> this_c =
          Utilities::split_string_list(const_list[i], '=');
        AssertThrow(this_c.size() == 2, ExcMessage("Invalid format"));
        double tmp;
        AssertThrow( std::sscanf(this_c[1].c_str(), "%lf", &tmp),
                     ExcMessage("Double number?"));
        constants[this_c[0]] = tmp;
      }

    // set pi and Pi as synonyms for the corresponding value. note that
    // this overrides any value a user may have given
    constants["pi"] = numbers::PI;
    constants["Pi"] = numbers::PI;

    const unsigned int nn = (Utilities::split_string_list(vnames)).size();
    switch (nn)
      {
      case dim:
        // Time independent function
        function_object.initialize(vnames, expression, constants);
        break;
      case dim+1:
        // Time dependent function
        function_object.initialize(vnames, expression, constants, true);
        break;
      default:
        AssertThrow(false,
                    ExcMessage("The list of variables specified is <" + vnames
                               + "> which is a list of length "
                               + Utilities::int_to_string(nn)
                               + " but it has to be a list of length equal to"
                               + " either dim (for a time-independent function)"
                               + " or dim+1 (for a time-dependent function)."));
      }
  }



  template <int dim>
  void ParsedFunction<dim>::vector_value (const Point<dim> &p,
                                          Vector<double>   &values) const
  {
    function_object.vector_value(p, values);
  }



  template <int dim>
  double ParsedFunction<dim>::value (const Point<dim>   &p,
                                     unsigned int comp) const
  {
    return function_object.value(p, comp);
  }



  template <int dim>
  void ParsedFunction<dim>::set_time (const double newtime)
  {
    function_object.set_time(newtime);
    AutoDerivativeFunction<dim>::set_time(newtime);
  }


// Explicit instantiations
  template class ParsedFunction<1>;
  template class ParsedFunction<2>;
  template class ParsedFunction<3>;
}
DEAL_II_NAMESPACE_CLOSE
