//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <cstdio>

DEAL_II_NAMESPACE_OPEN

namespace Functions {
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
		      "function, separated by ','.");

				     // The expression of the function
    std::string expr = "0";
    for (unsigned int i=1; i<n_components; ++i)
      expr += "; 0";

    prm.declare_entry("Function expression", expr, Patterns::Anything(),
		      "Separate vector valued expressions by ';' as ',' "
		      "is used internally by the function parser.");
    prm.declare_entry("Function constants", "", Patterns::Anything(),
		      "Any constant used inside the function which is not "
		      "a variable name.");
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
    for(unsigned int i = 0; i < const_list.size(); ++i)
      {
	std::vector<std::string> this_c = 
	  Utilities::split_string_list(const_list[i], '=');
	AssertThrow(this_c.size() == 2, ExcMessage("Invalid format"));
	double tmp;
	AssertThrow( std::sscanf(this_c[1].c_str(), "%lf", &tmp),
		     ExcMessage("Double number?"));
	constants[this_c[0]] = tmp;
      }
    
    constants["pi"] = numbers::PI;
    constants["Pi"] = numbers::PI;

    unsigned int nn = (Utilities::split_string_list(vnames)).size();
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
	      AssertThrow(false, ExcMessage("Not the correct size. Check your code."));
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
