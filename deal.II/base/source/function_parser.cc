//----------------------------  function_parser.cc ---------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_parser.cc ---------------------

#include <base/function_parser.h>
#include <base/point.h>
#include <lac/vector.h>

template <int dim>
FunctionParser<dim>::FunctionParser(const unsigned int n_components,
				    const double       initial_time)
  : Function<dim>(n_components, initial_time)
{ 
  fp = new fparser::FunctionParser[n_components];
}

template <int dim>
FunctionParser<dim>::~FunctionParser() {
  delete fp;
}



template <int dim>
void FunctionParser<dim>::initialize(const std::string variables,
				     const std::vector<std::string> expressions,
				     std::map<std::string, double> constants,
				     bool time_dependent,
				     bool use_degrees)
{
  // We distinguish between two cases: Time dependent problems, and
  // not time dependent problems. In the first case the number of
  // variables is given by the number of components plus one. In the
  // other case, the number of variables is equal to the number of
  // components. If none of this is the case, then an exception is
  // thrown.
  if (time_dependent) 
    n_vars = this->n_components+1;
  else 
    n_vars = this->n_components;

  AssertThrow(this->n_components == expressions.size(),
	      ExcInvalidExpressionSize(this->n_components,
				       expressions.size()) );

  for (unsigned int i=0; i<this->n_components; ++i) {
    // Add the variuous constants to the parser.
    std::map< std::string, double >::iterator
      constant = constants.begin();
    std::map< std::string, double >::const_iterator
      endc  = constants.end();
    for(;constant != endc; ++constant) {
      AssertThrow( fp[i].AddConstant(constant->first, constant->second),
	      ExcMessage("Invalid Constant Name"));
    }


    int ret_value = fp[i].Parse(expressions[i], 
				variables, 
				use_degrees);
    AssertThrow(ret_value == -1, 
		ExcParseError(ret_value, fp[i].ErrorMsg()));
    
    // The fact that the parser did not throw an error does not mean
    // that everything went ok... we can still have problems with the
    // number of variables...

  } 
  
  // Let's check if the number of variables is correct...
  AssertThrow(n_vars == fp[0].NVars(), 
	      ExcDimensionMismatch(n_vars,fp[0].NVars()));
    
  // Now set the initialization bit.
  initialized = true;
}
template <int dim>
void FunctionParser<dim>::initialize(const std::string variables,
				     const std::string expression,
				     std::map<std::string, double> constants,
				     bool time_dependent,
				     bool use_degrees)
{
  std::vector<std::string> expressions(1,expression);
  initialize(variables, expressions, constants, time_dependent, use_degrees); 
}

template <int dim>
inline
double FunctionParser<dim>::value (const Point<dim> &p,
				   unsigned int component) const 
{
  AssertThrow(initialized==true, ExcNotInitialized());
  AssertThrow(component < this->n_components,
	      ExcIndexRange(component, 0, this->n_components));
  // Statically allocates dim+1 double variables.
  double vars[dim+1];
  
  for(unsigned int i=0; i<dim; ++i)
    vars[i] = p(i);
  
  // We need the time variable only if the number of variables is
  // different from the dimension. In this case it can only be dim+1,
  // otherwise an exception would have already been thrown
  if(dim != n_vars)
    vars[dim] = this->get_time();
  
  double my_value = fp[component].Eval(vars);
  
  AssertThrow(fp[component].EvalError() == 0,
	      ExcMessage(fp[component].ErrorMsg()));
  
  return my_value;
}


template <int dim>
inline
void FunctionParser<dim>::vector_value (const Point<dim> &p,
					   Vector<double>   &values) const 
{
  AssertThrow(initialized==true, ExcNotInitialized());
  AssertThrow(values.size() == this->n_components,
	      ExcDimensionMismatch (values.size(), this->n_components));
  // Statically allocates dim+1 double variables.
  double vars[dim+1];
  
  for(unsigned int i=0; i<dim; ++i)
    vars[i] = p(i);
  
  // We need the time variable only if the number of variables is
  // different from the dimension. In this case it can only be dim+1,
  // otherwise an exception would have already been thrown
  if(dim != n_vars)
    vars[dim] = this->get_time();
  
  for(unsigned int component = 0; component < this->n_components;
      ++component) {
    values(component) = fp[component].Eval(vars);
    Assert(fp[component].EvalError() == 0,
	   ExcMessage(fp[component].ErrorMsg()));
  }
}

// Explicit Instantiations.

template class FunctionParser<1>;
template class FunctionParser<2>;
template class FunctionParser<3>;
