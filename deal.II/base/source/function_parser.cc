//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------

#include <base/function_parser.h>
#include <base/point.h>
#include <lac/vector.h>

// Utility to split a string given a delimiter.

unsigned int SplitString(const std::string& input, 
			 const std::string& delimiter, 
			 std::vector<std::string>& results)
{
  // Size of the delimiter. Use it to calculate the offsets between
  // different pieces of the string
  unsigned int sizeS2 = delimiter.size();
  // Counter for Current Position
  unsigned int iPos = 0;
  // Counter for Position of next delimiter
  std::string::size_type newPos = 0;  
  // Counter for Number of strings
  unsigned int numFound = 0;
  
  while( newPos != std::string::npos )
    { 
      // Look for a delimiter starting with the correct offset. Note
      // that if it is the beginning of the cycle the offset is zero.
      newPos = input.find (delimiter, iPos);
      // We either found a delimiter inside the string or the
      // delimiter was not found. Increment the found counter anyway,
      // as if there are no delimiters, we return one vector
      // containing exacly one string anyway
      numFound++;
      // Add the string we just found to the output vector
      Assert(iPos <= input.size(), ExcInternalError());
      results.push_back(input.substr(iPos, newPos-iPos));
      // Update the current position with the correct offset.
      iPos = newPos+sizeS2;
    }
  // At least we found one string and no delimiters.
  Assert(numFound > 0, ExcInternalError());
  return numFound;
}


template <int dim>
FunctionParser<dim>::FunctionParser(const unsigned int n_components,
				    const double       initial_time)
                :
                Function<dim>(n_components, initial_time),
                fp (0)
{ 
  fp = new fparser::FunctionParser[n_components];
}

template <int dim>
FunctionParser<dim>::~FunctionParser() {
  delete[] fp;
}

#ifndef DEAL_II_DISABLE_PARSER

template <int dim>
void FunctionParser<dim>::initialize(const std::string variables,
				     const std::vector<std::string> expressions,
				     std::map<std::string, double> constants,
				     bool time_dependent,
				     bool use_degrees)
{
  // We check that the number of components of this function matches
  // the number of components passed in as a vector of strings.
  AssertThrow(this->n_components == expressions.size(),
	      ExcInvalidExpressionSize(this->n_components,
				       expressions.size()) );
  
  for (unsigned int i=0; i<this->n_components; ++i) {
    // Add the variuous constants to the parser.
    std::map< std::string, double >::const_iterator
      constant = constants.begin(),
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
    
  // Now we define how many variables we expect to read in.  We
  // distinguish between two cases: Time dependent problems, and not
  // time dependent problems. In the first case the number of
  // variables is given by the dimension plus one. In the other case,
  // the number of variables is equal to the dimension. Once we parsed
  // the variables string, if none of this is the case, then an
  // exception is thrown.
  if (time_dependent) 
    n_vars = dim+1;
  else 
    n_vars = dim;
    
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
  // Start up with an empty vector of expressions.
  std::vector<std::string> expressions;
  
  // In this case a vector function is built with the number of
  // components found between the separators ';' 
  SplitString(expression, ";", expressions);
  // Now initialize with the things we got.
  initialize
    (variables, expressions, constants, time_dependent, use_degrees);  
}

template <int dim>
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

#else


template <int dim>
void FunctionParser<dim>::initialize(
  const std::string, const std::vector<std::string>,
  std::map<std::string, double>, bool, bool)
{
  Assert(false, ExcDisabled("parser"));
}


template <int dim>
void FunctionParser<dim>::initialize(
  const std::string, const std::string,
  std::map<std::string, double>, bool, bool)
{
  Assert(false, ExcDisabled("parser"));
}


template <int dim>
double FunctionParser<dim>::value (
  const Point<dim> &, unsigned int) const 
{
  Assert(false, ExcDisabled("parser"));
  return 0.;
}


template <int dim>
void FunctionParser<dim>::vector_value (
  const Point<dim>&, Vector<double>&) const 
{
  Assert(false, ExcDisabled("parser"));
}


#endif

// Explicit Instantiations.

template class FunctionParser<1>;
template class FunctionParser<2>;
template class FunctionParser<3>;
