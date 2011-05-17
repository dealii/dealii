/*
 * Immersed Boundary Problem:
 * 
 * Header Files
 *
 * Author: 
 * Luca Heltai <heltai@dimat.unipv.it>
 * =============
 * License: GPL.
 * =============
 * $Id: ibm_rhs.cc,v 1.34 2005/04/05 14:46:49 luca Exp $
 *
 */
#include "../include/parsed_symmetric_tensor_function.h"
#include <deal.II/base/utilities.h>

template <int rank, int dim>
  void ParsedSymmetricTensorFunction<rank, dim>::value_list(const std::vector< Point<dim> > &points,
							    std::vector< SymmetricTensor<rank, dim> > &values)
{

  for (unsigned int n=0; n<points.size(); ++n) 
    values[n] = this->operator()(points[n]);
      
}

template <int rank, int dim>
  unsigned int ParsedSymmetricTensorFunction<rank, dim>::get_dim_triangular() 
{
  if (rank == 2) 
    return dim;
  else if (rank == 4)
    return (dim*(dim+1)/2);
  else {
    AssertThrow(false, ExcInternalError());
    return 0;
  }
}


template <int rank, int dim>
void ParsedSymmetricTensorFunction<rank, dim>::set_time (double t) {
  f.set_time(t);
}


template <int rank, int dim>
ParsedSymmetricTensorFunction<rank, dim>::ParsedSymmetricTensorFunction () :
  // Number of vector functions.
  f(get_dim_triangular()*(get_dim_triangular()+1)/2)
{
  Assert(rank == 2 || rank == 4,
	 ExcMessage("Rank has to be even, either 2 or 4."));
}

template <int rank, int dim>
void ParsedSymmetricTensorFunction<rank, dim>::declare_parameters(ParameterHandler &prm) 
{
  std::string vnames;
  switch (dim) {
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
		    "The name of the variables as they will be used in the function, separated by ','.");
  // The expressions of the function
  std::vector<std::string> expr(get_dim_triangular());
  for(unsigned int i=0; i<get_dim_triangular(); ++i) {
    for(unsigned int j=i; j<get_dim_triangular(); ++j)
      expr[i] += (i==j) ? "1" : "; 0";
    char tmp[100];
    sprintf(tmp, "Row %d", i+1);
    // Now we have the function expressions
    prm.declare_entry(tmp, expr[i], Patterns::Anything(), (i == 0)  ?
		      "Separate different components expressions by ';' "
		      "as ',' is used internally by the function parser." 
		      : "");
  }
  prm.declare_entry("Function constants", "", Patterns::Anything(),
		    "Any constant used inside the functions which is not a variable name.");
}

template <int rank, int dim>
  void ParsedSymmetricTensorFunction<rank, dim>::parse_parameters(ParameterHandler &prm) 
{
  std::string vnames = prm.get("Variable names");
  std::string expr;
  for(unsigned int i=0; i < get_dim_triangular(); ++i) {
    char tmp[100];
    sprintf(tmp, "Row %d", i+1);
    expr += prm.get(tmp);
    if(i+1 < get_dim_triangular())
      expr += "; "; 
  }
  std::string constants_list = prm.get("Function constants");
  
  std::vector<std::string> const_list = 
    Utilities::split_string_list(constants_list, ',');
  std::map<std::string, double> constants;
  for(unsigned int i = 0; i < const_list.size(); ++i) {
    std::vector<std::string> this_c = 
      Utilities::split_string_list(const_list[i], '=');
    AssertThrow(this_c.size() == 2, ExcMessage("Invalid format"));
    double tmp;
    AssertThrow( sscanf(this_c[1].c_str(), "%lf", &tmp), ExcMessage("Double number?"));
    constants[this_c[0]] = tmp;
  }
  
  constants["pi"] = M_PI;
  constants["Pi"] = M_PI;
  
  unsigned int nn = (Utilities::split_string_list(vnames)).size();
  switch (nn) {
  case dim:
    // Time independent function
    f.initialize(vnames, expr, constants); 
    break;
  case dim+1:
    // Time dependent function
    f.initialize(vnames, expr, constants, true);
    break;
  default:
    AssertThrow(false, ExcMessage("Not the correct size. Check your code."));
  }
}
