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
 * $Id: ibm_rhs.h,v 1.15 2005/04/05 14:46:49 luca Exp $
 *
 */
#ifndef PARSED_SYMMETRIC_TENSOR_FUNCTION_H
#define PARSED_SYMMETRIC_TENSOR_FUNCTION_H

#include <base/parameter_handler.h>
#include <base/function.h>
#include <base/function_parser.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/symmetric_tensor.h>

#include <lac/vector.h>

using namespace dealii;

template <int rank, int dim>
class ParsedSymmetricTensorFunction
{
 public:
  /** Constructor. It initializes the function parser.*/
  ParsedSymmetricTensorFunction();
  
  /** Declare parameters needed by this class. */
  static void declare_parameters(ParameterHandler &prm);

  /** Parse parameters needed by this class. */
  void parse_parameters(ParameterHandler &prm);
  
  /** Get ONE value at the given point. Evaluate the parsed function
      at the given point and return the symmetric tensor. */
  const SymmetricTensor<rank, dim> & operator() (const Point<dim> &p) const;
  
  /** Get a vector of the values at a number of specified points */
  void value_list(const std::vector< Point<dim> > &points,
		  std::vector< SymmetricTensor<rank, dim> > &values);
    
  /** Get the tensor currently stored. */
  inline const SymmetricTensor<rank, dim> & operator() () const {
    return t;
  };
  
  /** Get the size of the equivalent triangular matrix. If rank = 2,
      this number is equal to dim, otherwise it is equal to dim*(dim+1)/2. */
  static unsigned int get_dim_triangular();
  
  /** Set internal time of the function. */
  void set_time(double);
 private:
  FunctionParser<dim> f;
  mutable SymmetricTensor<rank, dim> t;
};


template<int dim>
inline 
Tensor<2,dim> operator*(const SymmetricTensor<4,dim> & C, 
			const Tensor<2,dim> &gradv) {
  Tensor<2,dim> T;
  for(unsigned int i=0; i<dim; ++i)
    for(unsigned int j=0; j<dim; ++j)
      for(unsigned int k=0; k<dim; ++k)
	for(unsigned int l=0; l<dim; ++l)
	  T[i][j] += C[i][j][k][l] * gradv[k][l];
  return T;
}


template<int dim>
inline 
double operator*(const Tensor<2,dim> &gradv,
		 const SymmetricTensor<2,dim> & C) {
  double T = 0;
  for(unsigned int i=0; i<dim; ++i)
    for(unsigned int j=0; j<dim; ++j)
      T += C[i][j] * gradv[i][j];
  return T;
}


template<int dim>
inline 
double double_contract(const SymmetricTensor<2,dim> & C,
		       const Tensor<2,dim> &gradv) {
  double T = 0;
  for(unsigned int i=0; i<dim; ++i)
    for(unsigned int j=0; j<dim; ++j)
      T += C[i][j] * gradv[i][j];
  return T;
}



template<int dim>
inline 
Tensor<1,dim> operator*(const SymmetricTensor<2,dim> & C,
			const Tensor<1,dim> &gradt) {
  Tensor<1,dim> T;
  for(unsigned int j=0; j<dim; ++j)
    for(unsigned int i=0; i<dim; ++i)
      T[j] += C[j][i] * gradt[i];
  return T;
}


/* template<int dim> */
/* inline  */
/* double double_contract(const Tensor<2,dim> &gradv, */
/* 		       const Tensor<2,dim> & C) { */
/*   double T = 0; */
/*   for(unsigned int i=0; i<dim; ++i) */
/*     for(unsigned int j=0; j<dim; ++j) */
/*       T += C[i][j] * gradv[i][j]; */
/*   return T; */
/* } */



/* template<int dim> */
/* inline  */
/* Tensor<2,dim> operator*(const SymmetricTensor<4,dim> & C,  */
/* 			const Tensor<2,dim> &gradv) { */
/*   Tensor<2,dim> T; */
/*   for(unsigned int i=0; i<dim; ++i) */
/*     for(unsigned int j=0; j<dim; ++j) */
/*       for(unsigned int k=0; k<dim; ++k) */
/* 	for(unsigned int l=0; l<dim; ++l) */
/* 	  T[i][j] += C[i][j][k][l] * gradv[k][l]; */
/*   return T; */
/* } */


/* template<int dim> */
/* inline  */
/* double operator*(const Tensor<2,dim> &gradv, */
/* 		 const SymmetricTensor<2,dim> & C) { */
/*   double T = 0; */
/*   for(unsigned int i=0; i<dim; ++i) */
/*     for(unsigned int j=0; j<dim; ++j) */
/*       T += C[i][j] * gradv[i][j]; */
/*   return T; */
/* } */


/* template<int dim> */
/* inline  */
/* double double_contract(const SymmetricTensor<2,dim> & C, */
/* 		       const Tensor<2,dim> &gradv) { */
/*   double T = 0; */
/*   for(unsigned int i=0; i<dim; ++i) */
/*     for(unsigned int j=0; j<dim; ++j) */
/*       T += C[i][j] * gradv[i][j]; */
/*   return T; */
/* } */



/* template<int dim> */
/* inline  */
/* Tensor<1,dim> operator*(const SymmetricTensor<2,dim> & C, */
/* 			const Tensor<1,dim> &gradt) { */
/*   Tensor<1,dim> T; */
/*   for(unsigned int j=0; j<dim; ++j) */
/*     for(unsigned int i=0; i<dim; ++i) */
/*       T[j] += C[j][i] * gradt[i]; */
/*   return T; */
/* } */


template<int dim>
inline 
double double_contract(const Tensor<2,dim> &gradv,
		       const Tensor<2,dim> & C) {
  double T = 0;
  for(unsigned int i=0; i<dim; ++i)
    for(unsigned int j=0; j<dim; ++j)
      T += C[i][j] * gradv[i][j];
  return T;
}


#endif
