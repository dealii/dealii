//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------
#ifndef __deal2__function_derivative_h
#define __deal2__function_derivative_h

#include <base/exceptions.h>
#include <base/function.h>

/**
 * Names of difference formulas.
 */
enum DifferenceFormula
{
  Euler,
  UpwindEuler,
  FourthOrder
};


/**
 * Derivative of a function object.  The value access functions of
 * this class return the directional derivative of a function with
 * respect to a direction provided on construction. If @p{b} is the
 * vector, the derivative @p{b . grad f} is computed. This derivative
 * is evaluated directly, not by computing the gradient of @p{f} and
 * its scalar product with @p{b}.
 *
 * The derivative is computed numerically, using one of the provided
 * difference formulas (see @p{set_formula} for available
 * schemes). Experimenting with @p{h} and the difference scheme may be
 * necessary to obtain sufficient results.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FunctionDerivative : public Function<dim>
{
public:
				   /**
				    * Constructor. Provided are the
				    * function to compute derivatives
				    * of and the direction vector of
				    * the differentiation.
				    */
  FunctionDerivative (const Function<dim>& f,
		      const Point<dim>& direction);

				   /**
				    * Set step size of the difference
				    * formula. This is set to the
				    * default in the constructor.
				    */
  void set_h (double h_new = 1.e-6);
  
				   /**
				    * Choose the difference formula.
				    * This is set to the default in
				    * the constructor.
				    *
				    * Formulas implemented right now
				    * are first order backward Euler
				    * (@p{UpwindEuler}), second order
				    * symmetric Euler (@p{Euler}) and
				    * a symmetric fourth order formula
				    * (@p{FourthOrder}).
				    */
  void set_formula (DifferenceFormula formula = Euler);
  
				   /**
				    * Function value at one point.
				    */
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;
  
				   /**
				    * Function values at multiple points.
				    */
  virtual void value_list (const vector<Point<dim> > &points,
			   vector<double>            &values,
			   const unsigned int         component = 0) const;


  DeclException0(ExcInvalidFormula);
  
private:
				   /**
				    * Function for differentiation.
				    */
  const Function<dim>& f;
				   /**
				    * Differentiation vector.
				    */
  Point<dim> direction;
  
				   /**
				    * Step size of the difference formula.
				    */
  double h;
				   /**
				    * Difference formula.
				    */
  DifferenceFormula formula;
				   /**
				    * Helper object. Contains the
				    * increment vector for the
				    * formula.
				    */
  Point<dim> incr;
};

#endif
