//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
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
				      * Names of difference formulas.
				      */
    enum DifferenceFormula
    {
	  Euler,
	  UpwindEuler,
	  FourthOrder
    };


				     /**
				      * Constructor. Provided are the
				      * function to compute
				      * derivatives of, the direction
				      * vector of the differentiation
				      * and the step size @p{h} of the
				      * difference formula.
				      */
    FunctionDerivative (const Function<dim>& f,
			const Point<dim>& direction,
			const double h = 1.e-6);
    
				     /**
				      * Constructor. Provided are the
				      * function to compute
				      * derivatives of and the
				      * direction vector of the
				      * differentiation in each
				      * quadrature point and the
				      * difference step size.
				      *
				      * This is the constructor for a
				      * variable velocity field. Most
				      * probably, a new object of
				      * @p{FunctionDerivative} has to
				      * be constructed for each set of
				      * quadrature points.
				      *
				      * The number of quadrature point
				      * must still be the same, when
				      * values are accessed.
				      */
    FunctionDerivative (const Function<dim>& f,
			const vector<Point<dim> >& direction,
			const double h = 1.e-6);
    
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
    virtual void value_list (const typename std::vector<Point<dim> > &points,
			     std::vector<double>                     &values,
			     const unsigned int                       component = 0) const;    

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL @p{std::map} type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    unsigned int memory_consumption () const;
    
    DeclException0(ExcInvalidFormula);
  private:
				     /**
				      * Function for differentiation.
				      */
    const Function<dim>& f;

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
    vector<Point<dim> > incr;
};

#endif
