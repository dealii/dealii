//----------------------------  auto_derivative_function.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  auto_derivative_function.h  ---------------------------
#ifndef __deal2__auto_derivative_function_h
#define __deal2__auto_derivative_function_h


#include <base/exceptions.h>
#include <base/function.h>


/**
 * This class automatically computes the gradient of a function by
 * employing numerical difference quotients. This only, if the user
 * function does not provide the gradient function himself.
 *
 * @sect3{Usage}
 * The following example of an user defined function overloads and
 * implements only the @p{value} function but not the @p{gradient}
 * function. If the @p{gradient} function is invoked then the gradient
 * function implemented by the @p{AutoDerivativeFunction} is called,
 * where the latter function imployes numerical difference quotients.
 *
 * @begin{verbatim}
 * class UserFunction: public AutoDerivativeFunction
 * {               // access to one component at one point
 *   double value (const Point<dim> &p, const
 *                 unsigned int component = 0) const
 *          { // Implementation ....  };
 * } user_function;
 *
 *            // gradient by employing difference quotients.
 * Tensor<1,dim> grad=user_function.gradient(some_point);
 * @end{verbatim}
 * 
 * If the user overloads and implements also the gradient function,
 * then, of course, the users gradient function is called.
 *
 * Note, that the usage of the @p{value} and @p{gradient} functions
 * explained above, also applies to the @p{value_list} and
 * @p{gradient_list} functions as well as to the vector valued
 * versions of these functions, see e.g. @p{vector_value},
 * @p{vector_gradient}, @p{vector_value_list} and
 * @p{vector_gradient_list}.
 *
 * The @p{gradient} and @p{gradient_list} functions make use of the
 * @p{value} function. The @p{vector_gradient} and
 * @p{vector_gradient_list} make use of the @p{vector_value}
 * function. Make sure that the user defined function implements the
 * @p{value} function and the @p{vector_value} function, respectively.
 *
 * Furthermore note, that an object of this class does @em{not} represent
 * the derivative of a function, like @ref{FunctionDerivative}, that
 * gives a directional derivate by calling the @p{value} function. In
 * fact, this class (the @p{AutoDerivativeFunction} class) can
 * substitute the @p{Function} class as base class for user defined
 * classes. This class implements the @p{gradient} functions for
 * automatic computation of numerical difference quotients and serves
 * as intermediate class between the base @p{Function} class and the
 * user defined function class.
 *
 * @author Ralf Hartmann, 2001
 */
template <int dim>
class AutoDerivativeFunction : public Function<dim>
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
				      * Constructor. Takes the
				      * difference step size
				      * @p{h}. It's within the user's
				      * responsibility to choose an
				      * appropriate value here. @p{h}
				      * should be chosen taking into
				      * account the absolute value as
				      * well as the amount of local
				      * variation of the function.
				      * Setting @p{h=1e-6} might be a
				      * good choice for functions with
				      * an absolute value of about 1,
				      * that furthermore does not vary
				      * to much.
				      *
				      * @p{h} can be changed later
				      * using the @p{set_h} function.
				      *
				      * Sets @p{DifferenceFormula}
				      * @p{formula} to the default
				      * @p{Euler} formula of the
				      * @p{set_formula}
				      * function. Change this preset
				      * formula by calling the
				      * @p{set_formula} function.
				      */
    AutoDerivativeFunction (const double h,
			    const unsigned int n_components = 1,
			    const double       initial_time = 0.0);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~AutoDerivativeFunction ();
    
				     /**
				      * Choose the difference formula.
				      *
				      * Formulas implemented right now
				      * are first order backward Euler
				      * (@p{UpwindEuler}), second
				      * order symmetric Euler
				      * (@p{Euler}) and a symmetric
				      * fourth order formula
				      * (@p{FourthOrder}).
				      */
    void set_formula (const DifferenceFormula formula = Euler);

				     /**
				      * Takes the difference step size
				      * @p{h}. It's within the user's
				      * responsibility to choose an
				      * appropriate value here. @p{h}
				      * should be chosen taking into
				      * account the absolute value of
				      * as well as the amount of local
				      * variation of the function.
				      * Setting @p{h=1e-6} might be a
				      * good choice for functions with
				      * an absolute value of about 1,
				      * that furthermore does not vary
				      * to much.
				      */
    void set_h (const double h);
    
				     /**
				      * Return the gradient of the
				      * specified component of the
				      * function at the given point.
				      *
				      * Imployes numerical difference
				      * quotients using the preset
				      * @p{DifferenceFormula}
				      * @p{formula}.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Return the gradient of all
				      * components of the
				      * function at the given point.
				      *
				      * Imployes numerical difference
				      * quotients using the preset
				      * @p{DifferenceFormula}
				      * @p{formula}.
				      */
    virtual void vector_gradient (const Point<dim>            &p,
				  typename std::vector<Tensor<1,dim> > &gradients) const;
    
				     /**
				      * Set @p{gradients} to the
				      * gradients of the specified
				      * component of the function at
				      * the @p{points}.  It is assumed
				      * that @p{gradients} already has the
				      * right size, i.e.  the same
				      * size as the @p{points} array.
				      *
				      * Imployes numerical difference
				      * quotients using the preset
				      * @p{DifferenceFormula}
				      * @p{formula}.
				      */
    virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				typename std::vector<Tensor<1,dim> >    &gradients,
				const unsigned int              component = 0) const;
    
				     /**
				      * Set @p{gradients} to the gradients of
				      * the function at the @p{points},
				      * for all components.
				      * It is assumed that @p{gradients} 
				      * already has the right size, i.e.
				      * the same size as the @p{points} array.
				      *
				      * The outer loop over
				      * @p{gradients} is over the points
				      * in the list, the inner loop
				      * over the different components
				      * of the function.
				      *
				      * Imployes numerical difference
				      * quotients using the preset
				      * @p{DifferenceFormula}
				      * @p{formula}.
				      */
    virtual void vector_gradient_list (const typename std::vector<Point<dim> > &points,
				       typename std::vector<typename std::vector<Tensor<1,dim> > > &gradients) const;

				     /**
				      * Returns a
				      * @p{DifferenceFormula} of the
				      * order @p{ord} at minimum.
				      */
    static DifferenceFormula get_formula_of_order(const unsigned int ord);

				     /**
				      * Exception.
				      */
    DeclException0(ExcInvalidFormula);

  private:
    
				     /**
				      * Step size of the difference
				      * formula. Set by the @p{set_h}
				      * function.
				      */
    double h;

				     /**
				      * Includes the unit vectors
				      * scaled by @p{h}.
				      */
    typename std::vector<Tensor<1,dim> > ht;
    
				     /**
				      * Difference formula. Set by the
				      * @p{set_formula} function.
				      */
    DifferenceFormula formula;
};



#endif
