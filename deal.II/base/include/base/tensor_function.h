//----------------------------  tensor_function.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_function.h  ---------------------------
#ifndef __deal2__tensor_function_h
#define __deal2__tensor_function_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <base/function.h>
#include <base/point.h>
#include <base/function_time.h>

#include <vector>


/**
 *  This class is a model for a tensor valued function. The interface
 *  of the class is mostly the same as that for the @ref{Function}
 *  class, with the exception that it does not support vector-valued
 *  functions with several components, but that the return type is
 *  always tensor-valued. The returned values of the evaluation of
 *  objects of this type are always whole tensors, while for the
 *  @p{Function} class, one can ask for a specific component only, or
 *  use the @p{vector_value} function, which however does not return
 *  the value, but rather writes it into the address provided by its
 *  second argument. The reason for the different behaviour of the
 *  classes is that in the case if tensor valued functions, the size
 *  of the argument is known to the compiler a priori, such that the
 *  correct amount of memory can be allocated on the stack for the
 *  return value; on the other hand, for the vector valued functions,
 *  the size is not known to the compiler, so memory has to be
 *  allocated on the heap, resulting in relatively expensive copy
 *  operations. One can therefore consider this class a specialization
 *  of the @p{Function} class for which the size is known. An
 *  additional benefit is that tensors of arbitrary rank can be
 *  returned, not only vectors, as for them the size can be determined
 *  similarly simply.
 *
 *  @author Guido Kanschat, 1999 
 */
template <int rank, int dim>
class TensorFunction : public FunctionTime,
		       public Subscriptor
{
  public:
				     /**
				      * Define typedefs for the return
				      * types of the @p{value}
				      * functions.
				      */
    typedef Tensor<rank,dim> value_type;

#ifdef DEAL_II_LOCAL_TYPEDEF_COMP_WORKAROUND
				     /**
				      * Define a typedef for the
				      * return value of the gradient
				      * function.
				      *
				      * The construct here is only
				      * used in case we hit a certain
				      * bug in Sun's Forte
				      * compiler. See the respective
				      * macro in the aclocal.m4 file
				      * for a full description of the
				      * bug, or the documentation to
				      * the quadrature class.
				      *
				      * For better readability we
				      * later typedef this so-created
				      * type to one in the enclosing
				      * class.
				      */
    template <int dim2>
    struct GradientTypeHelper
    {
	typedef Tensor<rank+1,dim> type;
    };

				     /**
				      * Typedef the kludge declared
				      * above to a type in the class
				      * in which we would like to use
				      * it.
				      *
				      * This typedef is only used if
				      * the respective bug in the
				      * compiler is encountered,
				      * otherwise the proper typedef
				      * below is used.
				      */
    typedef typename GradientTypeHelper<dim>::type gradient_type;
#else
    typedef Tensor<rank+1,dim> gradient_type;
#endif    
    
				     /**
				      * Constructor. May take an
				      * initial value for the time
				      * variable, which defaults to
				      * zero.  
				      */
    TensorFunction (const double initial_time = 0.0);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case, as
				      * classes are usually not used
				      * by their true type, but rather
				      * through pointers to this base
				      * class.  
				      */
    virtual ~TensorFunction ();
    
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual value_type value (const Point<dim> &p) const;

				     /**
				      * Set @p{values} to the point
				      * values of the function at the
				      * @p{points}.  It is assumed
				      * that @p{values} already has
				      * the right size, i.e.  the same
				      * size as the @p{points} array.  
				      */
    virtual void value_list (const typename std::vector<Point<dim> > &points,
			     typename std::vector<value_type> &values) const;

				     /**
				      * Return the gradient of the
				      * function at the given point.
				      */
    virtual gradient_type gradient (const Point<dim> &p) const;

				     /**
				      * Set @p{gradients} to the
				      * gradients of the function at
				      * the @p{points}.  It is assumed
				      * that @p{values} already has
				      * the right size, i.e.  the same
				      * size as the @p{points} array.  
				      */
    virtual void gradient_list (const typename std::vector<Point<dim> >   &points,
				typename std::vector<gradient_type> &gradients) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcPureFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException2 (ExcVectorHasWrongSize,
		    int, int,
		    << "The vector has size " << arg1 << " but should have "
		    << arg2 << " elements.");
    
};


#endif
