//----------------------------  tensor_function.h  ---------------------------
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
//----------------------------  tensor_function.h  ---------------------------
#ifndef __deal2__tensor_function_h
#define __deal2__tensor_function_h


#include <base/exceptions.h>
#include <vector>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <base/function.h>
#include <base/point.h>
#include <base/function_time.h>
#include <base/forward_declarations.h>


/**
 *  This class is a model for a tensor valued function.
 *  It returns the value
 *  at a given point through the #operator ()# member functions,
 *  which are virtual. It also has a function to return a whole list of
 *  values at different points to reduce the overhead of the virtual function
 *  calls; this function is preset to successively call the function returning
 *  one value at a time.
 *
 *  There are other functions return the gradient of the function at one or
 *  several points. You only have to overload those functions you need; the
 *  functions returning several values at a time will call those returning
 *  only one value, while those ones will throw an exception when called but
 *  not overloaded.
 *
 *  Usually, efficiency of your program increases if you overload the
 *  complex virtual functions, too.
 *
 *  @author Guido Kanschat, 1999
 */
template <int rank, int dim>
class TensorFunction : public FunctionTime,
		       public Subscriptor
{
  public:
				     /**
				      * Constructor. May take an initial vakue
				      * for the time variable, which defaults
				      * to zero.
				      */
    TensorFunction (const double initial_time = 0.0);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~TensorFunction ();
    
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual Tensor<rank, dim> value (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<Tensor<rank,dim> > &values) const;

				     /**
				      * Return the gradient of the
				      * function at the given point.
				      */
    virtual Tensor<rank+1,dim> gradient (const Point<dim> &p) const;

				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#.
				      * It is assumed that #values# 
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<rank+1,dim> > &gradients) const;

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

