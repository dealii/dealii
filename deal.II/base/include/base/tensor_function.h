/*----------------------------   tensorfunction.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright G. Kanschat, University of Heidelberg, 1999 */
#ifndef __tensorfunction_H
#define __tensorfunction_H
/*----------------------------   tensorfunction.h     ---------------------------*/


#include <base/exceptions.h>
#include <vector>
#include <base/point.h>
#include <base/functiontime.h>
#include <base/tensorindex.h>


template<int rank_, int dim>
class Tensor;

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
template <int rank_, int dim>
class TensorFunction :
  public FunctionTime
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
    virtual Tensor<rank_, dim> operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<Tensor<rank_,dim> > &values) const;

//				   /**
//				    * Return one component of the value.
//				    */
//  virtual double operator() (TensorIndex<rank_> i, const Point<dim>& p) const;
  

				     /**
				      * Return the gradient of the function
				      * at the given point.
				      */
    virtual Tensor<rank_+1,dim> gradient (const Point<dim> &p) const;

				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#.
				      * It is assumed that #values# 
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<rank_+1,dim> > &gradients) const;
    
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






/*---------------------------- tensorfunction.h ---------------------------*/
/* end of #ifndef __tensorfunction_H */
#endif
/*---------------------------- tensorfunction.h ---------------------------*/
