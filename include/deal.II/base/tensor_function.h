// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__tensor_function_h
#define __deal2__tensor_function_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function_time.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 *  This class is a model for a tensor valued function. The interface
 *  of the class is mostly the same as that for the Function
 *  class, with the exception that it does not support vector-valued
 *  functions with several components, but that the return type is
 *  always tensor-valued. The returned values of the evaluation of
 *  objects of this type are always whole tensors, while for the
 *  <tt>Function</tt> class, one can ask for a specific component only, or
 *  use the <tt>vector_value</tt> function, which however does not return
 *  the value, but rather writes it into the address provided by its
 *  second argument. The reason for the different behaviour of the
 *  classes is that in the case of tensor valued functions, the size
 *  of the argument is known to the compiler a priori, such that the
 *  correct amount of memory can be allocated on the stack for the
 *  return value; on the other hand, for the vector valued functions,
 *  the size is not known to the compiler, so memory has to be
 *  allocated on the heap, resulting in relatively expensive copy
 *  operations. One can therefore consider this class a specialization
 *  of the <tt>Function</tt> class for which the size is known. An
 *  additional benefit is that tensors of arbitrary rank can be
 *  returned, not only vectors, as for them the size can be determined
 *  similarly simply.
 *
 *  @ingroup functions
 *  @author Guido Kanschat, 1999
 */
template <int rank, int dim, typename Number=double>
class TensorFunction : public FunctionTime<Number>,
  public Subscriptor
{
public:
  /**
   * Define typedefs for the return types of the <tt>value</tt> functions.
   */
  typedef Tensor<rank,dim,Number> value_type;

  typedef Tensor<rank+1,dim,Number> gradient_type;

  /**
   * Constructor. May take an initial value for the time variable, which
   * defaults to zero.
   */
  TensorFunction (const Number initial_time = Number(0.0));

  /**
   * Virtual destructor; absolutely necessary in this case, as classes are
   * usually not used by their true type, but rather through pointers to
   * this base class.
   */
  virtual ~TensorFunction ();

  /**
   * Return the value of the function at the given point.
   */
  virtual value_type value (const Point<dim, Number> &p) const;

  /**
   * Set <tt>values</tt> to the point values of the function at the
   * <tt>points</tt>.  It is assumed that <tt>values</tt> already has the
   * right size, i.e.  the same size as the <tt>points</tt> array.
   */
  virtual void value_list (const std::vector<Point<dim, Number> > &points,
                           std::vector<value_type> &values) const;

  /**
   * Return the gradient of the function at the given point.
   */
  virtual gradient_type gradient (const Point<dim, Number> &p) const;

  /**
   * Set <tt>gradients</tt> to the gradients of the function at the
   * <tt>points</tt>.  It is assumed that <tt>values</tt> already has the
   * right size, i.e.  the same size as the <tt>points</tt> array.
   */
  virtual void gradient_list (const std::vector<Point<dim, Number> >   &points,
                              std::vector<gradient_type> &gradients) const;

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



/**
 * Provide a tensor valued function which always returns a constant tensor
 * value. Obviously, all derivates of this function are zero.
 *
 * @ingroup functions
 * @author Matthias Maier, 2013
 */
template <int rank, int dim, typename Number=double>
class ConstantTensorFunction : public TensorFunction<rank, dim, Number>
{
public:
  /**
   * Constructor; takes the constant tensor value as an argument.
   * The reference value is copied internally.
   *
   * An initial value for the time variable may be specified, otherwise it
   * defaults to zero.
   */
  ConstantTensorFunction (const dealii::Tensor<rank, dim, Number> &value,
                          const Number initial_time = 0.0);

  virtual ~ConstantTensorFunction ();

  virtual typename dealii::TensorFunction<rank, dim, Number>::value_type value (const Point<dim, Number> &p) const;

  virtual void value_list (const std::vector<Point<dim, Number> > &points,
                           std::vector<typename dealii::TensorFunction<rank, dim, Number>::value_type> &values) const;

  virtual typename dealii::TensorFunction<rank, dim, Number>::gradient_type gradient (const Point<dim, Number> &p) const;

  virtual void gradient_list (const std::vector<Point<dim, Number> > &points,
                              std::vector<typename dealii::TensorFunction<rank, dim, Number>::gradient_type> &gradients) const;

private:
  const dealii::Tensor<rank, dim, Number> _value;
};



/**
 * Provide a tensor valued function which always returns zero.
 * Obviously, all derivates of this function are zero.
 *
 * @ingroup functions
 * @author Matthias Maier, 2013
 */
template <int rank, int dim, typename Number=double>
class ZeroTensorFunction : public ConstantTensorFunction<rank, dim, Number>
{
public:
  /**
   * Constructor.
   *
   * An initial value for the time variable may be specified, otherwise it
   * defaults to zero.
   */
  ZeroTensorFunction (const Number initial_time = 0.0);
};


DEAL_II_NAMESPACE_CLOSE

#endif
