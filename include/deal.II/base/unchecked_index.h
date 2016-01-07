// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#ifndef dealii__unchecked_index_h
#define dealii__unchecked_index_h


#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN


/**
 * A class that is used in place of regular indices in places where you
 * want to avoid index checks even in debug mode because you are sure
 * that the indices are correct. This class is supported in indexing
 * operators, i.e., <code>operator[]</code>, by the Tensor and
 * SymmetricTensor classes.
 *
 * In general, it is a bad idea to elide index checks because index
 * checks are among the most frequent bugs in codes. Furthermore, the
 * index checks deal.II performs are only enabled in debug, but not
 * in release mode, and consequently won't slow down large computations
 * after you have already debugged your program and switched to
 * release mode. However, there are cases in the innermost loops where
 * you are certain that indices are correct and where index checks
 * would be performed so many times that they slow down programs in
 * debug mode to an undue degree. An example would be if you
 * implemented the assembly loop of a Laplace solver in the following
 * way (a small variation of the loop in step-4, where the innermost
 * loop is in fact hidden behind the <code>operator*</code> of the
 * Tensor class):
 * @code
 *  for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
 *    for (unsigned int i=0; i<dofs_per_cell; ++i)
 *      for (unsigned int j=0; j<dofs_per_cell; ++j)
 *        {
 *          const Tensor<1,dim> grad_phi_i = fe_values.shape_grad (i, q_index);
 *          const Tensor<1,dim> grad_phi_j = fe_values.shape_grad (j, q_index);
 *
 *          double grad_i_times_grad_j = 0;
 *          for (unsigned int d=0; d<dim; ++d)
 *            grad_i_times_grad_j += grad_phi_i[d] * grad_phi_j[d];  // ***
 *
 *          cell_matrix(i,j) += (grad_i_times_grad_j *
 *                               fe_values.JxW (q_index));
 *      }
 * @endcode
 * Here, we compute the dot product $\nabla\varphi_i \cdot \nabla\varphi_j$
 * using a loop over the <code>dim</code> elements of the tensors, both
 * of which are of type Tensor<1,dim> and thus have exactly <code>dim</code>
 * elements. In other words, we are not accessing <i>any</i> tensor element
 * in the marked line, but elements whose indices we very carefully control
 * and whose values are well known from the loop in the line immediately
 * above. Consequently, there is no risk that the indices could be invalid,
 * and the repeated index checks performed by <code>operator[]</code> (of which
 * there are
 * <code>n_q_points * dofs_per_cell * dofs_per_cell * dim</code>) are not
 * necessary after the code has been debugged; given their frequency, they
 * do however slow down the overall assembly loop in significant ways.
 *
 * This class avoids the dilemma by offering the following alternative syntax:
 * @code
 *          double grad_i_times_grad_j = 0;
 *          for (unsigned int d=0; d<dim; ++d)
 *            grad_i_times_grad_j += grad_phi_i[unchecked_index(d)] *
 *                                   grad_phi_j[unchecked_index(d)];
 *      }
 * @endcode
 * Here, the function unchecked_index() returns an object of type
 * UncheckedIndex, and there is an overload of Tensor::operator[]
 * that accepts such objects, returning the requested element without
 * performing the index check.
 *
 * Obviously, such unchecked index accesses are almost never a good idea
 * and should only be used in very closely controlled contexts. In
 * practice, this means that they should only be used in a very small
 * number of primitives inside the library itself.
 *
 *
 * @author Wolfgang Bangerth, 2015
 */
template <typename IndexType>
class UncheckedIndex
{
public:
  /**
   * Construct an UncheckedIndex object that wraps the given index.
   */
  UncheckedIndex (const IndexType index);

  /**
   * Return the wrapped index set by the constructor.
   */
  IndexType value() const;

private:
  /**
   * The value of the wrapped index.
   */
  const IndexType index;
};



/**
 * Create an UncheckedIndex object with the same index type as the
 * argument given to this function.
 *
 * @param index The desired index.
 * @return An UncheckedIndex object that simply wraps @p index.
 */
template <typename IndexType>
UncheckedIndex<IndexType>
unchecked_index (const IndexType index)
{
  return UncheckedIndex<IndexType>(index);
}



// ----------- inline functions ---------------

template <typename IndexType>
inline
UncheckedIndex<IndexType>::UncheckedIndex (const IndexType index)
  :
  index (index)
{}



template <typename IndexType>
inline
IndexType
UncheckedIndex<IndexType>::value () const
{
  return index;
}




DEAL_II_NAMESPACE_CLOSE

#endif
