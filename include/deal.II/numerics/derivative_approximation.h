// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#ifndef __deal2__derivative_approximation_h
#define __deal2__derivative_approximation_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/filtered_iterator.h>

#include <utility>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class DoFHandler;
namespace hp
{
  template <int dim, int spacedim> class DoFHandler;
}



/**
 * This namespace provides functions that compute a cell-wise approximation of the norm of a
 * derivative of a finite element field by taking difference quotients
 * between neighboring cells. This is a rather simple but efficient
 * form to get an error indicator, since it can be computed with
 * relatively little numerical effort and yet gives a reasonable
 * approximation.
 *
 * The way the difference quotients are computed on cell $K$ is the
 * following (here described for the approximation of the gradient of
 * a finite element field, but see below for higher derivatives): let
 * $K'$ be a neighboring cell, and let $y_{K'}=x_{K'}-x_K$ be the
 * distance vector between the centers of the two cells, then
 *   $ \frac{u_h(x_{K'}) - u_h(x_K)}{ \|y_{K'}\| }$
 * is an approximation of the directional derivative
 *   $ \nabla u(x_K) \cdot \frac{y_{K'}}{ \|y_{K'}\| }.$
 * By multiplying both terms by $\frac{y_{K'}}{ \|y_{K'}\| }$ from the
 * left and summing over all neighbors $K'$, we obtain
 *   $ \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|}
 *                      \frac{y_{K'}^T}{ \|y_{K'}\| } \right) \nabla u(x_K)
 *     \approx
 *     \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|}
 *                      \frac{u_h(x_{K'}) - u_h(x_K)}{ \|y_{K'}\| }  \right).$
 *
 * Thus, if the matrix
 *   $ Y =  \sum_{K'} \left( \frac{y_{K'}}{\|y_{K'}\|}
 *                           \frac{y_{K'}^T}{ \|y_{K'}\| } \right)$ is
 * regular (which is the case when the vectors $y_{K'}$ to all neighbors span
 * the whole space), we can obtain an approximation to the true gradient by
 *   $ \nabla u(x_K)
 *     \approx
 *     Y^{-1} \sum_{K'} \left( \frac{y_{K'}}{\|y_{K'}\|}
 *                             \frac{u_h(x_{K'}) - u_h(x_K)}{ \|y_{K'}\| }
 *                      \right).$
 * This is a quantity that is easily computed. The value returned for
 * each cell when calling the @p approximate_gradient function of
 * this class is the $l_2$ norm of this approximation to the
 * gradient. To make this a useful quantity, you may want to scale
 * each element by the correct power of the respective cell size.
 *
 * The computation of this quantity must fail if a cell has only
 * neighbors for which the direction vectors $y_K$ do not span the
 * whole space, since then the matrix $Y$ is no longer invertible. If
 * this happens, you will get an error similar to this one:
 * @code
 * --------------------------------------------------------
 * An error occurred in line <749> of file <source/numerics/derivative_approximation.cc> in function
 *     void DerivativeApproximation::approximate(const Mapping<dim,spacedim>&, const DH<dim,spacedim>&, const InputVector&, unsigned int, const
 *  std::pair<unsigned int, unsigned int>&, Vector<float>&) [with DerivativeDescription = DerivativeApproximation::Gradient<3>, int
 * dim = 3, DH = DoFHandler, InputVector = Vector<double>]
 * The violated condition was:
 *     determinant(Y) != 0
 * The name and call sequence of the exception was:
 *     ExcInsufficientDirections()
 * Additional Information:
 * (none)
 * --------------------------------------------------------
 * @endcode
 * As can easily be verified, this can only happen on very
 * coarse grids, when some cells and all their neighbors have not been
 * refined even once. You should therefore only call the functions of
 * this class if all cells are at least once refined. In practice this
 * is not much of a restriction.
 *
 *
 * <h3>Approximation of higher derivatives</h3>
 *
 * Similar to the reasoning above, approximations to higher
 * derivatives can be computed in a similar fashion. For example, the
 * tensor of second derivatives is approximated by the formula
 *   $ \nabla^2 u(x_K)
 *     \approx
 *     Y^{-1}
 *     \sum_{K'}
 *        \left(
 *           \frac{y_{K'}}{\|y_{K'}\|} \otimes
 *           \frac{\nabla u_h(x_{K'}) - \nabla u_h(x_K)}{ \|y_{K'}\| }
 *        \right),
 *   $
 * where $\otimes$ denotes the outer product of two vectors. Note that
 * unlike the true tensor of second derivatives, its approximation is
 * not necessarily symmetric. This is due to the fact that in the
 * derivation, it is not clear whether we shall consider as projected
 * second derivative the term $\nabla^2 u y_{KK'}$ or $y_{KK'}^T
 * \nabla^2 u$. Depending on which choice we take, we obtain one
 * approximation of the tensor of second derivatives or its
 * transpose. To avoid this ambiguity, as result we take the
 * symmetrized form, which is the mean value of the approximation and
 * its transpose.
 *
 * The returned value on each cell is the spectral norm of the
 * approximated tensor of second derivatives, i.e. the largest
 * eigenvalue by absolute value. This equals the largest curvature of
 * the finite element field at each cell, and the spectral norm is the
 * matrix norm associated to the $l_2$ vector norm.
 *
 * Even higher than the second derivative can be obtained along the
 * same lines as exposed above.
 *
 *
 * <h3>Refinement indicators based on the derivatives</h3>
 *
 * If you would like to base a refinement criterion upon these
 * approximation of the derivatives, you will have to scale the results
 * of this class by an appropriate power of the mesh width. For
 * example, since
 * $\|u-u_h\|^2_{L_2} \le C h^2 \|\nabla u\|^2_{L_2}$, it might be the
 * right thing to scale the indicators as $\eta_K = h \|\nabla u\|_K$,
 * i.e. $\eta_K = h^{1+d/2} \|\nabla u\|_{\infty;K}$, i.e. the right
 * power is $1+d/2$.
 *
 * Likewise, for the second derivative, one should choose a power of
 * the mesh size $h$ one higher than for the gradient.
 *
 *
 * <h3>Implementation</h3>
 *
 * The formulae for the computation of approximations to the gradient
 * and to the tensor of second derivatives shown above are very much
 * alike. The basic difference is that in one case the finite
 * difference quotient is a scalar, while in the other case it is a
 * vector. For higher derivatives, this would be a tensor of even
 * higher rank. We then have to form the outer product of this
 * difference quotient with the distance vector $y_{KK'}$, symmetrize
 * it, contract it with the matrix $Y^{-1}$ and compute its norm. To
 * make the implementation simpler and to allow for code reuse, all
 * these operations that are dependent on the actual order of the
 * derivatives to be approximated, as well as the computation of the
 * quantities entering the difference quotient, have been separated
 * into auxiliary nested classes (names @p Gradient and
 * @p SecondDerivative) and the main algorithm is simply passed one
 * or the other data types and asks them to perform the order
 * dependent operations. The main framework that is independent of
 * this, such as finding all active neighbors, or setting up the
 * matrix $Y$ is done in the main function @p approximate.
 *
 * Due to this way of operation, the class may be easily extended for
 * higher oder derivatives than are presently implemented. Basically,
 * only an additional class along the lines of the derivative
 * descriptor classes @p Gradient and @p SecondDerivative has to be
 * implemented, with the respective typedefs and functions replaced by
 * the appropriate analogues for the derivative that is to be
 * approximated.
 *
 * @ingroup numerics
 * @author Wolfgang Bangerth, 2000
 */
namespace DerivativeApproximation
{
  /**
   * This function is used to obtain an approximation of the gradient. Pass it
   * the DoF handler object that describes the finite element field, a nodal
   * value vector, and receive the cell-wise Euclidian norm of the
   * approximated gradient.
   *
   * The last parameter denotes the solution component, for which the gradient
   * is to be computed. It defaults to the first component. For scalar
   * elements, this is the only valid choice; for vector-valued ones, any
   * component between zero and the number of vector components can be given
   * here.
   *
   * In a parallel computation the @p solution vector needs to contain the
   * locally relevant unknowns.
   */
  template <int dim, template <int, int> class DH, class InputVector, int spacedim>
  void
  approximate_gradient (const Mapping<dim,spacedim>    &mapping,
                        const DH<dim,spacedim>         &dof,
                        const InputVector     &solution,
                        Vector<float>         &derivative_norm,
                        const unsigned int     component = 0);

  /**
   * Calls the @p interpolate function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, template <int, int> class DH, class InputVector, int spacedim>
  void
  approximate_gradient (const DH<dim,spacedim>         &dof,
                        const InputVector     &solution,
                        Vector<float>         &derivative_norm,
                        const unsigned int     component = 0);

  /**
   * This function is the analogue to the one above, computing finite
   * difference approximations of the tensor of second derivatives. Pass it
   * the DoF handler object that describes the finite element field, a nodal
   * value vector, and receive the cell-wise spectral norm of the approximated
   * tensor of second derivatives. The spectral norm is the matrix norm
   * associated to the $l_2$ vector norm.
   *
   * The last parameter denotes the solution component, for which the gradient
   * is to be computed. It defaults to the first component. For scalar
   * elements, this is the only valid choice; for vector-valued ones, any
   * component between zero and the number of vector components can be given
   * here.
   *
   * In a parallel computation the @p solution vector needs to contain the
   * locally relevant unknowns.
   */
  template <int dim, template <int, int> class DH, class InputVector, int spacedim>
  void
  approximate_second_derivative (const Mapping<dim,spacedim>    &mapping,
                                 const DH<dim,spacedim>         &dof,
                                 const InputVector     &solution,
                                 Vector<float>         &derivative_norm,
                                 const unsigned int     component = 0);

  /**
   * Calls the @p interpolate function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, template <int, int> class DH, class InputVector, int spacedim>
  void
  approximate_second_derivative (const DH<dim,spacedim>         &dof,
                                 const InputVector     &solution,
                                 Vector<float>         &derivative_norm,
                                 const unsigned int     component = 0);

  /**
   * This function calculates the <tt>order</tt>-th order approximate
   * derivative and returns the full tensor for a single cell.
   *
   * The last parameter denotes the solution component, for which the gradient
   * is to be computed. It defaults to the first component. For scalar
   * elements, this is the only valid choice; for vector-valued ones, any
   * component between zero and the number of vector components can be given
   * here.
   *
   * In a parallel computation the @p solution vector needs to contain the
   * locally relevant unknowns.
   */
  template <class DH, class InputVector, int order>
  void
  approximate_derivative_tensor (const Mapping<DH::dimension,DH::space_dimension> &mapping,
                                 const DH                                     &dof,
                                 const InputVector                            &solution,
                                 const typename DH::active_cell_iterator      &cell,
                                 Tensor<order,DH::dimension>                  &derivative,
                                 const unsigned int                            component = 0);

  /**
   * Same as above, with <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <class DH, class InputVector, int order>
  void
  approximate_derivative_tensor (const DH                                     &dof,
                                 const InputVector                            &solution,
                                 const typename DH::active_cell_iterator      &cell,
                                 Tensor<order,DH::dimension>                  &derivative,
                                 const unsigned int                            component = 0);

  /**
   * Return the norm of the derivative.
   */
  template <int dim, int order>
  double
  derivative_norm (const Tensor<order,dim> &derivative);

  /**
   * Exception
   */
  DeclException2 (ExcInvalidVectorLength,
                  int, int,
                  << "Vector has length " << arg1 << ", but should have "
                  << arg2);
  /**
   * Exception
   */
  DeclException0 (ExcInsufficientDirections);
}



DEAL_II_NAMESPACE_CLOSE

#endif
