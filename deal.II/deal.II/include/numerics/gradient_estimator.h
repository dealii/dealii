//----------------------------  gradient_estimator.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  gradient_estimator.h  ---------------------------
#ifndef __deal2__gradient_estimator_h
#define __deal2__gradient_estimator_h


#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>
#include <base/exceptions.h>

#include <utility>



/**
 * This class computes a cell-wise estimate of the gradient by taking
 * difference quotients between neighboring cells. This is a rather
 * simple but efficient form to get an error indicator, since it can
 * be computed with relatively little numerical effort and yet gives a
 * reasonable approximation.
 *
 * The way the difference quotients are computed on cell $K$ is the
 * following: let $K'$ be a neighboring cell, and let
 * $y_{K'}=x_{K'}-x_K$ be the distance vector between the centers of
 * the two cells, then
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
 *   $ Y =  \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|} 
 *                           \frac{y_{K'}^T}{ \|y_{K'}\| } \right)$ is
 * regular (which is the case when the vectors $y_{K'}$ to all neighbors span
 * the whole space), we can obtain an approximation to the true gradient by
 *   $ \nabla u(x_K)
 *     \approx
 *     Y^{-1} \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|} 
 *                             \frac{u_h(x_{K'}) - u_h(x_K)}{ \|y_{K'}\| }  \right).$
 * This is a quantity that is easily computed. The value returned for
 * each cell when calling the main function of this class is the $l_2$
 * norm of this approximation to the gradient. To make this a useful
 * quantity, you may want to scale each element by the correct power
 * of the respective cell size.
 *
 * The computation of this quantity must fail if a cell has only
 * neighbors for which the direction vectors do not span the whole
 * space. As can easily be verified, this can only happen on very
 * coarse grids, when some cells and all their neighbors have not been
 * refined even once. You should therefore only call the functions of
 * this class if all cells are at least once refined. In practice this
 * is not much of a restriction. If for some cells, the neighbors do
 * not span the whole space, an exception is thrown.
 *
 * Note that for the computation of the quantities of this class, only
 * the values of the finite element field at the centers of the cells
 * are taken. It might therefore only be useful to use this class for
 * discontinuous, piecewise constant elements (i.e. using the
 * #FEDG_Q0# class), since all other finite elements can approximate
 * gradients themselves.
 *
 *
 * \section{Refinement indicators based on the gradients}
 *
 * If you would like to base a refinement criterion upon this
 * approximation of the gradient, you will have to scale the results
 * of this class by an appropriate power of the mesh width. For
 * example, since
 * $\|u-u_h\|^2_{L_2} \le C h^2 \|\nabla u\|^2_{L_2}$, it might be the
 * right thing to scale the indicators as $\eta_K = h \|\nabla u\|_K$,
 * i.e. $\eta_K = h^{1+d/2} \|\nabla u\|_{\infty;K}$, i.e. the right
 * power is $1+d/2$.
 *
 * @author Wolfgang Bangerth, 2000
 */
class GradientEstimator 
{
  public:
				     /**
				      * This is the main function that
				      * does what is announced in the
				      * general documentation of this
				      * class. Pass it the DoF handler
				      * object that describes the
				      * finite element field, a nodal
				      * value vector, and receive the
				      * cell-wise norm of the
				      * approximated gradient.
				      */
    template <int dim>
    static void estimate (const DoFHandler<dim> &dof,
			  const Vector<double>  &solution,
			  Vector<float>         &error_per_cell);

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

  private:
				     /**
				      * Convenience typedef denoting
				      * the range of indices on which
				      * a certain thread shall
				      * operate.
				      */
    typedef pair<unsigned int,unsigned int> IndexInterval;

				     /**
				      * Compute the error estimator on
				      * the cells in the range given
				      * by the third parameter.
				      */
    template <int dim>
    static void estimate_threaded (const DoFHandler<dim> &dof,
				   const Vector<double>  &solution,
				   const IndexInterval   &index_interval,
				   Vector<float>         &error_per_cell);    
};


#endif


