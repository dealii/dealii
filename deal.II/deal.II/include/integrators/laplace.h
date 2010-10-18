//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_laplace_h
#define __deal2__integrators_laplace_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <fe/mapping.h>
#include <fe/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
/**
 * @brief Local integrators related to the Laplacian and its DG formulations
 *
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Laplace
  {
/**
 * Laplacian in weak form, namely on the cell <i>Z</i> the matrix
 * \f[
 * \int_Z \nu \nabla u \cdot \nabla v \, dx.
 * \f]
 *
 * The FiniteElement in <tt>fe</tt> may be scalar or vector valued. In
 * the latter case, the Laplacian is applied to each component
 * separately.
 *
 * @author Guido Kanschat
 * @date 2008, 2009, 2010
 */
    template<int dim>
      void cell_matrix (
	FullMatrix<double>& M,
	const FEValuesBase<dim>& fe,
	const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();
      
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = fe.JxW(k) * factor;
	  for (unsigned i=0;i<n_dofs;++i)
	    {
	      for (unsigned j=0;j<n_dofs;++j)
		for (unsigned int d=0;d<n_components;++d)
		  M(i,j) += dx *
		    (fe.shape_grad_component(j,k,d) * fe.shape_grad_component(i,k,d));
	    }
	}
    }
/**
 * Weak boundary condition of Nitsche type for the Laplacian, namely on the face <i>F<//i> the matrix
 * @f[
 * \int_F \Bigl(\gamma u v - \partial_n u v - u \partial_n v\Bigr)\;ds.
 * @f]
 *
 * Here, &gamma; is the <tt>penalty</tt> parameter suitably computed
 * with compute_penalty().
 */
      template <int dim>
      void nitsche_matrix (
	FullMatrix<double>& M,
	const FEValuesBase<dim>& fe,
	double penalty,
	double factor = 1.)
      {
	const unsigned int n_dofs = fe.dofs_per_cell;
	const unsigned int n_comp = fe.get_fe().n_components();
	
	Assert (M.m() == n_dofs, ExcDimensionMismatch(M.m(), n_dofs));
	Assert (M.n() == n_dofs, ExcDimensionMismatch(M.n(), n_dofs));
	
	for (unsigned k=0;k<fe.n_quadrature_points;++k)
	  {
	    const double dx = fe.JxW(k) * factor;
	    const Point<dim>& n = fe.normal_vector(k);
	    for (unsigned i=0;i<n_dofs;++i)
	      for (unsigned j=0;j<n_dofs;++j)
		for (unsigned int d=0;d<n_comp;++d)
		  M(i,j) += dx *
			    (fe.shape_value_component(i,k,d) * penalty * fe.shape_value_component(j,k,d)
			     - (n * fe.shape_grad_component(i,k,d)) * fe.shape_value_component(j,k,d)
			     - (n * fe.shape_grad_component(j,k,d)) * fe.shape_value_component(i,k,d));
	  }
      }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
