//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_elasticity_h
#define __deal2__integrators_elasticity_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
/**
 * @brief Local integrators related to elasticity problems.
 *
 * @ingroup Integrators
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Elasticity
  {
/**
 * Scalar product of symmetric gradients.
 *
 * \f[
 * (\varepsilon(u), \varepsilon(v))
 * \f]				      
 */
    template <int dim>
    inline void cell_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      
      for (unsigned int k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = factor * fe.JxW(k);
	  for (unsigned int i=0;i<n_dofs;++i)
	    for (unsigned int j=0;j<n_dofs;++j)
	      for (unsigned int d1=0;d1<dim;++d1)
		for (unsigned int d2=0;d2<dim;++d2)
		  M(i,j) += dx * .25 *
			    (fe.shape_grad_component(j,k,d1)[d2] + fe.shape_grad_component(j,k,d2)[d1]) *
			    (fe.shape_grad_component(i,k,d1)[d2] + fe.shape_grad_component(i,k,d2)[d1]);
	}
    }


/**
 * The weak boundary condition
 * of Nitsche type for
 * symmetric gradients.
 */
    template <int dim>
    inline void nitsche_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      
      for (unsigned int k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = factor * fe.JxW(k);
	  const Point<dim>& n = fe.normal_vector(k);
	  for (unsigned int i=0;i<n_dofs;++i)
	    for (unsigned int j=0;j<n_dofs;++j)
	      for (unsigned int d1=0;d1<dim;++d1)
		{
		  const double u = fe.shape_value_component(j,k,d1);
		  const double v = fe.shape_value_component(i,k,d1);
		  M(i,j) += dx * penalty * u * v;
		  for (unsigned int d2=0;d2<dim;++d2)
		    {
						       // v . nabla u n
		      M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d1)[d2] *n(d2)* v;
						       // v (nabla u)^T n		      
		      M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d2)[d1] *n(d2)* v;
						       // u  nabla v n
		      M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d1)[d2] *n(d2)* u;
						       // u (nabla v)^T n		      
		      M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d2)[d1] *n(d2)* u;
		    }
		}
	}
    }
    
				     /**
				      * The interior penalty flux
				      * for symmetric gradients.
				      */
    template <int dim>
    inline void ip_matrix (
      FullMatrix<double>& M11,
      FullMatrix<double>& M12,
      FullMatrix<double>& M21,
      FullMatrix<double>& M22,
      const FEValuesBase<dim>& fe1,
      const FEValuesBase<dim>& fe2,
      const double pen,
      const double int_factor = 1.,
      const double ext_factor = -1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;
	
      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
      AssertDimension(M11.m(), n_dofs);
      AssertDimension(M11.n(), n_dofs);
      AssertDimension(M12.m(), n_dofs);
      AssertDimension(M12.n(), n_dofs);
      AssertDimension(M21.m(), n_dofs);
      AssertDimension(M21.n(), n_dofs);
      AssertDimension(M22.m(), n_dofs);
      AssertDimension(M22.n(), n_dofs);
	
      const double nu1 = int_factor;
      const double nu2 = (ext_factor < 0) ? int_factor : ext_factor;
      const double penalty = .5 * pen * (nu1 + nu2);
	
      for (unsigned int k=0;k<fe1.n_quadrature_points;++k)
	{
	  const double dx = fe1.JxW(k);
	  const Point<dim>& n = fe1.normal_vector(k);
	  for (unsigned int i=0;i<n_dofs;++i)
	    for (unsigned int j=0;j<n_dofs;++j)
	      for (unsigned int d1=0;d1<dim;++d1)
		{
		  const double u1 = fe1.shape_value_component(j,k,d1);
		  const double u2 = fe2.shape_value_component(j,k,d1);
		  const double v1 = fe1.shape_value_component(i,k,d1);
		  const double v2 = fe2.shape_value_component(i,k,d1);
		    
		  M11(i,j) += dx * penalty * u1*v1;
		  M12(i,j) -= dx * penalty * u2*v1;
		  M21(i,j) -= dx * penalty * u1*v2;
		  M22(i,j) += dx * penalty * u2*v2;
		    
		  for (unsigned int d2=0;d2<dim;++d2)
		    {
						       // v . nabla u n
		      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(j,k,d1)[d2] * n(d2) * v1;
		      M12(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(j,k,d1)[d2] * n(d2) * v1;
		      M21(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(j,k,d1)[d2] * n(d2) * v2;
		      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(j,k,d1)[d2] * n(d2) * v2;
						       // v (nabla u)^T n		      
		      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(j,k,d2)[d1] * n(d2) * v1;
		      M12(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(j,k,d2)[d1] * n(d2) * v1;
		      M21(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(j,k,d2)[d1] * n(d2) * v2;
		      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(j,k,d2)[d1] * n(d2) * v2;
						       // u  nabla v n
		      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(i,k,d1)[d2] * n(d2) * u1;
		      M12(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(i,k,d1)[d2] * n(d2) * u2;
		      M21(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(i,k,d1)[d2] * n(d2) * u1;
		      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(i,k,d1)[d2] * n(d2) * u2;
						       // u (nabla v)^T n		      
		      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(i,k,d2)[d1] * n(d2) * u1;
		      M12(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(i,k,d2)[d1] * n(d2) * u2;
		      M21(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(i,k,d2)[d1] * n(d2) * u1;
		      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(i,k,d2)[d1] * n(d2) * u2;
		    }
		}
	}
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
