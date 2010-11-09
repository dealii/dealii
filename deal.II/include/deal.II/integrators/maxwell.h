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
#ifndef __deal2__integrators_maxwell_h
#define __deal2__integrators_maxwell_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <fe/mapping.h>
#include <fe/fe_values.h>
#include <numerics/mesh_worker_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
/**
 * @brief Local integrators related to curl operators and their traces.
 *
 * @ingroup Integrators
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Maxwell
  {
				     /**
				      * The curl-curl operator
				      * @f[
				      * \int_Z \nabla\!\times\! u \cdot
				      * \nabla\!\times\! v \,dx
				      * @f]
				      */
    template <int dim>
    void curl_curl_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      
				       // Depending on the dimension,
				       // the cross product is either
				       // a scalar (2d) or a vector
				       // (3d). Accordingly, in the
				       // latter case we have to sum
				       // up three bilinear forms, but
				       // in 2d, we don't. Thus, we
				       // need to adapt the loop over
				       // all dimensions
      const unsigned int d_max = (dim==2) ? 1 : dim;
      
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = factor * fe.JxW(k);
	  for (unsigned i=0;i<n_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<d_max;++d)
		{
		  const unsigned int d1 = (d+1)%dim;
		  const unsigned int d2 = (d+2)%dim;
		  
		  const double cv = fe.shape_grad_component(i,k,d1)[d2] - fe.shape_grad_component(i,k,d2)[d1];
		  const double cu = fe.shape_grad_component(j,k,d1)[d2] - fe.shape_grad_component(j,k,d2)[d1];
		  
		  M(i,j) += dx * cu * cv;
		}
	}
    }
    
    				     /**
				      * The curl operator
				      * @f[
				      * \int_Z \nabla\!\times\! u \cdot v \,dx.
				      * @f]
				      *
				      * This is the standard curl
				      * operator in 3D and the scalar
				      * curl in 2D.
				      */
    template <int dim>
    void curl_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const FEValuesBase<dim>& fetest,
      double factor = 1.)
    {
      unsigned int t_comp = (dim==3) ? dim : 1;
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(fetest.get_fe().n_components(), t_comp);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);
      
      const unsigned int d_max = (dim==2) ? 1 : dim;
      
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = fe.JxW(k) * factor;
	  for (unsigned i=0;i<t_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<d_max;++d)
		{
		  const unsigned int d1 = (d+1)%dim;
		  const unsigned int d2 = (d+2)%dim;
		  
		  const double vv = fetest.shape_value_component(i,k,d);
		  const double cu = fe.shape_grad_component(j,k,d1)[d2] - fe.shape_grad_component(j,k,d2)[d1];
		  M(i,j) += dx * cu * vv;
		}
	}
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
