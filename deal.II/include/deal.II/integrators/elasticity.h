//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_elasticity_h
#define __deal2__integrators_elasticity_h


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
 * @brief Local integrators related to elasticity problems.
 *
 * @ingroup Integrators
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Elasticity
  {
				     /**
				      * The weak form of the grad-div
				      * operator penalizing volume changes
				      * @f[
				      * \int_Z \nabla\!\cdot\!u
				      * \nabla\!\cdot\!v \,dx
				      * @f]
				      */
    template <int dim>
    void grad_div_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = factor * fe.JxW(k);
	  for (unsigned i=0;i<n_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      {
		double dv = 0.;
		double du = 0.;
		for (unsigned int d=0;d<dim;++d)
		  {
		    dv += fe.shape_grad_component(i,k,d)[d];
		    du += fe.shape_grad_component(j,k,d)[d];
		  }
		
		M(i,j) += dx * du * dv;
	      }
	}
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
