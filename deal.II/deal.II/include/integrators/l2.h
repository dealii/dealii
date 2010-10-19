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
#ifndef __deal2__integrators_l2_h
#define __deal2__integrators_l2_h


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
 * @brief Local integrators related to <i>L<sup>2</sup</i>-inner products.
 *
 * @author Guido Kanschat
 * @date 2010
 */
  namespace L2
  {
/**
 * The mass matrix.
 *
 * \f[
 * (a u,v)
 * \f]
 */
    template <int dim>
    void mass_matrix (
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
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<n_components;++d)
		M(i,j) += dx
			  * fe.shape_value_component(j,k,d)
			  * fe.shape_value_component(i,k,d);
	}
    }

/**
 * <i>L<sup>2</sup</i>-inner product for scalar functions.
 *
 * \f[
 * (f,v)
 * \f]
 */
    template <int dim>
    void L2 (
      Vector<double>& result,
      const FEValuesBase<dim>& fe,
      const std::vector<double>& input,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(result.size(), n_dofs);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(input.size(), fe.n_quadrature_points);
      
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	for (unsigned i=0;i<n_dofs;++i)
	  result(i) += fe.JxW(k) * factor * input[k] * fe.shape_value(i,k);
    }

/**
 * <i>L<sup>2</sup</i>-inner product for a slice of a vector valued
 * right hand side.
 *
 * \f[
 * \int_Z f\cdot v\,dx
 * \f]
 */
    template <int dim>
    void L2 (
      Vector<double>& result,
      const FEValuesBase<dim>& fe,
      const VectorSlice<const std::vector<std::vector<double> > >& input,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int fe_components = fe.get_fe().n_components();
      const unsigned int n_components = input.size();
	
      AssertDimension(result.size(), n_dofs);
      AssertDimension(input.size(), fe_components);
	
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	for (unsigned i=0;i<n_dofs;++i)
	  for (unsigned int d=0;d<n_components;++d)
	    result(i) += fe.JxW(k) * factor * fe.shape_value_component(i,k,d) * input[d][k];
    } 
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
