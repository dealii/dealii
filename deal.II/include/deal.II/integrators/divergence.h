//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_divergence_h
#define __deal2__integrators_divergence_h


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
 * @brief Local integrators related to the divergence operator and its trace.
 *
 * @ingroup Integrators
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Divergence
  {
/**
 * Auxiliary function. Computes the grad-dic-operator from a set of
 * Hessians.
 *
 * @note The third tensor argument is not used in two dimensions and
 * can for instance duplicate one of the previous.
 *
 * @author Guido Kanschat
 * @date 2011
 */
    template <int dim>
    Tensor<1,dim>
    grad_div(
      const Tensor<2,dim>& h0,
      const Tensor<2,dim>& h1,
      const Tensor<2,dim>& h2)
    {
      Tensor<1,dim> result;
      for (unsigned int d=0;d<dim;++d)
	{
	  result[d] += h0[d][0];
	  if (dim >=2) result[d] += h1[d][1];
	  if (dim >=3) result[d] += h2[d][2];
	}
      return result;
    }
    
    
/**
 * Cell matrix for divergence. The derivative is on the trial function.
 *
 * \f[
 * \int_Z v\nabla \cdot \mathbf u \,dx
 * \f]
 */
    template <int dim>
    void cell_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const FEValuesBase<dim>& fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);
	
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = fe.JxW(k) * factor;
	  for (unsigned i=0;i<t_dofs;++i)
	    {
	      const double vv = fetest.shape_value(i,k);
	      for (unsigned int d=0;d<dim;++d)
		for (unsigned j=0;j<n_dofs;++j)
		  {
		    const double du = fe.shape_grad_component(j,k,d)[d];
		    M(i,j) += dx * du * vv;
		  }
	    }
	}
    }
/**
 * Cell matrix for divergence. The derivative is on the test function.
 *
 * \f[
 * \int_Z \nabla u \cdot \mathbf v\,dx
 * \f]
 */
    template <int dim>
    void gradient_matrix(
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const FEValuesBase<dim>& fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fetest.get_fe().n_components();
      const unsigned int t_dofs = fetest.dofs_per_cell;
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fecomp, dim);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);
	
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const double dx = fe.JxW(k) * factor;
	  for (unsigned int d=0;d<dim;++d)
	    for (unsigned i=0;i<t_dofs;++i)
	      {
		const double vv = fetest.shape_value_component(i,k,d);
		for (unsigned j=0;j<n_dofs;++j)
		  {
		    const Tensor<1,dim>& Du = fe.shape_grad(j,k);
		    M(i,j) += dx * vv * Du[d];
		  }
	      }
	}
    }

/**
 * The trace of the divergence
 * operator, namely the product
 * of the normal component of the
 * vector valued trial space and
 * the test space.
 * @f[
 * \int_F (\mathbf u\cdot \mathbf n) v \,ds
 * @f]
 */
    template<int dim>
    void
    u_dot_n_matrix (
      FullMatrix<double>& M,
      const FEValuesBase<dim>& fe,
      const FEValuesBase<dim>& fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
	
      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);
	
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const Tensor<1,dim> ndx = factor * fe.JxW(k) * fe.normal_vector(k);
	  for (unsigned i=0;i<t_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<dim;++d)
		M(i,j) += ndx[d] * fe.shape_value_component(j,k,d)
			  * fetest.shape_value(i,k);
	}
    }

/**
 * The trace of the divergence
 * operator, namely the product
 * of the normal component of the
 * vector valued trial space and
 * the test space.
 * @f[
 * \int_F (\mathbf f\cdot \mathbf n) v \,ds
 * @f]
 */
    template<int dim>
    void
    u_dot_n_residual (
      Vector<double>& result,
      const FEValuesBase<dim>& fe,
      const FEValuesBase<dim>& fetest,
      const VectorSlice<const std::vector<std::vector<double> > >& data,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int t_dofs = fetest.dofs_per_cell;
	
      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(result.size(), t_dofs);	
      AssertVectorVectorDimension (data, dim, fe.n_quadrature_points);	
	
      for (unsigned k=0;k<fe.n_quadrature_points;++k)
	{
	  const Tensor<1,dim> ndx = factor * fe.normal_vector(k) * fe.JxW(k);      
	    
	  for (unsigned i=0;i<t_dofs;++i)
	    for (unsigned int d=0;d<dim;++d)
	      result(i) += ndx[d] * fetest.shape_value(i,k) * data[d][k];
	}
    }
/**
 * The trace of the divergence
 * operator, namely the product
 * of the jump of the normal component of the
 * vector valued trial function and
 * the mean value of the test function.
 * @f[
 * \int_F (\mathbf u_1\cdot \mathbf n_1 + \mathbf u_2 \cdot \mathbf n_2) \frac{v_1+v_2}{2} \,ds
 * @f]
 */
    template<int dim>
    void
    u_dot_n_matrix (
      FullMatrix<double>& M11,
      FullMatrix<double>& M12,
      FullMatrix<double>& M21,
      FullMatrix<double>& M22,
      const FEValuesBase<dim>& fe1,
      const FEValuesBase<dim>& fe2,
      const FEValuesBase<dim>& fetest1,
      const FEValuesBase<dim>& fetest2,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;
      const unsigned int t_dofs = fetest1.dofs_per_cell;

      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
      AssertDimension(fetest1.get_fe().n_components(), 1);
      AssertDimension(fetest2.get_fe().n_components(), 1);
      AssertDimension(M11.m(), t_dofs);
      AssertDimension(M11.n(), n_dofs);
      AssertDimension(M12.m(), t_dofs);
      AssertDimension(M12.n(), n_dofs);
      AssertDimension(M21.m(), t_dofs);
      AssertDimension(M21.n(), n_dofs);
      AssertDimension(M22.m(), t_dofs);
      AssertDimension(M22.n(), n_dofs);
	
      for (unsigned k=0;k<fe1.n_quadrature_points;++k)
	{
	  const double dx = factor * fe1.JxW(k);
	  for (unsigned i=0;i<t_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<dim;++d)
		{
		  const double un1 = fe1.shape_value_component(j,k,d) * fe1.normal_vector(k)[d];
		  const double un2 =-fe2.shape_value_component(j,k,d) * fe1.normal_vector(k)[d];
		  const double v1 = fetest1.shape_value(i,k);
		  const double v2 = fetest2.shape_value(i,k);
		  
		  M11(i,j) += .5 * dx * un1 * v1;
		  M12(i,j) += .5 * dx * un2 * v1;
		  M21(i,j) += .5 * dx * un1 * v2;
		  M22(i,j) += .5 * dx * un2 * v2;
		}
	}
    }

/**
 * The weak form of the grad-div operator penalizing volume changes
 * @f[
 *  \int_Z \nabla\!\cdot\!u \nabla\!\cdot\!v \,dx
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

/**
 * The jump of the normal component
 * @f[
 * \int_F
 *  (\mathbf u_1\cdot \mathbf n_1 + \mathbf u_2 \cdot \mathbf n_2)
 *  (\mathbf v_1\cdot \mathbf n_1 + \mathbf v_2 \cdot \mathbf n_2)
 * \,ds
 * @f]
 */
    template<int dim>
    void
    u_dot_n_jump_matrix (
      FullMatrix<double>& M11,
      FullMatrix<double>& M12,
      FullMatrix<double>& M21,
      FullMatrix<double>& M22,
      const FEValuesBase<dim>& fe1,
      const FEValuesBase<dim>& fe2,
      double factor = 1.)
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
	
      for (unsigned k=0;k<fe1.n_quadrature_points;++k)
	{
	  const double dx = factor * fe1.JxW(k);
	  for (unsigned i=0;i<n_dofs;++i)
	    for (unsigned j=0;j<n_dofs;++j)
	      for (unsigned int d=0;d<dim;++d)
		{
		  const double un1 = fe1.shape_value_component(j,k,d) * fe1.normal_vector(k)[d];
		  const double un2 =-fe2.shape_value_component(j,k,d) * fe1.normal_vector(k)[d];
		  const double vn1 = fe1.shape_value_component(i,k,d) * fe1.normal_vector(k)[d];
		  const double vn2 =-fe2.shape_value_component(i,k,d) * fe1.normal_vector(k)[d];
		  
		  M11(i,j) += dx * un1 * vn1;
		  M12(i,j) += dx * un2 * vn1;
		  M21(i,j) += dx * un1 * vn2;
		  M22(i,j) += dx * un2 * vn2;
		}
	}
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
