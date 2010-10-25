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
#include <numerics/mesh_worker_info.h>

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
 * Weak boundary condition of Nitsche type for the Laplacian, namely on the face <i>F</i> the matrix
 * @f[
 * \int_F \Bigl(\gamma u v - \partial_n u v - u \partial_n v\Bigr)\;ds.
 * @f]
 *
 * Here, $\gamma$ is the <tt>penalty</tt> parameter suitably computed
 * with compute_penalty().
 *
 * @author Guido Kanschat
 * @date 2008, 2009, 2010
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
			  (2. * fe.shape_value_component(i,k,d) * penalty * fe.shape_value_component(j,k,d)
			   - (n * fe.shape_grad_component(i,k,d)) * fe.shape_value_component(j,k,d)
			   - (n * fe.shape_grad_component(j,k,d)) * fe.shape_value_component(i,k,d));
	}
    }

/**
 * Weak boundary condition for the Laplace operator by Nitsche, namely on the face <i>F</i>
 * the vector
 * @f[
 * \int_F \Bigl(\gamma (u-g) v - \partial_n u v - (u-g) \partial_n v\Bigr)\;ds.
 * @f]
 *
 * Here, <i>u</i> is the finite element function whose values and
 * gradient are given in the arguments <tt>input</tt> and
 * <tt>Dinput</tt>, respectively. <i>g</i> is the inhomogeneous
 * boundary value in the argument <tt>data</tt>. $\gamma$ is the usual
 * penalty parameter.
 */
      template <int dim>
      void nitsche_residual (
	Vector<double>& result,
	const FEValuesBase<dim>& fe,
	const VectorSlice<const std::vector<std::vector<double> > >& input,
	const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > >& Dinput,
	const VectorSlice<const std::vector<std::vector<double> > >& data,
	double penalty,
	double factor = 1.)
      {
	const unsigned int n_dofs = fe.dofs_per_cell;

	const unsigned int n_comp = fe.get_fe().n_components();
	AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
	AssertVectorVectorDimension(Dinput, n_comp, fe.n_quadrature_points);
	AssertVectorVectorDimension(data, n_comp, fe.n_quadrature_points);

	for (unsigned k=0;k<fe.n_quadrature_points;++k)
	  {
	    const double dx = factor * fe.JxW(k);
	    const Point<dim>& n = fe.normal_vector(k);
	    for (unsigned i=0;i<n_dofs;++i)
	      for (unsigned int d=0;d<n_comp;++d)
		{
		  const double dnv = fe.shape_grad_component(i,k,d) * n;
		  const double dnu = Dinput[d][k] * n;
		  const double v= fe.shape_value_component(i,k,d);
		  const double u= input[d][k];
		  const double g= data[d][k];

		  result(i) += dx*(2.*penalty*(u-g)*v - dnv*(u-g) - dnu*v);
		}
	  }
      }

/**
 * Flux for the interior penalty method for the Laplacian, namely on
 * the face <i>F</i> the matrices associated with the bilinear form
 * @f[
 * \int_F \Bigl( \gamma [u][v] - \{\nabla u\}[v\mathbf n] - [u\mathbf
 * n]\{\nabla v\} \Bigr) \; ds.
 * @f]
 *
 * The penalty parameter should always be the mean value of the
 * penalties needed for stability on each side. In the case of
 * constant coefficients, it can be computed using compute_penalty().
 *
 * If <tt>factor2</tt> is missing or negative, the factor is assumed
 * the same on both sides. If factors differ, note that the penalty
 * parameter has to be computed accordingly.
 *
 * @author Guido Kanschat
 * @date 2008, 2009, 2010
 */
    template <int dim>
    void ip_matrix (
      FullMatrix<double>& M11,
      FullMatrix<double>& M12,
      FullMatrix<double>& M21,
      FullMatrix<double>& M22,
      const FEValuesBase<dim>& fe1,
      const FEValuesBase<dim>& fe2,
      double penalty,
      double factor1 = 1.,
      double factor2 = -1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;
      AssertDimension(M11.n(), n_dofs);
      AssertDimension(M11.m(), n_dofs);
      AssertDimension(M12.n(), n_dofs);
      AssertDimension(M12.m(), n_dofs);
      AssertDimension(M21.n(), n_dofs);
      AssertDimension(M21.m(), n_dofs);
      AssertDimension(M22.n(), n_dofs);
      AssertDimension(M22.m(), n_dofs);

      const double nui = factor1;
      const double nue = (factor2 < 0) ? factor1 : factor2;

      for (unsigned k=0;k<fe1.n_quadrature_points;++k)
	{
	  const double dx = fe1.JxW(k);
	  const Point<dim>& n = fe1.normal_vector(k);
	  for (unsigned i=0;i<n_dofs;++i)
	    {
	      for (unsigned j=0;j<n_dofs;++j)
		{
		  if (fe1.get_fe().n_components() == 1)
		    {
		      const double vi = fe1.shape_value(i,k);
		      const double dnvi = n * fe1.shape_grad(i,k);
		      const double ve = fe2.shape_value(i,k);
		      const double dnve = n * fe2.shape_grad(i,k);
		      const double ui = fe1.shape_value(j,k);
		      const double dnui = n * fe1.shape_grad(j,k);
		      const double ue = fe2.shape_value(j,k);
		      const double dnue = n * fe2.shape_grad(j,k);

		      M11(i,j) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+penalty*ui*vi);
		      M12(i,j) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-penalty*vi*ue);
		      M21(i,j) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-penalty*ui*ve);
		      M22(i,j) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+penalty*ue*ve);
		    }
		  else
		    for (unsigned int d=0;d<dim;++d)
		      {
			const double vi = fe1.shape_value_component(i,k,d);
			const double dnvi = n * fe1.shape_grad_component(i,k,d);
			const double ve = fe2.shape_value_component(i,k,d);
			const double dnve = n * fe2.shape_grad_component(i,k,d);
			const double ui = fe1.shape_value_component(j,k,d);
			const double dnui = n * fe1.shape_grad_component(j,k,d);
			const double ue = fe2.shape_value_component(j,k,d);
			const double dnue = n * fe2.shape_grad_component(j,k,d);

			M11(i,j) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+penalty*ui*vi);
			M12(i,j) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-penalty*vi*ue);
			M21(i,j) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-penalty*ui*ve);
			M22(i,j) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+penalty*ue*ve);
		      }
		}
	    }
	}
    }

/**
 * Auxiliary function computing the penalty parameter for interior
 * penalty methods on rectangles.
 *
 * Computation is done in two steps: first, we compute on each cell
 * <i>Z<sub>i</sub></i> the value <i>P<sub>i</sub> =
 * p<sub>i</sub>(p<sub>i</sub>+1)/h<sub>i</sub></i>, where <i>p<sub>i</sub></i> is
 * the polynomial degree on cell <i>Z<sub>i</sub></i> and
 * <i>h<sub>i</sub></i> is the length of <i>Z<sub>i</sub></i>
 * orthogonal to the current face.
 *
 * @author Guido Kanschat
 * @date 2010
 */
    template <int dim>
    double compute_penalty(
      const MeshWorker::DoFInfo<dim>& dinfo1,
      const MeshWorker::DoFInfo<dim>& dinfo2,
      unsigned int deg1,
      unsigned int deg2)
    {
      const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[dinfo1.face_number];
      const unsigned int normal2 = GeometryInfo<dim>::unit_normal_direction[dinfo2.face_number];
      double penalty1 = deg1 * (deg1+1) / dinfo1.cell->extent_in_direction(normal1);
      double penalty2 = deg2 * (deg2+1) / dinfo2.cell->extent_in_direction(normal2);
      if (dinfo1.cell->has_children() ^ dinfo2.cell->has_children())
	{
	  Assert (dinfo1.face == dinfo2.face, ExcInternalError());
	  Assert (dinfo1.face->has_children(), ExcInternalError());
	  penalty1 *= 2;
	}
      const double penalty = 0.5*(penalty1 + penalty2);
      return penalty;
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
