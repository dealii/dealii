// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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

#ifndef __deal2__integrators_advection_h
#define __deal2__integrators_advection_h


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
   * @brief Local integrators related to advection along a vector field and its DG formulations
   *
   * All advection operators depend on an advection velocity denoted by
   * <b>w</b> in the formulas below. It is denoted as <tt>velocity</tt>
   * in the parameter lists.
   *
   * The functions cell_matrix() and both upwind_value_matrix() are
   * taking the equation in weak form, that is, the directional
   * derivative is on the test function.
   *
   * @ingroup Integrators
   * @author Guido Kanschat
   * @date 2012
   */
  namespace Advection
  {
    /**
     * Advection along the direction <b>w</b> in weak form
     * with derivative on the test function
     * \f[
     * m_{ij} = \int_Z u_j\,(\mathbf w \cdot \nabla) v_i \, dx.
     * \f]
     *
     * The FiniteElement in <tt>fe</tt> may be scalar or vector valued. In
     * the latter case, the advection operator is applied to each component
     * separately.
     *
     * @param M: The advection matrix obtained as result
     * @param fe: The FEValues object describing the local trial function
     * space. #update_values and #update_gradients, and #update_JxW_values
     * must be set.
     * @param fetest: The FEValues object describing the local test
     * function space. #update_values and #update_gradients must be set.
     * @param velocity: The advection velocity, a vector of dimension
     * <tt>dim</tt>. Each component may either contain a vector of length
     * one, in which case a constant velocity is assumed, or a vector with
     * as many entries as quadrature points if the velocity is not constant.
     * @param factor is an optional multiplication factor for the result.
     *
     * @author Guido Kanschat
     * @date 2012
     */
    template<int dim>
    void cell_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<double> > > &velocity,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();

      AssertDimension(velocity.size(), dim);
      // If the size of the
      // velocity vectors is one,
      // then do not increment
      // between quadrature points.
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;

      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      AssertDimension(M.n(), n_dofs);
      AssertDimension(M.m(), t_dofs);

      for (unsigned k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const unsigned int vindex = k * v_increment;

          for (unsigned j=0; j<n_dofs; ++j)
            for (unsigned i=0; i<t_dofs; ++i)
              for (unsigned int c=0; c<n_components; ++c)
                {
                  double wgradv = velocity[0][vindex]
                                  * fe.shape_grad_component(i,k,c)[0];
                  for (unsigned int d=1; d<dim; ++d)
                    wgradv += velocity[d][vindex]
                              * fe.shape_grad_component(i,k,c)[d];
                  M(i,j) -= dx * wgradv * fe.shape_value_component(j,k,c);
                }
        }
    }

    /**
     * Advection residual operator in strong form
     *
     * \f[
     * r_i = \int_Z  (\mathbf w \cdot \nabla)u\, v_i \, dx.
     * \f]
     */
    template <int dim>
    inline void
    cell_residual  (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const std::vector<Tensor<1,dim> > &input,
      const VectorSlice<const std::vector<std::vector<double> > > &velocity,
      double factor = 1.)
    {
      const unsigned int nq = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      Assert(input.size() == nq, ExcDimensionMismatch(input.size(), nq));
      Assert(result.size() == n_dofs, ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k=0; k<nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i=0; i<n_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) += dx * input[k][d]
                           * fe.shape_value(i,k) * velocity[d][k * v_increment];
        }
    }


    /**
     * Vector-valued advection residual operator in strong form
     *
     *
     * \f[
     * r_i = \int_Z \bigl((\mathbf w \cdot \nabla) \mathbf u\bigr) \cdot\mathbf v_i \, dx.
     * \f]
     */
    template <int dim>
    inline void
    cell_residual  (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &input,
      const VectorSlice<const std::vector<std::vector<double> > > &velocity,
      double factor = 1.)
    {
      const unsigned int nq = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      Assert(result.size() == n_dofs, ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k=0; k<nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i=0; i<n_dofs; ++i)
            for (unsigned int c=0; c<n_comp; ++c)
              for (unsigned int d=0; d<dim; ++d)
                result(i) += dx * input[c][k][d]
                             * fe.shape_value_component(i,k,c) * velocity[d][k * v_increment];
        }
    }

    /**
     * Upwind flux at the boundary
     * for weak advection
     * operator. This is the value of
     * the trial function at the
     * outflow boundary and zero
     * else:
     * @f[
     * a_{ij} = \int_{\partial\Omega}
     * [\mathbf w\cdot\mathbf n]_+
     * u_i v_j \, ds
     * @f]
     *
     * The <tt>velocity</tt> is
     * provided as a VectorSlice,
     * having <tt>dim</tt> vectors,
     * one for each velocity
     * component. Each of the
     * vectors must either have only
     * a single entry, if t he
     * advection velocity is
     * constant, or have an entry
     * for each quadrature point.
     *
     * The finite element can have
     * several components, in which
     * case each component is
     * advected by the same velocity.
     */
    template <int dim>
    void upwind_value_matrix(
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<double> > > &velocity,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      unsigned int n_components = fe.get_fe().n_components();
      AssertDimension (M.m(), n_dofs);
      AssertDimension (M.n(), n_dofs);

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);

          double nv = 0.;
          for (unsigned int d=0; d<dim; ++d)
            nv += fe.normal_vector(k)[d] * velocity[d][k * v_increment];

          if (nv > 0)
            {
              for (unsigned i=0; i<t_dofs; ++i)
                for (unsigned j=0; j<n_dofs; ++j)
                  {
                    if (fe.get_fe().is_primitive())
                      M(i,j) += dx * nv * fe.shape_value(i,k) * fe.shape_value(j,k);
                    else
                      for (unsigned int c=0; c<n_components; ++c)
                        M(i,j) += dx * nv * fetest.shape_value_component(i,k,c)
                                  * fe.shape_value_component(j,k,c);
                  }
            }
        }
    }

    /**
    * Upwind flux in the interior
    * for weak advection
    * operator. Matrix entries
    * correspond to the upwind value
    * of the trial function, multiplied
    * by the jump of the test
    * functions
    * @f[
    * a_{ij} = \int_F \left|\mathbf w
    * \cdot \mathbf n\right|
    * u^\uparrow
    * (v^\uparrow-v^\downarrow)
    * \,ds
    * @f]
    *
    * The <tt>velocity</tt> is
    * provided as a VectorSlice,
    * having <tt>dim</tt> vectors,
    * one for each velocity
    * component. Each of the
    * vectors must either have only
    * a single entry, if t he
    * advection velocity is
    * constant, or have an entry
    * for each quadrature point.
    *
    * The finite element can have
    * several components, in which
    * case each component is
    * advected the same way.
    */
    template <int dim>
    void upwind_value_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      const FEValuesBase<dim> &fetest1,
      const FEValuesBase<dim> &fetest2,
      const VectorSlice<const std::vector<std::vector<double> > > &velocity,
      const double factor = 1.)
    {
      const unsigned int n1 = fe1.dofs_per_cell;
      // Multiply the quadrature point
      // index below with this factor to
      // have simpler data for constant
      // velocities.
      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe1.n_quadrature_points);
        }

      for (unsigned k=0; k<fe1.n_quadrature_points; ++k)
        {
          double nbeta = fe1.normal_vector(k)[0] * velocity[0][k * v_increment];
          for (unsigned int d=1; d<dim; ++d)
            nbeta += fe1.normal_vector(k)[d] * velocity[d][k * v_increment];
          const double dx_nbeta = factor * nbeta * fe1.JxW(k);

          for (unsigned i=0; i<n1; ++i)
            for (unsigned j=0; j<n1; ++j)
              if (fe1.get_fe().is_primitive())
                {
                  if (nbeta > 0)
                    {
                      M11(i,j) += dx_nbeta
                                  * fe1.shape_value(j,k)
                                  * fetest1.shape_value(i,k);
                      M21(i,j) -= dx_nbeta
                                  * fe1.shape_value(j,k)
                                  * fetest2.shape_value(i,k);
                    }
                  else
                    {
                      M22(i,j) -= dx_nbeta
                                  * fe2.shape_value(j,k)
                                  * fetest2.shape_value(i,k);
                      M12(i,j) += dx_nbeta
                                  * fe2.shape_value(j,k)
                                  * fetest1.shape_value(i,k);
                    }
                }
              else
                {
                  for (unsigned int d=0; d<fe1.get_fe().n_components(); ++d)
                    if (nbeta > 0)
                      {
                        M11(i,j) += dx_nbeta
                                    * fe1.shape_value_component(j,k,d)
                                    * fetest1.shape_value_component(i,k,d);
                        M21(i,j) -= dx_nbeta
                                    * fe1.shape_value_component(j,k,d)
                                    * fetest2.shape_value_component(i,k,d);
                      }
                    else
                      {
                        M22(i,j) -= dx_nbeta
                                    * fe2.shape_value_component(j,k,d)
                                    * fetest2.shape_value_component(i,k,d);
                        M12(i,j) += dx_nbeta
                                    * fe2.shape_value_component(j,k,d)
                                    * fetest1.shape_value_component(i,k,d);
                      }
                }
        }
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
