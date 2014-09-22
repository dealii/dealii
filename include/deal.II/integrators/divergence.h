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
     * Auxiliary function. Computes the grad-div-operator from a set of
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
      const Tensor<2,dim> &h0,
      const Tensor<2,dim> &h1,
      const Tensor<2,dim> &h2)
    {
      Tensor<1,dim> result;
      for (unsigned int d=0; d<dim; ++d)
        {
          result[d] += h0[d][0];
          if (dim >=2) result[d] += h1[d][1];
          if (dim >=3) result[d] += h2[d][2];
        }
      return result;
    }


    /**
     * Cell matrix for divergence. The derivative is on the trial
     * function.
     * \f[
     * \int_Z v\nabla \cdot \mathbf u \,dx
     * \f]
     * This is the strong divergence operator and the trial
     * space should be at least <b>H</b><sup>div</sup>. The test functions
     * may be discontinuous.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void cell_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i=0; i<t_dofs; ++i)
            {
              const double vv = fetest.shape_value(i,k);
              for (unsigned int d=0; d<dim; ++d)
                for (unsigned int j=0; j<n_dofs; ++j)
                  {
                    const double du = fe.shape_grad_component(j,k,d)[d];
                    M(i,j) += dx * du * vv;
                  }
            }
        }
    }

    /**
     * The residual of the divergence operator in strong form.
     * \f[
     * \int_Z v\nabla \cdot \mathbf u \,dx
     * \f]
     * This is the strong divergence operator and the trial
     * space should be at least <b>H</b><sup>div</sup>. The test functions
     * may be discontinuous.
     *
     * The function cell_matrix() is the Frechet derivative of this function with respect
     * to the test functions.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim, typename number>
    void cell_residual(
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &input,
      const double factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert (result.size() == t_dofs, ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) += dx * input[d][k][d] * fetest.shape_value(i,k);
        }
    }


    /**
     * The residual of the divergence operator in weak form.
     * \f[
     * - \int_Z \nabla v \cdot \mathbf u \,dx
     * \f]
     * This is the weak divergence operator and the test
     * space should be at least <b>H</b><sup>1</sup>. The trial functions
     * may be discontinuous.
     *
     * @todo Verify: The function cell_matrix() is the Frechet derivative of this function with respect
     * to the test functions.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template <int dim, typename number>
    void cell_residual(
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<double> > > &input,
      const double factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert (result.size() == t_dofs, ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) -= dx * input[d][k] * fetest.shape_grad(i,k)[d];
        }
    }


    /**
     * Cell matrix for gradient. The derivative is on the trial function.
     * \f[
     * \int_Z \nabla u \cdot \mathbf v\,dx
     * \f]
     *
     * This is the strong gradient and the trial space should be at least
     * in <i>H</i><sup>1</sup>. The test functions can be discontinuous.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void gradient_matrix(
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fetest.get_fe().n_components();
      const unsigned int t_dofs = fetest.dofs_per_cell;
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fecomp, dim);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int i=0; i<t_dofs; ++i)
              {
                const double vv = fetest.shape_value_component(i,k,d);
                for (unsigned int j=0; j<n_dofs; ++j)
                  {
                    const Tensor<1,dim> &Du = fe.shape_grad(j,k);
                    M(i,j) += dx * vv * Du[d];
                  }
              }
        }
    }

    /**
     * The residual of the gradient operator in strong form.
     * \f[
     * \int_Z \mathbf v\cdot\nabla u \,dx
     * \f]
     * This is the strong gradient operator and the trial
     * space should be at least <b>H</b><sup>1</sup>. The test functions
     * may be discontinuous.
     *
     * The function gradient_matrix() is the Frechet derivative of this function with respect to the test functions.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim, typename number>
    void gradient_residual(
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const std::vector<Tensor<1,dim> > &input,
      const double factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(input.size(), fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert (result.size() == t_dofs, ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) += dx * input[k][d] * fetest.shape_value_component(i,k,d);
        }
    }

    /**
     * The residual of the gradient operator in weak form.
     * \f[
     * -\int_Z \nabla\cdot \mathbf v u \,dx
     * \f]
     * This is the weak gradient operator and the test
     * space should be at least <b>H</b><sup>div</sup>. The trial functions
     * may be discontinuous.
     *
     * @todo Verify: The function gradient_matrix() is the Frechet derivative of this function with respect to the test functions.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template <int dim, typename number>
    void gradient_residual(
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const std::vector<double> &input,
      const double factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(input.size(), fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert (result.size() == t_dofs, ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) -= dx * input[k] * fetest.shape_grad_component(i,k,d)[d];
        }
    }

    /**
     * The trace of the divergence operator, namely the product of the
     * normal component of the vector valued trial space and the test
     * space.
     * @f[ \int_F (\mathbf u\cdot \mathbf n) v \,ds @f]
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template<int dim>
    void
    u_dot_n_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const Tensor<1,dim> ndx = factor * fe.JxW(k) * fe.normal_vector(k);
          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<dim; ++d)
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
     * \int_F (\mathbf u\cdot \mathbf n) v \,ds
     * @f]
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template<int dim, typename number>
    void
    u_dot_n_residual (
      Vector<number> &result,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<double> > > &data,
      double factor = 1.)
    {
      unsigned int fecomp = fe.get_fe().n_components();
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fecomp, dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(result.size(), t_dofs);
      AssertVectorVectorDimension (data, dim, fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const Tensor<1,dim> ndx = factor * fe.normal_vector(k) * fe.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) += ndx[d] * fetest.shape_value(i,k) * data[d][k];
        }
    }

    /**
     * The trace of the gradient
     * operator, namely the product
     * of the normal component of the
     * vector valued test space and
     * the trial space.
     * @f[
     * \int_F u (\mathbf v\cdot \mathbf n) \,ds
     * @f]
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template<int dim, typename number>
    void
    u_times_n_residual (
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const std::vector<double> &data,
      double factor = 1.)
    {
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(result.size(), t_dofs);
      AssertDimension(data.size(), fetest.n_quadrature_points);

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const Tensor<1,dim> ndx = factor * fetest.normal_vector(k) * fetest.JxW(k);

          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int d=0; d<dim; ++d)
              result(i) += ndx[d] * fetest.shape_value_component(i,k,d) * data[k];
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
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template<int dim>
    void
    u_dot_n_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      const FEValuesBase<dim> &fetest1,
      const FEValuesBase<dim> &fetest2,
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

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = factor * fe1.JxW(k);
          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<dim; ++d)
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
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void grad_div_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              {
                double dv = 0.;
                double du = 0.;
                for (unsigned int d=0; d<dim; ++d)
                  {
                    dv += fe.shape_grad_component(i,k,d)[d];
                    du += fe.shape_grad_component(j,k,d)[d];
                  }

                M(i,j) += dx * du * dv;
              }
        }
    }

    /**
     * The weak form of the grad-div residual
     * @f[
     *  \int_Z \nabla\!\cdot\!u \nabla\!\cdot\!v \,dx
     * @f]
     *
     * @author Guido Kanschat
     * @date 2014
     */
    template <int dim, typename number>
    void grad_div_residual (
      Vector<number> &result,
      const FEValuesBase<dim> &fetest,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &input,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fetest.dofs_per_cell;

      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);

      for (unsigned int k=0; k<fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            {
              double dv = 0.;
              double du = 0.;
              for (unsigned int d=0; d<dim; ++d)
                {
                  dv += fetest.shape_grad_component(i,k,d)[d];
                  du += input[d][k][d];
                }

              result(i) += dx * du * dv;
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
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template<int dim>
    void
    u_dot_n_jump_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
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

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = factor * fe1.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<dim; ++d)
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

    /**
     * The <i>L</i><sup>2</sup>-norm of the divergence over the
     * quadrature set determined by the FEValuesBase object.
     *
     * The vector is expected to consist of dim vectors of length
     * equal to the number of quadrature points. The number of
     * components of the finite element has to be equal to the space
     * dimension.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template <int dim>
    double norm(const FEValuesBase<dim> &fe,
                const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Du)
    {
      unsigned int fecomp = fe.get_fe().n_components();

      AssertDimension(fecomp, dim);
      AssertVectorVectorDimension (Du, dim, fe.n_quadrature_points);

      double result = 0;
      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          double div = Du[0][k][0];
          for (unsigned int d=1; d<dim; ++d)
            div += Du[d][k][d];
          result += div*div*fe.JxW(k);
        }
      return result;
    }

  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
