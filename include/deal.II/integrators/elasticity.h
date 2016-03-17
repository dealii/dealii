// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

#ifndef dealii__integrators_elasticity_h
#define dealii__integrators_elasticity_h


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
     * The linear elasticity operator in weak form, namely double contraction
     * of symmetric gradients.
     *
     * \f[ \int_Z \varepsilon(u): \varepsilon(v)\,dx \f]
     */
    template <int dim>
    inline void cell_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const double factor = 1.)
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
              for (unsigned int d1=0; d1<dim; ++d1)
                for (unsigned int d2=0; d2<dim; ++d2)
                  M(i,j) += dx * .25 *
                            (fe.shape_grad_component(j,k,d1)[d2] + fe.shape_grad_component(j,k,d2)[d1]) *
                            (fe.shape_grad_component(i,k,d1)[d2] + fe.shape_grad_component(i,k,d2)[d1]);
        }
    }


    /**
     * Vector-valued residual operator for linear elasticity in weak form
     *
     * \f[ - \int_Z \varepsilon(u): \varepsilon(v) \,dx \f]
     */
    template <int dim, typename number>
    inline void
    cell_residual  (
      Vector<number> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &input,
      double factor = 1.)
    {
      const unsigned int nq = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(fe.get_fe().n_components(), dim);

      AssertVectorVectorDimension(input, dim, fe.n_quadrature_points);
      Assert(result.size() == n_dofs, ExcDimensionMismatch(result.size(), n_dofs));

      for (unsigned int k=0; k<nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int d1=0; d1<dim; ++d1)
              for (unsigned int d2=0; d2<dim; ++d2)
                {
                  result(i) += dx * .25 *
                               (input[d1][k][d2] + input[d2][k][d1]) *
                               (fe.shape_grad_component(i,k,d1)[d2] + fe.shape_grad_component(i,k,d2)[d1]);
                }
        }
    }


    /**
     * The matrix for the weak boundary condition of Nitsche type for linear elasticity:
     * @f[
     * \int_F \Bigl(\gamma u \cdot v - n^T \epsilon(u) v - u \epsilon(v) n\Bigr)\;ds.
     * @f]
     */
    template <int dim>
    inline void nitsche_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d1=0; d1<dim; ++d1)
                {
                  const double u = fe.shape_value_component(j,k,d1);
                  const double v = fe.shape_value_component(i,k,d1);
                  M(i,j) += dx * 2. * penalty * u * v;
                  for (unsigned int d2=0; d2<dim; ++d2)
                    {
                      // v . nabla u n
                      M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d1)[d2] *n[d2]* v;
                      // v (nabla u)^T n
                      M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d2)[d1] *n[d2]* v;
                      // u  nabla v n
                      M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d1)[d2] *n[d2]* u;
                      // u (nabla v)^T n
                      M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d2)[d1] *n[d2]* u;
                    }
                }
        }
    }

    /**
     * The matrix for the weak boundary condition of Nitsche type for the tangential displacement in linear elasticity:
     * @f[
     * \int_F \Bigl(\gamma u_\tau \cdot v_\tau - n^T \epsilon(u_\tau) v_\tau - u_\tau^T \epsilon(v_\tau) n\Bigr)\;ds.
     * @f]
     */
    template <int dim>
    inline void nitsche_tangential_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              {
                double udotn = 0.;
                double vdotn = 0.;
                double ngradun = 0.;
                double ngradvn = 0.;

                for (unsigned int d=0; d<dim; ++d)
                  {
                    udotn += n[d]*fe.shape_value_component(j,k,d);
                    vdotn += n[d]*fe.shape_value_component(i,k,d);
                    ngradun += n*fe.shape_grad_component(j,k,d)*n[d];
                    ngradvn += n*fe.shape_grad_component(i,k,d)*n[d];
                  }
                for (unsigned int d1=0; d1<dim; ++d1)
                  {
                    const double u = fe.shape_value_component(j,k,d1) - udotn * n[d1];
                    const double v = fe.shape_value_component(i,k,d1) - vdotn * n[d1];
                    M(i,j) += dx * 2. * penalty * u * v;
                    // Correct the gradients below and subtract normal component
                    M(i,j) += dx * (ngradun * v + ngradvn * u);
                    for (unsigned int d2=0; d2<dim; ++d2)
                      {
                        // v . nabla u n
                        M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d1)[d2] *n[d2]* v;
                        // v (nabla u)^T n
                        M(i,j) -= .5*dx* fe.shape_grad_component(j,k,d2)[d1] *n[d2]* v;
                        // u  nabla v n
                        M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d1)[d2] *n[d2]* u;
                        // u (nabla v)^T n
                        M(i,j) -= .5*dx* fe.shape_grad_component(i,k,d2)[d1] *n[d2]* u;
                      }
                  }
              }
        }
    }

    /**
     * Weak boundary condition for the elasticity operator by Nitsche, namely
     * on the face <i>F</i> the vector
     * @f[
     * \int_F \Bigl(\gamma (u-g) \cdot v - n^T \epsilon(u) v - (u-g) \epsilon(v) n^T\Bigr)\;ds.
     * @f]
     *
     * Here, <i>u</i> is the finite element function whose values and gradient
     * are given in the arguments <tt>input</tt> and <tt>Dinput</tt>,
     * respectively. <i>g</i> is the inhomogeneous boundary value in the
     * argument <tt>data</tt>. $n$ is the outer normal vector and $\gamma$ is
     * the usual penalty parameter.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template <int dim, typename number>
    void nitsche_residual (
      Vector<number> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<double> > > &input,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput,
      const VectorSlice<const std::vector<std::vector<double> > > &data,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertVectorVectorDimension(input, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(Dinput, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(data, dim, fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int d1=0; d1<dim; ++d1)
              {
                const double u= input[d1][k];
                const double v= fe.shape_value_component(i,k,d1);
                const double g= data[d1][k];
                result(i) += dx * 2.*penalty * (u-g) * v;

                for (unsigned int d2=0; d2<dim; ++d2)
                  {
                    // v . nabla u n
                    result(i) -= .5*dx* v * Dinput[d1][k][d2] * n[d2];
                    // v . (nabla u)^T n
                    result(i) -= .5*dx* v * Dinput[d2][k][d1] * n[d2];
                    // u  nabla v n
                    result(i) -= .5*dx * (u-g) * fe.shape_grad_component(i,k,d1)[d2] * n[d2];
                    // u  (nabla v)^T n
                    result(i) -= .5*dx * (u-g) * fe.shape_grad_component(i,k,d2)[d1] * n[d2];
                  }
              }
        }
    }

    /**
     * The weak boundary condition of Nitsche type for the tangential displacement in linear elasticity:
     * @f[
     * \int_F \Bigl(\gamma (u_\tau-g_\tau) \cdot v_\tau - n^T \epsilon(u_\tau) v - (u_\tau-g_\tau) \epsilon(v_\tau) n\Bigr)\;ds.
     * @f]
     */
    template <int dim, typename number>
    inline void nitsche_tangential_residual (
      Vector<number> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<double> > > &input,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput,
      const VectorSlice<const std::vector<std::vector<double> > > &data,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertVectorVectorDimension(input, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(Dinput, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(data, dim, fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            {
              double udotn = 0.;
              double gdotn = 0.;
              double vdotn = 0.;
              double ngradun = 0.;
              double ngradvn = 0.;

              for (unsigned int d=0; d<dim; ++d)
                {
                  udotn += n[d]*input[d][k];
                  gdotn += n[d]*data[d][k];
                  vdotn += n[d]*fe.shape_value_component(i,k,d);
                  ngradun += n*Dinput[d][k]*n[d];
                  ngradvn += n*fe.shape_grad_component(i,k,d)*n[d];
                }
              for (unsigned int d1=0; d1<dim; ++d1)
                {
                  const double u= input[d1][k] - udotn*n[d1];
                  const double v= fe.shape_value_component(i,k,d1) - vdotn*n[d1];
                  const double g= data[d1][k] - gdotn*n[d1];
                  result(i) += dx * 2.*penalty * (u-g) * v;
                  // Correct the gradients below and subtract normal component
                  result(i) += dx * (ngradun * v + ngradvn * (u-g));
                  for (unsigned int d2=0; d2<dim; ++d2)
                    {

                      // v . nabla u n
                      result(i) -= .5*dx* Dinput[d1][k][d2] *n[d2]* v;
                      // v (nabla u)^T n
                      result(i) -= .5*dx* Dinput[d2][k][d1] *n[d2]* v;
                      // u  nabla v n
                      result(i) -= .5*dx* (u-g) * fe.shape_grad_component(i,k,d1)[d2] *n[d2];
                      // u (nabla v)^T n
                      result(i) -= .5*dx* (u-g) * fe.shape_grad_component(i,k,d2)[d1] *n[d2];
                    }
                }
            }
        }
    }

    /**
     * Homogeneous weak boundary condition for the elasticity operator by
     * Nitsche, namely on the face <i>F</i> the vector
     * @f[
     * \int_F \Bigl(\gamma u \cdot v - n^T \epsilon(u) v - u \epsilon(v) n^T\Bigr)\;ds.
     * @f]
     *
     * Here, <i>u</i> is the finite element function whose values and gradient
     * are given in the arguments <tt>input</tt> and <tt>Dinput</tt>,
     * respectively. $n$ is the outer normal vector and $\gamma$ is the usual
     * penalty parameter.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template <int dim, typename number>
    void nitsche_residual_homogeneous (
      Vector<number> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<double> > > &input,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertVectorVectorDimension(input, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(Dinput, dim, fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int d1=0; d1<dim; ++d1)
              {
                const double u= input[d1][k];
                const double v= fe.shape_value_component(i,k,d1);
                result(i) += dx * 2.*penalty * u * v;

                for (unsigned int d2=0; d2<dim; ++d2)
                  {
                    // v . nabla u n
                    result(i) -= .5*dx* v * Dinput[d1][k][d2] * n[d2];
                    // v . (nabla u)^T n
                    result(i) -= .5*dx* v * Dinput[d2][k][d1] * n[d2];
                    // u  nabla v n
                    result(i) -= .5*dx * u * fe.shape_grad_component(i,k,d1)[d2] * n[d2];
                    // u  (nabla v)^T n
                    result(i) -= .5*dx * u * fe.shape_grad_component(i,k,d2)[d1] * n[d2];
                  }
              }
        }
    }

    /**
     * The interior penalty flux for symmetric gradients.
     */
    template <int dim>
    inline void ip_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
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

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Tensor<1,dim> n = fe1.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d1=0; d1<dim; ++d1)
                {
                  const double u1 = fe1.shape_value_component(j,k,d1);
                  const double u2 = fe2.shape_value_component(j,k,d1);
                  const double v1 = fe1.shape_value_component(i,k,d1);
                  const double v2 = fe2.shape_value_component(i,k,d1);

                  M11(i,j) += dx * penalty * u1*v1;
                  M12(i,j) -= dx * penalty * u2*v1;
                  M21(i,j) -= dx * penalty * u1*v2;
                  M22(i,j) += dx * penalty * u2*v2;

                  for (unsigned int d2=0; d2<dim; ++d2)
                    {
                      // v . nabla u n
                      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(j,k,d1)[d2] * n[d2] * v1;
                      M12(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(j,k,d1)[d2] * n[d2] * v1;
                      M21(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(j,k,d1)[d2] * n[d2] * v2;
                      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(j,k,d1)[d2] * n[d2] * v2;
                      // v (nabla u)^T n
                      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(j,k,d2)[d1] * n[d2] * v1;
                      M12(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(j,k,d2)[d1] * n[d2] * v1;
                      M21(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(j,k,d2)[d1] * n[d2] * v2;
                      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(j,k,d2)[d1] * n[d2] * v2;
                      // u  nabla v n
                      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(i,k,d1)[d2] * n[d2] * u1;
                      M12(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(i,k,d1)[d2] * n[d2] * u2;
                      M21(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(i,k,d1)[d2] * n[d2] * u1;
                      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(i,k,d1)[d2] * n[d2] * u2;
                      // u (nabla v)^T n
                      M11(i,j) -= .25 * dx * nu1 * fe1.shape_grad_component(i,k,d2)[d1] * n[d2] * u1;
                      M12(i,j) += .25 * dx * nu1 * fe1.shape_grad_component(i,k,d2)[d1] * n[d2] * u2;
                      M21(i,j) -= .25 * dx * nu2 * fe2.shape_grad_component(i,k,d2)[d1] * n[d2] * u1;
                      M22(i,j) += .25 * dx * nu2 * fe2.shape_grad_component(i,k,d2)[d1] * n[d2] * u2;
                    }
                }
        }
    }
    /**
     * Elasticity residual term for the symmetric interior penalty method.
     *
     * @author Guido Kanschat
     * @date 2013
     */
    template<int dim, typename number>
    void
    ip_residual(
      Vector<number> &result1,
      Vector<number> &result2,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      const VectorSlice<const std::vector<std::vector<double> > > &input1,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput1,
      const VectorSlice<const std::vector<std::vector<double> > > &input2,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput2,
      double pen,
      double int_factor = 1.,
      double ext_factor = -1.)
    {
      const unsigned int n1 = fe1.dofs_per_cell;

      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
      AssertVectorVectorDimension(input1, dim, fe1.n_quadrature_points);
      AssertVectorVectorDimension(Dinput1, dim, fe1.n_quadrature_points);
      AssertVectorVectorDimension(input2, dim, fe2.n_quadrature_points);
      AssertVectorVectorDimension(Dinput2, dim, fe2.n_quadrature_points);

      const double nu1 = int_factor;
      const double nu2 = (ext_factor < 0) ? int_factor : ext_factor;
      const double penalty = .5 * pen * (nu1 + nu2);


      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Tensor<1,dim> n = fe1.normal_vector(k);

          for (unsigned int i=0; i<n1; ++i)
            for (unsigned int d1=0; d1<dim; ++d1)
              {
                const double v1 = fe1.shape_value_component(i,k,d1);
                const double v2 = fe2.shape_value_component(i,k,d1);
                const double u1 = input1[d1][k];
                const double u2 = input2[d1][k];

                result1(i) += dx * penalty * u1*v1;
                result1(i) -= dx * penalty * u2*v1;
                result2(i) -= dx * penalty * u1*v2;
                result2(i) += dx * penalty * u2*v2;

                for (unsigned int d2=0; d2<dim; ++d2)
                  {
                    // v . nabla u n
                    result1(i) -= .25*dx* (nu1*Dinput1[d1][k][d2]+nu2*Dinput2[d1][k][d2]) * n[d2] * v1;
                    result2(i) += .25*dx* (nu1*Dinput1[d1][k][d2]+nu2*Dinput2[d1][k][d2]) * n[d2] * v2;
                    // v . (nabla u)^T n
                    result1(i) -= .25*dx* (nu1*Dinput1[d2][k][d1]+nu2*Dinput2[d2][k][d1]) * n[d2] * v1;
                    result2(i) += .25*dx* (nu1*Dinput1[d2][k][d1]+nu2*Dinput2[d2][k][d1]) * n[d2] * v2;
                    // u  nabla v n
                    result1(i) -= .25*dx* nu1*fe1.shape_grad_component(i,k,d1)[d2] * n[d2] * (u1-u2);
                    result2(i) -= .25*dx* nu2*fe2.shape_grad_component(i,k,d1)[d2] * n[d2] * (u1-u2);
                    // u  (nabla v)^T n
                    result1(i) -= .25*dx* nu1*fe1.shape_grad_component(i,k,d2)[d1] * n[d2] * (u1-u2);
                    result2(i) -= .25*dx* nu2*fe2.shape_grad_component(i,k,d2)[d1] * n[d2] * (u1-u2);
                  }
              }
        }
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
