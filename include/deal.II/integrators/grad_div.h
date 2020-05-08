// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_integrators_grad_div_h
#define dealii_integrators_grad_div_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
  /**
   * Local integrators related to the grad-div operator and its boundary
   * traces
   *
   * @ingroup Integrators
   * @author Guido Kanschat, Timo Heister
   * @date 2016
   */
  namespace GradDiv
  {
    /**
     * The weak form of the grad-div operator penalizing volume changes
     * @f[
     *  \int_Z \nabla\cdot u \nabla \cdot v \,dx
     * @f]
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void
    cell_matrix(FullMatrix<double> &     M,
                const FEValuesBase<dim> &fe,
                double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                const double divu =
                  fe[FEValuesExtractors::Vector(0)].divergence(j, k);
                const double divv =
                  fe[FEValuesExtractors::Vector(0)].divergence(i, k);

                M(i, j) += dx * divu * divv;
              }
        }
    }

    /**
     * The weak form of the grad-div residual
     * @f[
     *  \int_Z \nabla\cdot u \nabla \cdot v \,dx
     * @f]
     *
     * @author Guido Kanschat
     * @date 2014
     */
    template <int dim, typename number>
    void
    cell_residual(Vector<number> &                                    result,
                  const FEValuesBase<dim> &                           fetest,
                  const ArrayView<const std::vector<Tensor<1, dim>>> &input,
                  const double factor = 1.)
    {
      const unsigned int n_dofs = fetest.dofs_per_cell;

      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double divv =
                fetest[FEValuesExtractors::Vector(0)].divergence(i, k);
              double du = 0.;
              for (unsigned int d = 0; d < dim; ++d)
                du += input[d][k][d];

              result(i) += dx * du * divv;
            }
        }
    }

    /**
     * The matrix for the weak boundary condition of Nitsche type for linear
     * elasticity:
     * @f[
     * \int_F \Bigl(\gamma (u \cdot n)(v \cdot n)  - \nabla\cdot u
     * v\cdot n - u \cdot n \nabla \cdot v \Bigr)\;ds.
     * @f]
     */
    template <int dim>
    inline void
    nitsche_matrix(FullMatrix<double> &     M,
                   const FEValuesBase<dim> &fe,
                   double                   penalty,
                   double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double         dx = factor * fe.JxW(k);
          const Tensor<1, dim> n  = fe.normal_vector(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                const double divu =
                  fe[FEValuesExtractors::Vector(0)].divergence(j, k);
                const double divv =
                  fe[FEValuesExtractors::Vector(0)].divergence(i, k);
                double un = 0., vn = 0.;
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    un += fe.shape_value_component(j, k, d) * n[d];
                    vn += fe.shape_value_component(i, k, d) * n[d];
                  }

                M(i, j) += dx * 2. * penalty * un * vn;
                M(i, j) -= dx * (divu * vn + divv * un);
              }
        }
    }

    /**
     * Weak boundary condition for the Laplace operator by Nitsche, vector
     * valued version, namely on the face <i>F</i> the vector
     * @f[
     * \int_F \Bigl(\gamma (\mathbf u \cdot \mathbf n- \mathbf g \cdot
     * \mathbf n) (\mathbf v \cdot \mathbf n)
     * - \nabla \cdot \mathbf u (\mathbf v \cdot \mathbf n)
     * - (\mathbf u-\mathbf g) \cdot \mathbf n \nabla \cdot v\Bigr)\;ds.
     * @f]
     *
     * Here, <i>u</i> is the finite element function whose values and gradient
     * are given in the arguments <tt>input</tt> and <tt>Dinput</tt>,
     * respectively. <i>g</i> is the inhomogeneous boundary value in the
     * argument <tt>data</tt>. $\gamma$ is the usual penalty parameter.
     *
     * @author Guido Kanschat
     * @date 2008, 2009, 2010
     */
    template <int dim>
    void
    nitsche_residual(Vector<double> &                                    result,
                     const FEValuesBase<dim> &                           fe,
                     const ArrayView<const std::vector<double>> &        input,
                     const ArrayView<const std::vector<Tensor<1, dim>>> &Dinput,
                     const ArrayView<const std::vector<double>> &        data,
                     double penalty,
                     double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(fe.get_fe().n_components(), dim)
        AssertVectorVectorDimension(input, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(Dinput, dim, fe.n_quadrature_points);
      AssertVectorVectorDimension(data, dim, fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double         dx = factor * fe.JxW(k);
          const Tensor<1, dim> n  = fe.normal_vector(k);

          double umgn = 0.;
          double divu = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            {
              umgn += (input[d][k] - data[d][k]) * n[d];
              divu += Dinput[d][k][d];
            }

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              double       vn = 0.;
              const double divv =
                fe[FEValuesExtractors::Vector(0)].divergence(i, k);
              for (unsigned int d = 0; d < dim; ++d)
                vn += fe.shape_value_component(i, k, d) * n[d];

              result(i) +=
                dx * (2. * penalty * umgn * vn - divv * umgn - divu * vn);
            }
        }
    }

    /**
     * The interior penalty flux for the grad-div operator. See
     * ip_residual() for details.
     *
     * @author Guido Kanschat
     * @date 2016
     */

    template <int dim>
    void
    ip_matrix(FullMatrix<double> &     M11,
              FullMatrix<double> &     M12,
              FullMatrix<double> &     M21,
              FullMatrix<double> &     M22,
              const FEValuesBase<dim> &fe1,
              const FEValuesBase<dim> &fe2,
              double                   penalty,
              double                   factor1 = 1.,
              double                   factor2 = -1.)
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

      const double fi = factor1;
      const double fe = (factor2 < 0) ? factor1 : factor2;
      const double f  = .5 * (fi + fe);

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double         dx = fe1.JxW(k);
          const Tensor<1, dim> n  = fe1.normal_vector(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                double       uni = 0.;
                double       une = 0.;
                double       vni = 0.;
                double       vne = 0.;
                const double divui =
                  fe1[FEValuesExtractors::Vector(0)].divergence(j, k);
                const double divue =
                  fe2[FEValuesExtractors::Vector(0)].divergence(j, k);
                const double divvi =
                  fe1[FEValuesExtractors::Vector(0)].divergence(i, k);
                const double divve =
                  fe2[FEValuesExtractors::Vector(0)].divergence(i, k);

                for (unsigned int d = 0; d < dim; ++d)
                  {
                    uni += fe1.shape_value_component(j, k, d) * n[d];
                    une += fe2.shape_value_component(j, k, d) * n[d];
                    vni += fe1.shape_value_component(i, k, d) * n[d];
                    vne += fe2.shape_value_component(i, k, d) * n[d];
                  }
                M11(i, j) +=
                  dx * (-.5 * fi * divvi * uni - .5 * fi * divui * vni +
                        f * penalty * uni * vni);
                M12(i, j) +=
                  dx * (.5 * fi * divvi * une - .5 * fe * divue * vni -
                        f * penalty * vni * une);
                M21(i, j) +=
                  dx * (-.5 * fe * divve * uni + .5 * fi * divui * vne -
                        f * penalty * uni * vne);
                M22(i, j) +=
                  dx * (.5 * fe * divve * une + .5 * fe * divue * vne +
                        f * penalty * une * vne);
              }
        }
    }

    /**
     * Grad-div residual term for the symmetric interior penalty method:
     * @f[
     * \int_F \Bigl( \gamma [\mathbf u \cdot\mathbf n]
     * \cdot[\mathbf v \cdot \mathbf n]
     * - \{\nabla \cdot \mathbf u\}[\mathbf v\cdot \mathbf n]
     * - [\mathbf u\times \mathbf n]\{\nabla\cdot \mathbf v\} \Bigr) \; ds.
     * @f]
     *
     * See for instance Hansbo and Larson, 2002
     *
     * @author Guido Kanschat
     * @date 2016
     */
    template <int dim>
    void
    ip_residual(Vector<double> &                                    result1,
                Vector<double> &                                    result2,
                const FEValuesBase<dim> &                           fe1,
                const FEValuesBase<dim> &                           fe2,
                const ArrayView<const std::vector<double>> &        input1,
                const ArrayView<const std::vector<Tensor<1, dim>>> &Dinput1,
                const ArrayView<const std::vector<double>> &        input2,
                const ArrayView<const std::vector<Tensor<1, dim>>> &Dinput2,
                double                                              pen,
                double int_factor = 1.,
                double ext_factor = -1.)
    {
      const unsigned int n1 = fe1.dofs_per_cell;

      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertVectorVectorDimension(input1, dim, fe1.n_quadrature_points);
      AssertVectorVectorDimension(Dinput1, dim, fe1.n_quadrature_points);
      AssertVectorVectorDimension(input2, dim, fe2.n_quadrature_points);
      AssertVectorVectorDimension(Dinput2, dim, fe2.n_quadrature_points);

      const double fi      = int_factor;
      const double fe      = (ext_factor < 0) ? int_factor : ext_factor;
      const double penalty = .5 * pen * (fi + fe);


      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double         dx    = fe1.JxW(k);
          const Tensor<1, dim> n     = fe1.normal_vector(k);
          double               uni   = 0.;
          double               une   = 0.;
          double               divui = 0.;
          double               divue = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            {
              uni += input1[d][k] * n[d];
              une += input2[d][k] * n[d];
              divui += Dinput1[d][k][d];
              divue += Dinput2[d][k][d];
            }

          for (unsigned int i = 0; i < n1; ++i)
            {
              double       vni = 0.;
              double       vne = 0.;
              const double divvi =
                fe1[FEValuesExtractors::Vector(0)].divergence(i, k);
              const double divve =
                fe2[FEValuesExtractors::Vector(0)].divergence(i, k);
              for (unsigned int d = 0; d < dim; ++d)
                {
                  vni += fe1.shape_value_component(i, k, d) * n[d];
                  vne += fe2.shape_value_component(i, k, d) * n[d];
                }

              result1(i) += dx * (-.5 * fi * divvi * uni -
                                  .5 * fi * divui * vni + penalty * uni * vni);
              result1(i) += dx * (.5 * fi * divvi * une -
                                  .5 * fe * divue * vni - penalty * vni * une);
              result2(i) += dx * (-.5 * fe * divve * uni +
                                  .5 * fi * divui * vne - penalty * uni * vne);
              result2(i) += dx * (.5 * fe * divve * une +
                                  .5 * fe * divue * vne + penalty * une * vne);
            }
        }
    }
  } // namespace GradDiv
} // namespace LocalIntegrators

DEAL_II_NAMESPACE_CLOSE


#endif
