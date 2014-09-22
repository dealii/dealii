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

#ifndef __deal2__integrators_maxwell_h
#define __deal2__integrators_maxwell_h


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
   * @brief Local integrators related to curl operators and their
   * traces.
   *
   * We use the following conventions for curl
   * operators. First, in three space dimensions
   *
   * @f[
   * \nabla\times \mathbf u = \begin{pmatrix}
   *   \partial_3 u_2 - \partial_2 u_3 \\
   *   \partial_1 u_3 - \partial_3 u_1 \\
   *   \partial_2 u_1 - \partial_1 u_2
   * \end{pmatrix}.
   * @f]
   *
   * In two space dimensions, the curl is obtained by extending a vector
   * <b>u</b> to $(u_1, u_2, 0)^T$ and a scalar <i>p</i> to $(0,0,p)^T$.
   * Computing the nonzero components, we obtain the scalar
   * curl of a vector function and the vector curl of a scalar
   * function. The current implementation exchanges the sign and we have:
   * @f[
   *  \nabla \times \mathbf u = \partial_1 u_2 - \partial 2 u_1,
   *  \qquad
   *  \nabla \times p = \begin{pmatrix}
   *    \partial_2 p \\ -\partial_1 p
   *  \end{pmatrix}
   * @f]
   *
   * @ingroup Integrators
   * @author Guido Kanschat
   * @date 2010
   */
  namespace Maxwell
  {
    /**
     * Auxiliary function. Given the tensors of <tt>dim</tt> second derivatives,
     * compute the curl of the curl of a vector function. The result in
     * two and three dimensions is:
     * @f[
     * \nabla\times\nabla\times \mathbf u = \begin{pmatrix}
     * \partial_1\partial_2 u_2 - \partial_2^2 u_1 \\
     * \partial_1\partial_2 u_1 - \partial_1^2 u_2
     * \end{pmatrix}
     *
     * \nabla\times\nabla\times \mathbf u = \begin{pmatrix}
     * \partial_1\partial_2 u_2 + \partial_1\partial_3 u_3
     * - (\partial_2^2+\partial_3^2) u_1 \\
     * \partial_2\partial_3 u_3 + \partial_2\partial_1 u_1
     * - (\partial_3^2+\partial_1^2) u_2 \\
     * \partial_3\partial_1 u_1 + \partial_3\partial_2 u_2
     * - (\partial_1^2+\partial_2^2) u_3
     * \end{pmatrix}
     * @f]
     *
     * @note The third tensor argument is not used in two dimensions and
     * can for instance duplicate one of the previous.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    Tensor<1,dim>
    curl_curl (
      const Tensor<2,dim> &h0,
      const Tensor<2,dim> &h1,
      const Tensor<2,dim> &h2)
    {
      Tensor<1,dim> result;
      switch (dim)
        {
        case 2:
          result[0] = h1[0][1]-h0[1][1];
          result[1] = h0[0][1]-h1[0][0];
          break;
        case 3:
          result[0] = h1[0][1]+h2[0][2]-h0[1][1]-h0[2][2];
          result[1] = h2[1][2]+h0[1][0]-h1[2][2]-h1[0][0];
          result[2] = h0[2][0]+h1[2][1]-h2[0][0]-h2[1][1];
          break;
        default:
          Assert(false, ExcNotImplemented());
        }
      return result;
    }

    /**
     * Auxiliary function. Given <tt>dim</tt> tensors of first
     * derivatives and a normal vector, compute the tangential curl
     * @f[
     * \mathbf n \times \nabla \times u.
     * @f]
     *
     * @note The third tensor argument is not used in two dimensions and
     * can for instance duplicate one of the previous.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    Tensor<1,dim>
    tangential_curl (
      const Tensor<1,dim> &g0,
      const Tensor<1,dim> &g1,
      const Tensor<1,dim> &g2,
      const Tensor<1,dim> &normal)
    {
      Tensor<1,dim> result;

      switch (dim)
        {
        case 2:
          result[0] = normal[1] * (g1[0]-g0[1]);
          result[1] =-normal[0] * (g1[0]-g0[1]);
          break;
        case 3:
          result[0] = normal[2]*(g2[1]-g0[2])+normal[1]*(g1[0]-g0[1]);
          result[1] = normal[0]*(g0[2]-g1[0])+normal[2]*(g2[1]-g1[2]);
          result[2] = normal[1]*(g1[0]-g2[1])+normal[0]*(g0[2]-g2[0]);
          break;
        default:
          Assert(false, ExcNotImplemented());
        }
      return result;
    }

    /**
     * The curl-curl operator
     * @f[
     * \int_Z \nabla\!\times\! u \cdot
     * \nabla\!\times\! v \,dx
     * @f]
     * in weak form.
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void curl_curl_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
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

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<d_max; ++d)
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
     * The matrix for the curl operator
     * @f[
     * \int_Z \nabla\!\times\! u \cdot v \,dx.
     * @f]
     *
     * This is the standard curl operator in 3D and the scalar curl in
     * 2D. The vector curl operator can be obtained by exchanging test and
     * trial functions.
     *
     * @author Guido Kanschat
     * @date 2011
    */
    template <int dim>
    void curl_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const FEValuesBase<dim> &fetest,
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

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i=0; i<t_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<d_max; ++d)
                {
                  const unsigned int d1 = (d+1)%dim;
                  const unsigned int d2 = (d+2)%dim;

                  const double vv = fetest.shape_value_component(i,k,d);
                  const double cu = fe.shape_grad_component(j,k,d1)[d2] - fe.shape_grad_component(j,k,d2)[d1];
                  M(i,j) += dx * cu * vv;
                }
        }
    }

    /**
     * The matrix for weak boundary
     * condition of Nitsche type for
     * the tangential component in
     * Maxwell systems.
     *
     * @f[
     * \int_F \biggl( 2\gamma
     * (u\times n) (v\times n) -
     * (u\times n)(\nu \nabla\times
     * v) - (v\times
     * n)(\nu \nabla\times u)
     * \biggr)
     * @f]
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void nitsche_curl_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim==2) ? 1 : dim;

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Point<dim> &n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<d_max; ++d)
                {
                  const unsigned int d1 = (d+1)%dim;
                  const unsigned int d2 = (d+2)%dim;

                  const double cv = fe.shape_grad_component(i,k,d1)[d2] - fe.shape_grad_component(i,k,d2)[d1];
                  const double cu = fe.shape_grad_component(j,k,d1)[d2] - fe.shape_grad_component(j,k,d2)[d1];
                  const double v= fe.shape_value_component(i,k,d1)*n(d2) - fe.shape_value_component(i,k,d2)*n(d1);
                  const double u= fe.shape_value_component(j,k,d1)*n(d2) - fe.shape_value_component(j,k,d2)*n(d1);

                  M(i,j) += dx*(2.*penalty*u*v - cv*u - cu*v);
                }
        }
    }
    /**
     * The product of two tangential
     * traces,
     * @f[
     * \int_F (u\times n)(v\times n)
     * \, ds.
     * @f]
     *
     * @author Guido Kanschat
     * @date 2011
     */
    template <int dim>
    void tangential_trace_matrix (
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim==2) ? 1 : dim;

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Point<dim> &n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<d_max; ++d)
                {
                  const unsigned int d1 = (d+1)%dim;
                  const unsigned int d2 = (d+2)%dim;

                  const double v= fe.shape_value_component(i,k,d1)*n(d2) - fe.shape_value_component(i,k,d2)*n(d1);
                  const double u= fe.shape_value_component(j,k,d1)*n(d2) - fe.shape_value_component(j,k,d2)*n(d1);

                  M(i,j) += dx*u*v;
                }
        }
    }

    /**
     * The interior penalty fluxes
     * for Maxwell systems.
     *
     * @f[
     * \int_F \biggl( \gamma
     * \{u\times n\}\{v\times n\} -
     * \{u\times n\}\{\nu \nabla\times
     * v\}- \{v\times
     * n\}\{\nu \nabla\times u\}
     * \biggr)\;dx
     * @f]
    *
    * @author Guido Kanschat
    * @date 2011
     */
    template <int dim>
    inline void ip_curl_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      const double pen,
      const double factor1 = 1.,
      const double factor2 = -1.)
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

      const double nu1 = factor1;
      const double nu2 = (factor2 < 0) ? factor1 : factor2;
      const double penalty = .5 * pen * (nu1 + nu2);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim==2) ? 1 : dim;

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Point<dim> &n = fe1.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<d_max; ++d)
                {
                  const unsigned int d1 = (d+1)%dim;
                  const unsigned int d2 = (d+2)%dim;
                  // curl u, curl v
                  const double cv1 = nu1*fe1.shape_grad_component(i,k,d1)[d2] - fe1.shape_grad_component(i,k,d2)[d1];
                  const double cv2 = nu2*fe2.shape_grad_component(i,k,d1)[d2] - fe2.shape_grad_component(i,k,d2)[d1];
                  const double cu1 = nu1*fe1.shape_grad_component(j,k,d1)[d2] - fe1.shape_grad_component(j,k,d2)[d1];
                  const double cu2 = nu2*fe2.shape_grad_component(j,k,d1)[d2] - fe2.shape_grad_component(j,k,d2)[d1];

                  // u x n, v x n
                  const double u1= fe1.shape_value_component(j,k,d1)*n(d2) - fe1.shape_value_component(j,k,d2)*n(d1);
                  const double u2=-fe2.shape_value_component(j,k,d1)*n(d2) + fe2.shape_value_component(j,k,d2)*n(d1);
                  const double v1= fe1.shape_value_component(i,k,d1)*n(d2) - fe1.shape_value_component(i,k,d2)*n(d1);
                  const double v2=-fe2.shape_value_component(i,k,d1)*n(d2) + fe2.shape_value_component(i,k,d2)*n(d1);

                  M11(i,j) += .5*dx*(2.*penalty*u1*v1 - cv1*u1 - cu1*v1);
                  M12(i,j) += .5*dx*(2.*penalty*v1*u2 - cv1*u2 - cu2*v1);
                  M21(i,j) += .5*dx*(2.*penalty*u1*v2 - cv2*u1 - cu1*v2);
                  M22(i,j) += .5*dx*(2.*penalty*u2*v2 - cv2*u2 - cu2*v2);
                }
        }
    }


  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
