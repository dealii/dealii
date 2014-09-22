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

#ifndef __deal2__integrators_laplace_h
#define __deal2__integrators_laplace_h


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
   * @brief Local integrators related to the Laplacian and its DG formulations
   *
   * @ingroup Integrators
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
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      const double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i=0; i<n_dofs; ++i)
            {
              for (unsigned int j=0; j<n_dofs; ++j)
                for (unsigned int d=0; d<n_components; ++d)
                  M(i,j) += dx *
                            (fe.shape_grad_component(j,k,d) * fe.shape_grad_component(i,k,d));
            }
        }
    }

    /**
     * Laplacian residual operator in weak form
     *
     * \f[
     * \int_Z \nu \nabla u \cdot \nabla v \, dx.
     * \f]
     */
    template <int dim>
    inline void
    cell_residual  (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const std::vector<Tensor<1,dim> > &input,
      double factor = 1.)
    {
      const unsigned int nq = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      Assert(input.size() == nq, ExcDimensionMismatch(input.size(), nq));
      Assert(result.size() == n_dofs, ExcDimensionMismatch(result.size(), n_dofs));

      for (unsigned int k=0; k<nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            result(i) += dx * (input[k] * fe.shape_grad(i,k));
        }
    }


    /**
     * Vector-valued Laplacian residual operator in weak form
     *
     * \f[
     * \int_Z \nu \nabla u : \nabla v \, dx.
     * \f]
     */
    template <int dim>
    inline void
    cell_residual  (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &input,
      double factor = 1.)
    {
      const unsigned int nq = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      Assert(result.size() == n_dofs, ExcDimensionMismatch(result.size(), n_dofs));

      for (unsigned int k=0; k<nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int d=0; d<n_comp; ++d)
              {

                result(i) += dx * (input[d][k] * fe.shape_grad_component(i,k,d));
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
      FullMatrix<double> &M,
      const FEValuesBase<dim> &fe,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      Assert (M.m() == n_dofs, ExcDimensionMismatch(M.m(), n_dofs));
      Assert (M.n() == n_dofs, ExcDimensionMismatch(M.n(), n_dofs));

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          const Point<dim> &n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              for (unsigned int d=0; d<n_comp; ++d)
                M(i,j) += dx *
                          (2. * fe.shape_value_component(i,k,d) * penalty * fe.shape_value_component(j,k,d)
                           - (n * fe.shape_grad_component(i,k,d)) * fe.shape_value_component(j,k,d)
                           - (n * fe.shape_grad_component(j,k,d)) * fe.shape_value_component(i,k,d));
        }
    }

    /**
     * Weak boundary condition for the Laplace operator by Nitsche, scalar
     * version, namely on the face <i>F</i> the vector
     * @f[
     * \int_F \Bigl(\gamma (u-g) v - \partial_n u v - (u-g) \partial_n v\Bigr)\;ds.
     * @f]
     *
     * Here, <i>u</i> is the finite element function whose values and
     * gradient are given in the arguments <tt>input</tt> and
     * <tt>Dinput</tt>, respectively. <i>g</i> is the inhomogeneous
     * boundary value in the argument <tt>data</tt>. $\gamma$ is the usual
     * penalty parameter.
     *
     * @author Guido Kanschat
     * @date 2008, 2009, 2010
     */
    template <int dim>
    void nitsche_residual (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const std::vector<double> &input,
      const std::vector<Tensor<1,dim> > &Dinput,
      const std::vector<double> &data,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(input.size(), fe.n_quadrature_points);
      AssertDimension(Dinput.size(), fe.n_quadrature_points);
      AssertDimension(data.size(), fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Point<dim> &n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            {
              const double dnv = fe.shape_grad(i,k) * n;
              const double dnu = Dinput[k] * n;
              const double v= fe.shape_value(i,k);
              const double u= input[k];
              const double g= data[k];

              result(i) += dx*(2.*penalty*(u-g)*v - dnv*(u-g) - dnu*v);
            }
        }
    }

    /**
     * Weak boundary condition for the Laplace operator by Nitsche, vector
     * valued version, namely on the face <i>F</i>
     * the vector
     * @f[
     * \int_F \Bigl(\gamma (\mathbf u- \mathbf g) \cdot \mathbf v
     - \partial_n \mathbf u \cdot \mathbf v
     - (\mathbf u-\mathbf g) \cdot \partial_n \mathbf v\Bigr)\;ds.
     * @f]
     *
     * Here, <i>u</i> is the finite element function whose values and
     * gradient are given in the arguments <tt>input</tt> and
     * <tt>Dinput</tt>, respectively. <i>g</i> is the inhomogeneous
     * boundary value in the argument <tt>data</tt>. $\gamma$ is the usual
     * penalty parameter.
     *
     * @author Guido Kanschat
     * @date 2008, 2009, 2010
     */
    template <int dim>
    void nitsche_residual (
      Vector<double> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<double> > > &input,
      const VectorSlice<const std::vector<std::vector<Tensor<1,dim> > > > &Dinput,
      const VectorSlice<const std::vector<std::vector<double> > > &data,
      double penalty,
      double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();
      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      AssertVectorVectorDimension(Dinput, n_comp, fe.n_quadrature_points);
      AssertVectorVectorDimension(data, n_comp, fe.n_quadrature_points);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          const Point<dim> &n = fe.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int d=0; d<n_comp; ++d)
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
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
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
      const double nu = .5*(nui+nue);

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Point<dim> &n = fe1.normal_vector(k);
          for (unsigned int d=0; d<fe1.get_fe().n_components(); ++d)
            {
              for (unsigned int i=0; i<n_dofs; ++i)
                {
                  for (unsigned int j=0; j<n_dofs; ++j)
                    {
                      const double vi = fe1.shape_value_component(i,k,d);
                      const double dnvi = n * fe1.shape_grad_component(i,k,d);
                      const double ve = fe2.shape_value_component(i,k,d);
                      const double dnve = n * fe2.shape_grad_component(i,k,d);
                      const double ui = fe1.shape_value_component(j,k,d);
                      const double dnui = n * fe1.shape_grad_component(j,k,d);
                      const double ue = fe2.shape_value_component(j,k,d);
                      const double dnue = n * fe2.shape_grad_component(j,k,d);
                      M11(i,j) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+nu*penalty*ui*vi);
                      M12(i,j) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-nu*penalty*vi*ue);
                      M21(i,j) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-nu*penalty*ui*ve);
                      M22(i,j) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+nu*penalty*ue*ve);
                    }
                }
            }
        }
    }

    /**
     * Flux for the interior penalty method for the Laplacian applied
     * to the tangential components of a vector field, namely on
     * the face <i>F</i> the matrices associated with the bilinear form
     * @f[
     * \int_F \Bigl( \gamma [u_\tau][v_\tau] - \{\nabla u_\tau\}[v_\tau\mathbf n] - [u_\tau\mathbf
     * n]\{\nabla v_\tau\} \Bigr) \; ds.
     * @f]
     *
     * @warning This function is still under development!
     *
     * @author Bärbel Janssen, Guido Kanschat
     * @date 2013
     */
    template <int dim>
    void ip_tangential_matrix (
      FullMatrix<double> &M11,
      FullMatrix<double> &M12,
      FullMatrix<double> &M21,
      FullMatrix<double> &M22,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      double penalty,
      double factor1 = 1.,
      double factor2 = -1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;
      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
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
      const double nu = .5*(nui+nue);

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Point<dim> &n = fe1.normal_vector(k);
          for (unsigned int i=0; i<n_dofs; ++i)
            {
              for (unsigned int j=0; j<n_dofs; ++j)
                {
                  double u1dotn = 0.;
                  double v1dotn = 0.;
                  double u2dotn = 0.;
                  double v2dotn = 0.;

                  double ngradu1n = 0.;
                  double ngradv1n = 0.;
                  double ngradu2n = 0.;
                  double ngradv2n = 0.;

                  for (unsigned int d=0; d<dim; ++d)
                    {
                      u1dotn += n(d)*fe1.shape_value_component(j,k,d);
                      v1dotn += n(d)*fe1.shape_value_component(i,k,d);
                      u2dotn += n(d)*fe2.shape_value_component(j,k,d);
                      v2dotn += n(d)*fe2.shape_value_component(i,k,d);

                      ngradu1n += n*fe1.shape_grad_component(j,k,d)*n(d);
                      ngradv1n += n*fe1.shape_grad_component(i,k,d)*n(d);
                      ngradu2n += n*fe2.shape_grad_component(j,k,d)*n(d);
                      ngradv2n += n*fe2.shape_grad_component(i,k,d)*n(d);
                    }

                  for (unsigned int d=0; d<fe1.get_fe().n_components(); ++d)
                    {
                      const double vi = fe1.shape_value_component(i,k,d)-v1dotn*n(d);
                      const double dnvi = n * fe1.shape_grad_component(i,k,d)-ngradv1n*n(d);

                      const double ve = fe2.shape_value_component(i,k,d)-v2dotn*n(d);
                      const double dnve = n * fe2.shape_grad_component(i,k,d)-ngradv2n*n(d);

                      const double ui = fe1.shape_value_component(j,k,d)-u1dotn*n(d);
                      const double dnui = n * fe1.shape_grad_component(j,k,d)-ngradu1n*n(d);

                      const double ue = fe2.shape_value_component(j,k,d)-u2dotn*n(d);
                      const double dnue = n * fe2.shape_grad_component(j,k,d)-ngradu2n*n(d);

                      M11(i,j) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+nu*penalty*ui*vi);
                      M12(i,j) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-nu*penalty*vi*ue);
                      M21(i,j) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-nu*penalty*ui*ve);
                      M22(i,j) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+nu*penalty*ue*ve);
                    }
                }
            }
        }
    }

    /**
     * Residual term for the symmetric interior penalty method:
     * @f[
     * \int_F \Bigl( \gamma [u][v] - \{\nabla u\}[v\mathbf n] - [u\mathbf
     * n]\{\nabla v\} \Bigr) \; ds.
     * @f]
     *
     * @author Guido Kanschat
     * @date 2012
     */
    template<int dim>
    void
    ip_residual(
      Vector<double> &result1,
      Vector<double> &result2,
      const FEValuesBase<dim> &fe1,
      const FEValuesBase<dim> &fe2,
      const std::vector<double> &input1,
      const std::vector<Tensor<1,dim> > &Dinput1,
      const std::vector<double> &input2,
      const std::vector<Tensor<1,dim> > &Dinput2,
      double pen,
      double int_factor = 1.,
      double ext_factor = -1.)
    {
      Assert(fe1.get_fe().n_components() == 1,
             ExcDimensionMismatch(fe1.get_fe().n_components(), 1));
      Assert(fe2.get_fe().n_components() == 1,
             ExcDimensionMismatch(fe2.get_fe().n_components(), 1));

      const double nui = int_factor;
      const double nue = (ext_factor < 0) ? int_factor : ext_factor;
      const double penalty = .5 * pen * (nui + nue);

      const unsigned int n_dofs = fe1.dofs_per_cell;

      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Point<dim> &n = fe1.normal_vector(k);

          for (unsigned int i=0; i<n_dofs; ++i)
            {
              const double vi = fe1.shape_value(i,k);
              const Tensor<1,dim> &Dvi = fe1.shape_grad(i,k);
              const double dnvi = Dvi * n;
              const double ve = fe2.shape_value(i,k);
              const Tensor<1,dim> &Dve = fe2.shape_grad(i,k);
              const double dnve = Dve * n;

              const double ui = input1[k];
              const Tensor<1,dim> &Dui = Dinput1[k];
              const double dnui = Dui * n;
              const double ue = input2[k];
              const Tensor<1,dim> &Due = Dinput2[k];
              const double dnue = Due * n;

              result1(i) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+penalty*ui*vi);
              result1(i) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-penalty*vi*ue);
              result2(i) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-penalty*ui*ve);
              result2(i) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+penalty*ue*ve);
            }
        }
    }


    /**
     * Vector-valued residual term for the symmetric interior penalty method:
     * @f[
     * \int_F \Bigl( \gamma [\mathbf u]\cdot[\mathbf v]
     - \{\nabla \mathbf u\}[\mathbf v\otimes \mathbf n]
     - [\mathbf u\otimes \mathbf n]\{\nabla \mathbf v\} \Bigr) \; ds.
     * @f]
     *
     * @author Guido Kanschat
     * @date 2012
     */
    template<int dim>
    void
    ip_residual(
      Vector<double> &result1,
      Vector<double> &result2,
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
      const unsigned int n_comp = fe1.get_fe().n_components();
      const unsigned int n1 = fe1.dofs_per_cell;

      AssertVectorVectorDimension(input1, n_comp, fe1.n_quadrature_points);
      AssertVectorVectorDimension(Dinput1, n_comp, fe1.n_quadrature_points);
      AssertVectorVectorDimension(input2, n_comp, fe2.n_quadrature_points);
      AssertVectorVectorDimension(Dinput2, n_comp, fe2.n_quadrature_points);

      const double nui = int_factor;
      const double nue = (ext_factor < 0) ? int_factor : ext_factor;
      const double penalty = .5 * pen * (nui + nue);


      for (unsigned int k=0; k<fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);
          const Point<dim> &n = fe1.normal_vector(k);

          for (unsigned int i=0; i<n1; ++i)
            for (unsigned int d=0; d<n_comp; ++d)
              {
                const double vi = fe1.shape_value_component(i,k,d);
                const Tensor<1,dim> &Dvi = fe1.shape_grad_component(i,k,d);
                const double dnvi = Dvi * n;
                const double ve = fe2.shape_value_component(i,k,d);
                const Tensor<1,dim> &Dve = fe2.shape_grad_component(i,k,d);
                const double dnve = Dve * n;

                const double ui = input1[d][k];
                const Tensor<1,dim> &Dui = Dinput1[d][k];
                const double dnui = Dui * n;
                const double ue = input2[d][k];
                const Tensor<1,dim> &Due = Dinput2[d][k];
                const double dnue = Due * n;

                result1(i) += dx*(-.5*nui*dnvi*ui-.5*nui*dnui*vi+penalty*ui*vi);
                result1(i) += dx*( .5*nui*dnvi*ue-.5*nue*dnue*vi-penalty*vi*ue);
                result2(i) += dx*(-.5*nue*dnve*ui+.5*nui*dnui*ve-penalty*ui*ve);
                result2(i) += dx*( .5*nue*dnve*ue+.5*nue*dnue*ve+penalty*ue*ve);
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
    template <int dim, int spacedim, typename number>
    double compute_penalty(
      const MeshWorker::DoFInfo<dim,spacedim,number> &dinfo1,
      const MeshWorker::DoFInfo<dim,spacedim,number> &dinfo2,
      unsigned int deg1,
      unsigned int deg2)
    {
      const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[dinfo1.face_number];
      const unsigned int normal2 = GeometryInfo<dim>::unit_normal_direction[dinfo2.face_number];
      const unsigned int deg1sq = (deg1 == 0) ? 1 : deg1 * (deg1+1);
      const unsigned int deg2sq = (deg2 == 0) ? 1 : deg2 * (deg2+1);

      double penalty1 = deg1sq / dinfo1.cell->extent_in_direction(normal1);
      double penalty2 = deg2sq / dinfo2.cell->extent_in_direction(normal2);
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
