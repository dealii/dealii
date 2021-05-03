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

#ifndef dealii_integrators_l2_h
#define dealii_integrators_l2_h


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
   * @brief Local integrators related to <i>L<sup>2</sup></i>-inner products.
   *
   * @ingroup Integrators
   */
  namespace L2
  {
    /**
     * The mass matrix for scalar or vector values finite elements. \f[ \int_Z
     * uv\,dx \quad \text{or} \quad \int_Z \mathbf u\cdot \mathbf v\,dx \f]
     *
     * Likewise, this term can be used on faces, where it computes  the
     * integrals \f[ \int_F uv\,ds \quad \text{or} \quad \int_F \mathbf u\cdot
     * \mathbf v\,ds \f]
     *
     * @param M The mass matrix obtained as result.
     * @param fe The FEValues object describing the local trial function
     * space. #update_values and #update_JxW_values must be set.
     * @param factor A constant that multiplies the mass matrix.
     */
    template <int dim>
    void
    mass_matrix(FullMatrix<double> &     M,
                const FEValuesBase<dim> &fe,
                const double             factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              double Mii = 0.0;
              for (unsigned int d = 0; d < n_components; ++d)
                Mii += dx * fe.shape_value_component(i, k, d) *
                       fe.shape_value_component(i, k, d);

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < n_dofs; ++j)
                {
                  double Mij = 0.0;
                  for (unsigned int d = 0; d < n_components; ++d)
                    Mij += dx * fe.shape_value_component(j, k, d) *
                           fe.shape_value_component(i, k, d);

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }

    /**
     * The weighted mass matrix for scalar or vector values finite elements.
     * \f[ \int_Z \omega(x) uv\,dx \quad \text{or} \quad \int_Z \omega(x)
     * \mathbf u\cdot \mathbf v\,dx \f]
     *
     * Likewise, this term can be used on faces, where it computes  the
     * integrals \f[ \int_F \omega(x) uv\,ds \quad \text{or} \quad \int_F
     * \omega(x) \mathbf u\cdot \mathbf v\,ds \f]
     *
     * @param M The weighted mass matrix obtained as result.
     * @param fe The FEValues object describing the local trial function
     * space. #update_values and #update_JxW_values must be set.
     * @param weights The weights, $\omega(x)$, evaluated at the quadrature
     * points in the finite element (size must be equal to the number of
     * quadrature points in the element).
     */
    template <int dim>
    void
    weighted_mass_matrix(FullMatrix<double> &       M,
                         const FEValuesBase<dim> &  fe,
                         const std::vector<double> &weights)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      AssertDimension(weights.size(), fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * weights[k];
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              double Mii = 0.0;
              for (unsigned int d = 0; d < n_components; ++d)
                Mii += dx * fe.shape_value_component(i, k, d) *
                       fe.shape_value_component(i, k, d);

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < n_dofs; ++j)
                {
                  double Mij = 0.0;
                  for (unsigned int d = 0; d < n_components; ++d)
                    Mij += dx * fe.shape_value_component(j, k, d) *
                           fe.shape_value_component(i, k, d);

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }

    /**
     * <i>L<sup>2</sup></i>-inner product for scalar functions.
     *
     * \f[ \int_Z fv\,dx \quad \text{or} \quad \int_F fv\,ds \f]
     *
     * @param result The vector obtained as result.
     * @param fe The FEValues object describing the local trial function
     * space. #update_values and #update_JxW_values must be set.
     * @param input The representation of $f$ evaluated at the quadrature
     * points in the finite element (size must be equal to the number of
     * quadrature points in the element).
     * @param factor A constant that multiplies the result.
     */
    template <int dim, typename number>
    void
    L2(Vector<number> &           result,
       const FEValuesBase<dim> &  fe,
       const std::vector<double> &input,
       const double               factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(result.size(), n_dofs);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(input.size(), fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < n_dofs; ++i)
          result(i) += fe.JxW(k) * factor * input[k] * fe.shape_value(i, k);
    }

    /**
     * <i>L<sup>2</sup></i>-inner product for a slice of a vector valued right
     * hand side. \f[ \int_Z \mathbf f\cdot \mathbf v\,dx \quad \text{or}
     * \quad \int_F \mathbf f\cdot \mathbf v\,ds \f]
     *
     * @param result The vector obtained as result.
     * @param fe The FEValues object describing the local trial function
     * space. #update_values and #update_JxW_values must be set.
     * @param input The vector valued representation of $\mathbf f$ evaluated
     * at the quadrature points in the finite element (size of each component
     * must be equal to the number of quadrature points in the element).
     * @param factor A constant that multiplies the result.
     */
    template <int dim, typename number>
    void
    L2(Vector<number> &                            result,
       const FEValuesBase<dim> &                   fe,
       const ArrayView<const std::vector<double>> &input,
       const double                                factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = input.size();

      AssertDimension(result.size(), n_dofs);
      AssertDimension(input.size(), fe.get_fe().n_components());

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int d = 0; d < n_components; ++d)
            result(i) += fe.JxW(k) * factor *
                         fe.shape_value_component(i, k, d) * input[d][k];
    }

    /**
     * The jump matrix between two cells for scalar or vector values finite
     * elements. Note that the factor $\gamma$ can be used to implement
     * weighted jumps. \f[ \int_F [\gamma u][\gamma v]\,ds \quad \text{or}
     * \int_F [\gamma \mathbf u]\cdot [\gamma \mathbf v]\,ds \f]
     *
     * Using appropriate weights, this term can be used to penalize violation
     * of conformity in <i>H<sup>1</sup></i>.
     *
     * Note that for the parameters that follow, the external matrix refers to
     * the flux between cells, while the internal matrix refers to entries
     * coupling inside the cell.
     *
     * @param M11 The internal matrix for the first cell obtained as result.
     * @param M12 The external matrix for the first cell obtained as result.
     * @param M21 The external matrix for the second cell obtained as result.
     * @param M22 The internal matrix for the second cell obtained as result.
     * @param fe1 The FEValues object describing the local trial function
     * space for the first cell. #update_values and #update_JxW_values must be
     * set.
     * @param fe2 The FEValues object describing the local trial function
     * space for the second cell. #update_values and #update_JxW_values must be
     * set.
     * @param factor1 A constant that multiplies the shape functions for the
     * first cell.
     * @param factor2 A constant that multiplies the shape functions for the
     * second cell.
     */
    template <int dim>
    void
    jump_matrix(FullMatrix<double> &     M11,
                FullMatrix<double> &     M12,
                FullMatrix<double> &     M21,
                FullMatrix<double> &     M22,
                const FEValuesBase<dim> &fe1,
                const FEValuesBase<dim> &fe2,
                const double             factor1 = 1.,
                const double             factor2 = 1.)
    {
      const unsigned int n1_dofs      = fe1.n_dofs_per_cell();
      const unsigned int n2_dofs      = fe2.n_dofs_per_cell();
      const unsigned int n_components = fe1.get_fe().n_components();

      Assert(n1_dofs == n2_dofs, ExcNotImplemented());
      (void)n2_dofs;
      AssertDimension(n_components, fe2.get_fe().n_components());
      AssertDimension(M11.m(), n1_dofs);
      AssertDimension(M12.m(), n1_dofs);
      AssertDimension(M21.m(), n2_dofs);
      AssertDimension(M22.m(), n2_dofs);
      AssertDimension(M11.n(), n1_dofs);
      AssertDimension(M12.n(), n2_dofs);
      AssertDimension(M21.n(), n1_dofs);
      AssertDimension(M22.n(), n2_dofs);

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);

          for (unsigned int i = 0; i < n1_dofs; ++i)
            for (unsigned int j = 0; j < n1_dofs; ++j)
              for (unsigned int d = 0; d < n_components; ++d)
                {
                  const double u1 =
                    factor1 * fe1.shape_value_component(j, k, d);
                  const double u2 =
                    -factor2 * fe2.shape_value_component(j, k, d);
                  const double v1 =
                    factor1 * fe1.shape_value_component(i, k, d);
                  const double v2 =
                    -factor2 * fe2.shape_value_component(i, k, d);

                  M11(i, j) += dx * u1 * v1;
                  M12(i, j) += dx * u2 * v1;
                  M21(i, j) += dx * u1 * v2;
                  M22(i, j) += dx * u2 * v2;
                }
        }
    }
  } // namespace L2
} // namespace LocalIntegrators

DEAL_II_NAMESPACE_CLOSE

#endif
