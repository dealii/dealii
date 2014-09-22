// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN


template<int dim>
PolynomialsNedelec<dim>::PolynomialsNedelec (const unsigned int k) :
  my_degree (k), polynomial_space (create_polynomials (k)),
  n_pols (compute_n_pols (k))
{
}

template<int dim>
std::vector<std::vector< Polynomials::Polynomial<double> > >
PolynomialsNedelec<dim>::create_polynomials (const unsigned int k)
{
  std::vector<std::vector< Polynomials::Polynomial<double> > > pols (dim);

  pols[0] = Polynomials::Legendre::generate_complete_basis (k);

  for (unsigned int i = 1; i < dim; ++i)
    pols[i] = Polynomials::Lobatto::generate_complete_basis (k + 1);

  return pols;
}


// Compute the values, gradients
// and double gradients of the
// polynomial at the given point.
template<int dim>
void
PolynomialsNedelec<dim>::compute (const Point<dim> &unit_point,
                                  std::vector<Tensor<1,dim> > &values,
                                  std::vector<Tensor<2,dim> > &grads,
                                  std::vector<Tensor<3,dim> > &grad_grads)
const
{
  Assert(values.size () == n_pols || values.size () == 0,
         ExcDimensionMismatch(values.size (), n_pols));
  Assert(grads.size () == n_pols || grads.size () == 0,
         ExcDimensionMismatch(grads.size (), n_pols));
  Assert(grad_grads.size () == n_pols || grad_grads.size () == 0,
         ExcDimensionMismatch(grad_grads.size (), n_pols));

  // Declare the values, derivatives
  // and second derivatives vectors of
  // <tt>polynomial_space</tt> at
  // <tt>unit_point</tt>
  const unsigned int &n_basis = polynomial_space.n ();
  std::vector<double> unit_point_values ((values.size () == 0) ? 0 : n_basis);
  std::vector<Tensor<1, dim> >
  unit_point_grads ((grads.size () == 0) ? 0 : n_basis);
  std::vector<Tensor<2, dim> >
  unit_point_grad_grads ((grad_grads.size () == 0) ? 0 : n_basis);

  switch (dim)
    {
    case 1:
    {
      polynomial_space.compute (unit_point, unit_point_values,
                                unit_point_grads, unit_point_grad_grads);

      // Assign the correct values to the
      // corresponding shape functions.
      if (values.size () > 0)
        for (unsigned int i = 0; i < unit_point_values.size (); ++i)
          values[i][0] = unit_point_values[i];

      if (grads.size () > 0)
        for (unsigned int i = 0; i < unit_point_grads.size (); ++i)
          grads[i][0][0] = unit_point_grads[i][0];

      if (grad_grads.size () > 0)
        for (unsigned int i = 0; i < unit_point_grad_grads.size (); ++i)
          grad_grads[i][0][0][0] = unit_point_grad_grads[i][0][0];

      break;
    }

    case 2:
    {
      polynomial_space.compute (unit_point, unit_point_values,
                                unit_point_grads, unit_point_grad_grads);

      // Declare the values, derivatives and
      // second derivatives vectors of
      // <tt>polynomial_space</tt> at
      // <tt>unit_point</tt> with coordinates
      // shifted one step in positive direction
      Point<dim> p;

      p (0) = unit_point (1);
      p (1) = unit_point (0);

      std::vector<double> p_values ((values.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<1, dim> >
      p_grads ((grads.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<2, dim> >
      p_grad_grads ((grad_grads.size () == 0) ? 0 : n_basis);

      polynomial_space.compute (p, p_values, p_grads, p_grad_grads);

      // Assign the correct values to the
      // corresponding shape functions.
      if (values.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                values[i + j * (my_degree + 1)][0] = 0.0;
                values[i + j * (my_degree + 1)][1]
                  = p_values[i + j * (my_degree + 1)];
                values[i + (j + 2) * (my_degree + 1)][0]
                  = unit_point_values[i + j * (my_degree + 1)];
                values[i + (j + 2) * (my_degree + 1)][1] = 0.0;
              }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  values[(i + GeometryInfo<dim>::lines_per_cell) * my_degree
                         + j + GeometryInfo<dim>::lines_per_cell][0]
                    = unit_point_values[i + (j + 2) * (my_degree + 1)];
                  values[(i + GeometryInfo<dim>::lines_per_cell) * my_degree
                         + j + GeometryInfo<dim>::lines_per_cell][1] = 0.0;
                  values[i + (j + my_degree
                              + GeometryInfo<dim>::lines_per_cell)
                         * (my_degree + 1)][0] = 0.0;
                  values[i + (j + my_degree
                              + GeometryInfo<dim>::lines_per_cell)
                         * (my_degree + 1)][1]
                    = p_values[i + (j + 2) * (my_degree + 1)];
                }
        }

      if (grads.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                for (unsigned int k = 0; k < dim; ++k)
                  {
                    grads[i + j * (my_degree + 1)][0][k] = 0.0;
                    grads[i + (j + 2) * (my_degree + 1)][0][k]
                      = unit_point_grads[i + j * (my_degree + 1)][k];
                    grads[i + (j + 2) * (my_degree + 1)][1][k] = 0.0;
                  }

                grads[i + j * (my_degree + 1)][1][0]
                  = p_grads[i + j * (my_degree + 1)][1];
                grads[i + j * (my_degree + 1)][1][1]
                  = p_grads[i + j * (my_degree + 1)][0];
              }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      grads[(i + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][0][k]
                        = unit_point_grads[i + (j + 2) * (my_degree + 1)]
                          [k];
                      grads[(i + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][1][k]
                        = 0.0;
                      grads[i + (j + my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][0][k] = 0.0;
                    }

                  grads[i + (j + my_degree
                             + GeometryInfo<dim>::lines_per_cell)
                        * (my_degree + 1)][1][0]
                    = p_grads[i + (j + 2) * (my_degree + 1)][1];
                  grads[i + (j + my_degree
                             + GeometryInfo<dim>::lines_per_cell)
                        * (my_degree + 1)][1][1]
                    = p_grads[i + (j + 2) * (my_degree + 1)][0];
                }
        }

      if (grad_grads.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                for (unsigned int k = 0; k < dim; ++k)
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      grad_grads[i + j * (my_degree + 1)][0][k][l] = 0.0;
                      grad_grads[i + (j + 2) * (my_degree + 1)][0][k][l]
                        = unit_point_grad_grads[i + j * (my_degree + 1)][k]
                          [l];
                      grad_grads[i + (j + 2) * (my_degree + 1)][1][k][l]
                        = 0.0;
                    }

                grad_grads[i + j * (my_degree + 1)][1][0][0]
                  = p_grad_grads[i + j * (my_degree + 1)][1][1];
                grad_grads[i + j * (my_degree + 1)][1][0][1]
                  = p_grad_grads[i + j * (my_degree + 1)][1][0];
                grad_grads[i + j * (my_degree + 1)][1][1][0]
                  = p_grad_grads[i + j * (my_degree + 1)][0][1];
                grad_grads[i + j * (my_degree + 1)][1][1][1]
                  = p_grad_grads[i + j * (my_degree + 1)][0][0];
              }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      {
                        grad_grads[(i + GeometryInfo<dim>::lines_per_cell)
                                   * my_degree + j
                                   + GeometryInfo<dim>::lines_per_cell][0]
                        [k][l]
                          = unit_point_grad_grads[i + (j + 2)
                                                  * (my_degree + 1)][k][l];
                        grad_grads[(i + GeometryInfo<dim>::lines_per_cell)
                                   * my_degree + j
                                   + GeometryInfo<dim>::lines_per_cell][1]
                        [k][l] = 0.0;
                        grad_grads[i + (j + my_degree
                                        + GeometryInfo<dim>::lines_per_cell)
                                   * (my_degree + 1)][0][k][l] = 0.0;
                      }

                  grad_grads[i + (j + my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][1][0][0]
                    = p_grad_grads[i + (j + 2) * (my_degree + 1)][1][1];
                  grad_grads[i + (j + my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][1][0][1]
                    = p_grad_grads[i + (j + 2) * (my_degree + 1)][1][0];
                  grad_grads[i + (j + my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][1][1][0]
                    = p_grad_grads[i + (j + 2) * (my_degree + 1)][0][1];
                  grad_grads[i + (j + my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][1][1][1]
                    = p_grad_grads[i + (j + 2) * (my_degree + 1)][0][0];
                }
        }

      break;
    }

    case 3:
    {
      polynomial_space.compute (unit_point, unit_point_values,
                                unit_point_grads, unit_point_grad_grads);

      // Declare the values, derivatives
      // and second derivatives vectors of
      // <tt>polynomial_space</tt> at
      // <tt>unit_point</tt> with coordinates
      // shifted two steps in positive
      // direction
      Point<dim> p1, p2;
      std::vector<double> p1_values ((values.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<1, dim> >
      p1_grads ((grads.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<2, dim> >
      p1_grad_grads ((grad_grads.size () == 0) ? 0 : n_basis);
      std::vector<double> p2_values ((values.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<1, dim> >
      p2_grads ((grads.size () == 0) ? 0 : n_basis);
      std::vector<Tensor<2, dim> >
      p2_grad_grads ((grad_grads.size () == 0) ? 0 : n_basis);

      p1 (0) = unit_point (1);
      p1 (1) = unit_point (2);
      p1 (2) = unit_point (0);
      polynomial_space.compute (p1, p1_values, p1_grads, p1_grad_grads);
      p2 (0) = unit_point (2);
      p2 (1) = unit_point (0);
      p2 (2) = unit_point (1);
      polynomial_space.compute (p2, p2_values, p2_grads, p2_grad_grads);

      // Assign the correct values to the
      // corresponding shape functions.
      if (values.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            {
              for (unsigned int j = 0; j < 2; ++j)
                {
                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < 2; ++l)
                        {
                          values[i + (j + 4 * k) * (my_degree + 1)][2 * l]
                            = 0.0;
                          values[i + (j + 4 * k + 2) * (my_degree + 1)]
                          [l + 1] = 0.0;
                          values[i + (j + 2 * (k + 4)) * (my_degree + 1)][l]
                            = 0.0;
                        }

                      values[i + (j + 4 * k + 2) * (my_degree + 1)][0]
                        = unit_point_values[i + (j + k * (my_degree + 2))
                                            * (my_degree + 1)];
                      values[i + (j + 2 * (k + 4)) * (my_degree + 1)][2]
                        = p2_values[i + (j + k * (my_degree + 2))
                                    * (my_degree + 1)];
                    }

                  values[i + j * (my_degree + 1)][1]
                    = p1_values[i + j * (my_degree + 1) * (my_degree + 2)];
                }

              values[i + 4 * (my_degree + 1)][1]
                = p1_values[i + my_degree + 1];
              values[i + 5 * (my_degree + 1)][1]
                = p1_values[i + (my_degree + 1) * (my_degree + 3)];
            }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  for (unsigned int k = 0; k < my_degree; ++k)
                    {
                      for (unsigned int l = 0; l < 2; ++l)
                        {
                          values[((i + 2
                                   * GeometryInfo<dim>::faces_per_cell)
                                  * my_degree + j
                                  + GeometryInfo<dim>::lines_per_cell + 2
                                  * GeometryInfo<dim>::faces_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][l + 1]
                            = 0.0;
                          values[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree)
                                  * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][2 * l]
                            = 0.0;
                          values[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree))
                                      * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][l] = 0.0;
                        }

                      values[((i + 2 * GeometryInfo<dim>::faces_per_cell)
                              * my_degree + j
                              + GeometryInfo<dim>::lines_per_cell + 2
                              * GeometryInfo<dim>::faces_per_cell)
                             * my_degree + k
                             + GeometryInfo<dim>::lines_per_cell][0]
                        = unit_point_values[i + (j + (k + 2) * (my_degree
                                                                + 2) + 2)
                                            * (my_degree + 1)];
                      values[(i + (j + 2 * GeometryInfo<dim>::faces_per_cell
                                   + my_degree) * (my_degree + 1)
                              + GeometryInfo<dim>::lines_per_cell)
                             * my_degree + k
                             + GeometryInfo<dim>::lines_per_cell][1]
                        = p1_values[i + ((j + 2) * (my_degree + 2) + k + 2)
                                    * (my_degree + 1)];
                      values[i + (j + (k + 2
                                       * (GeometryInfo<dim>::faces_per_cell
                                          + my_degree))
                                  * my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][2]
                        = p2_values[i + (j + (k + 2) * (my_degree + 2) + 2)
                                    * (my_degree + 1)];
                    }

                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < 2; ++l)
                        {
                          for (unsigned int m = 0; m < 2; ++m)
                            {
                              values[i + (j + (2 * (k + 2 * l) + 1)
                                          * my_degree
                                          + GeometryInfo<dim>::lines_per_cell)
                                     * (my_degree + 1)][m + l] = 0.0;
                              values[(i + 2 * (k + 2 * (l + 1)) * (my_degree
                                                                   + 1)
                                      + GeometryInfo<dim>::lines_per_cell)
                                     * my_degree + j
                                     + GeometryInfo<dim>::lines_per_cell]
                              [m + l] = 0.0;
                            }

                          values[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2 * l]
                            = 0.0;
                          values[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2 * l] = 0.0;
                        }

                      values[(i + 2 * k * (my_degree + 1)
                              + GeometryInfo<dim>::lines_per_cell)
                             * my_degree + j
                             + GeometryInfo<dim>::lines_per_cell][1]
                        = p1_values[i + (j + k * (my_degree + 2) + 2)
                                    * (my_degree + 1)];
                      values[i + (j + (2 * k + 1) * my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][2]
                        = p2_values[i + ((j + 2) * (my_degree + 2) + k)
                                    * (my_degree + 1)];
                      values[(i + 2 * (k + 2) * (my_degree + 1)
                              + GeometryInfo<dim>::lines_per_cell)
                             * my_degree + j
                             + GeometryInfo<dim>::lines_per_cell][2]
                        = p2_values[i + (j + k * (my_degree + 2) + 2)
                                    * (my_degree + 1)];
                      values[i + (j + (2 * k + 5) * my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][0]
                        = unit_point_values[i + ((j + 2) * (my_degree + 2) + k)
                                            * (my_degree + 1)];
                      values[(i + 2 * (k + 4) * (my_degree + 1)
                              + GeometryInfo<dim>::lines_per_cell)
                             * my_degree + j
                             + GeometryInfo<dim>::lines_per_cell][0]
                        = unit_point_values[i + (j + k * (my_degree + 2)
                                                 + 2) * (my_degree + 1)];
                      values[i + (j + (2 * k + 9) * my_degree
                                  + GeometryInfo<dim>::lines_per_cell)
                             * (my_degree + 1)][1]
                        = p1_values[i + ((j + 2) * (my_degree + 2) + k)
                                    * (my_degree + 1)];
                    }
                }
        }

      if (grads.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            {
              for (unsigned int j = 0; j < 2; ++j)
                {
                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < 2; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            grads[i + (j + 4 * k) * (my_degree + 1)][2 * l]
                            [m] = 0.0;
                            grads[i + (j + 4 * k + 2) * (my_degree + 1)]
                            [l + 1][m] = 0.0;
                            grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                            [l][m] = 0.0;
                          }

                      for (unsigned int l = 0; l < dim; ++l)
                        grads[i + (j + 4 * k + 2) * (my_degree + 1)][0][l]
                          = unit_point_grads[i + (j + k * (my_degree + 2))
                                             * (my_degree + 1)][l];

                      grads[i + (j + 2 * (k + 4)) * (my_degree + 1)][2][0]
                        = p2_grads[i + (j + k * (my_degree + 2))
                                   * (my_degree + 1)][1];
                      grads[i + (j + 2 * (k + 4)) * (my_degree + 1)][2][1]
                        = p2_grads[i + (j + k * (my_degree + 2))
                                   * (my_degree + 1)][2];
                      grads[i + (j + 2 * (k + 4)) * (my_degree + 1)][2][2]
                        = p2_grads[i + (j + k * (my_degree + 2))
                                   * (my_degree + 1)][0];
                    }

                  grads[i + j * (my_degree + 1)][1][0]
                    = p1_grads[i + j * (my_degree + 1) * (my_degree + 2)]
                      [2];
                  grads[i + j * (my_degree + 1)][1][1]
                    = p1_grads[i + j * (my_degree + 1) * (my_degree + 2)]
                      [0];
                  grads[i + j * (my_degree + 1)][1][2]
                    = p1_grads[i + j * (my_degree + 1) * (my_degree + 2)]
                      [1];
                }

              grads[i + 4 * (my_degree + 1)][1][0]
                = p1_grads[i + my_degree + 1][2];
              grads[i + 4 * (my_degree + 1)][1][1]
                = p1_grads[i + my_degree + 1][0];
              grads[i + 4 * (my_degree + 1)][1][2]
                = p1_grads[i + my_degree + 1][1];
              grads[i + 5 * (my_degree + 1)][1][0]
                = p1_grads[i + (my_degree + 1) * (my_degree + 3)][2];
              grads[i + 5 * (my_degree + 1)][1][1]
                = p1_grads[i + (my_degree + 1) * (my_degree + 3)][0];
              grads[i + 5 * (my_degree + 1)][1][2]
                = p1_grads[i + (my_degree + 1) * (my_degree + 3)][1];
            }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  for (unsigned int k = 0; k < my_degree; ++k)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        {
                          for (unsigned int m = 0; m < 2; ++m)
                            {
                              grads[((i + 2
                                      * GeometryInfo<dim>::faces_per_cell)
                                     * my_degree + j
                                     + GeometryInfo<dim>::lines_per_cell
                                     + 2
                                     * GeometryInfo<dim>::faces_per_cell)
                                    * my_degree + k
                                    + GeometryInfo<dim>::lines_per_cell]
                              [m + 1][l] = 0.0;
                              grads[(i + (j + 2
                                          * GeometryInfo<dim>::faces_per_cell
                                          + my_degree) * (my_degree + 1)
                                     + GeometryInfo<dim>::lines_per_cell)
                                    * my_degree + k
                                    + GeometryInfo<dim>::lines_per_cell]
                              [2 * m][l] = 0.0;
                              grads[i + (j + (k + 2
                                              * (GeometryInfo<dim>::faces_per_cell
                                                 + my_degree)) * my_degree
                                         + GeometryInfo<dim>::lines_per_cell)
                                    * (my_degree + 1)][m][l] = 0.0;
                            }

                          grads[((i + 2 * GeometryInfo<dim>::faces_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell + 2
                                 * GeometryInfo<dim>::faces_per_cell)
                                * my_degree + k
                                + GeometryInfo<dim>::lines_per_cell][0][l]
                            = unit_point_grads[i + (j + (k + 2)
                                                    * (my_degree + 2) + 2)
                                               * (my_degree + 1)][l];
                        }

                      grads[(i + (j + 2 * GeometryInfo<dim>::faces_per_cell
                                  + my_degree) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + k
                            + GeometryInfo<dim>::lines_per_cell][1][0]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k + 2)
                                   * (my_degree + 1)][2];
                      grads[(i + (j + 2 * GeometryInfo<dim>::faces_per_cell
                                  + my_degree) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + k
                            + GeometryInfo<dim>::lines_per_cell][1][1]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k + 2)
                                   * (my_degree + 1)][0];
                      grads[(i + (j + 2 * GeometryInfo<dim>::faces_per_cell
                                  + my_degree) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + k
                            + GeometryInfo<dim>::lines_per_cell][1][2]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k + 2)
                                   * (my_degree + 1)][1];
                      grads[i + (j + (k + 2
                                      * (GeometryInfo<dim>::faces_per_cell
                                         + my_degree)) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][0]
                        = p2_grads[i + (j + (k + 2) * (my_degree + 2) + 2)
                                   * (my_degree + 1)][1];
                      grads[i + (j + (k + 2
                                      * (GeometryInfo<dim>::faces_per_cell
                                         + my_degree)) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][1]
                        = p2_grads[i + (j + (k + 2) * (my_degree + 2) + 2)
                                   * (my_degree + 1)][2];
                      grads[i + (j + (k + 2
                                      * (GeometryInfo<dim>::faces_per_cell
                                         + my_degree))
                                 * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][2]
                        = p2_grads[i + (j + (k + 2) * (my_degree + 2) + 2)
                                   * (my_degree + 1)][0];
                    }

                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < 2; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            for (unsigned int n = 0; n < 2; ++n)
                              {
                                grads[i + (j + (2 * (k + 2 * l) + 1)
                                           * my_degree
                                           + GeometryInfo<dim>::lines_per_cell)
                                      * (my_degree + 1)][n + l][m] = 0.0;
                                grads[(i + 2 * (k + 2 * (l + 1))
                                       * (my_degree + 1)
                                       + GeometryInfo<dim>::lines_per_cell)
                                      * my_degree + j
                                      + GeometryInfo<dim>::lines_per_cell]
                                [n + l][m] = 0.0;
                              }

                            grads[(i + 2 * k * (my_degree + 1)
                                   + GeometryInfo<dim>::lines_per_cell)
                                  * my_degree + j
                                  + GeometryInfo<dim>::lines_per_cell]
                            [2 * l][m] = 0.0;
                            grads[i + (j + (2 * k + 9) * my_degree
                                       + GeometryInfo<dim>::lines_per_cell)
                                  * (my_degree + 1)][2 * l][m] = 0.0;
                          }

                      for (unsigned int l = 0; l < dim; ++l)
                        {
                          grads[i + (j + (2 * k + 5) * my_degree
                                     + GeometryInfo<dim>::lines_per_cell)
                                * (my_degree + 1)][0][l]
                            = unit_point_grads[i + ((j + 2) * (my_degree
                                                               + 2) + k)
                                               * (my_degree + 1)][l];
                          grads[(i + 2 * (k + 4) * (my_degree + 1)
                                 + GeometryInfo<dim>::lines_per_cell)
                                * my_degree + j
                                + GeometryInfo<dim>::lines_per_cell][0][l]
                            = unit_point_grads[i + (j + k * (my_degree + 2)
                                                    + 2) * (my_degree + 1)]
                              [l];
                        }

                      grads[(i + 2 * k * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][1][0]
                        = p1_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][2];
                      grads[(i + 2 * k * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][1][1]
                        = p1_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][0];
                      grads[(i + 2 * k * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][1][2]
                        = p1_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][1];
                      grads[i + (j + (2 * k + 1) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][0]
                        = p2_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][1];
                      grads[i + (j + (2 * k + 1) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][1]
                        = p2_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][2];
                      grads[i + (j + (2 * k + 1) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][2][2]
                        = p2_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][0];
                      grads[(i + 2 * (k + 2) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][2][0]
                        = p2_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][1];
                      grads[(i + 2 * (k + 2) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][2][1]
                        = p2_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][2];
                      grads[(i + 2 * (k + 2) * (my_degree + 1)
                             + GeometryInfo<dim>::lines_per_cell)
                            * my_degree + j
                            + GeometryInfo<dim>::lines_per_cell][2][2]
                        = p2_grads[i + (j + k * (my_degree + 2) + 2)
                                   * (my_degree + 1)][0];
                      grads[i + (j + (2 * k + 9) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][1][0]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][2];
                      grads[i + (j + (2 * k + 9) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][1][1]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][0];
                      grads[i + (j + (2 * k + 9) * my_degree
                                 + GeometryInfo<dim>::lines_per_cell)
                            * (my_degree + 1)][1][2]
                        = p1_grads[i + ((j + 2) * (my_degree + 2) + k)
                                   * (my_degree + 1)][1];
                    }
                }
        }

      if (grad_grads.size () > 0)
        {
          for (unsigned int i = 0; i <= my_degree; ++i)
            {
              for (unsigned int j = 0; j < 2; ++j)
                {
                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            for (unsigned int n = 0; n < 2; ++n)
                              {
                                grad_grads[i + (j + 4 * k) * (my_degree
                                                              + 1)][2 * n]
                                [l][m] = 0.0;
                                grad_grads[i + (j + 4 * k + 2) * (my_degree
                                                                  + 1)]
                                [n + 1][l][m] = 0.0;
                                grad_grads[i + (j + 2 * (k + 4))
                                           * (my_degree + 1)][n][l][m]
                                  = 0.0;
                              }

                            grad_grads[i + (j + 4 * k + 2) * (my_degree
                                                              + 1)][0][l][m]
                              = unit_point_grad_grads[i + (j + k
                                                           * (my_degree
                                                              + 2))
                                                      * (my_degree + 1)][l]
                                [m];
                          }

                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][0][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][1][1];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][0][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][1][2];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][0][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][1][0];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][1][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][2][1];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][1][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][2][2];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][1][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][2][0];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][2][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][0][1];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][2][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][0][2];
                      grad_grads[i + (j + 2 * (k + 4)) * (my_degree + 1)]
                      [2][2][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2))
                                        * (my_degree + 1)][0][0];
                    }

                  grad_grads[i + j * (my_degree + 1)][1][0][0]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][2][2];
                  grad_grads[i + j * (my_degree + 1)][1][0][1]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][2][0];
                  grad_grads[i + j * (my_degree + 1)][1][0][2]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][2][1];
                  grad_grads[i + j * (my_degree + 1)][1][1][0]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][0][2];
                  grad_grads[i + j * (my_degree + 1)][1][1][1]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][0][0];
                  grad_grads[i + j * (my_degree + 1)][1][1][2]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][0][1];
                  grad_grads[i + j * (my_degree + 1)][1][2][0]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][1][2];
                  grad_grads[i + j * (my_degree + 1)][1][2][1]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][1][0];
                  grad_grads[i + j * (my_degree + 1)][1][2][2]
                    = p1_grad_grads[i + j * (my_degree + 1)
                                    * (my_degree + 2)][1][1];
                }

              grad_grads[i + 4 * (my_degree + 1)][1][0][0]
                = p1_grad_grads[i + my_degree + 1][2][2];
              grad_grads[i + 4 * (my_degree + 1)][1][0][1]
                = p1_grad_grads[i + my_degree + 1][2][0];
              grad_grads[i + 4 * (my_degree + 1)][1][0][2]
                = p1_grad_grads[i + my_degree + 1][2][1];
              grad_grads[i + 4 * (my_degree + 1)][1][1][0]
                = p1_grad_grads[i + my_degree + 1][0][2];
              grad_grads[i + 4 * (my_degree + 1)][1][1][1]
                = p1_grad_grads[i + my_degree + 1][0][0];
              grad_grads[i + 4 * (my_degree + 1)][1][1][2]
                = p1_grad_grads[i + my_degree + 1][0][1];
              grad_grads[i + 4 * (my_degree + 1)][1][2][0]
                = p1_grad_grads[i + my_degree + 1][1][2];
              grad_grads[i + 4 * (my_degree + 1)][1][2][1]
                = p1_grad_grads[i + my_degree + 1][1][0];
              grad_grads[i + 4 * (my_degree + 1)][1][2][2]
                = p1_grad_grads[i + my_degree + 1][1][1];
              grad_grads[i + 5 * (my_degree + 1)][1][0][0]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][2]
                  [2];
              grad_grads[i + 5 * (my_degree + 1)][1][0][1]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][2]
                  [0];
              grad_grads[i + 5 * (my_degree + 1)][1][0][2]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][2]
                  [1];
              grad_grads[i + 5 * (my_degree + 1)][1][1][0]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][0]
                  [2];
              grad_grads[i + 5 * (my_degree + 1)][1][1][1]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][0]
                  [0];
              grad_grads[i + 5 * (my_degree + 1)][1][1][2]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][0]
                  [1];
              grad_grads[i + 5 * (my_degree + 1)][1][2][0]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][1]
                  [2];
              grad_grads[i + 5 * (my_degree + 1)][1][2][1]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][1]
                  [0];
              grad_grads[i + 5 * (my_degree + 1)][1][2][2]
                = p1_grad_grads[i + (my_degree + 1) * (my_degree + 3)][1]
                  [1];
            }

          if (my_degree > 0)
            for (unsigned int i = 0; i <= my_degree; ++i)
              for (unsigned int j = 0; j < my_degree; ++j)
                {
                  for (unsigned int k = 0; k < my_degree; ++k)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            for (unsigned int n = 0; n < 2; ++n)
                              {
                                grad_grads[((i + 2
                                             * GeometryInfo<dim>::faces_per_cell)
                                            * my_degree + j
                                            + GeometryInfo<dim>::lines_per_cell
                                            + 2
                                            * GeometryInfo<dim>::faces_per_cell)
                                           * my_degree + k
                                           + GeometryInfo<dim>::lines_per_cell]
                                [n + 1][l][m] = 0.0;
                                grad_grads[(i + (j + 2
                                                 * GeometryInfo<dim>::faces_per_cell
                                                 + my_degree) * (my_degree
                                                                 + 1)
                                            + GeometryInfo<dim>::lines_per_cell)
                                           * my_degree + k
                                           + GeometryInfo<dim>::lines_per_cell]
                                [2 * n][l][m] = 0.0;
                                grad_grads[i + (j + (k + 2
                                                     * (GeometryInfo<dim>::faces_per_cell
                                                        + my_degree))
                                                * my_degree
                                                + GeometryInfo<dim>::lines_per_cell)
                                           * (my_degree + 1)][n][l][m]
                                  = 0.0;
                              }

                            grad_grads[((i + 2
                                         * GeometryInfo<dim>::faces_per_cell)
                                        * my_degree + j
                                        + GeometryInfo<dim>::lines_per_cell
                                        + 2
                                        * GeometryInfo<dim>::faces_per_cell)
                                       * my_degree + k
                                       + GeometryInfo<dim>::lines_per_cell]
                            [0][l][m]
                              = unit_point_grad_grads[i + (j + (k + 2)
                                                           * (my_degree + 2)
                                                           + 2) * (my_degree
                                                                   + 1)][l][m];
                          }

                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree)
                                  * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][2][2];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][2][0];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][2][1];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][0][2];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][0][0];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][0][1];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][1][2];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][1][0];
                      grad_grads[(i + (j + 2
                                       * GeometryInfo<dim>::faces_per_cell
                                       + my_degree) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + k
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k
                                             + 2) * (my_degree + 1)][1][1];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][0]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][1][1];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][1]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][1][2];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree))
                                      * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][2]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][1][0];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree))
                                      * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][0]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][2][1];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][1]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][2][2];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][2]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][2][0];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][0]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][0][1];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][1]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][0][2];
                      grad_grads[i + (j + (k + 2
                                           * (GeometryInfo<dim>::faces_per_cell
                                              + my_degree)) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][2]
                        = p2_grad_grads[i + (j + (k + 2) * (my_degree + 2)
                                             + 2) * (my_degree + 1)][0][0];
                    }

                  for (unsigned int k = 0; k < 2; ++k)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            for (unsigned int n = 0; n < 2; ++n)
                              {
                                for (unsigned int o = 0; o < 2; ++o)
                                  {
                                    grad_grads[i + (j + (2 * (k + 2 * n)
                                                         + 1) * my_degree
                                                    + GeometryInfo<dim>::lines_per_cell)
                                               * (my_degree + 1)][o + n][l][m]
                                      = 0.0;
                                    grad_grads[(i + 2 * (k + 2 * (n + 1))
                                                * (my_degree + 1)
                                                + GeometryInfo<dim>::lines_per_cell)
                                               * my_degree + j
                                               + GeometryInfo<dim>::lines_per_cell]
                                    [o + k][l][m] = 0.0;
                                  }

                                grad_grads[(i + 2 * k * (my_degree + 1)
                                            + GeometryInfo<dim>::lines_per_cell)
                                           * my_degree + j
                                           + GeometryInfo<dim>::lines_per_cell]
                                [2 * n][l][m] = 0.0;
                                grad_grads[i + (j + (2 * k + 9)
                                                * my_degree
                                                + GeometryInfo<dim>::lines_per_cell)
                                           * (my_degree + 1)][2 * n][l][m]
                                  = 0.0;
                              }

                            grad_grads[i + (j + (2 * k + 5) * my_degree
                                            + GeometryInfo<dim>::lines_per_cell)
                                       * (my_degree + 1)]
                            [0][l][m]
                              = unit_point_grad_grads[i + ((j + 2)
                                                           * (my_degree
                                                              + 2) + k)
                                                      * (my_degree + 1)][l]
                                [m];
                            grad_grads[(i + 2 * (k + 4) * (my_degree + 1)
                                        + GeometryInfo<dim>::lines_per_cell)
                                       * my_degree + j
                                       + GeometryInfo<dim>::lines_per_cell]
                            [0][l][m]
                              = unit_point_grad_grads[i + (j + k
                                                           * (my_degree
                                                              + 2) + 2)
                                                      * (my_degree + 1)][l]
                                [m];
                          }

                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [0]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][2];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [1]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][0];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][0]
                      [2]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][1];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [0]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][2];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [1]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][0];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][1]
                      [2]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][1];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [0]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][2];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [1]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][0];
                      grad_grads[(i + 2 * k * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][1][2]
                      [2]
                        = p1_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][1];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][0]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][1];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][1]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][2];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][0][2]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][0];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][0]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][1];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][1]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][2];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][1][2]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][0];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][0]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][1];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][1]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][2];
                      grad_grads[i + (j + (2 * k + 1) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][2][2][2]
                        = p2_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][0];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][0][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][1];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][0][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][2];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][0][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][1][0];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][1][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][1];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][1][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][2];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][1][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][2][0];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][2][0]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][1];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][2][1]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][2];
                      grad_grads[(i + 2 * (k + 2) * (my_degree + 1)
                                  + GeometryInfo<dim>::lines_per_cell)
                                 * my_degree + j
                                 + GeometryInfo<dim>::lines_per_cell][2][2][2]
                        = p2_grad_grads[i + (j + k * (my_degree + 2) + 2)
                                        * (my_degree + 1)][0][0];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][0][0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][2];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][0][1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][0];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][0][2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][2][1];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][1][0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][2];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][1][1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][0];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][1][2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][0][1];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][2][0]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][2];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][2][1]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][0];
                      grad_grads[i + (j + (2 * k + 9) * my_degree
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * (my_degree + 1)][1][2][2]
                        = p1_grad_grads[i + ((j + 2) * (my_degree + 2) + k)
                                        * (my_degree + 1)][1][1];
                    }
                }
        }

      break;
    }

    default:
      Assert (false, ExcNotImplemented ());
    }
}


template<int dim>
unsigned int
PolynomialsNedelec<dim>::compute_n_pols (unsigned int k)
{
  switch (dim)
    {
    case 1:
      return k + 1;

    case 2:
      return 2 * (k + 1) * (k + 2);

    case 3:
      return 3 * (k + 1) * (k + 2) * (k + 2);

    default:
    {
      Assert (false, ExcNotImplemented ());
      return 0;
    }
    }
}


template class PolynomialsNedelec<1>;
template class PolynomialsNedelec<2>;
template class PolynomialsNedelec<3>;


DEAL_II_NAMESPACE_CLOSE
