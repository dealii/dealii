// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_tensor_product_matrix_creator_h
#define dealii_tensor_product_matrix_creator_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <set>

DEAL_II_NAMESPACE_OPEN


/**
 * A namespace with functions that create input for certain standard matrices
 * for the classes TensorProductMatrixSymmetricSum and
 * TensorProductMatrixSymmetricSumCache.
 */
namespace TensorProductMatrixCreator
{
  /**
   * Boundary type that can be used in create_laplace_tensor_product_matrix();
   */
  enum LaplaceBoundaryType
  {
    dirichlet,
    neumann,
    internal_boundary,
  };

  /**
   * Create 1D mass matrix and 1D derivative matrix for a scalar
   * constant-coefficient
   * Laplacian for a @p dim dimensional Cartesian cell. Its boundary types
   * can be specified with @p boundary_ids. The cell extent (including the cell extent
   * of each neighbor) can be specified via @p cell_extent. With @p n_overlap, an
   * overlap with neighboring cells can be specified. The default value is one,
   * which correspond to all matrix entries restricted to the cell-local DoFs.
   */
  template <int dim, typename Number>
  std::pair<std::array<FullMatrix<Number>, dim>,
            std::array<FullMatrix<Number>, dim>>
  create_laplace_tensor_product_matrix(
    const FiniteElement<1> &                            fe,
    const Quadrature<1> &                               quadrature,
    const dealii::ndarray<LaplaceBoundaryType, dim, 2> &boundary_ids,
    const dealii::ndarray<double, dim, 3> &             cell_extent,
    const unsigned int                                  n_overlap = 1);

  /**
   * Same as above but the boundary IDs are extracted from the given @p cell
   * and are mapped to the boundary type via the sets @p dirichlet_boundaries and @p neumann_boundaries.
   */
  template <int dim, typename Number>
  std::pair<std::array<FullMatrix<Number>, dim>,
            std::array<FullMatrix<Number>, dim>>
  create_laplace_tensor_product_matrix(
    const typename Triangulation<dim>::cell_iterator &cell,
    const std::set<types::boundary_id> &              dirichlet_boundaries,
    const std::set<types::boundary_id> &              neumann_boundaries,
    const FiniteElement<1> &                          fe,
    const Quadrature<1> &                             quadrature,
    const dealii::ndarray<double, dim, 3> &           cell_extent,
    const unsigned int                                n_overlap = 1);

} // namespace TensorProductMatrixCreator



/*----------------------- Inline functions ----------------------------------*/


namespace TensorProductMatrixCreator
{
  namespace internal
  {
    template <typename Number>
    std::tuple<FullMatrix<Number>, FullMatrix<Number>, bool>
    create_reference_mass_and_stiffness_matrices(
      const FiniteElement<1> &fe,
      const Quadrature<1> &   quadrature)
    {
      Triangulation<1> tria;
      GridGenerator::hyper_cube(tria);

      DoFHandler<1> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      MappingQ1<1> mapping;

      const unsigned int n_dofs_1D = fe.n_dofs_per_cell();

      FullMatrix<Number> mass_matrix_reference(n_dofs_1D, n_dofs_1D);
      FullMatrix<Number> derivative_matrix_reference(n_dofs_1D, n_dofs_1D);

      FEValues<1> fe_values(mapping,
                            fe,
                            quadrature,
                            update_values | update_gradients |
                              update_JxW_values);

      fe_values.reinit(tria.begin());

      const auto lexicographic_to_hierarchic_numbering =
        Utilities::invert_permutation(
          FETools::hierarchic_to_lexicographic_numbering<1>(
            fe.tensor_degree()));

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          for (const unsigned int j : fe_values.dof_indices())
            {
              mass_matrix_reference(i, j) +=
                (fe_values.shape_value(lexicographic_to_hierarchic_numbering[i],
                                       q_index) *
                 fe_values.shape_value(lexicographic_to_hierarchic_numbering[j],
                                       q_index) *
                 fe_values.JxW(q_index));

              derivative_matrix_reference(i, j) +=
                (fe_values.shape_grad(lexicographic_to_hierarchic_numbering[i],
                                      q_index) *
                 fe_values.shape_grad(lexicographic_to_hierarchic_numbering[j],
                                      q_index) *
                 fe_values.JxW(q_index));
            }

      return {mass_matrix_reference, derivative_matrix_reference, false};
    }
  } // namespace internal



  template <int dim, typename Number>
  std::pair<std::array<FullMatrix<Number>, dim>,
            std::array<FullMatrix<Number>, dim>>
  create_laplace_tensor_product_matrix(
    const FiniteElement<1> &                            fe,
    const Quadrature<1> &                               quadrature,
    const dealii::ndarray<LaplaceBoundaryType, dim, 2> &boundary_ids,
    const dealii::ndarray<double, dim, 3> &             cell_extent,
    const unsigned int                                  n_overlap)
  {
    // 1) create element mass and siffness matrix (without overlap)
    const auto create_reference_mass_and_stiffness_matrices =
      internal::create_reference_mass_and_stiffness_matrices<Number>(
        fe, quadrature);

    const auto &M_ref =
      std::get<0>(create_reference_mass_and_stiffness_matrices);
    const auto &K_ref =
      std::get<1>(create_reference_mass_and_stiffness_matrices);
    const auto &is_dg =
      std::get<2>(create_reference_mass_and_stiffness_matrices);

    AssertIndexRange(n_overlap, M_ref.n());
    AssertIndexRange(0, n_overlap);
    AssertThrow(is_dg == false, ExcNotImplemented());

    // 2) loop over all dimensions and create 1D mass and stiffness
    // matrices so that boundary conditions and overlap are considered

    const unsigned int n_dofs_1D              = M_ref.n();
    const unsigned int n_dofs_1D_with_overlap = M_ref.n() - 2 + 2 * n_overlap;

    std::array<FullMatrix<Number>, dim> Ms;
    std::array<FullMatrix<Number>, dim> Ks;

    const auto clear_row_and_column = [&](const unsigned int n, auto &matrix) {
      for (unsigned int i = 0; i < n_dofs_1D_with_overlap; ++i)
        {
          matrix[i][n] = 0.0;
          matrix[n][i] = 0.0;
        }
    };

    for (unsigned int d = 0; d < dim; ++d)
      {
        Ms[d].reinit(n_dofs_1D_with_overlap, n_dofs_1D_with_overlap);
        Ks[d].reinit(n_dofs_1D_with_overlap, n_dofs_1D_with_overlap);

        // inner cell
        for (unsigned int i = 0; i < n_dofs_1D; ++i)
          for (unsigned int j = 0; j < n_dofs_1D; ++j)
            {
              const unsigned int i0 = i + n_overlap - 1;
              const unsigned int j0 = j + n_overlap - 1;
              Ms[d][i0][j0]         = M_ref[i][j] * cell_extent[d][1];
              Ks[d][i0][j0]         = K_ref[i][j] / cell_extent[d][1];
            }

        // left neighbor or left boundary
        if (boundary_ids[d][0] == LaplaceBoundaryType::internal_boundary)
          {
            // left neighbor
            Assert(cell_extent[d][0] > 0.0, ExcInternalError());

            for (unsigned int i = 0; i < n_overlap; ++i)
              for (unsigned int j = 0; j < n_overlap; ++j)
                {
                  const unsigned int i0 = n_dofs_1D - n_overlap + i;
                  const unsigned int j0 = n_dofs_1D - n_overlap + j;
                  Ms[d][i][j] += M_ref[i0][j0] * cell_extent[d][0];
                  Ks[d][i][j] += K_ref[i0][j0] / cell_extent[d][0];
                }
          }
        else
          {
            if (boundary_ids[d][0] == LaplaceBoundaryType::dirichlet)
              {
                // left DBC
                const unsigned i0 = n_overlap - 1;
                clear_row_and_column(i0, Ms[d]);
                clear_row_and_column(i0, Ks[d]);
              }
            else if (boundary_ids[d][0] == LaplaceBoundaryType::neumann)
              {
                // left NBC -> nothing to do
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
          }

        // right neighbor or right boundary
        if (boundary_ids[d][1] == LaplaceBoundaryType::internal_boundary)
          {
            Assert(cell_extent[d][2] > 0.0, ExcInternalError());

            for (unsigned int i = 0; i < n_overlap; ++i)
              for (unsigned int j = 0; j < n_overlap; ++j)
                {
                  const unsigned int i0 = n_overlap + n_dofs_1D + i - 2;
                  const unsigned int j0 = n_overlap + n_dofs_1D + j - 2;
                  Ms[d][i0][j0] += M_ref[i][j] * cell_extent[d][2];
                  Ks[d][i0][j0] += K_ref[i][j] / cell_extent[d][2];
                }
          }
        else
          {
            if (boundary_ids[d][1] == LaplaceBoundaryType::dirichlet)
              {
                // right DBC
                const unsigned i0 = n_overlap + n_dofs_1D - 2;
                clear_row_and_column(i0, Ms[d]);
                clear_row_and_column(i0, Ks[d]);
              }
            else if (boundary_ids[d][1] == LaplaceBoundaryType::neumann)
              {
                // right NBC -> nothing to do
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
          }
      }

    return {Ms, Ks};
  }


  template <int dim, typename Number>
  std::pair<std::array<FullMatrix<Number>, dim>,
            std::array<FullMatrix<Number>, dim>>
  create_laplace_tensor_product_matrix(
    const typename Triangulation<dim>::cell_iterator &cell,
    const std::set<types::boundary_id> &              dirichlet_boundaries,
    const std::set<types::boundary_id> &              neumann_boundaries,
    const FiniteElement<1> &                          fe,
    const Quadrature<1> &                             quadrature,
    const dealii::ndarray<double, dim, 3> &           cell_extent,
    const unsigned int                                n_overlap)
  {
    dealii::ndarray<LaplaceBoundaryType, dim, 2> boundary_ids;

    for (unsigned int d = 0; d < dim; ++d)
      {
        // left neighbor or left boundary
        if ((cell->at_boundary(2 * d) == false) ||
            cell->has_periodic_neighbor(2 * d))
          {
            // left neighbor
            Assert(cell_extent[d][0] > 0.0, ExcInternalError());

            boundary_ids[d][0] = LaplaceBoundaryType::internal_boundary;
          }
        else
          {
            const auto bid = cell->face(2 * d)->boundary_id();
            if (dirichlet_boundaries.find(bid) !=
                dirichlet_boundaries.end() /*DBC*/)
              {
                // left DBC
                boundary_ids[d][0] = LaplaceBoundaryType::dirichlet;
              }
            else if (neumann_boundaries.find(bid) !=
                     neumann_boundaries.end() /*NBC*/)
              {
                // left NBC
                boundary_ids[d][0] = LaplaceBoundaryType::neumann;
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
          }

        // right neighbor or right boundary
        if ((cell->at_boundary(2 * d + 1) == false) ||
            cell->has_periodic_neighbor(2 * d + 1))
          {
            Assert(cell_extent[d][2] > 0.0, ExcInternalError());

            boundary_ids[d][1] = LaplaceBoundaryType::internal_boundary;
          }
        else
          {
            const auto bid = cell->face(2 * d + 1)->boundary_id();
            if (dirichlet_boundaries.find(bid) !=
                dirichlet_boundaries.end() /*DBC*/)
              {
                // right DBC
                boundary_ids[d][1] = LaplaceBoundaryType::dirichlet;
              }
            else if (neumann_boundaries.find(bid) !=
                     neumann_boundaries.end() /*NBC*/)
              {
                // right NBC
                boundary_ids[d][1] = LaplaceBoundaryType::neumann;
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
          }
      }

    return create_laplace_tensor_product_matrix<dim, Number>(
      fe, quadrature, boundary_ids, cell_extent, n_overlap);
  }

} // namespace TensorProductMatrixCreator



DEAL_II_NAMESPACE_CLOSE

#endif
