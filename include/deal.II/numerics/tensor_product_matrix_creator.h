// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_product_matrix_creator_h
#define dealii_tensor_product_matrix_creator_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
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
 * TensorProductMatrixSymmetricSumCollection.
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
   * Create 1d @ref GlossMassMatrix "mass matrix" and 1d derivative matrix for a scalar
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
    const FiniteElement<1>                             &fe,
    const Quadrature<1>                                &quadrature,
    const dealii::ndarray<LaplaceBoundaryType, dim, 2> &boundary_ids,
    const dealii::ndarray<double, dim, 3>              &cell_extent,
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
    const std::set<types::boundary_id>               &dirichlet_boundaries,
    const std::set<types::boundary_id>               &neumann_boundaries,
    const FiniteElement<1>                           &fe,
    const Quadrature<1>                              &quadrature,
    const dealii::ndarray<double, dim, 3>            &cell_extent,
    const unsigned int                                n_overlap = 1);



  /**
   * Compute a 1D cell mass matrix for a given finite element. This function is
   * intended to provide a mass matrix for tensor product operators.
   *
   * @param fe The finite element object defining the shape functions.
   *
   * @param h  The cell size used to scale the integration (typically the
   * element size).
   *
   * @param include_endpoints A pair of booleans indicating whether to include
   * the first and last dof in the matrix. By excluding those DoFs, one can set
   * a Dirichlet BC on either side. Setting either one to false only makes sense
   * if the indexing of DoFs in the matrix is lexicographic, i.e., for FE_Q
   * numbering is required, while FE_DGQ can be processed as it is.
   *
   * @param numbering A vector of unsigned integers representing the
   * numbering of the degrees of freedom. If empty, the
   * numbering of the finite element is used.
   *
   * @return A FullMatrix<double> representing the assembled mass matrix.
   *
   */
  template <typename Number = double>
  FullMatrix<Number>
  create_1d_cell_mass_matrix(
    const FiniteElement<1>     &fe,
    const Number               &h,
    const std::pair<bool, bool> include_endpoints = {true, true},
    std::vector<unsigned int>   numbering = std::vector<unsigned int>());



  /**
   * Compute a 1D cell laplace matrix for a given finite element.
   *
   * @param fe The finite element object defining the shape functions.
   *
   * @param h The cell size used to scale the integration (typically the
   * element size).
   *
   * @param include_endpoints A pair of booleans indicating whether to include
   * the first and last dof in the matrix. By excluding those DoFs, one can set
   * a Dirichlet BC on either side. Setting either one to false only makes sense
   * if the indexing of DoFs in the matrix is lexicographic, i.e., for FE_Q
   * numbering is required, while FE_DGQ can be processed as it is.
   *
   * @param numbering A vector of unsigned integers representing the
   * numbering of the degrees of freedom. If empty, the
   * numbering of the finite element is used.
   */
  template <typename Number = double>
  FullMatrix<Number>
  create_1d_cell_laplace_matrix(
    const FiniteElement<1>     &fe,
    const Number               &h,
    const std::pair<bool, bool> include_endpoints = {true, true},
    std::vector<unsigned int>   numbering = std::vector<unsigned int>());


  /**
   * This function assembles a global matrix from a given cell matrix, assuming
   * a 1D discretization with a specified number of cells and overlap between
   * them.
   *
   * @param cell_matrix The local cell matrix to be assembled into the
   * global matrix.  It is assumed that the cell matrix has a size corresponding
   * to the number of degrees of freedom per cell. The cell matrix should be
   * ordered lexicographically, meaning that FE_Q should be renumbered
   *
   * @param n_cells The number of cells in the 1D discretization.
   *
   * @param overlap The number of degrees of freedom that overlap between
   * adjacent cells. For discretization with FE_Q elements, the overlap should
   * be set to 1, while for FE_DGQ elements, it should be set to 0.
   *
   * @param include_endpoints A pair of booleans indicating whether to
   * include the left and right most dofs of the domain in the constructed
   * matrix. The default value is {true, true}, which includes both endpoints.
   *
   * @return The assembled global matrix. The size of the matrix is determined by
   * the number of cells, the overlap, and whether the endpoints are included.
   *
   * @note The size of the cell matrix must be consistent with the overlap
   * parameter. Specifically, if the cell matrix is $n \times n$, then the
   * overlap should be less than $n$.
   *
   * @warning This function returns a full matrix,
   * if the number of cells is above 10 an exception will be thrown.
   *
   */
  template <typename Number = double>
  FullMatrix<Number>
  create_1D_discretization_matrix(
    FullMatrix<Number>         &cell_matrix,
    const unsigned int         &n_cells,
    const unsigned int         &overlap,
    const std::pair<bool, bool> include_endpoints = {true, true});


  /**
   * @brief Create a 1D ghost penalty matrix for a given finite element.
   * Ghost penalty is used for stabilization in cutFEM, see step-85.
   * Implemented only for FE_Q.
   *
   * @param fe The finite element space used for discretization.
   *
   * @param h The mesh size, representing the length of the element.
   *
   * @param coefficients A vector of coefficients for the ghost penalty terms.
   *                     If the vector is empty, default coefficients are used.
   *                     The size of the coefficient vector should match the
   *                     polynomial degree of the finite element space.
   *
   * @return A full matrix representing the 1D ghost penalty matrix.
   */
  template <typename Number = double>
  FullMatrix<Number>
  create_1d_ghost_penalty_matrix(
    const FiniteElement<1> &fe,
    const Number            h,
    std::vector<Number>     coefficients = std::vector<Number>());



  /**
   * @brief Create a 1D ghost penalty matrix.
   *
   * This function creates a single component of the sum over all derivatives in
   * the 1D ghost penalty matrix. The matrix is created using the
   * derivatives of the polynomial basis functions.
   *
   * @param polynomial_basis_derivative A vector of polynomial basis functions.
   * These should be the derivatives of the basis functions used to represent
   * the solution.
   *
   * @param overlap The number of overlapping dof to consider when
   * computing the ghost penalty. For continuous element it is 1,
   * for DG pick 0.
   *
   * @return A full matrix representing the 1D ghost penalty.
   */
  template <typename Number = double>
  FullMatrix<Number>
  create_1d_ghost_penalty_matrix(
    const std::vector<Polynomials::Polynomial<double>>
                      &polynomial_basis_derivative,
    const unsigned int overlap = 1);

} // namespace TensorProductMatrixCreator



/*----------------------- Inline functions ----------------------------------*/


namespace TensorProductMatrixCreator
{
  namespace internal
  {
    template <typename Number>
    void
    clear_row_and_column(const unsigned int  n_dofs_1D_with_overlap,
                         const unsigned int  n,
                         FullMatrix<Number> &matrix)
    {
      for (unsigned int i = 0; i < n_dofs_1D_with_overlap; ++i)
        {
          matrix[i][n] = 0.0;
          matrix[n][i] = 0.0;
        }
    }

    template <typename Number>
    std::tuple<FullMatrix<Number>, FullMatrix<Number>, bool>
    create_reference_mass_and_stiffness_matrices(
      const FiniteElement<1> &fe,
      const Quadrature<1>    &quadrature)
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

      return std::tuple<FullMatrix<Number>, FullMatrix<Number>, bool>{
        mass_matrix_reference, derivative_matrix_reference, false};
    }
  } // namespace internal



  template <int dim, typename Number>
  std::pair<std::array<FullMatrix<Number>, dim>,
            std::array<FullMatrix<Number>, dim>>
  create_laplace_tensor_product_matrix(
    const FiniteElement<1>                             &fe,
    const Quadrature<1>                                &quadrature,
    const dealii::ndarray<LaplaceBoundaryType, dim, 2> &boundary_ids,
    const dealii::ndarray<double, dim, 3>              &cell_extent,
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

    // 2) loop over all dimensions and create 1d mass and stiffness
    // matrices so that boundary conditions and overlap are considered

    const unsigned int n_dofs_1D              = M_ref.n();
    const unsigned int n_dofs_1D_with_overlap = M_ref.n() - 2 + 2 * n_overlap;

    std::array<FullMatrix<Number>, dim> Ms;
    std::array<FullMatrix<Number>, dim> Ks;

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
                internal::clear_row_and_column(n_dofs_1D_with_overlap,
                                               i0,
                                               Ms[d]);
                internal::clear_row_and_column(n_dofs_1D_with_overlap,
                                               i0,
                                               Ks[d]);
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
                internal::clear_row_and_column(n_dofs_1D_with_overlap,
                                               i0,
                                               Ms[d]);
                internal::clear_row_and_column(n_dofs_1D_with_overlap,
                                               i0,
                                               Ks[d]);
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
    const std::set<types::boundary_id>               &dirichlet_boundaries,
    const std::set<types::boundary_id>               &neumann_boundaries,
    const FiniteElement<1>                           &fe,
    const Quadrature<1>                              &quadrature,
    const dealii::ndarray<double, dim, 3>            &cell_extent,
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

  template <typename Number>
  FullMatrix<Number>
  create_1d_cell_mass_matrix(const FiniteElement<1>     &fe,
                             const Number               &h,
                             const std::pair<bool, bool> include_endpoints,
                             std::vector<unsigned int>   numbering)
  {
    if (dynamic_cast<const FE_DGQ<1> *>(&fe) == nullptr &&
        numbering.size() == 0)
      {
        Assert(
          include_endpoints.first == true && include_endpoints.second == true,
          ExcMessage(
            "You tried to generate a 1D mass matrix with excluding boundary "
            "dofs for a non-DGQ element without providing a numbering."));
      }

    if (numbering.size() == 0)
      {
        numbering.resize(fe.dofs_per_cell);
        std::iota(numbering.begin(), numbering.end(), 0);
      }

    const unsigned int degree          = fe.degree;
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
    QGauss<1>          quadrature(degree + 1);

    FullMatrix<Number> cell_matrix(n_dofs_per_cell, n_dofs_per_cell);
    cell_matrix = 0;

    unsigned int start_dof = include_endpoints.first ? 0 : 1;
    unsigned int end_dof =
      include_endpoints.second ? n_dofs_per_cell : n_dofs_per_cell - 1;
    const unsigned int shift = include_endpoints.first ? 0 : 1;

    for (unsigned int i = start_dof; i < end_dof; ++i)
      for (unsigned int j = start_dof; j < end_dof; ++j)
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          cell_matrix(i - shift, j - shift) +=
            (fe.shape_value(numbering[i], quadrature.point(q)) *
             fe.shape_value(numbering[j], quadrature.point(q))) *
            (h * quadrature.weight(q));

    return cell_matrix;
  }

  template <typename Number>
  FullMatrix<Number>
  create_1d_cell_laplace_matrix(const FiniteElement<1>     &fe,
                                const Number               &h,
                                const std::pair<bool, bool> include_endpoints,
                                std::vector<unsigned int>   numbering)
  {
    if (dynamic_cast<const FE_DGQ<1> *>(&fe) == nullptr &&
        numbering.size() == 0)
      {
        Assert(
          include_endpoints.first == true && include_endpoints.second == true,
          ExcMessage(
            "You tried to generate a 1D derivative matrix with excluding boundary "
            "dofs for a non-DGQ element without providing a numbering."));
      }

    if (numbering.size() == 0)
      {
        numbering.resize(fe.dofs_per_cell);
        std::iota(numbering.begin(), numbering.end(), 0);
      }

    const unsigned int degree          = fe.degree;
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
    const Number      &JxW             = h;
    QGauss<1>          quadrature(degree + 1);

    FullMatrix<Number> cell_matrix(n_dofs_per_cell, n_dofs_per_cell);
    cell_matrix = 0;

    unsigned int start_dof = include_endpoints.first ? 0 : 1;
    unsigned int end_dof =
      include_endpoints.second ? n_dofs_per_cell : n_dofs_per_cell - 1;
    const unsigned int shift = include_endpoints.first ? 0 : 1;

    for (unsigned int i = start_dof; i < end_dof; ++i)
      for (unsigned int j = start_dof; j < end_dof; ++j)
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          cell_matrix(i - shift, j - shift) +=
            (fe.shape_grad(numbering[i], quadrature.point(q)) / h *
             fe.shape_grad(numbering[j], quadrature.point(q))) /
            h * (h * quadrature.weight(q));

    return cell_matrix;
  }

  template <typename Number>
  FullMatrix<Number>
  create_1D_discretization_matrix(FullMatrix<Number>         &cell_matrix,
                                  const unsigned int         &n_cells,
                                  const unsigned int         &overlap,
                                  const std::pair<bool, bool> include_endpoints)
  {
    const unsigned int n_dofs_per_cell = cell_matrix.n();

    Assert(cell_matrix.m() == n_dofs_per_cell,
           ExcMessage(
             "The provided cell mass matrix must be a square matrix."));
    AssertThrow(
      n_cells <= 10,
      ExcMessage(
        "create_1D_discretization_matrix() returns a full matrix and is not meant to be used with a larger number of cells. "));
    Assert(n_cells > 0,
           ExcMessage("You are trying to get a mass matrix of zero cells."));
    Assert(overlap < n_dofs_per_cell,
           ExcMessage("The overlap must be smaller than the number of dofs."));

    unsigned int n_total_dofs =
      n_cells * n_dofs_per_cell - overlap * (n_cells - 1);

    if (!include_endpoints.first)
      n_total_dofs -= 1;
    if (!include_endpoints.second)
      n_total_dofs -= 1;

    FullMatrix<Number> result_matrix(n_total_dofs, n_total_dofs);
    result_matrix = 0;

    const unsigned int left_shift = include_endpoints.first ? 0 : 1;

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        const unsigned int dof_shift = cell * overlap + left_shift;

        const unsigned int start_dof =
          (cell == 0 && !include_endpoints.first) ? 1 : 0;

        const unsigned int end_dof =
          (cell == n_cells - 1 && !include_endpoints.second) ?
            n_dofs_per_cell - 1 :
            n_dofs_per_cell;
        for (unsigned int i = start_dof; i < end_dof; ++i)
          for (unsigned int j = start_dof; j < end_dof; ++j)
            {
              result_matrix(i + cell * n_dofs_per_cell - dof_shift,
                            j + cell * n_dofs_per_cell - dof_shift) +=
                cell_matrix(i, j);
            }
      }
    return result_matrix;
  }



  template <typename Number>
  FullMatrix<Number>
  create_1d_ghost_penalty_matrix(const FiniteElement<1> &fe,
                                 const Number            h,
                                 std::vector<Number>     coefficients)
  {
    Assert(dynamic_cast<const FE_Q<1> *>(&fe) != nullptr, ExcNotImplemented());
    Assert(h > 0, ExcMessage("Provided element size h is negative"));

    const unsigned int degree = fe.degree;
    Assert(degree > 0,
           ExcMessage("Provided element degree has to greater than 0"));


    Assert(coefficients.size() == 0 || coefficients.size() == degree,
           ExcMessage(
             "Provided coefficients vector has to be empty or the same size "
             "as the number of dofs"));

    if (coefficients.size() == 0)
      {
        coefficients.resize(degree);

        double inverse_factorial_square = 1.;
        coefficients[0]                 = 1.;
        for (unsigned int k = 2; k <= degree; ++k)
          {
            inverse_factorial_square /= (k * k);
            coefficients[k - 1] = inverse_factorial_square;
          }
      }

    std::vector<std::vector<Polynomials::Polynomial<double>>> polynomial_basis;

    polynomial_basis.resize(degree + 1);

    auto support_points = fe.get_unit_support_points();
    std::sort(support_points.begin(),
              support_points.end(),
              [](const Point<1> &p, const Point<1> &q) -> bool {
                return p(0) < q(0);
              });

    polynomial_basis[0] =
      Polynomials::generate_complete_Lagrange_basis(support_points);

    for (unsigned int k = 1; k < degree + 1; ++k)
      {
        polynomial_basis[k].reserve(degree + 1);
        for (unsigned int i = 0; i < degree + 1; ++i)
          polynomial_basis[k].push_back(
            polynomial_basis[k - 1][i].derivative());
      }


    FullMatrix<Number> penalty_matrix =
      create_1d_ghost_penalty_matrix(polynomial_basis[1]);
    penalty_matrix *= coefficients[0];

    for (unsigned int k = 2; k < degree + 1; ++k)
      {
        FullMatrix<Number> kth_matrix =
          create_1d_ghost_penalty_matrix(polynomial_basis[k]);
        penalty_matrix.add(coefficients[k - 1], kth_matrix);
      }

    penalty_matrix *= (1 / h);
    return penalty_matrix;
  }


  template <typename Number>
  FullMatrix<Number>
  create_1d_ghost_penalty_matrix(
    const std::vector<Polynomials::Polynomial<double>>
                      &polynomial_basis_derivative,
    const unsigned int overlap)
  {
    const unsigned int n_dofs_per_cell = polynomial_basis_derivative.size();
    const unsigned int n_total_dofs    = 2 * n_dofs_per_cell - overlap;
    const unsigned int shift           = n_dofs_per_cell - overlap;

    FullMatrix<Number> penalty_matrix(n_total_dofs, n_total_dofs);

    std::vector<double> values_left(n_dofs_per_cell);
    std::vector<double> values_right(n_dofs_per_cell);

    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      {
        values_left[i]  = polynomial_basis_derivative[i].value(0);
        values_right[i] = polynomial_basis_derivative[i].value(1);
      }

    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
        penalty_matrix(i, j) += values_right[i] * values_right[j];


    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
        penalty_matrix(i + shift, j) -= values_left[i] * values_right[j];

    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
        penalty_matrix(i, j + shift) -= values_right[i] * values_left[j];

    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
        penalty_matrix(i + shift, j + shift) += values_left[i] * values_left[j];

    return penalty_matrix;
  }
} // namespace TensorProductMatrixCreator



DEAL_II_NAMESPACE_CLOSE

#endif
