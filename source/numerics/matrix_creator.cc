// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/numerics/matrix_creator.templates.h>


DEAL_II_NAMESPACE_OPEN

// explicit instantiations
#define SPLIT_INSTANTIATIONS_COUNT 3
#ifndef SPLIT_INSTANTIATIONS_INDEX
#  define SPLIT_INSTANTIATIONS_INDEX 0
#endif
#include "numerics/matrix_creator.inst"

#if defined(SPLIT_INSTANTIATIONS_INDEX) && SPLIT_INSTANTIATIONS_INDEX == 0
namespace MatrixCreator
{

  FullMatrix<double>
  create_1d_cell_mass_matrix(const FiniteElement<1>     &fe,
                             const double               &h,
                             const std::pair<bool, bool> include_endpoints,
                             std::vector<unsigned int>   numbering)
  {
    if (dynamic_cast<const FE_DGQ<1> *>(&fe) == nullptr &&
        numbering.size() == 0)
      {
        Assert(
          include_endpoints.first == true && include_endpoints.second == true,
          ExcMessage(
            "You tried to genereate a 1D mass matrix with excluding boundary "
            "dofs for a non-DGQ element without providing a numbering."));
      }


    if (numbering.size() == 0)
      {
        numbering.resize(fe.dofs_per_cell);
        std::iota(numbering.begin(), numbering.end(), 0);
      }
    const unsigned int degree          = fe.degree;
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
    const double      &JxW             = h;
    QGauss<1>          quadrature(degree + 1);

    FullMatrix<double> cell_mass_matrix(n_dofs_per_cell, n_dofs_per_cell);
    cell_mass_matrix = 0;

    unsigned int start_dof = include_endpoints.first ? 0 : 1;
    unsigned int end_dof =
      include_endpoints.second ? n_dofs_per_cell : n_dofs_per_cell - 1;
    const unsigned int shift = include_endpoints.first ? 0 : 1;

    for (unsigned int i = start_dof; i < end_dof; ++i)
      for (unsigned int j = start_dof; j < end_dof; ++j)
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          cell_mass_matrix(i - shift, j - shift) +=
            (fe.shape_value(numbering[i], quadrature.point(q)) *
             fe.shape_value(numbering[j], quadrature.point(q))) *
            JxW * quadrature.weight(q);

    return cell_mass_matrix;
  }


  FullMatrix<double>
  create_1d_cell_derivative_matrix(
    const FiniteElement<1>     &fe,
    const double               &h,
    const std::pair<bool, bool> include_endpoints,
    std::vector<unsigned int>   numbering)
  {
    if (numbering.size() == 0)
      {
        numbering.resize(fe.dofs_per_cell);
        std::iota(numbering.begin(), numbering.end(), 0);
      }

    const unsigned int degree          = fe.degree;
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
    const double      &JxW             = h;
    QGauss<1>          quadrature(degree + 1);

    FullMatrix<double> cell_matrix(n_dofs_per_cell, n_dofs_per_cell);
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
            h * JxW * quadrature.weight(q);

    return cell_matrix;
  }

  FullMatrix<double>
  create_1D_discretization_matrix(FullMatrix<double>         &cell_matrix,
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

    FullMatrix<double> result_matrix(n_total_dofs, n_total_dofs);
    result_matrix = 0;

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        const unsigned int dof_shift =
          cell * overlap + !include_endpoints.first;

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
} // namespace MatrixCreator

#endif // SPLIT_INSTANTIATIONS_INDEX == 0

DEAL_II_NAMESPACE_CLOSE
