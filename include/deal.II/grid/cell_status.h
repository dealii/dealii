// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cell_status_h
#define dealii_cell_status_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * The elements of this `enum` are used to inform functions how a specific cell
 * is going to change. This is used in the course of transferring data from one
 * mesh to a refined or coarsened version of the mesh, for example. Note that
 * this may me different than the refine_flag() and coarsen_flag() set on a
 * cell, for example in parallel calculations, because of refinement constraints
 * that an individual machine does not see.
 */
enum class CellStatus : unsigned int
{
  /**
   * The cell will not be refined or coarsened and might or might not move
   * to a different processor.
   */
  cell_will_persist,
  /**
   * The cell will be or was refined.
   */
  cell_will_be_refined,
  /**
   * The children of this cell will be or were coarsened into this cell.
   */
  children_will_be_coarsened,
  /**
   * Invalid status. Will not occur for the user.
   */
  cell_invalid,

  /**
   * @deprecated Use CellStatus::cell_will_persist instead
   */
  CELL_PERSIST DEAL_II_DEPRECATED = cell_will_persist,

  /**
   * @deprecated Use CellStatus::cell_will_be_refined instead
   */
  CELL_REFINE DEAL_II_DEPRECATED = cell_will_be_refined,

  /**
   * @deprecated Use CellStatus::cell_will_be_coarsened instead
   */
  CELL_COARSEN DEAL_II_DEPRECATED = children_will_be_coarsened,

  /**
   * @deprecated Use CellStatus::cell_invalid instead
   */
  CELL_INVALID DEAL_II_DEPRECATED = cell_invalid,
};

DEAL_II_NAMESPACE_CLOSE

#endif
