// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2023 by the deal.II authors
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
  cell_invalid
};

DEAL_II_NAMESPACE_CLOSE

#endif
