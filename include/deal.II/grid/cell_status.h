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
 * Used to inform functions in derived classes how the cell with the given
 * cell_iterator is going to change. Note that this may me different than
 * the refine_flag() and coarsen_flag() in the cell_iterator in parallel
 * calculations because of refinement constraints that this machine does not
 * see.
 */
enum CellStatus
{
  /**
   * The cell will not be refined or coarsened and might or might not move
   * to a different processor.
   */
  CELL_PERSIST,
  /**
   * The cell will be or was refined.
   */
  CELL_REFINE,
  /**
   * The children of this cell will be or were coarsened into this cell.
   */
  CELL_COARSEN,
  /**
   * Invalid status. Will not occur for the user.
   */
  CELL_INVALID
};

DEAL_II_NAMESPACE_CLOSE

#endif
