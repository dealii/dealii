// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_trilinos_index_access_h
#define dealii_trilinos_index_access_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/types.h>

#  include <Epetra_BlockMap.h>
#  include <Epetra_CrsGraph.h>
#  include <Epetra_CrsMatrix.h>
#  include <Epetra_MultiVector.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  /**
   * A helper function that queries the size of an Epetra_BlockMap object
   * and calls either the 32 or 64 bit function to get the number of global
   * elements in the map.
   */
  inline TrilinosWrappers::types::int_type
  n_global_elements(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.NumGlobalElements64();
#  else
    return map.NumGlobalElements();
#  endif
  }

  /**
   * A helper function that finds the minimum global index value on the
   * calling processor by calling either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  min_my_gid(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MinMyGID64();
#  else
    return map.MinMyGID();
#  endif
  }

  /**
   * A helper function that finds the maximum global index value on the
   * calling processor by calling either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  max_my_gid(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MaxMyGID64();
#  else
    return map.MaxMyGID();
#  endif
  }

  /**
   * A helper function that converts a local index to a global one calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  global_index(const Epetra_BlockMap &               map,
               const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.GID64(i);
#  else
    return map.GID(i);
#  endif
  }

  /**
   * A helper function that returns a pointer to the array containing the
   * global indices assigned to the current process by calling either the 32
   * or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type *
  my_global_elements(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MyGlobalElements64();
#  else
    return map.MyGlobalElements();
#  endif
  }

  /**
   * A helper function that finds the global number of rows by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  n_global_rows(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalRows64();
#  else
    return graph.NumGlobalRows();
#  endif
  }

  /**
   * A helper function that finds the global number of columns by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  n_global_cols(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalCols64();
#  else
    return graph.NumGlobalCols();
#  endif
  }

  /**
   * A helper function that finds the number of global entries by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  n_global_entries(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalEntries64();
#  else
    return graph.NumGlobalEntries();
#  endif
  }

  /**
   * A helper function that finds the global row index by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  global_row_index(const Epetra_CrsMatrix &              matrix,
                   const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.GRID64(i);
#  else
    return matrix.GRID(i);
#  endif
  }

  /**
   * A helper function that finds the global column index by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  global_column_index(const Epetra_CrsMatrix &              matrix,
                      const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.GCID64(i);
#  else
    return matrix.GCID(i);
#  endif
  }

  /**
   * A helper function that finds the global length of a vector by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  global_length(const Epetra_MultiVector &vector)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return vector.GlobalLength64();
#  else
    return vector.GlobalLength();
#  endif
  }

  /**
   * A helper function that finds the global number of rows by calling
   * either the 32 or 64 bit function.
   */
  inline TrilinosWrappers::types::int_type
  n_global_rows(const Epetra_RowMatrix &matrix)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.NumGlobalRows64();
#  else
    return matrix.NumGlobalRows();
#  endif
  }
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE
#endif // DEAL_II_WITH_TRILINOS
#endif // dealii_trilinos_index_access_h
