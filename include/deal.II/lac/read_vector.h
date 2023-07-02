// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef dealii_lac_read_vector_h
#define dealii_lac_read_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/types.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Vectors
 * @{
 */

/**
 * Base class for providing read-only access to vector elements.
 *
 * deal.II supports a large number of vector classes, including both its own
 * serial and parallel vectors as well as vector classes from external
 * libraries like PETSc and Trilinos. ReadVector is a common base class for
 * all vector classes and defines a minimal interface for efficiently
 * accessing vector elements.
 */
template <typename Number>
class ReadVector
{
public:
  using size_type = types::global_dof_index;

  /**
   * Return the size of the vector.
   */
  virtual size_type
  size() const = 0;

  /**
   * Extract a subset of the vector specified by @p indices into the output
   * array @p elements.
   */
  virtual void
  extract_subvector_to(const ArrayView<const types::global_dof_index> &indices,
                       ArrayView<Number> &elements) const = 0;
};

/** @} */

DEAL_II_NAMESPACE_CLOSE
#endif
