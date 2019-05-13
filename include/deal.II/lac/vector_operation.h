// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#ifndef dealii_lac_vector_operation_h
#define dealii_lac_vector_operation_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Vectors
 *@{
 */

/**
 * This enum keeps track of the current operation in parallel linear algebra
 * objects like Vectors and Matrices.
 *
 * It is used in the various compress() functions. They also exist in serial
 * codes for compatibility and are empty there.
 *
 * See
 * @ref GlossCompress "Compressing distributed objects"
 * for more information.
 */
struct VectorOperation
{
  enum values
  {
    /**
     * The current operation is unknown.
     */
    unknown,
    /**
     * The current operation is an insertion.
     */
    insert,
    /**
     * The current operation is an addition.
     */
    add,
    /**
     * The current operation is a minimization.
     */
    min,
    /**
     * The current operation is a maximization.
     */
    max
  };
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
