// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_lac_vector_operation_h
#define dealii_lac_vector_operation_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Vectors
 * @{
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

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
