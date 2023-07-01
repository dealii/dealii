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

#ifndef dealii_lac_orthogonalization_h
#define dealii_lac_orthogonalization_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  /**
   * Supported orthogonalization strategies within SolverGMRES and
   * SolverFGMRES.
   */
  enum class OrthogonalizationStrategy
  {
    /**
     * Use modified Gram-Schmidt algorithm.
     */
    modified_gram_schmidt,
    /**
     * Use classical Gram-Schmidt algorithm. Since this approach works on
     * multi-vectors with a single global reduction (of multiple elements), it
     * is more efficient than the modified Gram-Schmidt algorithm. However, it
     * is less stable in terms of roundoff error propagation, requiring
     * additional re-orthogonalization steps more frequently.
     */
    classical_gram_schmidt
  };
} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
