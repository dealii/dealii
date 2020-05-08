// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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

#ifndef dealii_parallel_block_vector_h
#define dealii_parallel_block_vector_h

#include <deal.II/base/config.h>

#include <deal.II/lac/la_parallel_block_vector.h>

DEAL_II_WARNING(
  "This file is deprecated. Use <deal.II/lac/la_parallel_block_vector.h> and LinearAlgebra::distributed::BlockVector instead.")

#include <cstring>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace distributed
  {
    /*! @addtogroup Vectors
     *@{
     */

    /**
     * An implementation of block vectors based on distributed deal.II
     * vectors. While the base class provides for most of the interface, this
     * class handles the actual allocation of vectors and provides functions
     * that are specific to the underlying vector type.
     *
     * @note Instantiations for this template are provided for <tt>@<float@>
     * and @<double@></tt>; others can be generated in application programs
     * (see the section on
     * @ref Instantiations
     * in the manual).
     *
     * @see
     * @ref GlossBlockLA "Block (linear algebra)"
     * @author Katharina Kormann, Martin Kronbichler, 2011
     *
     * @deprecated Use LinearAlgebra::distributed::BlockVector instead.
     */
    template <typename Number>
    using BlockVector DEAL_II_DEPRECATED =
      LinearAlgebra::distributed::BlockVector<Number>;

    /*@}*/
  } // namespace distributed
} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
