// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__parallel_block_vector_h
#define dealii__parallel_block_vector_h

#include <deal.II/base/config.h>
#include <deal.II/lac/la_parallel_block_vector.h>

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
     */
    using LinearAlgebra::distributed::BlockVector;

    /*@}*/
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
