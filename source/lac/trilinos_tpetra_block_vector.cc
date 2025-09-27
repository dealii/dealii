// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_tpetra_block_vector.templates.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    // Instantiate these vectors types for specific scalar types.
    //
    // While there:
    // Check that the class we declare here satisfies the
    // vector-space-vector concept. If we catch it here,
    // any mistake in the vector class declaration would
    // show up in uses of this class later on as well.

#  ifdef HAVE_TPETRA_INST_FLOAT
#    ifdef DEAL_II_HAVE_CXX20
    static_assert(
      concepts::is_vector_space_vector<BlockVector<float, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<float, MemorySpace::Default>>);
#    endif
    template class BlockVector<float, MemorySpace::Host>;
    template class BlockVector<float, MemorySpace::Default>;
#  endif
#  ifdef HAVE_TPETRA_INST_DOUBLE
#    ifdef DEAL_II_HAVE_CXX20
    static_assert(
      concepts::is_vector_space_vector<BlockVector<double, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<double, MemorySpace::Default>>);
#    endif
    template class BlockVector<double, MemorySpace::Host>;
    template class BlockVector<double, MemorySpace::Default>;
#  endif
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<std::complex<float>, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<std::complex<float>, MemorySpace::Default>>);
#      endif
    template class BlockVector<std::complex<float>, MemorySpace::Host>;
    template class BlockVector<std::complex<float>, MemorySpace::Default>;
#    endif
#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<std::complex<double>, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  BlockVector<std::complex<double>, MemorySpace::Default>>);
#      endif
    template class BlockVector<std::complex<double>, MemorySpace::Host>;
    template class BlockVector<std::complex<double>, MemorySpace::Default>;
#    endif
#  endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
