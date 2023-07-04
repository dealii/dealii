// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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

#include <deal.II/lac/trilinos_tpetra_vector.templates.h>

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
    static_assert(concepts::is_vector_space_vector<Vector<float>>);
#    endif
    template class Vector<float>;
#  endif
#  ifdef HAVE_TPETRA_INST_DOUBLE
#    ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<Vector<double>>);
#    endif
    template class Vector<double>;
#  endif
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector <
                  Vector<std::complex<float>>);
#      endif
    template class Vector<std::complex<float>>;
#    endif
#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector <
                  Vector<std::complex<double>>);
#      endif
    template class Vector<std::complex<double>>;
#    endif
#  endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
