// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
    template Vector<float> &
    Vector<float>::operator=<float>(const dealii::Vector<float> &);
    namespace internal
    {
      template class VectorReference<float>;
    }
#  endif

#  ifdef HAVE_TPETRA_INST_DOUBLE
#    ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<Vector<double>>);
#    endif
    template class Vector<double>;
    template Vector<double> &
    Vector<double>::operator=<double>(const dealii::Vector<double> &);
    namespace internal
    {
      template class VectorReference<double>;
    }
#  endif

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(
      concepts::is_vector_space_vector<Vector<std::complex<float>>>);
#      endif
    template class Vector<std::complex<float>>;
    template Vector<std::complex<float>> &
    Vector<std::complex<float>>::operator=
      <std::complex<float>>(const dealii::Vector<std::complex<float>> &);
    namespace internal
    {
      template class VectorReference<std::complex<float>>;
    }
#    endif

#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(
      concepts::is_vector_space_vector<Vector<std::complex<double>>>);
#      endif
    template class Vector<std::complex<double>>;
    template Vector<std::complex<double>> &
    Vector<std::complex<double>>::operator=
      <std::complex<double>>(const dealii::Vector<std::complex<double>> &);
    namespace internal
    {
      template class VectorReference<std::complex<double>>;
    }
#    endif
#  endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
