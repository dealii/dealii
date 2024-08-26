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
    static_assert(
      concepts::is_vector_space_vector<Vector<float, MemorySpace::Host>>);
    static_assert(
      concepts::is_vector_space_vector<Vector<float, MemorySpace::Default>>);
#    endif
    template class Vector<float, MemorySpace::Host>;
    template class Vector<float, MemorySpace::Default>;
    template Vector<float, MemorySpace::Host> &
    Vector<float, MemorySpace::Host>::operator=
      <float>(const dealii::Vector<float> &);
    template Vector<float, MemorySpace::Default> &
    Vector<float, MemorySpace::Default>::operator=
      <float>(const dealii::Vector<float> &);
    namespace internal
    {
      template class VectorReference<float, MemorySpace::Host>;
      template class VectorReference<float, MemorySpace::Default>;
    } // namespace internal
#  endif

#  ifdef HAVE_TPETRA_INST_DOUBLE
#    ifdef DEAL_II_HAVE_CXX20
    static_assert(
      concepts::is_vector_space_vector<Vector<double, MemorySpace::Host>>);
    static_assert(
      concepts::is_vector_space_vector<Vector<double, MemorySpace::Default>>);
#    endif
    template class Vector<double, MemorySpace::Host>;
    template class Vector<double, MemorySpace::Default>;
    template Vector<double, MemorySpace::Host> &
    Vector<double, MemorySpace::Host>::operator=
      <double>(const dealii::Vector<double> &);
    template Vector<double, MemorySpace::Default> &
    Vector<double, MemorySpace::Default>::operator=
      <double>(const dealii::Vector<double> &);
    namespace internal
    {
      template class VectorReference<double, MemorySpace::Host>;
      template class VectorReference<double, MemorySpace::Default>;
    } // namespace internal
#  endif

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<
                  Vector<std::complex<float>, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  Vector<std::complex<float>, MemorySpace::Default>>);
#      endif
    template class Vector<std::complex<float>, MemorySpace::Host>;
    template class Vector<std::complex<float>, MemorySpace::Default>;
    template Vector<std::complex<float>, MemorySpace::Host> &
    Vector<std::complex<float>, MemorySpace::Host>::operator=
      <std::complex<float>>(const dealii::Vector<std::complex<float>> &);
    template Vector<std::complex<float>, MemorySpace::Default> &
    Vector<std::complex<float>, MemorySpace::Default>::operator=
      <std::complex<float>>(const dealii::Vector<std::complex<float>> &);
    namespace internal
    {
      template class VectorReference<std::complex<float>, MemorySpace::Host>;
      template class VectorReference<std::complex<float>, MemorySpace::Default>;
    } // namespace internal
#    endif

#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#      ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<
                  Vector<std::complex<double>, MemorySpace::Host>>);
    static_assert(concepts::is_vector_space_vector<
                  Vector<std::complex<double>, MemorySpace::Default>>);
#      endif
    template class Vector<std::complex<double>, MemorySpace::Host>;
    template class Vector<std::complex<double>, MemorySpace::Default>;
    template Vector<std::complex<double>, MemorySpace::Host> &
    Vector<std::complex<double>, MemorySpace::Host>::operator=
      <std::complex<double>>(const dealii::Vector<std::complex<double>> &);
    template Vector<std::complex<double>, MemorySpace::Default> &
    Vector<std::complex<double>, MemorySpace::Default>::operator=
      <std::complex<double>>(const dealii::Vector<std::complex<double>> &);
    namespace internal
    {
      template class VectorReference<std::complex<double>, MemorySpace::Host>;
      template class VectorReference<std::complex<double>,
                                     MemorySpace::Default>;
    } // namespace internal
#    endif
#  endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
