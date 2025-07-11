// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/read_write_vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "lac/read_write_vector.inst"

namespace LinearAlgebra
{
  // do a few functions that currently don't fit the scheme because they have
  // two template arguments that need to be different (the case of same
  // arguments is covered by the default copy constructor and copy operator that
  // is declared separately)

#define TEMPL_COPY_CONSTRUCTOR(S1, S2)                         \
  template ReadWriteVector<S1> &ReadWriteVector<S1>::operator= \
    <S2>(const ReadWriteVector<S2> &)

  TEMPL_COPY_CONSTRUCTOR(double, float);
  TEMPL_COPY_CONSTRUCTOR(float, double);
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  TEMPL_COPY_CONSTRUCTOR(std::complex<double>, std::complex<float>);
  TEMPL_COPY_CONSTRUCTOR(std::complex<float>, std::complex<double>);
#endif

#undef TEMPL_COPY_CONSTRUCTOR

#ifndef DOXYGEN
  template void
  ReadWriteVector<float>::import_elements(
    const distributed::Vector<float, ::dealii::MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<float>::import_elements(
    const distributed::Vector<float, ::dealii::MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);

  template void
  ReadWriteVector<double>::import_elements(
    const distributed::Vector<double, ::dealii::MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<double>::import_elements(
    const distributed::Vector<double, ::dealii::MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
  template void
  ReadWriteVector<std::complex<float>>::import_elements(
    const distributed::Vector<std::complex<float>, ::dealii::MemorySpace::Host>
      &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);

  template void
  ReadWriteVector<std::complex<double>>::import_elements(
    const distributed::Vector<std::complex<double>, ::dealii::MemorySpace::Host>
      &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#  endif



#  ifdef HAVE_TPETRA_INST_FLOAT
  template void
  ReadWriteVector<float>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<float>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#  endif
#  ifdef HAVE_TPETRA_INST_DOUBLE
  template void
  ReadWriteVector<double>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<double>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#  endif
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  template void
  ReadWriteVector<std::complex<float>>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                                MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<std::complex<float>>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<std::complex<float>,
                                                MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#    endif
#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  template void
  ReadWriteVector<std::complex<double>>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                                MemorySpace::Host> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<std::complex<double>>::import_elements(
    const LinearAlgebra::TpetraWrappers::Vector<std::complex<double>,
                                                MemorySpace::Default> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
#    endif
#  endif

#endif
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE
