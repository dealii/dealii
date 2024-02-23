// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
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

#ifndef DOXYGEN

namespace LinearAlgebra
{
  template void
  ReadWriteVector<float>::import_elements(
    const CUDAWrappers::Vector<float> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<float>::import_elements(
    const distributed::Vector<float, ::dealii::MemorySpace::CUDA> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);

  template void
  ReadWriteVector<double>::import_elements(
    const CUDAWrappers::Vector<double> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<double>::import_elements(
    const distributed::Vector<double, ::dealii::MemorySpace::CUDA> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
} // namespace LinearAlgebra

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE
