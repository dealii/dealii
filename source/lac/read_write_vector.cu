// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/read_write_vector.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  template void
  ReadWriteVector<float>::import(
    const CUDAWrappers::Vector<float> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<float>::import(
    const distributed::Vector<float, ::dealii::MemorySpace::CUDA> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);

  template void
  ReadWriteVector<double>::import(
    const CUDAWrappers::Vector<double> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
  template void
  ReadWriteVector<double>::import(
    const distributed::Vector<double, ::dealii::MemorySpace::CUDA> &,
    VectorOperation::values,
    const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &);
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE
