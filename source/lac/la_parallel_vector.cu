// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace distributed
  {
    template class Vector<float, ::dealii::MemorySpace::CUDA>;
    template class Vector<double, ::dealii::MemorySpace::CUDA>;
    template void
    Vector<float, ::dealii::MemorySpace::Host>::import<
      ::dealii::MemorySpace::CUDA>(
      const Vector<float, ::dealii::MemorySpace::CUDA> &,
      VectorOperation::values);
    template void
    Vector<double, ::dealii::MemorySpace::Host>::import<
      ::dealii::MemorySpace::CUDA>(
      const Vector<double, ::dealii::MemorySpace::CUDA> &,
      VectorOperation::values);

    template void
    Vector<float, ::dealii::MemorySpace::CUDA>::import<
      ::dealii::MemorySpace::Host>(
      const Vector<float, ::dealii::MemorySpace::Host> &,
      VectorOperation::values);
    template void
    Vector<double, ::dealii::MemorySpace::CUDA>::import<
      ::dealii::MemorySpace::Host>(
      const Vector<double, ::dealii::MemorySpace::Host> &,
      VectorOperation::values);

    template void
    Vector<float, ::dealii::MemorySpace::CUDA>::import<
      ::dealii::MemorySpace::CUDA>(
      const Vector<float, ::dealii::MemorySpace::CUDA> &,
      VectorOperation::values);
    template void
    Vector<double, ::dealii::MemorySpace::CUDA>::import<
      ::dealii::MemorySpace::CUDA>(
      const Vector<double, ::dealii::MemorySpace::CUDA> &,
      VectorOperation::values);
  } // namespace distributed
} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE
