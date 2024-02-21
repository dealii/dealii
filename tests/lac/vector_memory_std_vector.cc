// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/vector_memory.templates.h>

#include <vector>

#include "../tests.h"

// Check that we can create VectorMemory objects with vector-like classes in
// the standard library.

template <typename VectorType>
void
test_std_vector_pointer()
{
  GrowingVectorMemory<VectorType> mem;

  std::vector<typename VectorMemory<VectorType>::Pointer> va;
  va.push_back(typename VectorMemory<VectorType>::Pointer(mem));
  va.emplace_back(mem);
}



int
main()
{
  initlog();
  PrimitiveVectorMemory<std::vector<double>> primitive_memory;
  test_std_vector_pointer<std::vector<double>>();

  deallog << "OK" << std::endl;
}
