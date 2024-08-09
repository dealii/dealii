// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II Authors
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

// Tests the GinkgoWrappers::Csr interface

#include "deal.II/lac/dynamic_sparsity_pattern.h"

#include "../tests.h"

#include "test_macros.h"

// all include files you need here
#include <deal.II/lac/ginkgo_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>


auto exec = gko::ReferenceExecutor::create();


TEST(can_create_coo)
{
  GinkgoWrappers::Coo<double> m(exec, 3, 5);

  m.compress();

  TEST_ASSERT(m.m() == 3);
  TEST_ASSERT(m.n() == 5);
}


TEST(can_create_ell)
{
  GinkgoWrappers::Ell<double> m(exec, 3, 5);

  m.compress();

  TEST_ASSERT(m.m() == 3);
  TEST_ASSERT(m.n() == 5);
}


TEST(can_create_hybrid)
{
  GinkgoWrappers::Hybrid<double> m(exec, 3, 5);

  m.compress();

  TEST_ASSERT(m.m() == 3);
  TEST_ASSERT(m.n() == 5);
}


TEST(can_create_sellp)
{
  GinkgoWrappers::Sellp<double> m(exec, 3, 5);

  m.compress();

  TEST_ASSERT(m.m() == 3);
  TEST_ASSERT(m.n() == 5);
}



int
main()
{
  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  can_create_coo();
  can_create_ell();
  can_create_hybrid();
  can_create_sellp();

  return 0;
}
