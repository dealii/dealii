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

// Tests the GinkgoWrappers::Csr manipulation (set/add) interface

#include "../tests.h"

#include "test_macros.h"

// all include files you need here
#include <deal.II/lac/ginkgo_sparse_matrix.h>


auto exec       = gko::ReferenceExecutor::create();
using size_type = GinkgoWrappers::Csr<double>::size_type;


TEST(can_set_single_value)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);

  m.set(2, 4, 4.4);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 4.4);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 1);
}

TEST(can_set_full_matrix)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;
  std::vector<size_type> indices{0, 2};

  m.set(indices, k, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 4);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 2);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 0);
  TEST_ASSERT(obj->get_const_values()[3] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 4);
}

TEST(can_set_full_matrix_elide)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;
  std::vector<size_type> indices{0, 2};

  m.set(indices, k, true);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 2);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 3);
}

TEST(can_set_full_matrix_different_row_col_idxs)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;
  std::vector<size_type> row_indices{0, 2};
  std::vector<size_type> col_indices{3, 4};

  m.set(row_indices, col_indices, k, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 4);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 0);
  TEST_ASSERT(obj->get_const_values()[3] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 4);
}

TEST(can_set_row)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  std::vector<double>         values{1, 2};
  std::vector<size_type>      col_indices{3, 4};

  m.set(1, col_indices, values, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 2);
}

TEST(can_set_row_raw)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  std::vector<double>         values{1, 2};
  std::vector<size_type>      col_indices{3, 4};

  m.set(1, col_indices.size(), col_indices.data(), values.data(), false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 2);
}

TEST(can_add_single_value)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);

  m.add(2, 4, 1.1);
  m.add(2, 4, 4.4);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 5.5);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 1);
}

TEST(can_add_full_matrix)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;

  m.add({0, 2}, k, false);
  m.add({1, 2}, k, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 7);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[4] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[5] == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[6] == 2);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 1);
  TEST_ASSERT(obj->get_const_values()[3] == 2);
  TEST_ASSERT(obj->get_const_values()[4] == 0);
  TEST_ASSERT(obj->get_const_values()[5] == 0);
  TEST_ASSERT(obj->get_const_values()[6] == 8);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 7);
}

TEST(can_add_full_matrix_elide)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;

  m.add({0, 2}, k, true);
  m.add({1, 2}, k, true);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 5);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[4] == 2);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 1);
  TEST_ASSERT(obj->get_const_values()[3] == 2);
  TEST_ASSERT(obj->get_const_values()[4] == 8);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 5);
}

TEST(can_add_full_matrix_different_row_col_idxs)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  FullMatrix<double>          k(2, 2);
  k(0, 0) = 1;
  k(0, 1) = 2;
  k(1, 0) = 0;
  k(1, 1) = 4;
  std::vector<size_type> row_indices{0, 2};
  std::vector<size_type> col_indices{3, 4};

  m.add(row_indices, col_indices, k, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 4);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 0);
  TEST_ASSERT(obj->get_const_values()[3] == 4);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 4);
}

TEST(can_add_row)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  std::vector<double>         values{1, 2};
  std::vector<size_type>      col_indices{3, 4};

  m.add(1, col_indices, values, false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 2);
}

TEST(can_add_row_raw)
{
  GinkgoWrappers::Csr<double> m(exec, 3, 5);
  std::vector<double>         values{1, 2};
  std::vector<size_type>      col_indices{3, 4};

  m.add(1, col_indices.size(), col_indices.data(), values.data(), false);
  m.compress();

  auto obj = m.get_gko_object();
  TEST_ASSERT(obj->get_num_stored_elements() == 2);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 4);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 2);
}

int
main()
{
  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  can_set_single_value();
  can_set_full_matrix();
  can_set_full_matrix_elide();
  can_set_full_matrix_different_row_col_idxs();
  can_set_row();
  can_set_row_raw();
  can_add_single_value();
  can_add_full_matrix();
  can_add_full_matrix_elide();
  can_add_full_matrix_different_row_col_idxs();
  can_add_row();
  can_add_row_raw();

  return 0;
}
