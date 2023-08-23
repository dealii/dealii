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

// Tests the GinkgoWrappers::Vector interface for vector space operations

#include "../tests.h"

#include "test_macros.h"

// all include files you need here
#include <deal.II/lac/ginkgo_vector.h>

auto exec = gko::ReferenceExecutor::create();



TEST(can_set_to_number)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});

  v = 3.3;

  TEST_ASSERT_NEAR(v[0], 3.3, tol<double>);
  TEST_ASSERT_NEAR(v[1], 3.3, tol<double>);
  TEST_ASSERT_NEAR(v[2], 3.3, tol<double>);
  TEST_ASSERT_NEAR(v[3], 3.3, tol<double>);
}

TEST(can_scale_with_number)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});

  v *= 3.3;

  TEST_ASSERT_NEAR(v[0], 3.3, tol<double>);
  TEST_ASSERT_NEAR(v[1], 6.6, tol<double>);
  TEST_ASSERT_NEAR(v[2], 9.9, tol<double>);
  TEST_ASSERT_NEAR(v[3], 13.2, tol<double>);
}

TEST(can_divide_by_number)
{
  GinkgoWrappers::Vector<double> v(exec, {3.3, 6.6, 9.9, 13.2});

  v /= 3.3;

  TEST_ASSERT_NEAR(v[0], 1.0, tol<double>);
  TEST_ASSERT_NEAR(v[1], 2.0, tol<double>);
  TEST_ASSERT_NEAR(v[2], 3.0, tol<double>);
  TEST_ASSERT_NEAR(v[3], 4.0, tol<double>);
}

TEST(can_add_vector)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  v += w;

  for (size_t i = 0; i < v.size(); ++i)
    {
      TEST_ASSERT_NEAR(v[i], 5.0, tol<double>);
    }
}

TEST(can_substract_vector)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  v -= w;

  TEST_ASSERT_NEAR(v[0], -3.0, tol<double>);
  TEST_ASSERT_NEAR(v[1], -1.0, tol<double>);
  TEST_ASSERT_NEAR(v[2], 1.0, tol<double>);
  TEST_ASSERT_NEAR(v[3], 3.0, tol<double>);
}

TEST(can_add_number)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});

  v.add(3.3);

  TEST_ASSERT_NEAR(v[0], 4.3, tol<double>);
  TEST_ASSERT_NEAR(v[1], 5.3, tol<double>);
  TEST_ASSERT_NEAR(v[2], 6.3, tol<double>);
  TEST_ASSERT_NEAR(v[3], 7.3, tol<double>);
}

TEST(can_add_number_vector)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  v.add(3.3, w);

  TEST_ASSERT_NEAR(v[0], 14.2, tol<double>);
  TEST_ASSERT_NEAR(v[1], 11.9, tol<double>);
  TEST_ASSERT_NEAR(v[2], 9.6, tol<double>);
  TEST_ASSERT_NEAR(v[3], 7.3, tol<double>);
}

TEST(can_add_number_vector_number_vector)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});
  GinkgoWrappers::Vector<double> u(exec, {2, 1, 4, 3});

  v.add(1.1, w, 2.2, u);

  TEST_ASSERT_NEAR(v[0], 9.8, tol<double>);
  TEST_ASSERT_NEAR(v[1], 7.5, tol<double>);
  TEST_ASSERT_NEAR(v[2], 14.0, tol<double>);
  TEST_ASSERT_NEAR(v[3], 11.7, tol<double>);
}

TEST(can_sadd_number_number_vector)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  v.sadd(1.1, 2.2, w);

  TEST_ASSERT_NEAR(v[0], 9.9, tol<double>);
  TEST_ASSERT_NEAR(v[1], 8.8, tol<double>);
  TEST_ASSERT_NEAR(v[2], 7.7, tol<double>);
  TEST_ASSERT_NEAR(v[3], 6.6, tol<double>);
}

TEST(can_equ)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  v.equ(2.2, w);

  TEST_ASSERT_NEAR(v[0], 8.8, tol<double>);
  TEST_ASSERT_NEAR(v[1], 6.6, tol<double>);
  TEST_ASSERT_NEAR(v[2], 4.4, tol<double>);
  TEST_ASSERT_NEAR(v[3], 2.2, tol<double>);
}

int
main()
{
  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  can_set_to_number();
  can_scale_with_number();
  can_divide_by_number();
  can_add_vector();
  can_substract_vector();
  can_add_number();
  can_add_number_vector();
  can_add_number_vector_number_vector();
  can_sadd_number_number_vector();
  can_equ();

  return 0;
}
