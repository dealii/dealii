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

// Tests the GinkgoWrappers::Vector interface for vector space operations that
// rely on reductions

#include "../tests.h"

#include "test_macros.h"

// all include files you need here
#include <deal.II/lac/ginkgo_vector.h>

auto exec = gko::ReferenceExecutor::create();

TEST(can_compute_scalar_product)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});

  auto sp = v * w;

  TEST_ASSERT_NEAR(sp, 20.0, tol<double>);
}

TEST(can_compute_l1_norm)
{
  GinkgoWrappers::Vector<double> v(exec, {1, -2, 3, -4});

  auto norm = v.l1_norm();

  TEST_ASSERT_NEAR(norm, 10.0, tol<double>);
}

TEST(can_compute_l2_norm)
{
  GinkgoWrappers::Vector<double> v(exec, {1, -2, 3, -4});

  auto norm = v.l2_norm();

  TEST_ASSERT_NEAR(norm, std::sqrt(30), tol<double>);
}

TEST(throws_on_linfty_norm)
{
  GinkgoWrappers::Vector<double> v(exec);

  TEST_ASSERT_THROW(v.linfty_norm(), ExcNotImplemented);
}

TEST(can_compute_add_and_dot)
{
  GinkgoWrappers::Vector<double> v(exec, {1, 2, 3, 4});
  GinkgoWrappers::Vector<double> w(exec, {4, 3, 2, 1});
  GinkgoWrappers::Vector<double> u(exec, {2, 1, 4, 3});

  auto sp = v.add_and_dot(3.3, w, u);

  TEST_ASSERT_NEAR(sp, 100.6, tol<double>);
}

TEST(can_compute_all_zero)
{
  {
    GinkgoWrappers::Vector<double> v(exec, {1, 5, 3, 9});

    TEST_ASSERT(v.all_zero() == false);
  }
  {
    GinkgoWrappers::Vector<double> v(exec, {0, 0, 0});

    TEST_ASSERT(v.all_zero() == true);
  }
}

TEST(throws_on_mean_value)
{
  GinkgoWrappers::Vector<double> v(exec);

  TEST_ASSERT_THROW(v.mean_value(), ExcNotImplemented);
}


int
main()
{
  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  can_compute_scalar_product();
  can_compute_l1_norm();
  can_compute_l2_norm();
  throws_on_linfty_norm();
  can_compute_add_and_dot();
  can_compute_all_zero();
  throws_on_mean_value();

  return 0;
}
