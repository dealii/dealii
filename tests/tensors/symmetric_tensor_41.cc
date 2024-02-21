// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that SymmetricTensor eigenector calculations are correct

#include <deal.II/base/symmetric_tensor.h>

#include <array>

#include "../tests.h"

void
check_value(const int    dim,
            const int    index,
            const double expected,
            const double actual,
            const double tol = 1e-12)
{
  const double rel_error = std::abs(expected - actual) / std::abs(actual);
  if (rel_error > tol)
    {
      deallog << "Incorrect eigenvalue calculated: "
              << "Dim " << dim << ", index " << index << ". Expected "
              << expected << ", actual: " << actual
              << ", relative error: " << rel_error << std::endl;
    }
};

template <int dim>
bool
is_unit_vector(const Tensor<1, dim> &v)
{
  return std::abs(v.norm() - 1.0) < 1e-9;
}

template <int dim>
bool
check_orientation(Tensor<1, dim> v1, Tensor<1, dim> v2, const double tol = 1e-9)
{
  v1 /= v1.norm();
  v2 /= v2.norm();
  return std::abs(std::abs(v1 * v2) - 1.0) < tol;
}

template <int dim>
void
check_vector(const int            index,
             const Tensor<1, dim> expected,
             const Tensor<1, dim> actual,
             const double         tol = 1e-12)
{
  const bool orientation = check_orientation(expected, actual);
  const bool unit_vec    = is_unit_vector(actual);
  if (!(orientation & unit_vec))
    {
      deallog << "Incorrect eigenvector calculated: "
              << "Dim " << dim << ", index " << index << ". Expected "
              << expected << ", actual: " << actual
              << ", orientation check: " << orientation
              << ", unit vector check: " << unit_vec << std::endl;
    }
};

void
test_dim_1(const enum SymmetricTensorEigenvectorMethod method,
           const double                                e1,
           const double                                tol = 1e-12)
{
  const unsigned int      dim = 1;
  SymmetricTensor<2, dim> T;
  T[0][0]                  = e1;
  const auto eig_vals_vecs = eigenvectors(T, method);

  check_value(dim, 0, e1, eig_vals_vecs[0].first, tol);
  check_vector(0, Tensor<1, dim>({1}), eig_vals_vecs[0].second);
}

void
test_dim_2(const enum SymmetricTensorEigenvectorMethod method,
           const double                                e1,
           Tensor<1, 2>                                v1,
           const double                                e2,
           const double                                tol = 1e-12)
{
  const unsigned int dim = 2;
  v1 /= v1.norm();
  const Tensor<1, dim> v2 = cross_product_2d(v1);

  Assert(is_unit_vector(v1), ExcMessage("Vector is not of unit length."));
  Assert(is_unit_vector(v2), ExcMessage("Vector is not of unit length."));
  Assert(e1 >= e2, ExcMessage("Input eigenvalue ordering is not correct."));

  const SymmetricTensor<2, dim> T = e1 * symmetrize(outer_product(v1, v1)) +
                                    e2 * symmetrize(outer_product(v2, v2));

  const auto eig_vals_vecs = eigenvectors(T, method);

  check_value(dim, 0, e1, eig_vals_vecs[0].first, tol);
  check_value(dim, 1, e2, eig_vals_vecs[1].first, tol);
  check_vector(0, v1, eig_vals_vecs[0].second);
  check_vector(1, v2, eig_vals_vecs[1].second);
}

void
test_dim_3(const enum SymmetricTensorEigenvectorMethod method,
           const double                                e1,
           Tensor<1, 3>                                v1,
           const double                                e2,
           Tensor<1, 3>                                v2,
           const double                                e3,
           const double                                tol = 1e-12)
{
  const unsigned int dim = 3;

  // Note: We do not necessarily expect the seed directors v1,v2 to be
  // orthogonal
  v1 /= v1.norm();
  v2 /= v2.norm();
  Tensor<1, dim> v3 = cross_product_3d(v1, v2);
  v3 /= v3.norm();
  v1 = cross_product_3d(v2, v3);

  Assert(is_unit_vector(v1), ExcMessage("Vector is not of unit length."));
  Assert(is_unit_vector(v2), ExcMessage("Vector is not of unit length."));
  Assert(is_unit_vector(v3), ExcMessage("Vector is not of unit length."));
  Assert(check_orientation(v2, cross_product_3d(v3, v1)),
         ExcMessage("Vectors are not orthogonal."));
  Assert(e1 >= e2, ExcMessage("Input eigenvalue ordering is not correct."));
  Assert(e2 >= e3, ExcMessage("Input eigenvalue ordering is not correct."));

  const SymmetricTensor<2, dim> T = e1 * symmetrize(outer_product(v1, v1)) +
                                    e2 * symmetrize(outer_product(v2, v2)) +
                                    e3 * symmetrize(outer_product(v3, v3));

  const auto eig_vals_vecs = eigenvectors(T, method);

  check_value(dim, 0, e1, eig_vals_vecs[0].first, tol);
  check_value(dim, 1, e2, eig_vals_vecs[1].first, tol);
  check_value(dim, 2, e3, eig_vals_vecs[2].first, tol);
  check_vector(0, v1, eig_vals_vecs[0].second);
  check_vector(1, v2, eig_vals_vecs[1].second);
  check_vector(2, v3, eig_vals_vecs[2].second);
}


void
run_tests(const enum SymmetricTensorEigenvectorMethod method)
{
  // Dim = 1
  {
    deallog.push("Test 1");
    {
      test_dim_1(method, 3.6);
    }
    deallog.pop();
  }

  // Dim = 2
  {
    // Diagonal
    deallog.push("Test 2a");
    {
      test_dim_2(method, 3.6, Tensor<1, 2>({1, 0}), 2.4);
    }
    deallog.pop();

    // Diagonal (large difference)
    deallog.push("Test 2b");
    {
      test_dim_2(method, 1.2e7, Tensor<1, 2>({1, 0}), -0.2e-8);
    }
    deallog.pop();

    // Diagonal (equal)
    deallog.push("Test 2c");
    {
      test_dim_2(method, 16.7, Tensor<1, 2>({1, 0}), 16.7);
    }
    deallog.pop();

    // Non-diagonal
    deallog.push("Test 2d");
    {
      test_dim_2(method, 115.7, Tensor<1, 2>({1, 1}), 13.6);
    }
    deallog.pop();

    // Non-diagonal (large difference)
    deallog.push("Test 2e");
    {
      const double tol = 1e-10;
      test_dim_2(method, 7.2956e8, Tensor<1, 2>({3, 2}), -5.284e3, tol);
    }
    deallog.pop();
  }

  // Dim = 3
  {
    // Diagonal
    deallog.push("Test 3a");
    {
      test_dim_3(method,
                 3.6,
                 Tensor<1, 3>({1, 0, 0}),
                 2.4,
                 Tensor<1, 3>({0, 1, 0}),
                 1.2);
    }
    deallog.pop();

    // Diagonal (large difference)
    deallog.push("Test 3b");
    {
      test_dim_3(method,
                 1.2e7,
                 Tensor<1, 3>({1, 0, 0}),
                 -0.2e-8,
                 Tensor<1, 3>({0, 1, 0}),
                 -6.5e8);
    }
    deallog.pop();

    // Diagonal (2 equal)
    deallog.push("Test 3c");
    {
      test_dim_3(method,
                 16.7,
                 Tensor<1, 3>({1, 0, 0}),
                 16.7,
                 Tensor<1, 3>({0, 1, 0}),
                 1e-6);
    }
    deallog.pop();

    // Diagonal (3 equal)
    deallog.push("Test 3d");
    {
      test_dim_3(method,
                 4.2,
                 Tensor<1, 3>({1, 0, 0}),
                 4.2,
                 Tensor<1, 3>({0, 1, 0}),
                 4.2);
    }
    deallog.pop();

    // Non-diagonal
    deallog.push("Test 3e");
    {
      test_dim_3(method,
                 115.7,
                 Tensor<1, 3>({1, 1, 1}),
                 13.6,
                 Tensor<1, 3>({-1, 1, -1}),
                 -45.2);
    }
    deallog.pop();

    // Non-diagonal (1 large difference)
    deallog.push("Test 3f");
    {
      const double tol =
        (method == SymmetricTensorEigenvectorMethod::hybrid ?
           1e-9 :
           (method == SymmetricTensorEigenvectorMethod::ql_implicit_shifts ?
              1e-10 :
              5e-11));
      test_dim_3(method,
                 7.2956e8,
                 Tensor<1, 3>({3, 2, 5}),
                 -4.856e3,
                 Tensor<1, 3>({-0.2, 3, 1}),
                 -5.284e3,
                 tol);
    }
    deallog.pop();

    // Non-diagonal (2 large difference)
    deallog.push("Test 3g");
    {
      const double tol = 1.e-7;
      test_dim_3(method,
                 9.274e7,
                 Tensor<1, 3>({2, -0.7, 1.4}),
                 2.59343,
                 Tensor<1, 3>({0.5, -0.22, -1.42}),
                 -5.292e8,
                 tol);
    }
    deallog.pop();
  }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  deallog.push("Hybrid");
  run_tests(SymmetricTensorEigenvectorMethod::hybrid);
  deallog.pop();

  deallog.push("QL");
  run_tests(SymmetricTensorEigenvectorMethod::ql_implicit_shifts);
  deallog.pop();
  //
  deallog.push("Jacobi");
  run_tests(SymmetricTensorEigenvectorMethod::jacobi);
  deallog.pop();
}
