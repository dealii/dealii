// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test operator[] and norm_square of cuda_tensor.

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int rank, int dim, typename Number>
struct InitFunctor;

template <int dim, typename Number>
struct InitFunctor<0, dim, Number>
{
  Kokkos::View<Tensor<0, dim, Number>, MemorySpace::Default::kokkos_space> t;

  KOKKOS_FUNCTION void
  operator()(int) const
  {
    t() = 1.;
  }
};

template <int dim, typename Number>
struct InitFunctor<1, dim, Number>
{
  Kokkos::View<Tensor<1, dim, Number>, MemorySpace::Default::kokkos_space> t;

  KOKKOS_FUNCTION void
  operator()(int k) const
  {
    t()[k] = k + 1.;
  }
};

template <int dim, typename Number>
struct InitFunctor<2, dim, Number>
{
  Kokkos::View<Tensor<2, dim, Number>, MemorySpace::Default::kokkos_space> t;

  KOKKOS_FUNCTION void
  operator()(int k) const
  {
    int i     = k / dim;
    int j     = k % dim;
    t()[i][j] = k + 1.;
  }
};


template <int rank, int dim, typename Number>
void
test_gpu()
{
  const double tolerance = 1.e-8;

  Kokkos::View<Tensor<rank, dim, Number>, MemorySpace::Default::kokkos_space>
    t_dev("t_dev");
  Kokkos::View<Tensor<rank, dim, Number>, MemorySpace::Default::kokkos_space>
    t1_dev("t1_dev");
  Kokkos::View<Tensor<rank, dim, Number>, MemorySpace::Default::kokkos_space>
    t2_dev("t2_dev");

  Tensor<rank, dim, Number> t_host;
  Tensor<rank, dim, Number> t1_host;
  Tensor<rank, dim, Number> t2_host;

  Tensor<rank, dim, Number> reference_host;

  using ExecutionSpace = MemorySpace::Default::kokkos_space::execution_space;
  ExecutionSpace exec;

  // Initialize
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(
                         exec, 0, Utilities::fixed_power<dim>(rank)),
                       InitFunctor<rank, dim, Number>{t_dev});
  Kokkos::deep_copy(reference_host, t_dev);

  // Test multiplication.
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
      t1_dev() = t_dev() * Number(2.);
      t2_dev() = Number(2.) * t_dev();
      t_dev() *= 2.;
    });
  Kokkos::deep_copy(t_host, t_dev);
  Kokkos::deep_copy(t1_host, t1_dev);
  Kokkos::deep_copy(t2_host, t2_dev);

  reference_host *= 2;
  AssertThrow((t_host - reference_host).norm() < tolerance, ExcInternalError());
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());

  deallog << "multiplication OK" << std::endl;

  // Test division.
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
      t1_dev() = t_dev() / Number(2.);
      t_dev() /= 2.;
      t2_dev() = t1_dev();
    });
  Kokkos::deep_copy(t_host, t_dev);
  Kokkos::deep_copy(t1_host, t1_dev);

  reference_host /= 2.;
  AssertThrow((t_host - reference_host).norm() < tolerance, ExcInternalError());
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());

  deallog << "division OK" << std::endl;

  // Test summation
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
      t2_dev() += t_dev();
      t1_dev() = t1_dev() + t_dev();
    });
  Kokkos::deep_copy(t1_host, t1_dev);
  Kokkos::deep_copy(t2_host, t2_dev);

  reference_host *= 2.;
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());


  // Test subtraction
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
      t2_dev() -= t_dev();
      t1_dev() = t1_dev() - t_dev();
    });
  Kokkos::deep_copy(t1_host, t1_dev);
  Kokkos::deep_copy(t2_host, t2_dev);

  reference_host /= 2.;
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());

  // Miscellaneous
  {
    Kokkos::View<Number, MemorySpace::Default::kokkos_space> check_1("check_1");
    Kokkos::View<Number, MemorySpace::Default::kokkos_space> check_2("check_2");
    Kokkos::View<Number, MemorySpace::Default::kokkos_space> check_3("check_3");
    Kokkos::View<Number, MemorySpace::Default::kokkos_space> check_4("check_4");
    Kokkos::View<Number, MemorySpace::Default::kokkos_space> check_5("check_5");

    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
        // constructors
        typename Tensor<rank, dim, Number>::array_type array{};
        Tensor<rank, dim, Number>                      dummy_1(array);
        check_1() = dummy_1.norm_square();
        Tensor<rank, dim, Number> dummy_2;
        check_2()                         = dummy_2.norm_square();
        Tensor<rank, dim, Number> dummy_3 = dummy_2;
        check_3()                         = dummy_3.norm_square();

        // access
        Tensor<rank + 1, dim, Number>   initializer_1;
        const Tensor<rank, dim, Number> dummy_5 = initializer_1[0];
        check_4()                               = dummy_5.norm_square();

        // assignment
        dummy_2   = dummy_3;
        check_5() = dummy_2.norm_square();
      });

    Number check_1_host, check_2_host, check_3_host, check_4_host, check_5_host;

    Kokkos::deep_copy(check_1_host, check_1);
    Kokkos::deep_copy(check_2_host, check_2);
    Kokkos::deep_copy(check_3_host, check_3);
    Kokkos::deep_copy(check_4_host, check_4);
    Kokkos::deep_copy(check_5_host, check_5);

    AssertThrow(std::abs(check_1_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_2_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_3_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_4_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_5_host) < tolerance, ExcInternalError());
  }
}

int
main()
{
  initlog();

  Kokkos::initialize();

  test_gpu<0, 3, double>();
  test_gpu<1, 3, double>();
  test_gpu<2, 3, double>();
  test_gpu<0, 3, float>();
  test_gpu<1, 3, float>();
  test_gpu<2, 3, float>();

  Kokkos::finalize();
}
