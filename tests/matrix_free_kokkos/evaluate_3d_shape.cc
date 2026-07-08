// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// this function tests the correctness of the 3d evaluation functions used in
// Portable::FEEvaluation. These functions are marked 'internal' but it is
// much easier to check their correctness directly rather than from the results
// in dependent functions

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/read_write_vector.h>

#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

// FIXME Inlining the functor creates "invalid device code" errors
template <int M, int N, int type, bool add, bool dof_to_quad>
struct EvaluateTensorProduct
{
  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    Kokkos::View<
      double *,
      MemorySpace::Default::kokkos_space::execution_space::scratch_memory_space,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      dummy_scratch(team_member.team_shmem(), 0);

    Portable::internal::EvaluatorTensorProduct<
      Portable::internal::evaluate_general,
      3,
      M,
      N,
      double>
      evaluator(team_member,
                shape_values,
                shape_gradients,
                co_shape_gradients,
                dummy_scratch);

    // Direction order matches portable_evaluation_kernels.h: 0,1,2 for
    // dof-to-quadrature (evaluate) and 2,1,0 for quadrature-to-dof (integrate).
    constexpr int d0 = dof_to_quad ? 0 : 2;
    constexpr int d1 = 1;
    constexpr int d2 = dof_to_quad ? 2 : 0;

    if (type == 0)
      {
        evaluator.template values<d0, dof_to_quad, false, false>(src, tmp1);
        team_member.team_barrier();
        evaluator.template values<d1, dof_to_quad, false, false>(tmp1, tmp2);
        team_member.team_barrier();
        evaluator.template values<d2, dof_to_quad, add, false>(tmp2, dst);
      }
    if (type == 1)
      {
        evaluator.template gradients<d0, dof_to_quad, false, false>(src, tmp1);
        team_member.team_barrier();
        evaluator.template gradients<d1, dof_to_quad, false, false>(tmp1, tmp2);
        team_member.team_barrier();
        evaluator.template gradients<d2, dof_to_quad, add, false>(tmp2, dst);
      }
  }

  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_values;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_gradients;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> co_shape_gradients;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> dst;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> src;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> tmp1;
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> tmp2;
};

template <int M, int N, int type, bool add>
void
test()
{
  deallog << "Test " << M << " x " << N << std::endl;
  unsigned int                           size_shape_values = M * N;
  LinearAlgebra::ReadWriteVector<double> shape_host(size_shape_values);
  for (unsigned int i = 0; i < (M + 1) / 2; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        shape_host[i * N + j] =
          -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;
        if (type == 1)
          shape_host[(M - 1 - i) * N + N - 1 - j] = -shape_host[i * N + j];
        else
          shape_host[(M - 1 - i) * N + N - 1 - j] = shape_host[i * N + j];
      }
  if (type == 0 && M % 2 == 1 && N % 2 == 1)
    {
      for (unsigned int i = 0; i < M; ++i)
        shape_host[i * N + N / 2] = 0.;
      shape_host[M / 2 * N + N / 2] = 1.;
    }
  if (type == 1 && M % 2 == 1 && N % 2 == 1)
    shape_host[M / 2 * N + N / 2] = 0.;

  constexpr int                          M_3d = M * M * M;
  constexpr int                          N_3d = N * N * N;
  constexpr int                          MN   = M * N;
  constexpr int                          M2N  = M * M * N;
  LinearAlgebra::ReadWriteVector<double> x_ref(N_3d), y_ref(M_3d);
  Kokkos::View<double[N_3d], MemorySpace::Default::kokkos_space> x_dev(
    Kokkos::view_alloc("x_dev", Kokkos::WithoutInitializing));
  Kokkos::View<double[M_3d], MemorySpace::Default::kokkos_space> y_dev(
    Kokkos::view_alloc("y_dev", Kokkos::WithoutInitializing));
  Kokkos::View<double[M2N], MemorySpace::Default::kokkos_space> tmp1_dev(
    Kokkos::view_alloc("tmp1_dev", Kokkos::WithoutInitializing));
  Kokkos::View<double[MN * N], MemorySpace::Default::kokkos_space> tmp2_dev(
    Kokkos::view_alloc("tmp2_dev", Kokkos::WithoutInitializing));
  auto x_host = Kokkos::create_mirror_view(x_dev);
  auto y_host = Kokkos::create_mirror_view(y_dev);

  for (unsigned int i = 0; i < N_3d; ++i)
    x_host(i) = static_cast<double>(Testing::rand()) / RAND_MAX;

  // The 3d operator is the Kronecker product S (x) S (x) S of the 1d shape
  // matrix S with itself.
  FullMatrix<double> shape_3d(M_3d, N_3d);
  for (unsigned int m0 = 0; m0 < M; ++m0)
    for (unsigned int m1 = 0; m1 < M; ++m1)
      for (unsigned int m2 = 0; m2 < M; ++m2)
        for (unsigned int n0 = 0; n0 < N; ++n0)
          for (unsigned int n1 = 0; n1 < N; ++n1)
            for (unsigned int n2 = 0; n2 < N; ++n2)
              shape_3d(m0 * M * M + m1 * M + m2, n0 * N * N + n1 * N + n2) =
                shape_host[m0 * N + n0] * shape_host[m1 * N + n1] *
                shape_host[m2 * N + n2];

  // Compute reference
  for (unsigned int i = 0; i < M_3d; ++i)
    {
      y_host(i) = 1.;
      y_ref[i]  = add ? y_host(i) : 0.;
      for (unsigned int j = 0; j < N_3d; ++j)
        y_ref[i] += shape_3d(i, j) * x_host(j);
    }

  // Copy data to the GPU.
  Kokkos::deep_copy(x_dev, x_host);
  Kokkos::deep_copy(y_dev, y_host);

  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_values(
    Kokkos::view_alloc("shape_values", Kokkos::WithoutInitializing),
    size_shape_values);
  Kokkos::View<double *,
               MemorySpace::Host::kokkos_space,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    shape_host_view(shape_host.begin(), size_shape_values);
  Kokkos::deep_copy(shape_values, shape_host_view);

  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_gradients(
    Kokkos::view_alloc("shape_gradients", Kokkos::WithoutInitializing),
    size_shape_values);
  Kokkos::deep_copy(shape_gradients, shape_host_view);

  Kokkos::View<double *, MemorySpace::Default::kokkos_space> co_shape_gradients(
    Kokkos::view_alloc("co_shape_gradients", Kokkos::WithoutInitializing),
    size_shape_values);
  Kokkos::deep_copy(co_shape_gradients, shape_host_view);

  // Launch the kernel
  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, 1, Kokkos::AUTO);
  EvaluateTensorProduct<M, N, type, add, false> functor_to_dof{
    shape_values,
    shape_gradients,
    co_shape_gradients,
    y_dev,
    x_dev,
    tmp2_dev,
    tmp1_dev};
  Kokkos::parallel_for(team_policy, functor_to_dof);

  // Check the results on the host
  Kokkos::deep_copy(y_host, y_dev);
  deallog << "Errors no transpose: ";

  for (unsigned int i = 0; i < M_3d; ++i)
    deallog << y_host(i) - y_ref[i] << " ";
  deallog << std::endl;

  for (unsigned int i = 0; i < M_3d; ++i)
    y_host(i) = static_cast<double>(Testing::rand()) / RAND_MAX;

  // Copy y_host to the device
  Kokkos::deep_copy(y_dev, y_host);

  // Compute reference
  for (unsigned int i = 0; i < N_3d; ++i)
    {
      x_host(i) = 2.;
      x_ref[i]  = add ? x_host(i) : 0.;
      for (unsigned int j = 0; j < M_3d; ++j)
        x_ref[i] += shape_3d(j, i) * y_host(j);
    }

  // Copy x_host to the device
  Kokkos::deep_copy(x_dev, x_host);

  // Launch the kernel
  EvaluateTensorProduct<M, N, type, add, true> functor_to_quad{
    shape_values,
    shape_gradients,
    co_shape_gradients,
    x_dev,
    y_dev,
    tmp1_dev,
    tmp2_dev};
  Kokkos::parallel_for(team_policy, functor_to_quad);

  // Check the results on the host
  Kokkos::deep_copy(x_host, x_dev);
  deallog << "Errors transpose:    ";
  for (unsigned int i = 0; i < N_3d; ++i)
    deallog << x_host(i) - x_ref[i] << " ";
  deallog << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Kokkos::initialize();

  deallog.push("values");
  test<4, 4, 0, false>();
  test<3, 3, 0, false>();
  test<3, 4, 0, false>();
  test<3, 5, 0, false>();
  deallog.pop();

  deallog.push("gradients");
  test<4, 4, 1, false>();
  test<3, 3, 1, false>();
  test<3, 4, 1, false>();
  test<3, 5, 1, false>();
  deallog.pop();

  deallog.push("add");

  deallog.push("values");
  test<4, 4, 0, true>();
  test<3, 3, 0, true>();
  test<3, 4, 0, true>();
  test<3, 5, 0, true>();
  deallog.pop();

  deallog.push("gradients");
  test<4, 4, 1, true>();
  test<3, 3, 1, true>();
  test<3, 4, 1, true>();
  test<3, 5, 1, true>();
  deallog.pop();

  deallog.pop();

  Kokkos::finalize();

  return 0;
}
