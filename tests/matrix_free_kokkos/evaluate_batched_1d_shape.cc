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



// this function tests the correctness of the 1d evaluation functions used in
// Portable::internal::batched::EvaluatorTensorProduct -- the batched
// counterpart of Portable::internal::EvaluatorTensorProduct already checked
// in evaluate_1d_shape.cc, where n_cells share one team and are
// all transformed by the same shape matrix in a single apply() call, laid
// out contiguously per element.

#include <deal.II/base/memory_space.h>

#include <deal.II/lac/read_write_vector.h>

#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

template <int M, int N, int type, bool add, bool dof_to_quad>
DEAL_II_HOST_DEVICE void
evaluate_tensor_product(
  const TeamHandle                                          &team_member,
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_values,
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> shape_gradients,
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> co_shape_gradients,
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> dst,
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> src,
  const int                                                  n_cells)
{
  using ScratchView = Kokkos::View<
    double *,
    MemorySpace::Default::kokkos_space::execution_space::scratch_memory_space,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Copy the shape matrices into team scratch memory
  ScratchView s_shape_values(team_member.team_shmem(), M * N);
  ScratchView s_shape_gradients(team_member.team_shmem(), M * N);
  ScratchView s_co_shape_gradients(team_member.team_shmem(), M * N);

  for (int tid = team_member.team_rank(); tid < M * N;
       tid += team_member.team_size())
    {
      s_shape_values(tid)       = shape_values(tid);
      s_shape_gradients(tid)    = shape_gradients(tid);
      s_co_shape_gradients(tid) = co_shape_gradients(tid);
    }
  team_member.team_barrier();

  ScratchView dummy_scratch(team_member.team_shmem(), 0);

  Portable::internal::batched::EvaluatorTensorProduct<
    Portable::internal::evaluate_general,
    1,
    M,
    N,
    double>
    evaluator(team_member,
              s_shape_values,
              s_shape_gradients,
              s_co_shape_gradients,
              dummy_scratch,
              n_cells);

  if (type == 0)
    evaluator.template values<0, dof_to_quad, add, false>(src, dst);
  if (type == 1)
    evaluator.template gradients<0, dof_to_quad, add, false>(src, dst);
}

template <int M, int N, int type, bool add>
void
test(const int n_cells)
{
  deallog << "Test " << M << " x " << N << ", n_cells = " << n_cells
          << std::endl;
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

  const unsigned int                     N_batch = n_cells * N;
  const unsigned int                     M_batch = n_cells * M;
  LinearAlgebra::ReadWriteVector<double> x_ref(N_batch), y_ref(M_batch);
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> x_dev(
    Kokkos::view_alloc("x_dev", Kokkos::WithoutInitializing), N_batch);
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> y_dev(
    Kokkos::view_alloc("y_dev", Kokkos::WithoutInitializing), M_batch);
  auto x_host = Kokkos::create_mirror_view(x_dev);
  auto y_host = Kokkos::create_mirror_view(y_dev);

  for (unsigned int i = 0; i < N_batch; ++i)
    x_host(i) = static_cast<double>(Testing::rand()) / RAND_MAX;

  // Compute reference: each of the n_cells is transformed
  // independently by the same 1d shape matrix.
  for (int e = 0; e < n_cells; ++e)
    for (unsigned int i = 0; i < M; ++i)
      {
        y_host(e * M + i) = 1.;
        y_ref[e * M + i]  = add ? y_host(e * M + i) : 0.;
        for (unsigned int j = 0; j < N; ++j)
          y_ref[e * M + i] += shape_host[i * N + j] * x_host(e * N + j);
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
  // Scratch for the 3 staged shape matrices (shape_values/shape_gradients/
  // co_shape_gradients)
  team_policy.set_scratch_size(0, Kokkos::PerTeam(3 * M * N * sizeof(double)));
  Kokkos::parallel_for(
    team_policy, KOKKOS_LAMBDA(const TeamHandle &team_member) {
      evaluate_tensor_product<M, N, type, add, false>(team_member,
                                                      shape_values,
                                                      shape_gradients,
                                                      co_shape_gradients,
                                                      y_dev,
                                                      x_dev,
                                                      n_cells);
    });

  // Check the results on the host
  Kokkos::deep_copy(y_host, y_dev);
  deallog << "Errors no transpose: ";
  for (unsigned int i = 0; i < M_batch; ++i)
    deallog << y_host(i) - y_ref[i] << " ";
  deallog << std::endl;

  for (unsigned int i = 0; i < M_batch; ++i)
    y_host(i) = static_cast<double>(Testing::rand()) / RAND_MAX;

  // Copy y_host to the device
  Kokkos::deep_copy(y_dev, y_host);

  // Compute reference
  for (int e = 0; e < n_cells; ++e)
    for (unsigned int i = 0; i < N; ++i)
      {
        x_host(e * N + i) = 2.;
        x_ref[e * N + i]  = add ? x_host(e * N + i) : 0.;
        for (unsigned int j = 0; j < M; ++j)
          x_ref[e * N + i] += shape_host[j * N + i] * y_host(e * M + j);
      }

  // Copy x_host to the device
  Kokkos::deep_copy(x_dev, x_host);

  // Launch the kernel
  Kokkos::parallel_for(
    team_policy, KOKKOS_LAMBDA(const TeamHandle &team_member) {
      evaluate_tensor_product<M, N, type, add, true>(team_member,
                                                     shape_values,
                                                     shape_gradients,
                                                     co_shape_gradients,
                                                     x_dev,
                                                     y_dev,
                                                     n_cells);
    });

  // Check the results on the host
  Kokkos::deep_copy(x_host, x_dev);
  deallog << "Errors transpose:    ";
  for (unsigned int i = 0; i < N_batch; ++i)
    deallog << x_host(i) - x_ref[i] << " ";
  deallog << std::endl;
}


template <int M, int N, int type, bool add>
void
test_n_cells()
{
  // test a batch with 1 cell
  test<M, N, type, add>(1);
  // test a batch with 3 cells
  test<M, N, type, add>(3);
}

int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Kokkos::initialize();

  deallog.push("values");
  test_n_cells<4, 4, 0, false>();
  test_n_cells<3, 3, 0, false>();
  test_n_cells<3, 4, 0, false>();
  test_n_cells<3, 5, 0, false>();
  deallog.pop();

  deallog.push("gradients");
  test_n_cells<4, 4, 1, false>();
  test_n_cells<3, 3, 1, false>();
  test_n_cells<3, 4, 1, false>();
  test_n_cells<3, 5, 1, false>();
  deallog.pop();

  deallog.push("add");

  deallog.push("values");
  test_n_cells<4, 4, 0, true>();
  test_n_cells<3, 3, 0, true>();
  test_n_cells<3, 4, 0, true>();
  test_n_cells<3, 5, 0, true>();
  deallog.pop();

  deallog.push("gradients");
  test_n_cells<4, 4, 1, true>();
  test_n_cells<3, 3, 1, true>();
  test_n_cells<3, 4, 1, true>();
  test_n_cells<3, 5, 1, true>();
  deallog.pop();

  deallog.pop();

  Kokkos::finalize();

  return 0;
}
