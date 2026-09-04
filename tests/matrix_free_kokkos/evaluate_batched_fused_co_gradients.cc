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



// this function tests the correctness of the fused, all-directions-at-once
// Portable::internal::batched::EvaluatorTensorProduct::co_gradients<transpose,
// add>(in, out) overload against the already-established per-direction
// co_gradients<direction, dof_to_quad, add, in_place>() chain -- the same
// sequence FEEvaluationImplTransformToCollocation::evaluate()/integrate()
// build the full gradient from, one direction at a time. transpose = false
// is the forward/evaluate_gradients() case (one scalar field -> dim gradient
// components, matches the dof_to_quad = true chain); transpose = true is the
// adjoint/integrate_gradients() case (dim gradient components -> one scalar
// field, matches the dof_to_quad = false chain, later directions
// accumulating onto the first). n_cells "elements" share one team, laid out
// contiguously per element.

#include <deal.II/base/memory_space.h>

#include <deal.II/lac/read_write_vector.h>

#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

using ScratchView = Kokkos::View<
  double *,
  MemorySpace::Default::kokkos_space::execution_space::scratch_memory_space,
  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

using DeviceView1d = Kokkos::View<double *, MemorySpace::Default::kokkos_space>;
using DeviceView2d =
  Kokkos::View<double **, MemorySpace::Default::kokkos_space>;

template <int dim, int N, bool add>
struct FusedCoGradientsForward
{
  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    // Copy the co_shape_gradients into team scratch memory
    ScratchView s_co_shape_gradients(team_member.team_shmem(), N * N);
    for (int tid = team_member.team_rank(); tid < N * N;
         tid += team_member.team_size())
      s_co_shape_gradients(tid) = co_shape_gradients(tid);
    team_member.team_barrier();

    ScratchView dummy_shape(team_member.team_shmem(), 0);
    ScratchView dummy_scratch(team_member.team_shmem(), 0);

    Portable::internal::batched::EvaluatorTensorProduct<
      Portable::internal::evaluate_general,
      dim,
      N,
      N,
      double>
      evaluator(team_member,
                dummy_shape,
                dummy_shape,
                s_co_shape_gradients,
                dummy_scratch,
                n_cells);

    // Reference: per-direction consecutive evaluation
    evaluator.template co_gradients<0, true, add, false>(
      u, Kokkos::subview(grad_ref, Kokkos::ALL, 0));
    if constexpr (dim > 1)
      evaluator.template co_gradients<1, true, add, false>(
        u, Kokkos::subview(grad_ref, Kokkos::ALL, 1));
    if constexpr (dim > 2)
      evaluator.template co_gradients<2, true, add, false>(
        u, Kokkos::subview(grad_ref, Kokkos::ALL, 2));

    // Fused
    evaluator.template co_gradients<false, add>(u, grad_fused);
  }

  DeviceView1d co_shape_gradients;
  DeviceView1d u;
  DeviceView2d grad_ref;
  DeviceView2d grad_fused;
  int          n_cells;
};

template <int dim, int N, bool add>
struct FusedCoGradientsBackward
{
  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    // Copy the co_shape_gradients into team scratch memory
    ScratchView s_co_shape_gradients(team_member.team_shmem(), N * N);
    for (int tid = team_member.team_rank(); tid < N * N;
         tid += team_member.team_size())
      s_co_shape_gradients(tid) = co_shape_gradients(tid);
    team_member.team_barrier();

    ScratchView dummy_shape(team_member.team_shmem(), 0);
    ScratchView dummy_scratch(team_member.team_shmem(), 0);

    Portable::internal::batched::EvaluatorTensorProduct<
      Portable::internal::evaluate_general,
      dim,
      N,
      N,
      double>
      evaluator(team_member,
                dummy_shape,
                dummy_shape,
                s_co_shape_gradients,
                dummy_scratch,
                n_cells);

    // Reference: per-direction consecutive evaluation
    if constexpr (dim == 1)
      {
        evaluator.template co_gradients<0, false, add, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 0), values_ref);
      }
    else if constexpr (dim == 2)
      {
        evaluator.template co_gradients<1, false, add, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 1), values_ref);
        evaluator.template co_gradients<0, false, true, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 0), values_ref);
      }
    else // dim == 3
      {
        evaluator.template co_gradients<2, false, add, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 2), values_ref);
        evaluator.template co_gradients<1, false, true, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 1), values_ref);
        evaluator.template co_gradients<0, false, true, false>(
          Kokkos::subview(grad_in, Kokkos::ALL, 0), values_ref);
      }

    // Fused.
    evaluator.template co_gradients<true, add>(grad_in, values_fused);
  }

  DeviceView1d co_shape_gradients;
  DeviceView2d grad_in;
  DeviceView1d values_ref;
  DeviceView1d values_fused;
  int          n_cells;
};

template <int dim, int N, bool add>
void
test(const int n_cells)
{
  constexpr int      n_q_points = Utilities::pow(N, dim);
  const unsigned int size_batch = n_cells * n_q_points;

  deallog << "Test dim = " << dim << ", N = " << N << ", n_cells = " << n_cells
          << std::endl;

  LinearAlgebra::ReadWriteVector<double> co_shape_host(N * N);
  for (unsigned int i = 0; i < N * N; ++i)
    co_shape_host[i] =
      -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;

  DeviceView1d co_shape_gradients(
    Kokkos::view_alloc("co_shape_gradients", Kokkos::WithoutInitializing),
    N * N);
  {
    Kokkos::View<double *,
                 MemorySpace::Host::kokkos_space,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      co_shape_host_view(co_shape_host.begin(), N * N);
    Kokkos::deep_copy(co_shape_gradients, co_shape_host_view);
  }

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, 1, Kokkos::AUTO);
  // Scratch for the staged co_shape_gradients; dummy_shape/
  // dummy_scratch above need exactly zero.
  team_policy.set_scratch_size(0, Kokkos::PerTeam(N * N * sizeof(double)));

  // ---- forward: evaluate_gradients (transpose = false) ----
  {
    DeviceView1d u(Kokkos::view_alloc("u", Kokkos::WithoutInitializing),
                   size_batch);
    auto         u_host = Kokkos::create_mirror_view(u);
    for (unsigned int i = 0; i < size_batch; ++i)
      u_host(i) = static_cast<double>(Testing::rand()) / RAND_MAX;
    Kokkos::deep_copy(u, u_host);

    DeviceView2d grad_ref(Kokkos::view_alloc("grad_ref",
                                             Kokkos::WithoutInitializing),
                          size_batch,
                          dim);
    DeviceView2d grad_fused(Kokkos::view_alloc("grad_fused",
                                               Kokkos::WithoutInitializing),
                            size_batch,
                            dim);
    auto         grad_ref_host   = Kokkos::create_mirror_view(grad_ref);
    auto         grad_fused_host = Kokkos::create_mirror_view(grad_fused);
    for (unsigned int i = 0; i < size_batch; ++i)
      for (int d = 0; d < dim; ++d)
        {
          const double val =
            add ? static_cast<double>(Testing::rand()) / RAND_MAX : 0.;
          grad_ref_host(i, d)   = val;
          grad_fused_host(i, d) = val;
        }
    Kokkos::deep_copy(grad_ref, grad_ref_host);
    Kokkos::deep_copy(grad_fused, grad_fused_host);

    FusedCoGradientsForward<dim, N, add> functor{
      co_shape_gradients, u, grad_ref, grad_fused, n_cells};
    Kokkos::parallel_for(team_policy, functor);

    Kokkos::deep_copy(grad_ref_host, grad_ref);
    Kokkos::deep_copy(grad_fused_host, grad_fused);

    deallog << "Errors evaluate_gradients: ";
    for (unsigned int i = 0; i < size_batch; ++i)
      for (int d = 0; d < dim; ++d)
        deallog << grad_fused_host(i, d) - grad_ref_host(i, d) << " ";
    deallog << std::endl;
  }

  // ---- backward: integrate_gradients (transpose = true) ----
  {
    DeviceView2d grad_in(Kokkos::view_alloc("grad_in",
                                            Kokkos::WithoutInitializing),
                         size_batch,
                         dim);
    auto         grad_in_host = Kokkos::create_mirror_view(grad_in);
    for (unsigned int i = 0; i < size_batch; ++i)
      for (int d = 0; d < dim; ++d)
        grad_in_host(i, d) = static_cast<double>(Testing::rand()) / RAND_MAX;
    Kokkos::deep_copy(grad_in, grad_in_host);

    DeviceView1d values_ref(Kokkos::view_alloc("values_ref",
                                               Kokkos::WithoutInitializing),
                            size_batch);
    DeviceView1d values_fused(Kokkos::view_alloc("values_fused",
                                                 Kokkos::WithoutInitializing),
                              size_batch);
    auto         values_ref_host   = Kokkos::create_mirror_view(values_ref);
    auto         values_fused_host = Kokkos::create_mirror_view(values_fused);
    for (unsigned int i = 0; i < size_batch; ++i)
      {
        const double val =
          add ? static_cast<double>(Testing::rand()) / RAND_MAX : 0.;
        values_ref_host(i)   = val;
        values_fused_host(i) = val;
      }
    Kokkos::deep_copy(values_ref, values_ref_host);
    Kokkos::deep_copy(values_fused, values_fused_host);

    FusedCoGradientsBackward<dim, N, add> functor{
      co_shape_gradients, grad_in, values_ref, values_fused, n_cells};
    Kokkos::parallel_for(team_policy, functor);

    Kokkos::deep_copy(values_ref_host, values_ref);
    Kokkos::deep_copy(values_fused_host, values_fused);

    deallog << "Errors integrate_gradients: ";
    for (unsigned int i = 0; i < size_batch; ++i)
      deallog << values_fused_host(i) - values_ref_host(i) << " ";
    deallog << std::endl;
  }
}


template <int dim, int N, bool add>
void
test_n_cells()
{
  // test a batch with 1 cell
  test<dim, N, add>(1);
  // test a batch with 3 cells
  test<dim, N, add>(3);
}

int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Kokkos::initialize();

  deallog.push("1d");
  test_n_cells<1, 3, false>();
  test_n_cells<1, 4, false>();
  test_n_cells<1, 3, true>();
  deallog.pop();

  deallog.push("2d");
  test_n_cells<2, 3, false>();
  test_n_cells<2, 4, false>();
  test_n_cells<2, 3, true>();
  deallog.pop();

  deallog.push("3d");
  test_n_cells<3, 3, false>();
  test_n_cells<3, 4, false>();
  test_n_cells<3, 3, true>();
  deallog.pop();

  Kokkos::finalize();

  return 0;
}
