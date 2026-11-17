// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_tests_evaluate_batched_fe_evaluation_common_h
#define dealii_tests_evaluate_batched_fe_evaluation_common_h

// Common harness for the tests of the batched evaluation kernels
// Portable::internal::batched::FEEvaluationImpl,
// FEEvaluationImplCollocation and FEEvaluationImplTransformToCollocation.
//
// One team processes a batch of n_cells cells: the 1d shape data is staged
// into shared memory through batched::SharedData::reinit(), the dof or
// quadrature point values of all cells are loaded into SharedData::values
// and ::gradients, one kernel is applied and the result is compared on the
// host against the plain tensor product (A_{dim-1} x ... x A_0) applied cell
// by cell with the same random 1d matrices.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>

#include <deal.II/matrix_free/portable_evaluation_kernels.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_matrix_free.templates.h>

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"


enum class BatchedKernel
{
  general,
  collocation,
  transform_to_collocation
};

using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

template <typename T>
using DeviceView = Kokkos::View<T, MemorySpace::Default::kokkos_space>;



// Random M x N 1d shape matrix, stored as in ShapeInfo as [i * N + j], with
// the symmetry of real shape functions about the midpoint of the unit
// interval as in evaluate_batched_*d_shape.cc: symmetric for shape values,
// antisymmetric for shape gradients, and for odd M and N the middle
// quadrature point coincides with the middle node.
inline std::vector<double>
random_shape_matrix(const unsigned int M,
                    const unsigned int N,
                    const bool         antisymmetric)
{
  std::vector<double> shape(M * N);
  for (unsigned int i = 0; i < (M + 1) / 2; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        shape[i * N + j] =
          -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;
        shape[(M - 1 - i) * N + N - 1 - j] =
          (antisymmetric ? -1. : 1.) * shape[i * N + j];
      }
  if (M % 2 == 1 && N % 2 == 1)
    {
      if (antisymmetric)
        shape[M / 2 * N + N / 2] = 0.;
      else
        {
          for (unsigned int i = 0; i < M; ++i)
            shape[i * N + N / 2] = 0.;
          shape[M / 2 * N + N / 2] = 1.;
        }
    }
  return shape;
}



// An element type for which Portable::FEEvaluation dispatches to the given
// kernel; it selects which shape data SharedData::reinit() stages. Note that
// Portable::MatrixFree supports element types up to
// tensor_symmetric_no_collocation only.
template <BatchedKernel kernel>
constexpr ::dealii::internal::MatrixFreeFunctions::ElementType
element_type()
{
  using ElementType = ::dealii::internal::MatrixFreeFunctions::ElementType;
  if (kernel == BatchedKernel::general)
    return ElementType::tensor_symmetric_no_collocation;
  else if (kernel == BatchedKernel::collocation)
    return ElementType::tensor_symmetric_collocation;
  else
    return ElementType::tensor_symmetric;
}



// Stage the shape data and the input of all cells in shared memory, apply one
// kernel and copy the result back to global memory.
template <int dim, int fe_degree, int n_q_points_1d, BatchedKernel kernel>
struct BatchedEvaluation
{
  static_assert(
    kernel != BatchedKernel::collocation || n_q_points_1d == fe_degree + 1,
    "The collocation kernel requires n_q_points_1d == fe_degree + 1");

  using Number = double;
  using PrecomputedData =
    typename Portable::MatrixFree<dim, Number>::PrecomputedData;
  using SharedData  = Portable::internal::batched::SharedData<dim, Number>;
  using BatchedData = Portable::internal::batched::BatchedData<dim, Number>;

  static constexpr int n_q_points = Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    Kokkos::Array<PrecomputedData, Portable::n_max_dof_handlers>
      precomputed_data;
    precomputed_data[0] = gpu_data;

    Kokkos::Array<SharedData, Portable::n_max_dof_handlers> shared_data;
    shared_data[0].reinit(team_member,
                          n_q_points,
                          n_cells,
                          precomputed_data[0]);
    // reinit() stages the shape data in a TeamVectorRange loop without a
    // trailing barrier
    team_member.team_barrier();

    // The kernels may only depend on the size of the cell range, not on
    // where it starts.
    const int         first_cell = 7;
    const BatchedData data{team_member,
                           n_q_points,
                           1,
                           Kokkos::pair<int, int>(first_cell,
                                                  first_cell + n_cells),
                           precomputed_data,
                           shared_data};

    const int n_points = n_cells * n_q_points;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team_member,
                                                 n_points * n_components),
                         [&](const int idx) {
                           const int i                 = idx % n_points;
                           const int c                 = idx / n_points;
                           shared_data[0].values(i, c) = src_values(i, c);
                           for (int d = 0; d < dim; ++d)
                             shared_data[0].gradients(i, d, c) =
                               src_gradients(i, d, c);
                         });
    team_member.team_barrier();

    if constexpr (kernel == BatchedKernel::general)
      apply<Portable::internal::batched::
              FEEvaluationImpl<dim, fe_degree, n_q_points_1d, Number>>(&data);
    else if constexpr (kernel == BatchedKernel::collocation)
      apply<Portable::internal::batched::
              FEEvaluationImplCollocation<dim, fe_degree, Number>>(&data);
    else
      apply<Portable::internal::batched::FEEvaluationImplTransformToCollocation<
        dim,
        fe_degree,
        n_q_points_1d,
        Number>>(&data);
    team_member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team_member,
                                                 n_points * n_components),
                         [&](const int idx) {
                           const int i      = idx % n_points;
                           const int c      = idx / n_points;
                           dst_values(i, c) = shared_data[0].values(i, c);
                           for (int d = 0; d < dim; ++d)
                             dst_gradients(i, d, c) =
                               shared_data[0].gradients(i, d, c);
                         });
  }

  template <typename Impl>
  DEAL_II_HOST_DEVICE void
  apply(const BatchedData *data) const
  {
    if (integrate)
      Impl::integrate(0, n_components, flags, data);
    else
      Impl::evaluate(0, n_components, flags, data);
  }

  PrecomputedData                  gpu_data;
  int                              n_cells;
  unsigned int                     n_components;
  EvaluationFlags::EvaluationFlags flags;
  bool                             integrate;
  DeviceView<Number **>            src_values;
  DeviceView<Number ***>           src_gradients;
  DeviceView<Number **>            dst_values;
  DeviceView<Number ***>           dst_gradients;
};



// out = (A[dim-1] x ... x A[0]) in for each of the n_cells blocks of `in`,
// direction 0 being the fastest running index. All A[d] are n_out x n_in.
template <int dim>
std::vector<double>
apply_tensor_product(const std::array<FullMatrix<double>, dim> &A,
                     const int                                  n_cells,
                     const std::vector<double>                 &in)
{
  const unsigned int n_in       = A[0].n();
  const unsigned int n_out      = A[0].m();
  const unsigned int n_in_cell  = Utilities::pow(n_in, dim);
  const unsigned int n_out_cell = Utilities::pow(n_out, dim);
  AssertDimension(in.size(), n_cells * n_in_cell);

  std::vector<double> out(n_cells * n_out_cell, 0.);
  for (int cell = 0; cell < n_cells; ++cell)
    for (unsigned int o = 0; o < n_out_cell; ++o)
      {
        double sum = 0.;
        for (unsigned int i = 0; i < n_in_cell; ++i)
          {
            double       prod = 1.;
            unsigned int oo = o, ii = i;
            for (int d = 0; d < dim; ++d)
              {
                prod *= A[d](oo % n_out, ii % n_in);
                oo /= n_out;
                ii /= n_in;
              }
            sum += prod * in[cell * n_in_cell + i];
          }
        out[cell * n_out_cell + o] = sum;
      }
  return out;
}



// The 1d matrices of a tensor product that uses `base` in every direction
// except `special` in direction `direction` (none if direction < 0), e.g.,
// the values S (x) S (x) S for directional(S, S, -1) and the derivative in
// y-direction S (x) D (x) S for directional(S, D, 1). Direction 0 is the
// fastest running index, see apply_tensor_product().
template <int dim>
std::array<FullMatrix<double>, dim>
directional(const FullMatrix<double> &base,
            const FullMatrix<double> &special,
            const int                 direction)
{
  std::array<FullMatrix<double>, dim> result;
  for (int d = 0; d < dim; ++d)
    result[d] = (d == direction) ? special : base;
  return result;
}



template <int dim, int fe_degree, int n_q_points_1d, BatchedKernel kernel>
void
test(const int                              n_cells,
     const unsigned int                     n_components,
     const EvaluationFlags::EvaluationFlags flags,
     const bool                             integrate)
{
  using Number     = double;
  using Functor    = BatchedEvaluation<dim, fe_degree, n_q_points_1d, kernel>;
  using SharedData = typename Functor::SharedData;

  constexpr int n_dofs_1d  = fe_degree + 1;
  constexpr int n_dofs     = Utilities::pow(n_dofs_1d, dim);
  constexpr int n_q_points = Utilities::pow(n_q_points_1d, dim);
  const int     n_points   = n_cells * n_q_points;

  const bool do_values =
    (flags & EvaluationFlags::values) != EvaluationFlags::nothing;
  const bool do_gradients =
    (flags & EvaluationFlags::gradients) != EvaluationFlags::nothing;

  deallog << "dim=" << dim << " degree=" << fe_degree
          << " n_q_points_1d=" << n_q_points_1d << " n_cells=" << n_cells
          << " n_components=" << n_components << " "
          << (do_values ? "values" : "")
          << (do_values && do_gradients ? "|" : "")
          << (do_gradients ? "gradients" : "")
          << (integrate ? " integrate:" : " evaluate:");

  const auto random_unit = []() {
    return static_cast<double>(Testing::rand()) / RAND_MAX;
  };

  // 1d shape data: values symmetric, gradients antisymmetric
  const std::vector<Number> shape_values =
    random_shape_matrix(n_dofs_1d, n_q_points_1d, false);
  const std::vector<Number> shape_gradients =
    random_shape_matrix(n_dofs_1d, n_q_points_1d, true);
  const std::vector<Number> co_shape_gradients =
    random_shape_matrix(n_q_points_1d, n_q_points_1d, true);

  const auto to_device = [](const std::vector<Number> &v,
                            const std::string         &name) {
    DeviceView<Number *> d(Kokkos::view_alloc(name,
                                              Kokkos::WithoutInitializing),
                           v.size());
    Kokkos::deep_copy(d,
                      Kokkos::View<const Number *,
                                   Kokkos::HostSpace,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                        v.data(), v.size()));
    return d;
  };

  // The part of PrecomputedData that SharedData::reinit() and the kernels use
  typename Functor::PrecomputedData gpu_data{};
  gpu_data.shape_values    = to_device(shape_values, "shape_values");
  gpu_data.shape_gradients = to_device(shape_gradients, "shape_gradients");
  gpu_data.co_shape_gradients =
    to_device(co_shape_gradients, "co_shape_gradients");
  gpu_data.n_components     = n_components;
  gpu_data.element_type     = element_type<kernel>();
  gpu_data.scratch_pad_size = Portable::internal::compute_scratch_pad_size(
    gpu_data.element_type, dim, fe_degree, n_q_points_1d);

  // Random input: dof values (evaluate) or quadrature point values and
  // gradients (integrate), stored contiguously cell after cell
  DeviceView<Number **>  src_values("src_values", n_points, n_components);
  DeviceView<Number ***> src_gradients("src_gradients",
                                       n_points,
                                       dim,
                                       n_components);
  DeviceView<Number **>  dst_values("dst_values", n_points, n_components);
  DeviceView<Number ***> dst_gradients("dst_gradients",
                                       n_points,
                                       dim,
                                       n_components);

  auto src_values_host    = Kokkos::create_mirror_view(src_values);
  auto src_gradients_host = Kokkos::create_mirror_view(src_gradients);
  for (unsigned int c = 0; c < n_components; ++c)
    for (int i = 0; i < n_points; ++i)
      {
        src_values_host(i, c) = random_unit();
        for (int d = 0; d < dim; ++d)
          src_gradients_host(i, d, c) = random_unit();
      }
  Kokkos::deep_copy(src_values, src_values_host);
  Kokkos::deep_copy(src_gradients, src_gradients_host);

  const std::size_t scratch_bytes =
    SharedData::shmem_size(n_q_points, n_cells, gpu_data);

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, 1, Kokkos::AUTO);
  team_policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));
  Kokkos::parallel_for(team_policy,
                       Functor{gpu_data,
                               n_cells,
                               n_components,
                               flags,
                               integrate,
                               src_values,
                               src_gradients,
                               dst_values,
                               dst_gradients});
  Kokkos::fence();

  auto dst_values_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dst_values);
  auto dst_gradients_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dst_gradients);

  // Reference. S and D map dofs to quadrature points, C maps quadrature
  // points to quadrature points (collocation derivative).
  FullMatrix<double> S(n_q_points_1d, n_dofs_1d);
  FullMatrix<double> D(n_q_points_1d, n_dofs_1d);
  FullMatrix<double> C(n_q_points_1d, n_q_points_1d);
  for (int q = 0; q < n_q_points_1d; ++q)
    {
      for (int i = 0; i < n_dofs_1d; ++i)
        {
          S(q, i) = shape_values[i * n_q_points_1d + q];
          D(q, i) = shape_gradients[i * n_q_points_1d + q];
        }
      for (int j = 0; j < n_q_points_1d; ++j)
        C(q, j) = co_shape_gradients[j * n_q_points_1d + q];
    }
  FullMatrix<double> ST(n_dofs_1d, n_q_points_1d);
  FullMatrix<double> DT(n_dofs_1d, n_q_points_1d);
  FullMatrix<double> CT(n_q_points_1d, n_q_points_1d);
  ST.copy_transposed(S);
  DT.copy_transposed(D);
  CT.copy_transposed(C);
  const IdentityMatrix     identity(n_q_points_1d);
  const FullMatrix<double> I(identity);

  const auto add_to = [](std::vector<double> &a, const std::vector<double> &b) {
    AssertDimension(a.size(), b.size());
    for (unsigned int i = 0; i < a.size(); ++i)
      a[i] += b[i];
  };

  // max |device - ref| and max |ref| over the compared entries
  const auto compare = [](const std::vector<double> &ref,
                          const auto                &device,
                          double                    &max_error,
                          double                    &max_reference) {
    for (unsigned int i = 0; i < ref.size(); ++i)
      {
        max_error     = std::max(max_error, std::abs(device(i) - ref[i]));
        max_reference = std::max(max_reference, std::abs(ref[i]));
      }
  };

  double error_values = 0., reference_values = 0.;
  double error_gradients = 0., reference_gradients = 0.;

  for (unsigned int c = 0; c < n_components; ++c)
    {
      if (!integrate)
        {
          std::vector<double> u(n_cells * n_dofs);
          for (unsigned int i = 0; i < u.size(); ++i)
            u[i] = src_values_host(i, c);

          std::vector<double>                  values;
          std::array<std::vector<double>, dim> gradients;
          if constexpr (kernel == BatchedKernel::general)
            {
              values = apply_tensor_product<dim>(directional<dim>(S, S, -1),
                                                 n_cells,
                                                 u);
              for (int d = 0; d < dim; ++d)
                gradients[d] =
                  apply_tensor_product<dim>(directional<dim>(S, D, d),
                                            n_cells,
                                            u);
            }
          else if constexpr (kernel == BatchedKernel::collocation)
            {
              values = u;
              for (int d = 0; d < dim; ++d)
                gradients[d] =
                  apply_tensor_product<dim>(directional<dim>(I, C, d),
                                            n_cells,
                                            u);
            }
          else
            {
              values = apply_tensor_product<dim>(directional<dim>(S, S, -1),
                                                 n_cells,
                                                 u);
              for (int d = 0; d < dim; ++d)
                gradients[d] =
                  apply_tensor_product<dim>(directional<dim>(I, C, d),
                                            n_cells,
                                            values);
            }

          if (do_values)
            compare(
              values,
              [&](const int i) { return dst_values_host(i, c); },
              error_values,
              reference_values);
          if (do_gradients)
            for (int d = 0; d < dim; ++d)
              compare(
                gradients[d],
                [&](const int i) { return dst_gradients_host(i, d, c); },
                error_gradients,
                reference_gradients);
        }
      else
        {
          std::vector<double>                  v(n_points);
          std::array<std::vector<double>, dim> g;
          for (int i = 0; i < n_points; ++i)
            v[i] = src_values_host(i, c);
          for (int d = 0; d < dim; ++d)
            {
              g[d].resize(n_points);
              for (int i = 0; i < n_points; ++i)
                g[d][i] = src_gradients_host(i, d, c);
            }

          std::vector<double> out;
          if constexpr (kernel == BatchedKernel::general)
            {
              out.assign(n_cells * n_dofs, 0.);
              if (do_values)
                out = apply_tensor_product<dim>(directional<dim>(ST, ST, -1),
                                                n_cells,
                                                v);
              if (do_gradients)
                for (int d = 0; d < dim; ++d)
                  add_to(out,
                         apply_tensor_product<dim>(directional<dim>(ST, DT, d),
                                                   n_cells,
                                                   g[d]));
            }
          else
            {
              // collocation part: quadrature points to quadrature points
              out.assign(n_points, 0.);
              if (do_values)
                out = v;
              if (do_gradients)
                for (int d = 0; d < dim; ++d)
                  add_to(out,
                         apply_tensor_product<dim>(directional<dim>(I, CT, d),
                                                   n_cells,
                                                   g[d]));
              if constexpr (kernel == BatchedKernel::transform_to_collocation)
                out = apply_tensor_product<dim>(directional<dim>(ST, ST, -1),
                                                n_cells,
                                                out);
            }

          compare(
            out,
            [&](const int i) { return dst_values_host(i, c); },
            error_values,
            reference_values);
        }
    }

  // Print the relative errors, rounded to zero below the tolerance
  const auto print =
    [](const char *name, const double error, const double reference) {
      const double relative = error / std::max(reference, 1e-300);
      deallog << ' ' << name << ' ' << (relative < 1e-12 ? 0. : relative);
    };
  if (integrate)
    print("error", error_values, reference_values);
  else
    {
      if (do_values)
        print("values", error_values, reference_values);
      if (do_gradients)
        print("gradients", error_gradients, reference_gradients);
    }
  deallog << std::endl;
}



// All combinations of batch size, evaluation flags and evaluate/integrate
// for one element, plus a two-component case.
template <int dim, int fe_degree, int n_q_points_1d, BatchedKernel kernel>
void
run()
{
  const EvaluationFlags::EvaluationFlags all_flags[] = {
    EvaluationFlags::values,
    EvaluationFlags::gradients,
    EvaluationFlags::values | EvaluationFlags::gradients};

  for (const int n_cells : {1, 3})
    for (const auto flags : all_flags)
      for (const bool integrate : {false, true})
        test<dim, fe_degree, n_q_points_1d, kernel>(n_cells,
                                                    1,
                                                    flags,
                                                    integrate);

  for (const bool integrate : {false, true})
    test<dim, fe_degree, n_q_points_1d, kernel>(
      3, 2, EvaluationFlags::values | EvaluationFlags::gradients, integrate);
}

#endif
