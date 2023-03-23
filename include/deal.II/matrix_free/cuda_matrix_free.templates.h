// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2022 by the deal.II authors
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


#ifndef dealii_cuda_matrix_free_templates_h
#define dealii_cuda_matrix_free_templates_h

#include <deal.II/base/config.h>

#include "deal.II/base/memory_space.h"

#include <deal.II/matrix_free/cuda_matrix_free.h>

#include <string>

#ifdef DEAL_II_WITH_CUDA

#  include <deal.II/base/cuda.h>
#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/graph_coloring.h>

#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/fe/fe_dgq.h>
#  include <deal.II/fe/fe_values.h>
#  include <deal.II/fe/mapping_q1.h>

#  include <deal.II/matrix_free/cuda_hanging_nodes_internal.h>
#  include <deal.II/matrix_free/shape_info.h>

#  include <Kokkos_Core.hpp>
#  include <cuda_runtime_api.h>

#  include <cmath>
#  include <functional>


DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  namespace internal
  {
    /**
     * Helper class to (re)initialize MatrixFree object.
     */
    template <int dim, typename Number>
    class ReinitHelper
    {
    public:
      ReinitHelper(
        MatrixFree<dim, Number> *      data,
        const Mapping<dim> &           mapping,
        const FiniteElement<dim, dim> &fe,
        const Quadrature<1> &          quad,
        const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
          &                    shape_info,
        const DoFHandler<dim> &dof_handler,
        const UpdateFlags &    update_flags);

      void
      setup_color_arrays(const unsigned int n_colors);

      void
      setup_cell_arrays(const unsigned int color);

      template <typename CellFilter>
      void
      get_cell_data(
        const CellFilter &                                        cell,
        const unsigned int                                        cell_id,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      void
      alloc_and_copy_arrays(const unsigned int cell);

    private:
      MatrixFree<dim, Number> *data;
      Kokkos::View<types::global_dof_index **,
                   MemorySpace::Default::kokkos_space>
        local_to_global;
      Kokkos::View<Point<dim, Number> **, MemorySpace::Default::kokkos_space>
                                                                  q_points;
      Kokkos::View<Number **, MemorySpace::Default::kokkos_space> JxW;
      Kokkos::View<Number **[dim][dim], MemorySpace::Default::kokkos_space>
        inv_jacobian;
      Kokkos::View<dealii::internal::MatrixFreeFunctions::ConstraintKinds *,
                   MemorySpace::Default::kokkos_space>
        constraint_mask;
      // Local buffer
      std::vector<types::global_dof_index> local_dof_indices;
      FEValues<dim>                        fe_values;
      // Convert the default dof numbering to a lexicographic one
      const std::vector<unsigned int> &    lexicographic_inv;
      std::vector<types::global_dof_index> lexicographic_dof_indices;
      const unsigned int                   fe_degree;
      const unsigned int                   dofs_per_cell;
      const unsigned int                   q_points_per_cell;
      const UpdateFlags &                  update_flags;
      const unsigned int                   padding_length;
      dealii::internal::MatrixFreeFunctions::HangingNodes<dim> hanging_nodes;
    };



    template <int dim, typename Number>
    ReinitHelper<dim, Number>::ReinitHelper(
      MatrixFree<dim, Number> * data,
      const Mapping<dim> &      mapping,
      const FiniteElement<dim> &fe,
      const Quadrature<1> &     quad,
      const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
        &                    shape_info,
      const DoFHandler<dim> &dof_handler,
      const UpdateFlags &    update_flags)
      : data(data)
      , fe_values(mapping,
                  fe,
                  Quadrature<dim>(quad),
                  update_inverse_jacobians | update_quadrature_points |
                    update_values | update_gradients | update_JxW_values)
      , lexicographic_inv(shape_info.lexicographic_numbering)
      , fe_degree(data->fe_degree)
      , dofs_per_cell(data->dofs_per_cell)
      , q_points_per_cell(data->q_points_per_cell)
      , update_flags(update_flags)
      , padding_length(data->get_padding_length())
      , hanging_nodes(dof_handler.get_triangulation())
    {
      local_dof_indices.resize(data->dofs_per_cell);
      lexicographic_dof_indices.resize(dofs_per_cell);
    }



    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::setup_color_arrays(const unsigned int n_colors)
    {
      // We need at least three colors when we are using CUDA-aware MPI and
      // overlapping the communication
      data->n_cells.resize(std::max(n_colors, 3U), 0);
      data->grid_dim.resize(n_colors);
      data->block_dim.resize(n_colors);
      data->local_to_global.resize(n_colors);
      data->constraint_mask.resize(n_colors);

      data->row_start.resize(n_colors);

      if (update_flags & update_quadrature_points)
        data->q_points.resize(n_colors);

      if (update_flags & update_JxW_values)
        data->JxW.resize(n_colors);

      if (update_flags & update_gradients)
        data->inv_jacobian.resize(n_colors);
    }



    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::setup_cell_arrays(const unsigned int color)
    {
      const unsigned int n_cells         = data->n_cells[color];
      const unsigned int cells_per_block = data->cells_per_block;

      // Setup kernel parameters
      const double apply_n_blocks = std::ceil(
        static_cast<double>(n_cells) / static_cast<double>(cells_per_block));
      const auto apply_x_n_blocks =
        static_cast<unsigned int>(std::round(std::sqrt(apply_n_blocks)));
      const auto apply_y_n_blocks = static_cast<unsigned int>(
        std::ceil(apply_n_blocks / static_cast<double>(apply_x_n_blocks)));

      data->grid_dim[color] = dim3(apply_x_n_blocks, apply_y_n_blocks);

      // TODO this should be a templated parameter.
      const unsigned int n_dofs_1d = fe_degree + 1;

      if (dim == 1)
        data->block_dim[color] = dim3(n_dofs_1d * cells_per_block);
      else if (dim == 2)
        data->block_dim[color] = dim3(n_dofs_1d * cells_per_block, n_dofs_1d);
      else
        data->block_dim[color] =
          dim3(n_dofs_1d * cells_per_block, n_dofs_1d, n_dofs_1d);


      local_to_global = Kokkos::View<types::global_dof_index **,
                                     MemorySpace::Default::kokkos_space>(
        Kokkos::view_alloc("local_to_global_" + std::to_string(color),
                           Kokkos::WithoutInitializing),
        n_cells,
        dofs_per_cell);

      if (update_flags & update_quadrature_points)
        q_points = Kokkos::View<Point<dim, Number> **,
                                MemorySpace::Default::kokkos_space>(
          Kokkos::view_alloc("q_points_" + std::to_string(color),
                             Kokkos::WithoutInitializing),
          n_cells,
          q_points_per_cell);

      if (update_flags & update_JxW_values)
        JxW = Kokkos::View<Number **, MemorySpace::Default::kokkos_space>(
          Kokkos::view_alloc("JxW_" + std::to_string(color),
                             Kokkos::WithoutInitializing),
          n_cells,
          dofs_per_cell);

      if (update_flags & update_gradients)
        inv_jacobian =
          Kokkos::View<Number **[dim][dim], MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("inv_jacobian_" + std::to_string(color),
                               Kokkos::WithoutInitializing),
            n_cells,
            dofs_per_cell);

      // Initialize to zero, i.e., unconstrained cell
      constraint_mask =
        Kokkos::View<dealii::internal::MatrixFreeFunctions::ConstraintKinds *,
                     MemorySpace::Default::kokkos_space>(
          "constraint_mask_" + std::to_string(color), n_cells);
    }



    template <int dim, typename Number>
    template <typename CellFilter>
    void
    ReinitHelper<dim, Number>::get_cell_data(
      const CellFilter &                                        cell,
      const unsigned int                                        cell_id,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
    {
      cell->get_dof_indices(local_dof_indices);
      // When using MPI, we need to transform the local_dof_indices, which
      // contains global dof indices, to get local (to the current MPI process)
      // dof indices.
      if (partitioner)
        for (auto &index : local_dof_indices)
          index = partitioner->global_to_local(index);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        lexicographic_dof_indices[i] = local_dof_indices[lexicographic_inv[i]];

      // FIXME too many deep_copy
      auto constraint_mask_host =
        Kokkos::create_mirror_view_and_copy(MemorySpace::Host::kokkos_space{},
                                            constraint_mask);
      const ArrayView<dealii::internal::MatrixFreeFunctions::ConstraintKinds>
        cell_id_view(constraint_mask_host[cell_id]);

      hanging_nodes.setup_constraints(cell,
                                      partitioner,
                                      {lexicographic_inv},
                                      lexicographic_dof_indices,
                                      cell_id_view);
      Kokkos::deep_copy(constraint_mask, constraint_mask_host);

      // FIXME too many deep_copy
      auto local_to_global_host =
        Kokkos::create_mirror_view_and_copy(MemorySpace::Host::kokkos_space{},
                                            local_to_global);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        local_to_global_host(cell_id, i) = lexicographic_dof_indices[i];
      Kokkos::deep_copy(local_to_global, local_to_global_host);

      fe_values.reinit(cell);

      // Quadrature points
      if (update_flags & update_quadrature_points)
        {
          // FIXME too many deep_copy
          auto q_points_host = Kokkos::create_mirror_view_and_copy(
            MemorySpace::Host::kokkos_space{}, q_points);
          const std::vector<Point<dim>> &q_points_vec =
            fe_values.get_quadrature_points();
          for (unsigned int i = 0; i < q_points_per_cell; ++i)
            q_points_host(cell_id, i) = q_points_vec[i];
          Kokkos::deep_copy(q_points, q_points_host);
        }

      if (update_flags & update_JxW_values)
        {
          // FIXME too many deep_copy
          auto JxW_host = Kokkos::create_mirror_view_and_copy(
            MemorySpace::Host::kokkos_space{}, JxW);
          std::vector<double> JxW_values_double = fe_values.get_JxW_values();
          for (unsigned int i = 0; i < q_points_per_cell; ++i)
            JxW_host(cell_id, i) = static_cast<Number>(JxW_values_double[i]);
          Kokkos::deep_copy(JxW, JxW_host);
        }

      if (update_flags & update_gradients)
        {
          // FIXME too many deep_copy
          auto inv_jacobian_host = Kokkos::create_mirror_view_and_copy(
            MemorySpace::Host::kokkos_space{}, inv_jacobian);
          const std::vector<DerivativeForm<1, dim, dim>> &inv_jacobians =
            fe_values.get_inverse_jacobians();
          for (unsigned int i = 0; i < q_points_per_cell; ++i)
            for (unsigned int j = 0; j < dim; ++j)
              for (unsigned int k = 0; k < dim; ++k)
                inv_jacobian_host(cell_id, i, j, k) = inv_jacobians[i][j][k];
          Kokkos::deep_copy(inv_jacobian, inv_jacobian_host);
        }
    }



    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::alloc_and_copy_arrays(const unsigned int color)
    {
      const unsigned int n_cells = data->n_cells[color];

      // Local-to-global mapping
      data->local_to_global[color] = local_to_global;

      // Quadrature points
      if (update_flags & update_quadrature_points)
        {
          data->q_points[color] = q_points;
        }

      // Jacobian determinants/quadrature weights
      if (update_flags & update_JxW_values)
        {
          data->JxW[color] = JxW;
        }

      // Inverse jacobians
      if (update_flags & update_gradients)
        {
          data->inv_jacobian[color] = inv_jacobian;
        }

      data->constraint_mask[color] = constraint_mask;
    }



    template <int dim, typename number>
    std::vector<types::global_dof_index>
    get_conflict_indices(
      const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
        &                              cell,
      const AffineConstraints<number> &constraints)
    {
      std::vector<types::global_dof_index> local_dof_indices(
        cell->get_fe().n_dofs_per_cell());
      cell->get_dof_indices(local_dof_indices);
      constraints.resolve_indices(local_dof_indices);

      return local_dof_indices;
    }



    template <int dim, typename Number, typename Functor>
    __global__ void
    apply_kernel_shmem(Functor                                      func,
                       const typename MatrixFree<dim, Number>::Data gpu_data,
                       const Number *                               src,
                       Number *                                     dst)
    {
      constexpr unsigned int cells_per_block =
        cells_per_block_shmem(dim, Functor::n_dofs_1d - 1);

      constexpr unsigned int n_dofs_per_block =
        cells_per_block * Functor::n_local_dofs;
      constexpr unsigned int n_q_points_per_block =
        cells_per_block * Functor::n_q_points;
      // TODO make use of dynamically allocated shared memory
      __shared__ Number values[n_dofs_per_block];
      __shared__ Number gradients[dim][n_q_points_per_block];

      const unsigned int local_cell = threadIdx.x / Functor::n_dofs_1d;
      const unsigned int cell =
        local_cell + cells_per_block * (blockIdx.x + gridDim.x * blockIdx.y);

      Number *gq[dim];
      for (unsigned int d = 0; d < dim; ++d)
        gq[d] = &gradients[d][local_cell * Functor::n_q_points];

      SharedData<dim, Number> shared_data(
        &values[local_cell * Functor::n_local_dofs], gq);

      if (cell < gpu_data.n_cells)
        func(cell, &gpu_data, &shared_data, src, dst);
    }



    template <int dim, typename Number, typename Functor>
    __global__ void
    evaluate_coeff(Functor                                      func,
                   const typename MatrixFree<dim, Number>::Data gpu_data)
    {
      constexpr unsigned int cells_per_block =
        cells_per_block_shmem(dim, Functor::n_dofs_1d - 1);

      const unsigned int local_cell = threadIdx.x / Functor::n_dofs_1d;
      const unsigned int cell =
        local_cell + cells_per_block * (blockIdx.x + gridDim.x * blockIdx.y);

      if (cell < gpu_data.n_cells)
        func(cell, &gpu_data);
    }



    template <typename VectorType>
    struct VectorLocalSize
    {
      static unsigned int
      get(const VectorType &vec)
      {
        return vec.locally_owned_size();
      }
    };

    template <>
    struct VectorLocalSize<LinearAlgebra::CUDAWrappers::Vector<double>>
    {
      static unsigned int
      get(const LinearAlgebra::CUDAWrappers::Vector<double> &vec)
      {
        return vec.size();
      }
    };

    template <>
    struct VectorLocalSize<LinearAlgebra::CUDAWrappers::Vector<float>>
    {
      static unsigned int
      get(const LinearAlgebra::CUDAWrappers::Vector<float> &vec)
      {
        return vec.size();
      }
    };
  } // namespace internal



  template <int dim, typename Number>
  MatrixFree<dim, Number>::MatrixFree()
    : my_id(-1)
    , n_dofs(0)
    , padding_length(0)
    , dof_handler(nullptr)
  {}



  template <int dim, typename Number>
  template <typename IteratorFiltersType>
  void
  MatrixFree<dim, Number>::reinit(const Mapping<dim> &             mapping,
                                  const DoFHandler<dim> &          dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1> &            quad,
                                  const IteratorFiltersType &iterator_filter,
                                  const AdditionalData &     additional_data)
  {
    const auto &triangulation = dof_handler.get_triangulation();
    if (const auto parallel_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim> *>(
            &triangulation))
      internal_reinit(mapping,
                      dof_handler,
                      constraints,
                      quad,
                      iterator_filter,
                      std::make_shared<const MPI_Comm>(
                        parallel_triangulation->get_communicator()),
                      additional_data);
    else
      internal_reinit(mapping,
                      dof_handler,
                      constraints,
                      quad,
                      iterator_filter,
                      nullptr,
                      additional_data);
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::reinit(const Mapping<dim> &             mapping,
                                  const DoFHandler<dim> &          dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1> &            quad,
                                  const AdditionalData &additional_data)
  {
    IteratorFilters::LocallyOwnedCell locally_owned_cell_filter;
    reinit(mapping,
           dof_handler,
           constraints,
           quad,
           locally_owned_cell_filter,
           additional_data);
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::reinit(const DoFHandler<dim> &          dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1> &            quad,
                                  const AdditionalData &additional_data)
  {
    reinit(StaticMappingQ1<dim>::mapping,
           dof_handler,
           constraints,
           quad,
           additional_data);
  }



  template <int dim, typename Number>
  typename MatrixFree<dim, Number>::Data
  MatrixFree<dim, Number>::get_data(unsigned int color) const
  {
    Data data_copy;
    if (q_points.size() > 0)
      data_copy.q_points = q_points[color];
    if (inv_jacobian.size() > 0)
      data_copy.inv_jacobian = inv_jacobian[color];
    if (JxW.size() > 0)
      data_copy.JxW = JxW[color];
    data_copy.local_to_global    = local_to_global[color];
    data_copy.constraint_mask    = constraint_mask[color];
    data_copy.shape_values       = shape_values;
    data_copy.shape_gradients    = shape_gradients;
    data_copy.co_shape_gradients = co_shape_gradients;
    data_copy.constraint_weights = constraint_weights;
    data_copy.n_cells            = n_cells[color];
    data_copy.padding_length     = padding_length;
    data_copy.row_start          = row_start[color];
    data_copy.use_coloring       = use_coloring;

    return data_copy;
  }



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::copy_constrained_values(const VectorType &src,
                                                   VectorType &      dst) const
  {
    static_assert(
      std::is_same<Number, typename VectorType::value_type>::value,
      "VectorType::value_type and Number should be of the same type.");
    Assert(src.size() == dst.size(),
           ExcMessage("src and dst vectors have different size."));
    // FIXME When using C++17, we can use KOKKOS_CLASS_LAMBDA and this
    // work-around can be removed.
    auto               constr_dofs = constrained_dofs;
    const unsigned int size = internal::VectorLocalSize<VectorType>::get(dst);
    const Number *     src_ptr = src.get_values();
    Number *           dst_ptr = dst.get_values();
    Kokkos::parallel_for(
      "copy_constrained_values",
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        0, n_constrained_dofs),
      KOKKOS_LAMBDA(int dof) {
        // When working with distributed vectors, the constrained dofs are
        // computed for ghosted vectors but we want to copy the values of the
        // constrained dofs of non-ghosted vectors.
        if (constr_dofs[dof] < size)
          dst_ptr[constr_dofs[dof]] = src_ptr[constr_dofs[dof]];
      });
  }



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::set_constrained_values(Number      val,
                                                  VectorType &dst) const
  {
    static_assert(
      std::is_same<Number, typename VectorType::value_type>::value,
      "VectorType::value_type and Number should be of the same type.");
    Number *dst_ptr = dst.get_values();
    // FIXME When using C++17, we can use KOKKOS_CLASS_LAMBDA and this
    // work-around can be removed.
    auto constr_dofs = constrained_dofs;
    // When working with distributed vectors, the constrained dofs are
    // computed for ghosted vectors but we want to set the values of the
    // constrained dofs of non-ghosted vectors.
    const unsigned int size =
      partitioner ? dst.locally_owned_size() : dst.size();
    Kokkos::parallel_for(
      "set_constrained_values",
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        0, n_constrained_dofs),
      KOKKOS_LAMBDA(int dof) {
        if (constr_dofs[dof] < size)
          dst_ptr[constr_dofs[dof]] = val;
      });
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::initialize_dof_vector(
    LinearAlgebra::CUDAWrappers::Vector<Number> &vec) const
  {
    vec.reinit(n_dofs);
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::initialize_dof_vector(
    LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &vec) const
  {
    if (partitioner)
      vec.reinit(partitioner);
    else
      vec.reinit(n_dofs);
  }



  template <int dim, typename Number>
  unsigned int
  MatrixFree<dim, Number>::get_padding_length() const
  {
    return padding_length;
  }



  template <int dim, typename Number>
  template <typename Functor, typename VectorType>
  void
  MatrixFree<dim, Number>::cell_loop(const Functor &   func,
                                     const VectorType &src,
                                     VectorType &      dst) const
  {
    if (partitioner)
      distributed_cell_loop(func, src, dst);
    else
      serial_cell_loop(func, src, dst);
  }



  template <int dim, typename Number>
  template <typename Functor>
  void
  MatrixFree<dim, Number>::evaluate_coefficients(Functor func) const
  {
    for (unsigned int i = 0; i < n_colors; ++i)
      if (n_cells[i] > 0)
        {
          internal::evaluate_coeff<dim, Number, Functor>
            <<<grid_dim[i], block_dim[i]>>>(func, get_data(i));
          AssertCudaKernel();
        }
  }



  template <int dim, typename Number>
  std::size_t
  MatrixFree<dim, Number>::memory_consumption() const
  {
    // First compute the size of n_cells, row_starts, kernel launch parameters,
    // and constrained_dofs
    std::size_t bytes = n_cells.size() * sizeof(unsigned int) * 2 +
                        2 * n_colors * sizeof(dim3) +
                        n_constrained_dofs * sizeof(unsigned int);

    // For each color, add local_to_global, inv_jacobian, JxW, and q_points.
    // FIXME
    for (unsigned int i = 0; i < n_colors; ++i)
      {
        bytes += n_cells[i] * padding_length * sizeof(unsigned int) +
                 n_cells[i] * padding_length * dim * dim * sizeof(Number) +
                 n_cells[i] * padding_length * sizeof(Number) +
                 n_cells[i] * padding_length * sizeof(point_type) +
                 n_cells[i] * sizeof(unsigned int);
      }

    return bytes;
  }



  template <int dim, typename Number>
  template <typename IteratorFiltersType>
  void
  MatrixFree<dim, Number>::internal_reinit(
    const Mapping<dim> &             mapping,
    const DoFHandler<dim> &          dof_handler_,
    const AffineConstraints<Number> &constraints,
    const Quadrature<1> &            quad,
    const IteratorFiltersType &      iterator_filter,
    std::shared_ptr<const MPI_Comm>  comm,
    const AdditionalData             additional_data)
  {
    dof_handler = &dof_handler_;

    if (typeid(Number) == typeid(double))
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    UpdateFlags update_flags = additional_data.mapping_update_flags;
    if (update_flags & update_gradients)
      update_flags |= update_JxW_values;

    this->use_coloring = additional_data.use_coloring;
    this->overlap_communication_computation =
      additional_data.overlap_communication_computation;

    n_dofs = dof_handler->n_dofs();

    const FiniteElement<dim> &fe = dof_handler->get_fe();

    fe_degree = fe.degree;
    // TODO this should be a templated parameter
    const unsigned int n_dofs_1d     = fe_degree + 1;
    const unsigned int n_q_points_1d = quad.size();

    Assert(n_dofs_1d == n_q_points_1d,
           ExcMessage("n_q_points_1d must be equal to fe_degree+1."));

    // Set padding length to the closest power of two larger than or equal to
    // the number of threads.
    padding_length = 1 << static_cast<unsigned int>(
                       std::ceil(dim * std::log2(fe_degree + 1.)));

    dofs_per_cell     = fe.n_dofs_per_cell();
    q_points_per_cell = Utilities::fixed_power<dim>(n_q_points_1d);

    const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info(
      quad, fe);

    unsigned int size_shape_values = n_dofs_1d * n_q_points_1d;

    FE_DGQArbitraryNodes<1> fe_quad_co(quad);
    const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
      shape_info_co(quad, fe_quad_co);

    shape_values = Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
      Kokkos::view_alloc("shape_values", Kokkos::WithoutInitializing),
      size_shape_values);
    auto shape_values_host = Kokkos::create_mirror_view(shape_values);
    for (unsigned int i = 0; i < size_shape_values; ++i)
      {
        shape_values_host[i] = shape_info.data.front().shape_values[i];
      }
    Kokkos::deep_copy(shape_values, shape_values_host);

    if (update_flags & update_gradients)
      {
        shape_gradients =
          Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("shape_gradients", Kokkos::WithoutInitializing),
            size_shape_values);
        auto shape_gradients_host = Kokkos::create_mirror_view(shape_gradients);
        for (unsigned int i = 0; i < size_shape_values; ++i)
          {
            shape_gradients_host[i] =
              shape_info.data.front().shape_gradients[i];
          }
        Kokkos::deep_copy(shape_gradients, shape_gradients_host);


        co_shape_gradients =
          Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("co_shape_gradients",
                               Kokkos::WithoutInitializing),
            size_shape_values);
        auto co_shape_gradients_host =
          Kokkos::create_mirror_view(co_shape_gradients);
        for (unsigned int i = 0; i < size_shape_values; ++i)
          {
            co_shape_gradients_host[i] =
              shape_info_co.data.front().shape_gradients[i];
          }
        Kokkos::deep_copy(co_shape_gradients, co_shape_gradients_host);
      }

    // Setup the number of cells per CUDA thread block
    cells_per_block = cells_per_block_shmem(dim, fe_degree);

    internal::ReinitHelper<dim, Number> helper(
      this, mapping, fe, quad, shape_info, *dof_handler, update_flags);

    const unsigned int constraint_weights_size =
      shape_info.data.front().subface_interpolation_matrices[0].size();
    constraint_weights =
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
        Kokkos::view_alloc("constraint_weights", Kokkos::WithoutInitializing),
        constraint_weights_size);
    auto constraint_weights_host =
      Kokkos::create_mirror_view(constraint_weights);
    for (unsigned int i = 0; i < constraint_weights_size; ++i)
      {
        constraint_weights_host[i] =
          shape_info.data.front().subface_interpolation_matrices[0][i];
      }
    Kokkos::deep_copy(constraint_weights, constraint_weights_host);

    // Create a graph coloring
    CellFilter begin(iterator_filter, dof_handler->begin_active());
    CellFilter end(iterator_filter, dof_handler->end());

    if (begin != end)
      {
        if (additional_data.use_coloring)
          {
            const auto fun = [&](const CellFilter &filter) {
              return internal::get_conflict_indices<dim, Number>(filter,
                                                                 constraints);
            };
            graph = GraphColoring::make_graph_coloring(begin, end, fun);
          }
        else
          {
            graph.clear();
            if (additional_data.overlap_communication_computation)
              {
                // We create one color (1) with the cells on the boundary of the
                // local domain and two colors (0 and 2) with the interior
                // cells.
                graph.resize(3, std::vector<CellFilter>());

                std::vector<bool> ghost_vertices(
                  dof_handler->get_triangulation().n_vertices(), false);

                for (const auto &cell :
                     dof_handler->get_triangulation().active_cell_iterators())
                  if (cell->is_ghost())
                    for (unsigned int i = 0;
                         i < GeometryInfo<dim>::vertices_per_cell;
                         i++)
                      ghost_vertices[cell->vertex_index(i)] = true;

                std::vector<dealii::FilteredIterator<dealii::TriaActiveIterator<
                  dealii::DoFCellAccessor<dim, dim, false>>>>
                  inner_cells;

                for (auto cell = begin; cell != end; ++cell)
                  {
                    bool ghost_vertex = false;

                    for (unsigned int i = 0;
                         i < GeometryInfo<dim>::vertices_per_cell;
                         i++)
                      if (ghost_vertices[cell->vertex_index(i)])
                        {
                          ghost_vertex = true;
                          break;
                        }

                    if (ghost_vertex)
                      graph[1].emplace_back(cell);
                    else
                      inner_cells.emplace_back(cell);
                  }
                for (unsigned i = 0; i < inner_cells.size(); ++i)
                  if (i < inner_cells.size() / 2)
                    graph[0].emplace_back(inner_cells[i]);
                  else
                    graph[2].emplace_back(inner_cells[i]);
              }
            else
              {
                // If we are not using coloring, all the cells belong to the
                // same color.
                graph.resize(1, std::vector<CellFilter>());
                for (auto cell = begin; cell != end; ++cell)
                  graph[0].emplace_back(cell);
              }
          }
      }
    n_colors = graph.size();

    helper.setup_color_arrays(n_colors);

    IndexSet locally_relevant_dofs;
    if (comm)
      {
        locally_relevant_dofs =
          DoFTools::extract_locally_relevant_dofs(*dof_handler);
        partitioner = std::make_shared<Utilities::MPI::Partitioner>(
          dof_handler->locally_owned_dofs(), locally_relevant_dofs, *comm);
      }
    for (unsigned int i = 0; i < n_colors; ++i)
      {
        n_cells[i] = graph[i].size();
        helper.setup_cell_arrays(i);
        typename std::vector<CellFilter>::iterator cell     = graph[i].begin(),
                                                   end_cell = graph[i].end();
        for (unsigned int cell_id = 0; cell != end_cell; ++cell, ++cell_id)
          helper.get_cell_data(*cell, cell_id, partitioner);

        helper.alloc_and_copy_arrays(i);
      }

    // Setup row starts
    if (n_colors > 0)
      row_start[0] = 0;
    for (unsigned int i = 1; i < n_colors; ++i)
      row_start[i] = row_start[i - 1] + n_cells[i - 1] * get_padding_length();

    // Constrained indices
    n_constrained_dofs = constraints.n_constraints();

    if (n_constrained_dofs != 0)
      {
        const auto constraint_n_blocks = static_cast<unsigned int>(
          std::ceil(static_cast<double>(n_constrained_dofs) /
                    static_cast<double>(block_size)));
        const auto constraint_x_n_blocks =
          static_cast<unsigned int>(std::round(std::sqrt(constraint_n_blocks)));
        const auto constraint_y_n_blocks = static_cast<unsigned int>(
          std::ceil(static_cast<double>(constraint_n_blocks) /
                    static_cast<double>(constraint_x_n_blocks)));

        constraint_grid_dim =
          dim3(constraint_x_n_blocks, constraint_y_n_blocks);
        constraint_block_dim = dim3(block_size);

        std::vector<dealii::types::global_dof_index> constrained_dofs_host(
          n_constrained_dofs);

        if (partitioner)
          {
            const unsigned int n_local_dofs =
              locally_relevant_dofs.n_elements();
            unsigned int i_constraint = 0;
            for (unsigned int i = 0; i < n_local_dofs; ++i)
              {
                // is_constrained uses a global dof id but
                // constrained_dofs_host works on the local id
                if (constraints.is_constrained(partitioner->local_to_global(i)))
                  {
                    constrained_dofs_host[i_constraint] = i;
                    ++i_constraint;
                  }
              }
          }
        else
          {
            const unsigned int n_local_dofs = dof_handler->n_dofs();
            unsigned int       i_constraint = 0;
            for (unsigned int i = 0; i < n_local_dofs; ++i)
              {
                if (constraints.is_constrained(i))
                  {
                    constrained_dofs_host[i_constraint] = i;
                    ++i_constraint;
                  }
              }
          }

        constrained_dofs = Kokkos::View<types::global_dof_index *,
                                        MemorySpace::Default::kokkos_space>(
          Kokkos::view_alloc("constrained_dofs", Kokkos::WithoutInitializing),
          n_constrained_dofs);

        Kokkos::View<types::global_dof_index *,
                     MemorySpace::Default::kokkos_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          constrained_dofs_host_view(constrained_dofs_host.data(),
                                     constrained_dofs_host.size());
        Kokkos::deep_copy(constrained_dofs, constrained_dofs_host_view);
      }
  }



  template <int dim, typename Number>
  template <typename Functor, typename VectorType>
  void
  MatrixFree<dim, Number>::serial_cell_loop(const Functor &   func,
                                            const VectorType &src,
                                            VectorType &      dst) const
  {
    // Execute the loop on the cells
    for (unsigned int i = 0; i < n_colors; ++i)
      if (n_cells[i] > 0)
        {
          internal::apply_kernel_shmem<dim, Number, Functor>
            <<<grid_dim[i], block_dim[i]>>>(func,
                                            get_data(i),
                                            src.get_values(),
                                            dst.get_values());
          AssertCudaKernel();
        }
  }



  template <int dim, typename Number>
  template <typename Functor>
  void
  MatrixFree<dim, Number>::distributed_cell_loop(
    const Functor &                                                      func,
    const LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &src,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const
  {
    // in case we have compatible partitioners, we can simply use the provided
    // vectors
    if (src.get_partitioner().get() == partitioner.get() &&
        dst.get_partitioner().get() == partitioner.get())
      {
        // This code is inspired to the code in TaskInfo::loop.
        if (overlap_communication_computation)
          {
            src.update_ghost_values_start(0);
            // In parallel, it's possible that some processors do not own any
            // cells.
            if (n_cells[0] > 0)
              {
                internal::apply_kernel_shmem<dim, Number, Functor>
                  <<<grid_dim[0], block_dim[0]>>>(func,
                                                  get_data(0),
                                                  src.get_values(),
                                                  dst.get_values());
                AssertCudaKernel();
              }
            src.update_ghost_values_finish();

            // In serial this color does not exist because there are no ghost
            // cells
            if (n_cells[1] > 0)
              {
                internal::apply_kernel_shmem<dim, Number, Functor>
                  <<<grid_dim[1], block_dim[1]>>>(func,
                                                  get_data(1),
                                                  src.get_values(),
                                                  dst.get_values());
                AssertCudaKernel();
                // We need a synchronization point because we don't want
                // CUDA-aware MPI to start the MPI communication until the
                // kernel is done.
                cudaDeviceSynchronize();
              }

            dst.compress_start(0, VectorOperation::add);
            // When the mesh is coarse it is possible that some processors do
            // not own any cells
            if (n_cells[2] > 0)
              {
                internal::apply_kernel_shmem<dim, Number, Functor>
                  <<<grid_dim[2], block_dim[2]>>>(func,
                                                  get_data(2),
                                                  src.get_values(),
                                                  dst.get_values());
                AssertCudaKernel();
              }
            dst.compress_finish(VectorOperation::add);
          }
        else
          {
            src.update_ghost_values();

            // Execute the loop on the cells
            for (unsigned int i = 0; i < n_colors; ++i)
              if (n_cells[i] > 0)
                {
                  internal::apply_kernel_shmem<dim, Number, Functor>
                    <<<grid_dim[i], block_dim[i]>>>(func,
                                                    get_data(i),
                                                    src.get_values(),
                                                    dst.get_values());
                }
            dst.compress(VectorOperation::add);
          }
        src.zero_out_ghost_values();
      }
    else
      {
        // Create the ghosted source and the ghosted destination
        LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
          ghosted_src(partitioner);
        LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
          ghosted_dst(ghosted_src);
        ghosted_src = src;
        ghosted_dst = dst;

        // Execute the loop on the cells
        for (unsigned int i = 0; i < n_colors; ++i)
          if (n_cells[i] > 0)
            {
              internal::apply_kernel_shmem<dim, Number, Functor>
                <<<grid_dim[i], block_dim[i]>>>(func,
                                                get_data(i),
                                                ghosted_src.get_values(),
                                                ghosted_dst.get_values());
              AssertCudaKernel();
            }

        // Add the ghosted values
        ghosted_dst.compress(VectorOperation::add);
        dst = ghosted_dst;
      }
  }



  template <int dim, typename Number>
  template <typename Functor>
  void
  MatrixFree<dim, Number>::distributed_cell_loop(
    const Functor &,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &,
    LinearAlgebra::CUDAWrappers::Vector<Number> &) const
  {
    Assert(false, ExcInternalError());
  }
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
