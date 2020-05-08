// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#include <deal.II/matrix_free/cuda_matrix_free.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/graph_coloring.h>

#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/fe/fe_dgq.h>
#  include <deal.II/fe/fe_values.h>

#  include <deal.II/grid/filtered_iterator.h>

#  include <deal.II/matrix_free/cuda_hanging_nodes_internal.h>
#  include <deal.II/matrix_free/shape_info.h>

#  include <cuda_runtime_api.h>

#  include <cmath>
#  include <functional>


DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  namespace internal
  {
    // These variables are stored in the device constant memory.
    constexpr unsigned int max_elem_degree = 10;

    template <typename NumberType>
    using DataArray = NumberType[(max_elem_degree + 1) * (max_elem_degree + 1)];

    __constant__ double
      global_shape_values_d[(max_elem_degree + 1) * (max_elem_degree + 1)];
    __constant__ float
      global_shape_values_f[(max_elem_degree + 1) * (max_elem_degree + 1)];

    template <typename Number>
    __host__ __device__ inline DataArray<Number> &
             get_global_shape_values();

    template <>
    __host__ __device__ inline DataArray<double> &
             get_global_shape_values<double>()
    {
      return global_shape_values_d;
    }

    template <>
    __host__ __device__ inline DataArray<float> &
             get_global_shape_values<float>()
    {
      return global_shape_values_f;
    }

    __constant__ double
      global_shape_gradients_d[(max_elem_degree + 1) * (max_elem_degree + 1)];
    __constant__ float
      global_shape_gradients_f[(max_elem_degree + 1) * (max_elem_degree + 1)];

    template <typename Number>
    __host__ __device__ inline DataArray<Number> &
             get_global_shape_gradients();

    template <>
    __host__ __device__ inline DataArray<double> &
             get_global_shape_gradients<double>()
    {
      return global_shape_gradients_d;
    }

    template <>
    __host__ __device__ inline DataArray<float> &
             get_global_shape_gradients<float>()
    {
      return global_shape_gradients_f;
    }

    // for collocation methods
    __constant__ double global_co_shape_gradients_d[(max_elem_degree + 1) *
                                                    (max_elem_degree + 1)];
    __constant__ float  global_co_shape_gradients_f[(max_elem_degree + 1) *
                                                   (max_elem_degree + 1)];

    template <typename Number>
    __host__ __device__ inline DataArray<Number> &
             get_global_co_shape_gradients();

    template <>
    __host__ __device__ inline DataArray<double> &
             get_global_co_shape_gradients<double>()
    {
      return global_co_shape_gradients_d;
    }

    template <>
    __host__ __device__ inline DataArray<float> &
             get_global_co_shape_gradients<float>()
    {
      return global_co_shape_gradients_f;
    }

    template <typename Number>
    using CUDAVector = ::dealii::LinearAlgebra::CUDAWrappers::Vector<Number>;

    /**
     * Transpose a N x M matrix stored in a one-dimensional array to a M x N
     * matrix stored in a one-dimensional array.
     */
    template <typename Number>
    void
    transpose(const unsigned int N,
              const unsigned     M,
              const Number *     src,
              Number *           dst)
    {
      // src is N X M
      // dst is M X N
      for (unsigned int i = 0; i < N; ++i)
        for (unsigned int j = 0; j < M; ++j)
          dst[j * N + i] = src[i * M + j];
    }



    /**
     * Same as above but the source and the destination are the same vector.
     */
    template <typename Number>
    void
    transpose_in_place(std::vector<Number> &array_host,
                       const unsigned int   n,
                       const unsigned int   m)
    {
      // convert to structure-of-array
      std::vector<Number> old(array_host.size());
      old.swap(array_host);

      transpose(n, m, old.data(), array_host.data());
    }



    /**
     * Allocate an array to the device and copy @p array_host to the device.
     */
    template <typename Number1>
    void
    alloc_and_copy(Number1 **array_device,
                   const ArrayView<const Number1, MemorySpace::Host> array_host,
                   const unsigned int                                n)
    {
      cudaError_t error_code = cudaMalloc(array_device, n * sizeof(Number1));
      AssertCuda(error_code);
      AssertDimension(array_host.size(), n);

      error_code = cudaMemcpy(*array_device,
                              array_host.data(),
                              n * sizeof(Number1),
                              cudaMemcpyHostToDevice);
      AssertCuda(error_code);
    }



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
      // Host data
      std::vector<types::global_dof_index> local_to_global_host;
      std::vector<Point<dim, Number>>      q_points_host;
      std::vector<Number>                  JxW_host;
      std::vector<Number>                  inv_jacobian_host;
      std::vector<unsigned int>            constraint_mask_host;
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
      HangingNodes<dim>                    hanging_nodes;
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
      , fe_degree(data->fe_degree)
      , dofs_per_cell(data->dofs_per_cell)
      , q_points_per_cell(data->q_points_per_cell)
      , fe_values(mapping,
                  fe,
                  Quadrature<dim>(quad),
                  update_inverse_jacobians | update_quadrature_points |
                    update_values | update_gradients | update_JxW_values)
      , lexicographic_inv(shape_info.lexicographic_numbering)
      , update_flags(update_flags)
      , padding_length(data->get_padding_length())
      , hanging_nodes(fe_degree, dof_handler, lexicographic_inv)
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
      const unsigned int apply_x_n_blocks =
        std::round(std::sqrt(apply_n_blocks));
      const unsigned int apply_y_n_blocks =
        std::ceil(apply_n_blocks / static_cast<double>(apply_x_n_blocks));

      data->grid_dim[color] = dim3(apply_x_n_blocks, apply_y_n_blocks);

      // TODO this should be a templated parameter.
      const unsigned int n_dofs_1d = fe_degree + 1;

      if (data->parallelization_scheme ==
          MatrixFree<dim, Number>::parallel_in_elem)
        {
          if (dim == 1)
            data->block_dim[color] = dim3(n_dofs_1d * cells_per_block);
          else if (dim == 2)
            data->block_dim[color] =
              dim3(n_dofs_1d * cells_per_block, n_dofs_1d);
          else
            data->block_dim[color] =
              dim3(n_dofs_1d * cells_per_block, n_dofs_1d, n_dofs_1d);
        }
      else
        data->block_dim[color] = dim3(cells_per_block);

      local_to_global_host.resize(n_cells * padding_length);

      if (update_flags & update_quadrature_points)
        q_points_host.resize(n_cells * padding_length);

      if (update_flags & update_JxW_values)
        JxW_host.resize(n_cells * padding_length);

      if (update_flags & update_gradients)
        inv_jacobian_host.resize(n_cells * padding_length * dim * dim);

      constraint_mask_host.resize(n_cells);
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

      hanging_nodes.setup_constraints(lexicographic_dof_indices,
                                      cell,
                                      partitioner,
                                      constraint_mask_host[cell_id]);

      memcpy(&local_to_global_host[cell_id * padding_length],
             lexicographic_dof_indices.data(),
             dofs_per_cell * sizeof(types::global_dof_index));

      fe_values.reinit(cell);

      // Quadrature points
      if (update_flags & update_quadrature_points)
        {
          const std::vector<Point<dim>> &q_points =
            fe_values.get_quadrature_points();
          std::copy(q_points.begin(),
                    q_points.end(),
                    &q_points_host[cell_id * padding_length]);
        }

      if (update_flags & update_JxW_values)
        {
          std::vector<double> JxW_values_double = fe_values.get_JxW_values();
          const unsigned int  offset            = cell_id * padding_length;
          for (unsigned int i = 0; i < q_points_per_cell; ++i)
            JxW_host[i + offset] = static_cast<Number>(JxW_values_double[i]);
        }

      if (update_flags & update_gradients)
        {
          const std::vector<DerivativeForm<1, dim, dim>> &inv_jacobians =
            fe_values.get_inverse_jacobians();
          std::copy(&inv_jacobians[0][0][0],
                    &inv_jacobians[0][0][0] +
                      q_points_per_cell * sizeof(DerivativeForm<1, dim, dim>) /
                        sizeof(double),
                    &inv_jacobian_host[cell_id * padding_length * dim * dim]);
        }
    }



    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::alloc_and_copy_arrays(const unsigned int color)
    {
      const unsigned int n_cells = data->n_cells[color];

      // Local-to-global mapping
      if (data->parallelization_scheme ==
          MatrixFree<dim, Number>::parallel_over_elem)
        transpose_in_place(local_to_global_host, n_cells, padding_length);

      alloc_and_copy(
        &data->local_to_global[color],
        ArrayView<const types::global_dof_index>(local_to_global_host.data(),
                                                 local_to_global_host.size()),
        n_cells * padding_length);

      // Quadrature points
      if (update_flags & update_quadrature_points)
        {
          if (data->parallelization_scheme ==
              MatrixFree<dim, Number>::parallel_over_elem)
            transpose_in_place(q_points_host, n_cells, padding_length);

          alloc_and_copy(&data->q_points[color],
                         ArrayView<const Point<dim, Number>>(
                           q_points_host.data(), q_points_host.size()),
                         n_cells * padding_length);
        }

      // Jacobian determinants/quadrature weights
      if (update_flags & update_JxW_values)
        {
          if (data->parallelization_scheme ==
              MatrixFree<dim, Number>::parallel_over_elem)
            transpose_in_place(JxW_host, n_cells, padding_length);

          alloc_and_copy(&data->JxW[color],
                         ArrayView<const Number>(JxW_host.data(),
                                                 JxW_host.size()),
                         n_cells * padding_length);
        }

      // Inverse jacobians
      if (update_flags & update_gradients)
        {
          // Reorder so that all J_11 elements are together, all J_12 elements
          // are together, etc., i.e., reorder indices from
          // cell_id*q_points_per_cell*dim*dim + q*dim*dim +i to
          // i*q_points_per_cell*n_cells + cell_id*q_points_per_cell+q
          transpose_in_place(inv_jacobian_host,
                             padding_length * n_cells,
                             dim * dim);

          // Transpose second time means we get the following index order:
          // q*n_cells*dim*dim + i*n_cells + cell_id which is good for an
          // element-level parallelization
          if (data->parallelization_scheme ==
              MatrixFree<dim, Number>::parallel_over_elem)
            transpose_in_place(inv_jacobian_host,
                               n_cells * dim * dim,
                               padding_length);

          alloc_and_copy(&data->inv_jacobian[color],
                         ArrayView<const Number>(inv_jacobian_host.data(),
                                                 inv_jacobian_host.size()),
                         n_cells * dim * dim * padding_length);
        }

      alloc_and_copy(&data->constraint_mask[color],
                     ArrayView<const unsigned int>(constraint_mask_host.data(),
                                                   constraint_mask_host.size()),
                     n_cells);
    }



    template <int dim, typename number>
    std::vector<types::global_dof_index>
    get_conflict_indices(
      const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
        &                              cell,
      const AffineConstraints<number> &constraints)
    {
      std::vector<types::global_dof_index> local_dof_indices(
        cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      constraints.resolve_indices(local_dof_indices);

      return local_dof_indices;
    }



    template <typename Number>
    __global__ void
    copy_constrained_dofs(
      const dealii::types::global_dof_index *constrained_dofs,
      const unsigned int                     n_constrained_dofs,
      const unsigned int                     size,
      const Number *                         src,
      Number *                               dst)
    {
      const unsigned int dof =
        threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
      // When working with distributed vectors, the constrained dofs are
      // computed for ghosted vectors but we want to copy the values of the
      // constrained dofs of non-ghosted vectors.
      if ((dof < n_constrained_dofs) && (constrained_dofs[dof] < size))
        dst[constrained_dofs[dof]] = src[constrained_dofs[dof]];
    }



    template <typename Number>
    __global__ void
    set_constrained_dofs(
      const dealii::types::global_dof_index *constrained_dofs,
      const unsigned int                     n_constrained_dofs,
      const unsigned int                     size,
      Number                                 val,
      Number *                               dst)
    {
      const unsigned int dof =
        threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
      // When working with distributed vectors, the constrained dofs are
      // computed for ghosted vectors but we want to set the values of the
      // constrained dofs of non-ghosted vectors.
      if ((dof < n_constrained_dofs) && (constrained_dofs[dof] < size))
        dst[constrained_dofs[dof]] = val;
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
      for (int d = 0; d < dim; ++d)
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
  } // namespace internal



  template <int dim, typename Number>
  MatrixFree<dim, Number>::MatrixFree()
    : n_dofs(0)
    , constrained_dofs(nullptr)
    , padding_length(0)
  {}



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::reinit(const Mapping<dim> &             mapping,
                                  const DoFHandler<dim> &          dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1> &            quad,
                                  const AdditionalData &additional_data)
  {
    const auto &triangulation = dof_handler.get_triangulation();
    if (const auto parallel_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim> *>(
            &triangulation))
      internal_reinit(mapping,
                      dof_handler,
                      constraints,
                      quad,
                      std::make_shared<const MPI_Comm>(
                        parallel_triangulation->get_communicator()),
                      additional_data);
    else
      internal_reinit(
        mapping, dof_handler, constraints, quad, nullptr, additional_data);
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
  MatrixFree<dim, Number>::Data
  MatrixFree<dim, Number>::get_data(unsigned int color) const
  {
    Data data_copy;
    data_copy.q_points        = q_points[color];
    data_copy.local_to_global = local_to_global[color];
    data_copy.inv_jacobian    = inv_jacobian[color];
    data_copy.JxW             = JxW[color];
    data_copy.constraint_mask = constraint_mask[color];
    data_copy.n_cells         = n_cells[color];
    data_copy.padding_length  = padding_length;
    data_copy.row_start       = row_start[color];
    data_copy.use_coloring    = use_coloring;

    return data_copy;
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::free()
  {
    for (auto &q_points_color_ptr : q_points)
      Utilities::CUDA::free(q_points_color_ptr);
    q_points.clear();

    for (auto &local_to_global_color_ptr : local_to_global)
      Utilities::CUDA::free(local_to_global_color_ptr);
    local_to_global.clear();

    for (auto &inv_jacobian_color_ptr : inv_jacobian)
      Utilities::CUDA::free(inv_jacobian_color_ptr);
    inv_jacobian.clear();

    for (auto &JxW_color_ptr : JxW)
      Utilities::CUDA::free(JxW_color_ptr);
    JxW.clear();

    for (auto &constraint_mask_color_ptr : constraint_mask)
      Utilities::CUDA::free(constraint_mask_color_ptr);
    constraint_mask.clear();

    Utilities::CUDA::free(constrained_dofs);
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
    if (partitioner)
      distributed_copy_constrained_values(src, dst);
    else
      serial_copy_constrained_values(src, dst);
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
    if (partitioner)
      distributed_set_constrained_values(val, dst);
    else
      serial_set_constrained_values(val, dst);
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
  void
  MatrixFree<dim, Number>::internal_reinit(
    const Mapping<dim> &             mapping,
    const DoFHandler<dim> &          dof_handler,
    const AffineConstraints<Number> &constraints,
    const Quadrature<1> &            quad,
    std::shared_ptr<const MPI_Comm>  comm,
    const AdditionalData             additional_data)
  {
    if (typeid(Number) == typeid(double))
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    const UpdateFlags &update_flags = additional_data.mapping_update_flags;

    if (additional_data.parallelization_scheme != parallel_over_elem &&
        additional_data.parallelization_scheme != parallel_in_elem)
      AssertThrow(false, ExcMessage("Invalid parallelization scheme."));

    this->parallelization_scheme = additional_data.parallelization_scheme;
    this->use_coloring           = additional_data.use_coloring;
    this->overlap_communication_computation =
      additional_data.overlap_communication_computation;

    // TODO: only free if we actually need arrays of different length
    free();

    n_dofs = dof_handler.n_dofs();

    const FiniteElement<dim> &fe = dof_handler.get_fe();

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

    dofs_per_cell     = fe.dofs_per_cell;
    q_points_per_cell = std::pow(n_q_points_1d, dim);

    const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info(
      quad, fe);

    unsigned int size_shape_values = n_dofs_1d * n_q_points_1d * sizeof(Number);

    FE_DGQArbitraryNodes<1> fe_quad_co(quad);
    const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
      shape_info_co(quad, fe_quad_co);

    unsigned int size_co_shape_values =
      n_q_points_1d * n_q_points_1d * sizeof(Number);

    cudaError_t cuda_error =
      cudaMemcpyToSymbol(internal::get_global_shape_values<Number>(),
                         shape_info.data.front().shape_values.data(),
                         size_shape_values,
                         0,
                         cudaMemcpyHostToDevice);
    AssertCuda(cuda_error);

    if (update_flags & update_gradients)
      {
        cuda_error =
          cudaMemcpyToSymbol(internal::get_global_shape_gradients<Number>(),
                             shape_info.data.front().shape_gradients.data(),
                             size_shape_values,
                             0,
                             cudaMemcpyHostToDevice);
        AssertCuda(cuda_error);

        cuda_error =
          cudaMemcpyToSymbol(internal::get_global_co_shape_gradients<Number>(),
                             shape_info_co.data.front().shape_gradients.data(),
                             size_co_shape_values,
                             0,
                             cudaMemcpyHostToDevice);
        AssertCuda(cuda_error);
      }

    // Setup the number of cells per CUDA thread block
    cells_per_block = cells_per_block_shmem(dim, fe_degree);

    internal::ReinitHelper<dim, Number> helper(
      this, mapping, fe, quad, shape_info, dof_handler, update_flags);

    // Create a graph coloring
    using CellFilter =
      FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;
    CellFilter begin(IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active());
    CellFilter end(IteratorFilters::LocallyOwnedCell(), dof_handler.end());
    std::vector<std::vector<CellFilter>> graph;

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
            if (additional_data.overlap_communication_computation)
              {
                // We create one color (1) with the cells on the boundary of the
                // local domain and two colors (0 and 2) with the interior
                // cells.
                graph.resize(3, std::vector<CellFilter>());

                std::vector<bool> ghost_vertices(
                  dof_handler.get_triangulation().n_vertices(), false);

                for (const auto cell :
                     dof_handler.get_triangulation().active_cell_iterators())
                  if (cell->is_ghost())
                    for (unsigned int i = 0;
                         i < GeometryInfo<dim>::vertices_per_cell;
                         i++)
                      ghost_vertices[cell->vertex_index(i)] = true;

                std::vector<dealii::FilteredIterator<dealii::TriaActiveIterator<
                  dealii::DoFCellAccessor<dealii::DoFHandler<dim, dim>,
                                          false>>>>
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
                for (unsigned i = 0; i < inner_cells.size(); i++)
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
        DoFTools::extract_locally_relevant_dofs(dof_handler,
                                                locally_relevant_dofs);
        partitioner = std::make_shared<Utilities::MPI::Partitioner>(
          dof_handler.locally_owned_dofs(), locally_relevant_dofs, *comm);
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
        const unsigned int constraint_n_blocks =
          std::ceil(static_cast<double>(n_constrained_dofs) /
                    static_cast<double>(block_size));
        const unsigned int constraint_x_n_blocks =
          std::round(std::sqrt(constraint_n_blocks));
        const unsigned int constraint_y_n_blocks =
          std::ceil(static_cast<double>(constraint_n_blocks) /
                    static_cast<double>(constraint_x_n_blocks));

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
            const unsigned int n_local_dofs = dof_handler.n_dofs();
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

        cuda_error = cudaMalloc(&constrained_dofs,
                                n_constrained_dofs *
                                  sizeof(dealii::types::global_dof_index));
        AssertCuda(cuda_error);

        cuda_error = cudaMemcpy(constrained_dofs,
                                constrained_dofs_host.data(),
                                n_constrained_dofs *
                                  sizeof(dealii::types::global_dof_index),
                                cudaMemcpyHostToDevice);
        AssertCuda(cuda_error);
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
        src.zero_out_ghosts();
      }
    else
      {
        // Create the ghosted source and the ghosted destination
        LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
          ghosted_src(partitioner);
        LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
          ghosted_dst(ghosted_src);
        ghosted_src = src;

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



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::serial_copy_constrained_values(const VectorType &src,
                                                          VectorType &dst) const
  {
    Assert(src.size() == dst.size(),
           ExcMessage("src and dst vectors have different size."));
    internal::copy_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(constrained_dofs,
                                                      n_constrained_dofs,
                                                      src.size(),
                                                      src.get_values(),
                                                      dst.get_values());
    AssertCudaKernel();
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::distributed_copy_constrained_values(
    const LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &src,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const
  {
    Assert(src.size() == dst.size(),
           ExcMessage("src and dst vectors have different local size."));
    internal::copy_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(constrained_dofs,
                                                      n_constrained_dofs,
                                                      src.local_size(),
                                                      src.get_values(),
                                                      dst.get_values());
    AssertCudaKernel();
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::distributed_copy_constrained_values(
    const LinearAlgebra::CUDAWrappers::Vector<Number> &,
    LinearAlgebra::CUDAWrappers::Vector<Number> &) const
  {
    Assert(false, ExcInternalError());
  }



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::serial_set_constrained_values(const Number val,
                                                         VectorType & dst) const
  {
    internal::set_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(constrained_dofs,
                                                      n_constrained_dofs,
                                                      dst.size(),
                                                      val,
                                                      dst.get_values());
    AssertCudaKernel();
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::distributed_set_constrained_values(
    const Number                                                   val,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const
  {
    internal::set_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(constrained_dofs,
                                                      n_constrained_dofs,
                                                      dst.local_size(),
                                                      val,
                                                      dst.get_values());
    AssertCudaKernel();
  }



  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::distributed_set_constrained_values(
    const Number,
    LinearAlgebra::CUDAWrappers::Vector<Number> &) const
  {
    Assert(false, ExcInternalError());
  }
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
