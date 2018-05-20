// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_cuda_matrix_free_templates_h
#define dealii_cuda_matrix_free_templates_h

#include <deal.II/matrix_free/cuda_matrix_free.h>

#ifdef DEAL_II_WITH_CUDA

#  include <cuda_runtime_api.h>
#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/graph_coloring.h>
#  include <deal.II/fe/fe_values.h>
#  include <deal.II/grid/filtered_iterator.h>
#  include <deal.II/matrix_free/shape_info.h>
#  include <functional>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  namespace internal
  {
    // These variables are stored in the device constant memory.
    // TODO: use a template parameter
    constexpr unsigned int max_elem_degree = 10;
    __constant__ double
      global_shape_values[(max_elem_degree + 1) * (max_elem_degree + 1)];
    __constant__ double
      global_shape_gradients[(max_elem_degree + 1) * (max_elem_degree + 1)];

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
              const Number*      src,
              Number*            dst)
    {
      // src is N X M
      // dst is M X N
      for(unsigned int i = 0; i < N; ++i)
        for(unsigned int j = 0; j < M; ++j)
          dst[j * N + i] = src[i * M + j];
    }

    /**
     * Same as above but the source and the destination are the same vector.
     */
    template <typename Number>
    void
    transpose_in_place(std::vector<Number>& array_host,
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
    template <typename Number1, typename Number2>
    void
    alloc_and_copy(Number1**             array_device,
                   std::vector<Number2>& array_host,
                   const unsigned int    n)
    {
      cudaError_t error_code = cudaMalloc(array_device, n * sizeof(Number1));
      AssertCuda(error_code);

      error_code = cudaMemcpy(*array_device,
                              array_host.data(),
                              n * sizeof(Number1),
                              cudaMemcpyHostToDevice);
      AssertCuda(error_code);
    }

    /**
     * Helper class to (re)initialize MatrixFree object.
     */
    //TODO for now does not support hanging_nodes
    template <int dim, typename Number>
    class ReinitHelper
    {
    public:
      ReinitHelper(
        MatrixFree<dim, Number>*       data,
        const Mapping<dim>&            mapping,
        const FiniteElement<dim, dim>& fe,
        const Quadrature<1>&           quad,
        const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>&
                           shape_info,
        const UpdateFlags& update_flags);

      void
      setup_color_arrays(const unsigned int n_colors);

      void
      setup_cell_arrays(const unsigned int color);

      template <typename CellFilter>
      void
      get_cell_data(const CellFilter& cell, const unsigned int cell_id);

      void
      alloc_and_copy_arrays(const unsigned int cell);

    private:
      MatrixFree<dim, Number>* data;
      // Host data
      std::vector<unsigned int> local_to_global_host;
      std::vector<Point<dim>>   q_points_host;
      std::vector<Number>       JxW_host;
      std::vector<Number>       inv_jacobian_host;
      std::vector<unsigned int> constraint_mask_host;
      // Local buffer
      std::vector<types::global_dof_index> local_dof_indices;
      FEValues<dim>                        fe_values;
      // Convert the default dof numbering to a lexicographic one
      const std::vector<unsigned int>& lexicographic_inv;
      std::vector<unsigned int>        lexicographic_dof_indices;
      const unsigned int               fe_degree;
      const unsigned int               dofs_per_cell;
      const unsigned int               q_points_per_cell;
      const UpdateFlags&               update_flags;
      const unsigned int               padding_length;
    };

    template <int dim, typename Number>
    ReinitHelper<dim, Number>::ReinitHelper(
      MatrixFree<dim, Number>*  data,
      const Mapping<dim>&       mapping,
      const FiniteElement<dim>& fe,
      const Quadrature<1>&      quad,
      const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>&
                         shape_info,
      const UpdateFlags& update_flags)
      : data(data),
        fe_degree(data->fe_degree),
        dofs_per_cell(data->dofs_per_cell),
        q_points_per_cell(data->q_points_per_cell),
        fe_values(mapping,
                  fe,
                  Quadrature<dim>(quad),
                  update_inverse_jacobians | update_quadrature_points
                    | update_values | update_gradients | update_JxW_values),
        lexicographic_inv(shape_info.lexicographic_numbering),
        update_flags(update_flags),
        padding_length(data->get_padding_length())
    {
      local_dof_indices.resize(data->dofs_per_cell);
      lexicographic_dof_indices.resize(dofs_per_cell);
    }

    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::setup_color_arrays(const unsigned int n_colors)
    {
      data->n_cells.resize(n_colors);
      data->grid_dim.resize(n_colors);
      data->block_dim.resize(n_colors);
      data->local_to_global.resize(n_colors);
      data->constraint_mask.resize(n_colors);

      data->row_start.resize(n_colors);

      if(update_flags & update_quadrature_points)
        data->q_points.resize(n_colors);

      if(update_flags & update_JxW_values)
        data->JxW.resize(n_colors);

      if(update_flags & update_gradients)
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
      const unsigned int apply_x_n_blocks
        = std::round(std::sqrt(apply_n_blocks));
      const unsigned int apply_y_n_blocks
        = std::ceil(apply_n_blocks / static_cast<double>(apply_x_n_blocks));

      data->grid_dim[color] = dim3(apply_x_n_blocks, apply_y_n_blocks);

      // TODO this should be a templated parameter.
      const unsigned int n_dofs_1d = fe_degree + 1;

      if(data->parallelization_scheme
         == MatrixFree<dim, Number>::parallel_in_elem)
        {
          if(dim == 1)
            data->block_dim[color] = dim3(n_dofs_1d * cells_per_block);
          else if(dim == 2)
            data->block_dim[color]
              = dim3(n_dofs_1d * cells_per_block, n_dofs_1d);
          else
            data->block_dim[color]
              = dim3(n_dofs_1d * cells_per_block, n_dofs_1d, n_dofs_1d);
        }
      else
        data->block_dim[color] = dim3(cells_per_block);

      local_to_global_host.resize(n_cells * padding_length);

      if(update_flags & update_quadrature_points)
        q_points_host.resize(n_cells * padding_length);

      if(update_flags & update_JxW_values)
        JxW_host.resize(n_cells * padding_length);

      if(update_flags & update_gradients)
        inv_jacobian_host.resize(n_cells * padding_length * dim * dim);

      constraint_mask_host.resize(n_cells);
    }

    template <int dim, typename Number>
    template <typename CellFilter>
    void
    ReinitHelper<dim, Number>::get_cell_data(const CellFilter&  cell,
                                             const unsigned int cell_id)
    {
      cell->get_dof_indices(local_dof_indices);

      for(unsigned int i = 0; i < dofs_per_cell; ++i)
        lexicographic_dof_indices[i] = local_dof_indices[lexicographic_inv[i]];

      memcpy(&local_to_global_host[cell_id * padding_length],
             lexicographic_dof_indices.data(),
             dofs_per_cell * sizeof(unsigned int));

      fe_values.reinit(cell);

      // Quadrature points
      if(update_flags & update_quadrature_points)
        {
          const std::vector<Point<dim>>& q_points
            = fe_values.get_quadrature_points();
          memcpy(&q_points_host[cell_id * padding_length],
                 q_points.data(),
                 q_points_per_cell * sizeof(Point<dim>));
        }

      if(update_flags & update_JxW_values)
        {
          std::vector<double> JxW_values_double = fe_values.get_JxW_values();
          const unsigned int  offset            = cell_id * padding_length;
          for(unsigned int i = 0; i < q_points_per_cell; ++i)
            JxW_host[i + offset] = static_cast<Number>(JxW_values_double[i]);
        }

      if(update_flags & update_gradients)
        {
          const std::vector<DerivativeForm<1, dim, dim>>& inv_jacobians
            = fe_values.get_inverse_jacobians();
          memcpy(&inv_jacobian_host[cell_id * padding_length * dim * dim],
                 inv_jacobians.data(),
                 q_points_per_cell * sizeof(DerivativeForm<1, dim, dim>));
        }
    }

    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::alloc_and_copy_arrays(const unsigned int color)
    {
      const unsigned int n_cells = data->n_cells[color];

      // Local-to-global mapping
      if(data->parallelization_scheme
         == MatrixFree<dim, Number>::parallel_over_elem)
        internal::transpose_in_place(
          local_to_global_host, n_cells, padding_length);

      alloc_and_copy(&data->local_to_global[color],
                     local_to_global_host,
                     n_cells * padding_length);

      // Quadrature points
      if(update_flags & update_quadrature_points)
        {
          if(data->parallelization_scheme
             == MatrixFree<dim, Number>::parallel_over_elem)
            internal::transpose_in_place(
              q_points_host, n_cells, padding_length);

          alloc_and_copy(
            &data->q_points[color], q_points_host, n_cells * padding_length);
        }

      // Jacobian determinants/quadrature weights
      if(update_flags & update_JxW_values)
        {
          if(data->parallelization_scheme
             == MatrixFree<dim, Number>::parallel_over_elem)
            internal::transpose_in_place(JxW_host, n_cells, padding_length);

          alloc_and_copy(&data->JxW[color], JxW_host, n_cells * padding_length);
        }

      // Inverse jacobians
      if(update_flags & update_gradients)
        {
          // Reorder so that all J_11 elements are together, all J_12 elements are
          // together, etc., i.e., reorder indices from
          // cell_id*q_points_per_cell*dim*dim + q*dim*dim +i to
          // i*q_points_per_cell*n_cells + cell_id*q_points_per_cell+q
          internal::transpose_in_place(
            inv_jacobian_host, padding_length * n_cells, dim * dim);

          // Transpose second time means we get the following index order:
          // q*n_cells*dim*dim + i*n_cells + cell_id which is good for an
          // element-level parallelization
          if(data->parallelization_scheme
             == MatrixFree<dim, Number>::parallel_over_elem)
            internal::transpose_in_place(
              inv_jacobian_host, n_cells * dim * dim, padding_length);

          alloc_and_copy(&data->inv_jacobian[color],
                         inv_jacobian_host,
                         n_cells * dim * dim * padding_length);
        }

      alloc_and_copy(
        &data->constraint_mask[color], constraint_mask_host, n_cells);
    }

    template <int dim>
    std::vector<types::global_dof_index>
    get_conflict_indices(
      const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>&
                              cell,
      const ConstraintMatrix& constraints)
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
      const dealii::types::global_dof_index* constrained_dofs,
      const unsigned int                     n_constrained_dofs,
      const Number*                          src,
      Number*                                dst)
    {
      const unsigned int dof
        = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
      if(dof < n_constrained_dofs)
        dst[constrained_dofs[dof]] = src[constrained_dofs[dof]];
    }

    template <typename Number>
    __global__ void
    set_constrained_dofs(
      const dealii::types::global_dof_index* constrained_dofs,
      const unsigned int                     n_constrained_dofs,
      Number                                 val,
      Number*                                dst)
    {
      const unsigned int dof
        = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
      if(dof < n_constrained_dofs)
        dst[constrained_dofs[dof]] = val;
    }

    template <int dim, typename Number, typename functor>
    __global__ void
    apply_kernel_shmem(const functor&                               func,
                       const typename MatrixFree<dim, Number>::Data gpu_data,
                       const Number*                                src,
                       Number*                                      dst)
    {
      const unsigned int cells_per_block
        = cells_per_block_shmem(dim, functor::n_dofs_1d - 1);

      // TODO make use of dynamically allocated shared memory
      __shared__ Number values[cells_per_block * functor::n_local_dofs];
      __shared__ Number gradients[dim][cells_per_block * functor::n_q_points];

      const unsigned int local_cell = threadIdx.x / functor::n_dofs_1d;
      const unsigned int cell
        = local_cell + cells_per_block * (blockIdx.x + gridDim.x * blockIdx.y);

      Number* gq[dim];
      for(int d = 0; d < dim; ++d)
        gq[d] = &gradients[d][local_cell * functor::n_q_points];

      SharedData<dim, Number> shared_data(
        &values[local_cell * functor::n_local_dofs], gq);

      if(cell < gpu_data.n_cells)
        func(cell, &gpu_data, &shared_data, src, dst);
    }
  } // namespace internal

  template <int dim, typename Number>
  MatrixFree<dim, Number>::MatrixFree()
    : constrained_dofs(nullptr), padding_length(0)
  {}

  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::reinit(const Mapping<dim>&     mapping,
                                  const DoFHandler<dim>&  dof_handler,
                                  const ConstraintMatrix& constraints,
                                  const Quadrature<1>&    quad,
                                  const AdditionalData    additional_data)
  {
    if(typeid(Number) == typeid(double))
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    const UpdateFlags& update_flags = additional_data.mapping_update_flags;

    if(additional_data.parallelization_scheme != parallel_over_elem
       && additional_data.parallelization_scheme != parallel_in_elem)
      AssertThrow(false, ExcMessage("Invalid parallelization scheme."));

    this->parallelization_scheme = additional_data.parallelization_scheme;

    //TODO: only free if we actually need arrays of different length
    free();

    const FiniteElement<dim>& fe = dof_handler.get_fe();

    fe_degree = fe.degree;
    //TODO this should be a templated parameter
    const unsigned int n_dofs_1d     = fe_degree + 1;
    const unsigned int n_q_points_1d = quad.size();

    Assert(n_dofs_1d == n_q_points_1d,
           ExcMessage("n_q_points_1d must be equal to fe_degree+1."));

    // Set padding length to the closest power of two larger than or equal to the
    // number of threads.
    padding_length = 1 << static_cast<unsigned int>(
                       std::ceil(dim * std::log2(fe_degree + 1.)));

    dofs_per_cell     = fe.dofs_per_cell;
    q_points_per_cell = std::pow(n_q_points_1d, dim);

    const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info(
      quad, fe);

    unsigned int size_shape_values = n_dofs_1d * n_q_points_1d * sizeof(Number);

    cudaError_t cuda_error = cudaMemcpyToSymbol(internal::global_shape_values,
                                                &shape_info.shape_values[0],
                                                size_shape_values,
                                                0,
                                                cudaMemcpyHostToDevice);
    AssertCuda(cuda_error);

    if(update_flags & update_gradients)
      {
        cuda_error = cudaMemcpyToSymbol(internal::global_shape_gradients,
                                        &shape_info.shape_gradients[0],
                                        size_shape_values,
                                        0,
                                        cudaMemcpyHostToDevice);
        AssertCuda(cuda_error);
      }

    // Setup the number of cells per CUDA thread block
    cells_per_block = cells_per_block_shmem(dim, fe_degree);

    internal::ReinitHelper<dim, Number> helper(
      this, mapping, fe, quad, shape_info, update_flags);

    // Create a graph coloring
    typedef FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
               CellFilter;
    CellFilter begin(IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active());
    CellFilter end(IteratorFilters::LocallyOwnedCell(), dof_handler.end());
    typedef std::function<std::vector<types::global_dof_index>(
      CellFilter const&)>
                    fun_type;
    const fun_type& fun
      = static_cast<fun_type>(std::bind(&internal::get_conflict_indices<dim>,
                                        std::placeholders::_1,
                                        constraints));

    std::vector<std::vector<CellFilter>> graph
      = GraphColoring::make_graph_coloring(begin, end, fun);
    n_colors = graph.size();

    helper.setup_color_arrays(n_colors);
    for(unsigned int i = 0; i < n_colors; ++i)
      {
        n_cells[i] = graph[i].size();
        helper.setup_cell_arrays(i);
        typename std::vector<CellFilter>::iterator cell     = graph[i].begin(),
                                                   end_cell = graph[i].end();
        for(unsigned int cell_id = 0; cell != end_cell; ++cell, ++cell_id)
          helper.get_cell_data(*cell, cell_id);

        helper.alloc_and_copy_arrays(i);
      }

    // Setup row starts
    row_start[0] = 0;
    for(unsigned int i = 0; i < n_colors - 1; ++i)
      row_start[i + 1] = row_start[i] + n_cells[i] * get_padding_length();

    // Constrained indices
    n_constrained_dofs = constraints.n_constraints();

    if(n_constrained_dofs != 0)
      {
        const unsigned int constraint_n_blocks
          = std::ceil(static_cast<double>(n_constrained_dofs)
                      / static_cast<double>(block_size));
        const unsigned int constraint_x_n_blocks
          = std::round(std::sqrt(constraint_n_blocks));
        const unsigned int constraint_y_n_blocks
          = std::ceil(static_cast<double>(constraint_n_blocks)
                      / static_cast<double>(constraint_x_n_blocks));

        constraint_grid_dim
          = dim3(constraint_x_n_blocks, constraint_y_n_blocks);
        constraint_block_dim = dim3(block_size);

        std::vector<dealii::types::global_dof_index> constrained_dofs_host(
          n_constrained_dofs);

        unsigned int       i_constraint = 0;
        const unsigned int n_dofs       = dof_handler.n_dofs();
        for(unsigned int i = 0; i < n_dofs; ++i)
          {
            if(constraints.is_constrained(i))
              {
                constrained_dofs_host[i_constraint] = i;
                ++i_constraint;
              }
          }

        cuda_error = cudaMalloc(&constrained_dofs,
                                n_constrained_dofs
                                  * sizeof(dealii::types::global_dof_index));
        AssertCuda(cuda_error);

        cuda_error = cudaMemcpy(constrained_dofs,
                                constrained_dofs_host.data(),
                                n_constrained_dofs
                                  * sizeof(dealii::types::global_dof_index),
                                cudaMemcpyHostToDevice);
        AssertCuda(cuda_error);
      }
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

    return data_copy;
  }

  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::free()
  {
    for(unsigned int i = 0; i < q_points.size(); ++i)
      {
        if(q_points[i] != nullptr)
          {
            cudaError_t cuda_error = cudaFree(q_points[i]);
            AssertCuda(cuda_error);
            q_points[i] = nullptr;
          }
      }

    for(unsigned int i = 0; i < local_to_global.size(); ++i)
      {
        if(local_to_global[i] != nullptr)
          {
            cudaError_t cuda_error = cudaFree(local_to_global[i]);
            AssertCuda(cuda_error);
            local_to_global[i] = nullptr;
          }
      }

    for(unsigned int i = 0; i < inv_jacobian.size(); ++i)
      {
        if(inv_jacobian[i] != nullptr)
          {
            cudaError_t cuda_error = cudaFree(inv_jacobian[i]);
            AssertCuda(cuda_error);
            inv_jacobian[i] = nullptr;
          }
      }

    for(unsigned int i = 0; i < JxW.size(); ++i)
      {
        if(JxW[i] != nullptr)
          {
            cudaError_t cuda_error = cudaFree(JxW[i]);
            AssertCuda(cuda_error);
            JxW[i] = nullptr;
          }
      }

    for(unsigned int i = 0; i < constraint_mask.size(); ++i)
      {
        if(constraint_mask[i] != nullptr)
          {
            cudaError_t cuda_error = cudaFree(constraint_mask[i]);
            AssertCuda(cuda_error);
            constraint_mask[i] = nullptr;
          }
      }

    q_points.clear();
    local_to_global.clear();
    inv_jacobian.clear();
    JxW.clear();
    constraint_mask.clear();

    if(constrained_dofs != nullptr)
      {
        cudaError_t cuda_error = cudaFree(constrained_dofs);
        AssertCuda(cuda_error);
        constrained_dofs = nullptr;
      }
  }

  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::copy_constrained_values(
    const CUDAVector<Number>& src,
    CUDAVector<Number>&       dst) const
  {
    internal::copy_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(constrained_dofs,
                                                      n_constrained_dofs,
                                                      src.get_values(),
                                                      dst.get_values());
  }

  template <int dim, typename Number>
  void
  MatrixFree<dim, Number>::set_constrained_values(Number              val,
                                                  CUDAVector<Number>& dst) const
  {
    internal::set_constrained_dofs<Number>
      <<<constraint_grid_dim, constraint_block_dim>>>(
        constrained_dofs, n_constrained_dofs, val, dst.get_values());
  }

  template <int dim, typename Number>
  unsigned int
  MatrixFree<dim, Number>::get_padding_length() const
  {
    return padding_length;
  }

  template <int dim, typename Number>
  template <typename functor>
  void
  MatrixFree<dim, Number>::cell_loop(const functor&            func,
                                     const CUDAVector<Number>& src,
                                     CUDAVector<Number>&       dst) const
  {
    for(unsigned int i = 0; i < n_colors; ++i)
      internal::apply_kernel_shmem<dim, Number, functor>
        <<<grid_dim[i], block_dim[i]>>>(
          func, get_data(i), src.get_values(), dst.get_values());
  }

  template <int dim, typename Number>
  std::size_t
  MatrixFree<dim, Number>::memory_consumption() const
  {
    // First compute the size of n_cells, row_starts, kernel launch parameters,
    // and constrained_dofs
    std::size_t bytes = n_cells.size() * sizeof(unsigned int) * 2
                        + 2 * n_colors * sizeof(dim3)
                        + n_constrained_dofs * sizeof(unsigned int);

    // For each color, add local_to_global, inv_jacobian, JxW, and q_points.
    for(unsigned int i = 0; i < n_colors; ++i)
      {
        bytes += n_cells[i] * padding_length * sizeof(unsigned int)
                 + n_cells[i] * padding_length * dim * dim * sizeof(Number)
                 + n_cells[i] * padding_length * sizeof(Number)
                 + n_cells[i] * padding_length * sizeof(point_type)
                 + n_cells[i] * sizeof(unsigned int);
      }

    return bytes;
  }
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
