// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_portable_matrix_free_templates_h
#define dealii_portable_matrix_free_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/graph_coloring.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/block_vector_base.h>

#include <deal.II/matrix_free/portable_hanging_nodes_internal.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>

#include <Kokkos_Core.hpp>

#include <cmath>
#include <functional>
#include <string>


DEAL_II_NAMESPACE_OPEN

namespace Portable
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
        MatrixFree<dim, Number>       *data,
        const Mapping<dim>            &mapping,
        const FiniteElement<dim, dim> &fe,
        const Quadrature<1>           &quad,
        const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
                              &shape_info,
        const DoFHandler<dim> &dof_handler,
        const UpdateFlags     &update_flags);

      void
      resize(const unsigned int n_colors);

      template <typename CellFilter>
      void
      fill_data(
        const unsigned int                                        color,
        const std::vector<CellFilter>                            &graph,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

    private:
      MatrixFree<dim, Number> *data;
      // Local buffer
      std::vector<types::global_dof_index> local_dof_indices;
      FEValues<dim>                        fe_values;
      // Convert the default dof numbering to a lexicographic one
      const std::vector<unsigned int>     &lexicographic_inv;
      std::vector<types::global_dof_index> lexicographic_dof_indices;
      const unsigned int                   fe_degree;
      const unsigned int                   n_components;
      const unsigned int                   scalar_dofs_per_cell;
      const unsigned int                   dofs_per_cell;
      const unsigned int                   q_points_per_cell;
      const UpdateFlags                   &update_flags;
      const unsigned int                   padding_length;
      dealii::internal::MatrixFreeFunctions::HangingNodes<dim> hanging_nodes;
    };



    template <int dim, typename Number>
    ReinitHelper<dim, Number>::ReinitHelper(
      MatrixFree<dim, Number>  *data,
      const Mapping<dim>       &mapping,
      const FiniteElement<dim> &fe,
      const Quadrature<1>      &quad,
      const ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
                            &shape_info,
      const DoFHandler<dim> &dof_handler,
      const UpdateFlags     &update_flags)
      : data(data)
      , fe_values(mapping,
                  fe,
                  Quadrature<dim>(quad),
                  update_inverse_jacobians | update_quadrature_points |
                    update_values | update_gradients | update_JxW_values)
      , lexicographic_inv(shape_info.lexicographic_numbering)
      , fe_degree(data->fe_degree)
      , n_components(data->n_components)
      , scalar_dofs_per_cell(data->scalar_dofs_per_cell)
      , dofs_per_cell(data->dofs_per_cell)
      , q_points_per_cell(data->q_points_per_cell)
      , update_flags(update_flags)
      , padding_length(data->get_padding_length())
      , hanging_nodes(dof_handler.get_triangulation())
    {
      local_dof_indices.resize(data->dofs_per_cell);
      lexicographic_dof_indices.resize(dofs_per_cell);
      fe_values.always_allow_check_for_cell_similarity(true);
    }



    template <int dim, typename Number>
    void
    ReinitHelper<dim, Number>::resize(const unsigned int n_colors)
    {
      // We need at least three colors when we are using device-aware MPI and
      // overlapping the communication
      data->n_cells.resize(std::max(n_colors, 3U), 0);
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
    template <typename CellFilter>
    void
    ReinitHelper<dim, Number>::fill_data(
      const unsigned int                                        color,
      const std::vector<CellFilter>                            &graph,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
    {
      const unsigned int n_cells = data->n_cells[color];

      // Create the Views
      data->local_to_global[color] =
        Kokkos::View<types::global_dof_index **,
                     MemorySpace::Default::kokkos_space>(
          Kokkos::view_alloc("local_to_global_" + std::to_string(color),
                             Kokkos::WithoutInitializing),
          dofs_per_cell,
          n_cells);

      if (update_flags & update_quadrature_points)
        data->q_points[color] =
          Kokkos::View<Point<dim, Number> **,
                       MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("q_points_" + std::to_string(color),
                               Kokkos::WithoutInitializing),
            q_points_per_cell,
            n_cells);

      if (update_flags & update_JxW_values)
        data->JxW[color] =
          Kokkos::View<Number **, MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("JxW_" + std::to_string(color),
                               Kokkos::WithoutInitializing),
            q_points_per_cell,
            n_cells);

      if (update_flags & update_gradients)
        data->inv_jacobian[color] =
          Kokkos::View<Number **[dim][dim], MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("inv_jacobian_" + std::to_string(color),
                               Kokkos::WithoutInitializing),
            q_points_per_cell,
            n_cells);

      // Initialize to zero, i.e., unconstrained cell
      data->constraint_mask[color] =
        Kokkos::View<dealii::internal::MatrixFreeFunctions::ConstraintKinds *,
                     MemorySpace::Default::kokkos_space>(
          "constraint_mask_" + std::to_string(color), n_cells * n_components);

      // Create the host mirrow Views and fill them
      auto constraint_mask_host =
        Kokkos::create_mirror_view(data->constraint_mask[color]);

      typename std::remove_reference_t<
        decltype(data->q_points[color])>::HostMirror q_points_host;
      typename std::remove_reference_t<decltype(data->JxW[color])>::HostMirror
        JxW_host;
      typename std::remove_reference_t<
        decltype(data->inv_jacobian[color])>::HostMirror inv_jacobian_host;
#if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
      auto local_to_global_host =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                   data->local_to_global[color]);
      if (update_flags & update_quadrature_points)
        q_points_host = Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                                   data->q_points[color]);
      if (update_flags & update_JxW_values)
        JxW_host = Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                              data->JxW[color]);
      if (update_flags & update_gradients)
        inv_jacobian_host =
          Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                     data->inv_jacobian[color]);
#else
      auto local_to_global_host =
        Kokkos::create_mirror_view(data->local_to_global[color]);
      if (update_flags & update_quadrature_points)
        q_points_host = Kokkos::create_mirror_view(data->q_points[color]);
      if (update_flags & update_JxW_values)
        JxW_host = Kokkos::create_mirror_view(data->JxW[color]);
      if (update_flags & update_gradients)
        inv_jacobian_host =
          Kokkos::create_mirror_view(data->inv_jacobian[color]);
#endif

      auto cell = graph.cbegin(), end_cell = graph.cend();
      for (unsigned int cell_id = 0; cell != end_cell; ++cell, ++cell_id)
        {
          (*cell)->get_dof_indices(local_dof_indices);
          // When using MPI, we need to transform the local_dof_indices, which
          // contain global numbers of dof indices in the MPI universe, to get
          // local (to the current MPI process) dof indices.
          if (partitioner)
            for (auto &index : local_dof_indices)
              index = partitioner->global_to_local(index);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            lexicographic_dof_indices[i] =
              local_dof_indices[lexicographic_inv[i]];

          const ArrayView<
            dealii::internal::MatrixFreeFunctions::ConstraintKinds>
            cell_id_view(&constraint_mask_host[cell_id * n_components],
                         n_components);

          hanging_nodes.setup_constraints(*cell,
                                          partitioner,
                                          {lexicographic_inv},
                                          lexicographic_dof_indices,
                                          cell_id_view);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            local_to_global_host(i, cell_id) = lexicographic_dof_indices[i];

          fe_values.reinit(*cell);

          // Quadrature points
          if (update_flags & update_quadrature_points)
            {
              for (unsigned int i = 0; i < q_points_per_cell; ++i)
                q_points_host(i, cell_id) = fe_values.quadrature_point(i);
            }

          if (update_flags & update_JxW_values)
            {
              for (unsigned int i = 0; i < q_points_per_cell; ++i)
                JxW_host(i, cell_id) = fe_values.JxW(i);
            }

          if (update_flags & update_gradients)
            {
              for (unsigned int i = 0; i < q_points_per_cell; ++i)
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    inv_jacobian_host(i, cell_id, d, e) =
                      fe_values.inverse_jacobian(i)[d][e];
            }
        }

      // Copy the data to the device
      Kokkos::deep_copy(data->constraint_mask[color], constraint_mask_host);
      Kokkos::deep_copy(data->local_to_global[color], local_to_global_host);
      if (update_flags & update_quadrature_points)
        Kokkos::deep_copy(data->q_points[color], q_points_host);
      if (update_flags & update_JxW_values)
        Kokkos::deep_copy(data->JxW[color], JxW_host);
      if (update_flags & update_gradients)
        Kokkos::deep_copy(data->inv_jacobian[color], inv_jacobian_host);
    }



    template <int dim, typename number>
    std::vector<types::global_dof_index>
    get_conflict_indices(
      const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
                                      &cell,
      const AffineConstraints<number> &constraints)
    {
      std::vector<types::global_dof_index> local_dof_indices(
        cell->get_fe().n_dofs_per_cell());
      cell->get_dof_indices(local_dof_indices);
      constraints.resolve_indices(local_dof_indices);

      return local_dof_indices;
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



    template <int dim, typename Number, typename Functor, bool IsBlock>
    struct ApplyKernel
    {
      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;
      using SharedViewValues =
        Kokkos::View<Number **,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      using SharedViewGradients =
        Kokkos::View<Number ***,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      using SharedViewScratchPad =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      ApplyKernel(
        Functor                                                 func,
        const typename MatrixFree<dim, Number>::PrecomputedData gpu_data,
        const LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
                                                                         &src,
        LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> &dst)
        : func(func)
        , gpu_data(gpu_data)
        , src(DeviceVector<Number>(src.get_values(), src.locally_owned_size()))
        , dst(DeviceVector<Number>(dst.get_values(), dst.locally_owned_size()))
      {}

      ApplyKernel(
        Functor                                                 func,
        const typename MatrixFree<dim, Number>::PrecomputedData gpu_data,
        const LinearAlgebra::distributed::BlockVector<Number,
                                                      MemorySpace::Default>
          &src,
        LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>
          &dst)
        : func(func)
        , gpu_data(gpu_data)
        , src(src)
        , dst(dst)
      {}

      Functor                                                 func;
      const typename MatrixFree<dim, Number>::PrecomputedData gpu_data;
      const DeviceBlockVector<Number>                         src;
      DeviceBlockVector<Number>                               dst;


      // Provide the shared memory capacity. This function takes the team_size
      // as an argument, which allows team_size dependent allocations.
      std::size_t
      team_shmem_size(int /*team_size*/) const
      {
        return SharedViewValues::shmem_size(Functor::n_q_points,
                                            gpu_data.n_components) +
               SharedViewGradients::shmem_size(Functor::n_q_points,
                                               dim,
                                               gpu_data.n_components) +
               SharedViewScratchPad::shmem_size(gpu_data.scratch_pad_size);
      }


      DEAL_II_HOST_DEVICE
      void
      operator()(const TeamHandle &team_member) const
      {
        // Get the scratch memory
        SharedViewValues     values(team_member.team_shmem(),
                                Functor::n_q_points,
                                gpu_data.n_components);
        SharedViewGradients  gradients(team_member.team_shmem(),
                                      Functor::n_q_points,
                                      dim,
                                      gpu_data.n_components);
        SharedViewScratchPad scratch_pad(team_member.team_shmem(),
                                         gpu_data.scratch_pad_size);

        SharedData<dim, Number> shared_data(values, gradients, scratch_pad);

        const int cell_index = team_member.league_rank();

        typename MatrixFree<dim, Number>::Data data{team_member,
                                                    /* n_dofhandler */ 1,
                                                    cell_index,
                                                    &gpu_data,
                                                    &shared_data};

        if constexpr (IsBlock)
          {
            DeviceBlockVector<Number> nonconstdst = dst;
            func(&data, src, nonconstdst);
          }
        else
          {
            DeviceVector<Number> nonconstdst = dst.block(0);
            func(&data, src.block(0), nonconstdst);
          }
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
  MatrixFree<dim, Number>::reinit(const Mapping<dim>              &mapping,
                                  const DoFHandler<dim>           &dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1>             &quad,
                                  const IteratorFiltersType &iterator_filter,
                                  const AdditionalData      &additional_data)
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
                        parallel_triangulation->get_mpi_communicator()),
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
  MatrixFree<dim, Number>::reinit(const Mapping<dim>              &mapping,
                                  const DoFHandler<dim>           &dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1>             &quad,
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
  MatrixFree<dim, Number>::reinit(const DoFHandler<dim>           &dof_handler,
                                  const AffineConstraints<Number> &constraints,
                                  const Quadrature<1>             &quad,
                                  const AdditionalData &additional_data)
  {
    reinit(StaticMappingQ1<dim>::mapping,
           dof_handler,
           constraints,
           quad,
           additional_data);
  }



  template <int dim, typename Number>
  typename MatrixFree<dim, Number>::PrecomputedData
  MatrixFree<dim, Number>::get_data(unsigned int color) const
  {
    PrecomputedData data_copy;
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
    data_copy.n_components       = n_components;
    data_copy.padding_length     = padding_length;
    data_copy.row_start          = row_start[color];
    data_copy.use_coloring       = use_coloring;
    data_copy.element_type       = element_type;
    data_copy.scratch_pad_size   = scratch_pad_size;

    return data_copy;
  }



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::copy_constrained_values(const VectorType &src,
                                                   VectorType       &dst) const
  {
    static_assert(
      std::is_same_v<Number, typename VectorType::value_type>,
      "VectorType::value_type and Number should be of the same type.");
    Assert(src.size() == dst.size(),
           ExcMessage("src and dst vectors have different size."));
    // FIXME When using C++17, we can use KOKKOS_CLASS_LAMBDA and this
    // work-around can be removed.
    auto               constr_dofs = constrained_dofs;
    const unsigned int size = internal::VectorLocalSize<VectorType>::get(dst);
    const Number      *src_ptr = src.get_values();
    Number            *dst_ptr = dst.get_values();
    Kokkos::parallel_for(
      "dealii::copy_constrained_values",
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        0, n_constrained_dofs),
      KOKKOS_LAMBDA(int dof) {
        // When working with distributed vectors, the constrained dofs are
        // computed for ghosted vectors but we want to copy the values of the
        // constrained dofs of non-ghosted vectors.
        const auto constrained_dof = constr_dofs[dof];
        if (constrained_dof < size)
          dst_ptr[constrained_dof] = src_ptr[constrained_dof];
      });
  }



  template <int dim, typename Number>
  template <typename VectorType>
  void
  MatrixFree<dim, Number>::set_constrained_values(Number      val,
                                                  VectorType &dst) const
  {
    static_assert(
      std::is_same_v<Number, typename VectorType::value_type>,
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
      "dealii::set_constrained_values",
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        0, n_constrained_dofs),
      KOKKOS_LAMBDA(int dof) {
        if (constr_dofs[dof] < size)
          dst_ptr[constr_dofs[dof]] = val;
      });
  }



  template <int dim, typename Number>
  template <typename MemorySpaceType>
  void
  MatrixFree<dim, Number>::initialize_dof_vector(
    LinearAlgebra::distributed::Vector<Number, MemorySpaceType> &vec) const
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
  MatrixFree<dim, Number>::cell_loop(const Functor    &func,
                                     const VectorType &src,
                                     VectorType       &dst) const
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
    const unsigned int n_q_points = Functor::n_q_points;

    for (unsigned int i = 0; i < n_colors; ++i)
      if (n_cells[i] > 0)
        {
          auto color_data = get_data(i);

          MemorySpace::Default::kokkos_space::execution_space exec;
          Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>
            team_policy(exec, n_cells[i], Kokkos::AUTO);

          Kokkos::parallel_for(
            "dealii::MatrixFree::evaluate_coeff_cell_loop",
            team_policy,
            KOKKOS_LAMBDA(const Kokkos::TeamPolicy<
                          MemorySpace::Default::kokkos_space::execution_space>::
                            member_type &team_member) {
              Kokkos::parallel_for(Kokkos::TeamVectorRange(team_member,
                                                           n_q_points),
                                   [&](const int q_point) {
                                     const int cell_index =
                                       team_member.league_rank();

                                     Data data{team_member,
                                               /* n_dofhandler */ 1,
                                               cell_index,
                                               &color_data,
                                               /* shared_data */ nullptr};

                                     func(&data, cell_index, q_point);
                                   });
            });
        }
  }



  template <int dim, typename Number>
  std::size_t
  MatrixFree<dim, Number>::memory_consumption() const
  {
    // First compute the size of n_cells, row_starts, kernel launch parameters,
    // and constrained_dofs
    std::size_t bytes = n_cells.size() * sizeof(unsigned int) * 2 +
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



  namespace internal
  {
    /**
     * Helper function for determining the scratch pad size.
     */
    inline unsigned int
    compute_scratch_pad_size(
      const ::dealii::internal::MatrixFreeFunctions::ElementType element_type,
      const int                                                  dim,
      const int                                                  fe_degree,
      const int                                                  n_q_points_1d)
    {
      using ElementType = ::dealii::internal::MatrixFreeFunctions::ElementType;

      if (fe_degree >= 0 && element_type <= ElementType::tensor_symmetric)
        {
          // evaluate/integrate with FEEvaluationImplCollocation or
          // FEEvaluationImplTransformToCollocation
          return Utilities::pow(n_q_points_1d, dim);
        }
      else if (fe_degree >= 0 &&
               element_type <= ElementType::tensor_symmetric_no_collocation)
        {
          // evaluate/integrate with FEEvaluationImpl
          if (dim == 1)
            return n_q_points_1d;
          else if (dim == 2)
            return (fe_degree + 1) * n_q_points_1d;
          else if (dim == 3)
            return (fe_degree + 1) * n_q_points_1d *
                   (fe_degree + 1 + n_q_points_1d);
          else
            AssertThrow(false, ExcNotImplemented());
        }
      else
        AssertThrow(false, ExcNotImplemented());

      return numbers::invalid_unsigned_int;
    }
  } // namespace internal



  template <int dim, typename Number>
  template <typename IteratorFiltersType>
  void
  MatrixFree<dim, Number>::internal_reinit(
    const Mapping<dim>                    &mapping,
    const DoFHandler<dim>                 &dof_handler_,
    const AffineConstraints<Number>       &constraints,
    const Quadrature<1>                   &quad,
    const IteratorFiltersType             &iterator_filter,
    const std::shared_ptr<const MPI_Comm> &comm,
    const AdditionalData                   additional_data)
  {
    dof_handler = &dof_handler_;

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

    // TODO remove the limitation in the future
    AssertThrow(n_dofs_1d <= n_q_points_1d,
                ExcMessage("n_q_points_1d must be greater than or equal to "
                           "fe_degree + 1."));

    // Set padding length to the closest power of two larger than or equal to
    // the number of threads.
    padding_length = 1 << static_cast<unsigned int>(std::ceil(
                       dim * std::log2(static_cast<float>(n_q_points_1d))));

    dofs_per_cell        = fe.n_dofs_per_cell();
    n_components         = fe.n_components();
    scalar_dofs_per_cell = dofs_per_cell / n_components;
    q_points_per_cell    = Utilities::fixed_power<dim>(n_q_points_1d);

    ::dealii::internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info(quad,
                                                                          fe);

    this->element_type     = shape_info.element_type;
    this->scratch_pad_size = internal::compute_scratch_pad_size(
      shape_info.element_type, dim, fe_degree, n_q_points_1d);

    unsigned int size_shape_values = n_dofs_1d * n_q_points_1d;

    shape_values = Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
      Kokkos::view_alloc("shape_values", Kokkos::WithoutInitializing),
      size_shape_values);
    Kokkos::deep_copy(shape_values,
                      Kokkos::View<Number *, Kokkos::HostSpace>(
                        shape_info.data.front().shape_values.data(),
                        size_shape_values));

    if (update_flags & update_gradients)
      {
        shape_gradients =
          Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("shape_gradients", Kokkos::WithoutInitializing),
            size_shape_values);
        Kokkos::deep_copy(shape_gradients,
                          Kokkos::View<Number *, Kokkos::HostSpace>(
                            shape_info.data.front().shape_gradients.data(),
                            size_shape_values));


        co_shape_gradients =
          Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
            Kokkos::view_alloc("co_shape_gradients",
                               Kokkos::WithoutInitializing),
            n_q_points_1d * n_q_points_1d);
        Kokkos::deep_copy(
          co_shape_gradients,
          Kokkos::View<Number *, Kokkos::HostSpace>(
            shape_info.data.front().shape_gradients_collocation.data(),
            n_q_points_1d * n_q_points_1d));
      }

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

    helper.resize(n_colors);

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
        helper.fill_data(i, graph[i], partitioner);
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
  MatrixFree<dim, Number>::serial_cell_loop(const Functor    &func,
                                            const VectorType &src,
                                            VectorType       &dst) const
  {
    // Execute the loop on the cells
    for (unsigned int color = 0; color < n_colors; ++color)
      if (n_cells[color] > 0)
        {
          MemorySpace::Default::kokkos_space::execution_space exec;
          Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>
            team_policy(exec, n_cells[color], Kokkos::AUTO);


          internal::
            ApplyKernel<dim, Number, Functor, IsBlockVector<VectorType>::value>
              apply_kernel(func, get_data(color), src, dst);

          Kokkos::parallel_for("dealii::MatrixFree::serial_cell_loop",
                               team_policy,
                               apply_kernel);
        }
    Kokkos::fence();
  }



  template <int dim, typename Number>
  template <typename Functor>
  void
  MatrixFree<dim, Number>::distributed_cell_loop(
    const Functor &,
    const LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>
      &,
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> &)
    const
  {
    Assert(false, ExcNotImplemented());
  }



  template <int dim, typename Number>
  template <typename Functor>
  void
  MatrixFree<dim, Number>::distributed_cell_loop(
    const Functor &func,
    const LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> &src,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> &dst) const
  {
    MemorySpace::Default::kokkos_space::execution_space exec;

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
                Kokkos::TeamPolicy<
                  MemorySpace::Default::kokkos_space::execution_space>
                  team_policy(exec, n_cells[0], Kokkos::AUTO);

                internal::ApplyKernel<dim, Number, Functor, false> apply_kernel(
                  func, get_data(0), src, dst);

                Kokkos::parallel_for(
                  "dealii::MatrixFree::distributed_cell_loop_0",
                  team_policy,
                  apply_kernel);
              }
            src.update_ghost_values_finish();

            // In serial this color does not exist because there are no ghost
            // cells
            if (n_cells[1] > 0)
              {
                Kokkos::TeamPolicy<
                  MemorySpace::Default::kokkos_space::execution_space>
                  team_policy(exec, n_cells[1], Kokkos::AUTO);

                internal::ApplyKernel<dim, Number, Functor, false> apply_kernel(
                  func, get_data(1), src, dst);

                Kokkos::parallel_for(
                  "dealii::MatrixFree::distributed_cell_loop_1",
                  team_policy,
                  apply_kernel);

                // We need a synchronization point because we don't want
                // device-aware MPI to start the MPI communication until the
                // kernel is done.
                Kokkos::fence();
              }

            dst.compress_start(0, VectorOperation::add);
            // When the mesh is coarse it is possible that some processors do
            // not own any cells
            if (n_cells[2] > 0)
              {
                Kokkos::TeamPolicy<
                  MemorySpace::Default::kokkos_space::execution_space>
                  team_policy(exec, n_cells[2], Kokkos::AUTO);

                internal::ApplyKernel<dim, Number, Functor, false> apply_kernel(
                  func, get_data(2), src, dst);

                Kokkos::parallel_for(
                  "dealii::MatrixFree::distributed_cell_loop_2",
                  team_policy,
                  apply_kernel);
              }
            dst.compress_finish(VectorOperation::add);
          }
        else
          {
            src.update_ghost_values();
            std::vector<
              Kokkos::View<Number *, MemorySpace::Default::kokkos_space>>
              values_colors(n_colors);
            std::vector<
              Kokkos::View<Number *[dim], MemorySpace::Default::kokkos_space>>
              gradients_colors(n_colors);

            // Execute the loop on the cells
            for (unsigned int i = 0; i < n_colors; ++i)
              if (n_cells[i] > 0)
                {
                  Kokkos::TeamPolicy<
                    MemorySpace::Default::kokkos_space::execution_space>
                    team_policy(exec, n_cells[i], Kokkos::AUTO);

                  internal::ApplyKernel<dim, Number, Functor, false>
                    apply_kernel(func, get_data(i), src, dst);

                  Kokkos::parallel_for(
                    "dealii::MatrixFree::distributed_cell_loop_" +
                      std::to_string(i),
                    team_policy,
                    apply_kernel);
                }
            dst.compress(VectorOperation::add);
          }
        src.zero_out_ghost_values();
      }
    else
      {
        // Create the ghosted source and the ghosted destination
        LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
          ghosted_src(partitioner);
        LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
          ghosted_dst(ghosted_src);
        ghosted_src = src;
        ghosted_dst = dst;
        ghosted_dst.zero_out_ghost_values();

        // Execute the loop on the cells
        for (unsigned int i = 0; i < n_colors; ++i)
          if (n_cells[i] > 0)
            {
              Kokkos::TeamPolicy<
                MemorySpace::Default::kokkos_space::execution_space>
                team_policy(exec, n_cells[i], Kokkos::AUTO);

              internal::ApplyKernel<dim, Number, Functor, false> apply_kernel(
                func, get_data(i), ghosted_src, ghosted_dst);

              Kokkos::parallel_for(
                "dealii::MatrixFree::distributed_cell_loop_" +
                  std::to_string(i),
                team_policy,
                apply_kernel);
            }

        // Add the ghosted values
        ghosted_dst.compress(VectorOperation::add);
        dst = ghosted_dst;
      }
  }

} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
