// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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


#ifndef dealii_mg_transfer_global_coarsening_templates_h
#define dealii_mg_transfer_global_coarsening_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/vector_access_internal.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Helper class to select the right templated implementation.
   *
   * @note This class is similar to internal::FEEvaluationFactory
   */
  class CellTransferFactory
  {
  public:
    static const unsigned int max_degree = 9;

    CellTransferFactory(const unsigned int degree_fine,
                        const unsigned int degree_coarse)
      : degree_fine(degree_fine)
      , degree_coarse(degree_coarse)
    {}

    template <typename Fu, unsigned int deg = 1>
    bool
    run(Fu &fu)
    {
      if ((degree_fine == (2 * deg)) && (degree_coarse == deg))
        fu.template run<2 * deg, deg>(); // h-MG (FE_Q)
      else if ((degree_fine == (2 * deg + 1)) && (degree_coarse == deg))
        fu.template run<2 * deg + 1, deg>(); // h-MG
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg / 2, 1u)))
        {
          constexpr unsigned int degree_coarse_used = std::max(deg / 2u, 1u);
          fu.template run<deg, degree_coarse_used>(); // p-MG: bisection
        }
      else if ((degree_fine == deg) && (degree_coarse == deg))
        fu.template run<deg, deg>(); // identity (nothing to do)
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg - 1, 1u)))
        {
          constexpr unsigned int degree_coarse_used = std::max(deg - 1u, 1u);
          fu.template run<deg, degree_coarse_used>(); // p-MG: --
        }
      else if ((degree_fine == deg) && (degree_coarse == 1))
        fu.template run<deg, 1>(); // p-MG: jump to 1
      else if (deg < max_degree)
        return run<Fu, std::min(deg + 1, max_degree)>(fu); // try next degree
      else
        {
          // no match -> slow path
          fu.template run<-1, -1>(degree_fine, degree_coarse);
          return false; // indicate that slow path has been taken
        }

      return true; // indicate that fast path has been taken
    }

  private:
    const unsigned int degree_fine;
    const unsigned int degree_coarse;
  };

  /**
   * Helper class containing the cell-wise prolongation operation.
   */
  template <int dim, typename Number>
  class CellProlongator
  {
  public:
    CellProlongator(const AlignedVector<Number> &prolongation_matrix,
                    const AlignedVector<Number> &prolongation_matrix_1d,
                    const Number *               evaluation_data_coarse,
                    Number *                     evaluation_data_fine)
      : prolongation_matrix(prolongation_matrix)
      , prolongation_matrix_1d(prolongation_matrix_1d)
      , evaluation_data_coarse(evaluation_data_coarse)
      , evaluation_data_fine(evaluation_data_fine)
    {}

    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int degree_fine_   = numbers::invalid_unsigned_int,
        const unsigned int degree_coarse_ = numbers::invalid_unsigned_int)
    {
      Assert(prolongation_matrix_1d.size() > 0, ExcNotImplemented());

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        dim,
        degree_coarse + 1,
        degree_fine + 1>::do_forward(1,
                                     prolongation_matrix_1d,
                                     evaluation_data_coarse,
                                     evaluation_data_fine,
                                     degree_coarse_ + 1,
                                     degree_fine_ + 1);
    }

    void
    run_full(const unsigned int n_dofs_fine, const unsigned int n_dofs_coarse)
    {
      AssertDimension(prolongation_matrix.size(), n_dofs_coarse * n_dofs_fine);

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        1,
        0,
        0>::do_forward(1,
                       prolongation_matrix,
                       evaluation_data_coarse,
                       evaluation_data_fine,
                       n_dofs_coarse,
                       n_dofs_fine);
    }

  private:
    const AlignedVector<Number> &prolongation_matrix;
    const AlignedVector<Number> &prolongation_matrix_1d;
    const Number *               evaluation_data_coarse;
    Number *                     evaluation_data_fine;
  };

  /**
   * Helper class containing the cell-wise restriction operation.
   */
  template <int dim, typename Number>
  class CellRestrictor
  {
  public:
    CellRestrictor(const AlignedVector<Number> &prolongation_matrix,
                   const AlignedVector<Number> &prolongation_matrix_1d,
                   Number *                     evaluation_data_fine,
                   Number *                     evaluation_data_coarse)
      : prolongation_matrix(prolongation_matrix)
      , prolongation_matrix_1d(prolongation_matrix_1d)
      , evaluation_data_fine(evaluation_data_fine)
      , evaluation_data_coarse(evaluation_data_coarse)
    {}

    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int degree_fine_   = numbers::invalid_unsigned_int,
        const unsigned int degree_coarse_ = numbers::invalid_unsigned_int)
    {
      Assert(prolongation_matrix_1d.size() > 0, ExcNotImplemented());

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        dim,
        degree_coarse == 0 ? -1 : (degree_coarse + 1),
        degree_fine == 0 ? -1 : (degree_fine + 1)>::
        do_backward(1,
                    prolongation_matrix_1d,
                    false,
                    evaluation_data_fine,
                    evaluation_data_coarse,
                    degree_coarse_ + 1,
                    degree_fine_ + 1);
    }

    void
    run_full(const unsigned int n_dofs_fine, const unsigned int n_dofs_coarse)
    {
      AssertDimension(prolongation_matrix.size(), n_dofs_coarse * n_dofs_fine);

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        1,
        0,
        0>::do_backward(1,
                        prolongation_matrix,
                        false,
                        evaluation_data_fine,
                        evaluation_data_coarse,
                        n_dofs_coarse,
                        n_dofs_fine);
    }

  private:
    const AlignedVector<Number> &prolongation_matrix;
    const AlignedVector<Number> &prolongation_matrix_1d;
    Number *                     evaluation_data_fine;
    Number *                     evaluation_data_coarse;
  };

  class CellProlongatorTest
  {
  public:
    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int = numbers::invalid_unsigned_int,
        const unsigned int = numbers::invalid_unsigned_int)
    {}
  };

} // namespace



namespace internal
{
  namespace
  {
    template <typename MeshType, typename OPType>
    DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
    void loop_over_active_or_level_cells(const MeshType &   tria,
                                         const unsigned int level,
                                         const OPType &     op)
    {
      if (level == numbers::invalid_unsigned_int)
        {
          for (const auto &cell : tria.active_cell_iterators())
            if (cell->is_locally_owned())
              op(cell);
        }
      else
        {
          for (const auto &cell : tria.cell_iterators_on_level(level))
            if (cell->is_locally_owned_on_level())
              op(cell);
        }
    }

    template <int dim>
    unsigned int
    compute_shift_within_children(const unsigned int child,
                                  const unsigned int fe_shift_1d,
                                  const unsigned int fe_degree)
    {
      // we put the degrees of freedom of all child cells in lexicographic
      // ordering
      unsigned int c_tensor_index[dim];
      unsigned int tmp = child;
      for (unsigned int d = 0; d < dim; ++d)
        {
          c_tensor_index[d] = tmp % 2;
          tmp /= 2;
        }
      const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
      unsigned int       factor          = 1;
      unsigned int       shift           = fe_shift_1d * c_tensor_index[0];
      for (unsigned int d = 1; d < dim; ++d)
        {
          factor *= n_child_dofs_1d;
          shift = shift + factor * fe_shift_1d * c_tensor_index[d];
        }
      return shift;
    }

    template <int dim>
    void
    get_child_offset(const unsigned int         child,
                     const unsigned int         fe_shift_1d,
                     const unsigned int         fe_degree,
                     std::vector<unsigned int> &local_dof_indices)
    {
      const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
      const unsigned int shift =
        compute_shift_within_children<dim>(child, fe_shift_1d, fe_degree);
      const unsigned int n_components =
        local_dof_indices.size() / Utilities::fixed_power<dim>(fe_degree + 1);
      const unsigned int n_scalar_cell_dofs =
        Utilities::fixed_power<dim>(n_child_dofs_1d);
      for (unsigned int c = 0, m = 0; c < n_components; ++c)
        for (unsigned int k = 0; k < (dim > 2 ? (fe_degree + 1) : 1); ++k)
          for (unsigned int j = 0; j < (dim > 1 ? (fe_degree + 1) : 1); ++j)
            for (unsigned int i = 0; i < (fe_degree + 1); ++i, ++m)
              local_dof_indices[m] = c * n_scalar_cell_dofs +
                                     k * n_child_dofs_1d * n_child_dofs_1d +
                                     j * n_child_dofs_1d + i + shift;
    }

    template <int dim>
    std::vector<std::vector<unsigned int>>
    get_child_offsets(const unsigned int n_dofs_per_cell_coarse,
                      const unsigned int fe_shift_1d,
                      const unsigned int fe_degree)
    {
      std::vector<std::vector<unsigned int>> cell_local_children_indices(
        GeometryInfo<dim>::max_children_per_cell,
        std::vector<unsigned int>(n_dofs_per_cell_coarse));
      for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
           c++)
        get_child_offset<dim>(c,
                              fe_shift_1d,
                              fe_degree,
                              cell_local_children_indices[c]);
      return cell_local_children_indices;
    }

    template <int dim>
    std::vector<std::vector<unsigned int>>
    get_child_offsets_general(const unsigned int n_dofs_per_cell_coarse)
    {
      std::vector<std::vector<unsigned int>> cell_local_children_indices(
        GeometryInfo<dim>::max_children_per_cell,
        std::vector<unsigned int>(n_dofs_per_cell_coarse));
      for (unsigned int c = 0, k = 0;
           c < GeometryInfo<dim>::max_children_per_cell;
           c++)
        for (unsigned int d = 0; d < n_dofs_per_cell_coarse; ++d, ++k)
          cell_local_children_indices[c][d] = k;
      return cell_local_children_indices;
    }

    template <int dim, int spacedim>
    std::unique_ptr<FiniteElement<1>>
    create_1D_fe(const FiniteElement<dim, spacedim> &fe)
    {
      std::string fe_name = fe.get_name();
      {
        const std::size_t template_starts = fe_name.find_first_of('<');
        Assert(fe_name[template_starts + 1] ==
                 (dim == 1 ? '1' : (dim == 2 ? '2' : '3')),
               ExcInternalError());
        fe_name[template_starts + 1] = '1';
      }
      return FETools::get_fe_by_name<1, 1>(fe_name);
    }

    template <int dim, int spacedim>
    FullMatrix<double>
    get_restriction_matrix(const FiniteElement<dim, spacedim> &fe,
                           const unsigned int                  child)
    {
      auto matrix = fe.get_restriction_matrix(child);

      for (unsigned int c_other = 0; c_other < child; ++c_other)
        {
          auto matrix_other = fe.get_restriction_matrix(c_other);
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            {
              if (fe.restriction_is_additive(i) == true)
                continue;

              bool do_zero = false;
              for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                if (matrix_other(i, j) != 0.)
                  do_zero = true;

              if (do_zero)
                for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                  matrix(i, j) = 0.0;
            }
        }
      return matrix;
    }

    template <int dim>
    bool
    use_fast_hanging_node_algorithm(const DoFHandler<dim> &dof_handler_coarse,
                                    const unsigned int     mg_level_coarse)
    {
      // algorithm is only needed on active levels
      bool use_fast_hanging_node_algorithm =
        mg_level_coarse == numbers::invalid_unsigned_int;

      // algorithm can be only used on meshes consisting of hypercube and
      // simplices
      if (use_fast_hanging_node_algorithm)
        {
          const auto &reference_cells =
            dof_handler_coarse.get_triangulation().get_reference_cells();
          use_fast_hanging_node_algorithm =
            std::all_of(reference_cells.begin(),
                        reference_cells.end(),
                        [](const auto &r) {
                          return r.is_hyper_cube() || r.is_simplex();
                        });
        }

      // local p-refinement is not supported
      if (use_fast_hanging_node_algorithm)
        {
          const auto &fes = dof_handler_coarse.get_fe_collection();

          use_fast_hanging_node_algorithm &=
            std::all_of(fes.begin(), fes.end(), [&fes](const auto &fe) {
              return fes[0].compare_for_domination(fe) ==
                     FiniteElementDomination::Domination::
                       either_element_can_dominate;
            });
        }

      // check that all components are either supported or not
      if (use_fast_hanging_node_algorithm)
        {
          const std::vector<std::vector<bool>> supported_components =
            internal::MatrixFreeFunctions::HangingNodes<
              dim>::compute_supported_components(dof_handler_coarse
                                                   .get_fe_collection());

          use_fast_hanging_node_algorithm &= std::any_of(
            supported_components.begin(),
            supported_components.end(),
            [](const auto &supported_components_per_fe) {
              return std::all_of(supported_components_per_fe.begin(),
                                 supported_components_per_fe.end(),
                                 [](const auto &a) { return a == true; });
            });

          use_fast_hanging_node_algorithm &= std::all_of(
            supported_components.begin(),
            supported_components.end(),
            [](const auto &supported_components_per_fe) {
              return std::all_of(supported_components_per_fe.begin(),
                                 supported_components_per_fe.end(),
                                 [&supported_components_per_fe](const auto &a) {
                                   return a == supported_components_per_fe[0];
                                 });
            });
        }

      return use_fast_hanging_node_algorithm;
    }

  } // namespace

  class FineDoFHandlerViewCell
  {
  public:
    FineDoFHandlerViewCell(
      const std::function<bool()> &has_children_function,
      const std::function<void(std::vector<types::global_dof_index> &)>
        &                                  get_dof_indices_function,
      const std::function<unsigned int()> &active_fe_index_function)
      : has_children_function(has_children_function)
      , get_dof_indices_function(get_dof_indices_function)
      , active_fe_index_function(active_fe_index_function)
    {}

    bool
    has_children() const
    {
      return has_children_function();
    }

    void
    get_dof_indices(std::vector<types::global_dof_index> &dof_indices) const
    {
      get_dof_indices_function(dof_indices);
    }

    unsigned int
    active_fe_index() const
    {
      return active_fe_index_function();
    }


  private:
    const std::function<bool()> has_children_function;
    const std::function<void(std::vector<types::global_dof_index> &)>
                                        get_dof_indices_function;
    const std::function<unsigned int()> active_fe_index_function;
  };



  template <int dim>
  class FineDoFHandlerView
  {
  public:
    FineDoFHandlerView(const DoFHandler<dim> &dof_handler_fine,
                       const DoFHandler<dim> &dof_handler_coarse,
                       const unsigned int     mg_level_fine)
      : dof_handler_fine(dof_handler_fine)
      , dof_handler_coarse(dof_handler_coarse)
      , mg_level_fine(mg_level_fine)
      , communicator(
          dof_handler_fine.get_communicator() /*TODO: fix for different comms*/)
      , cell_id_translator(
          dof_handler_fine.get_triangulation().n_global_coarse_cells(),
          dof_handler_fine.get_triangulation().n_global_levels())
    {
      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_coarse.get_triangulation().n_global_coarse_cells());
      AssertIndexRange(dof_handler_coarse.get_triangulation().n_global_levels(),
                       dof_handler_fine.get_triangulation().n_global_levels() +
                         1);
    }

    void
    reinit(const IndexSet &is_dst_locally_owned,
           const IndexSet &is_dst_remote_input,
           const IndexSet &is_src_locally_owned,
           const bool      check_if_elements_in_is_dst_remote_exist = false)
    {
      IndexSet is_dst_remote = is_dst_remote_input;

      if (check_if_elements_in_is_dst_remote_exist)
        {
          IndexSet is_dst_remote_potentially_relevant = is_dst_remote;
          is_dst_remote.clear();

          is_dst_remote_potentially_relevant.subtract_set(is_dst_locally_owned);

          std::vector<unsigned int> owning_ranks_of_ghosts(
            is_dst_remote_potentially_relevant.n_elements());

          {
            Utilities::MPI::internal::ComputeIndexOwner::
              ConsensusAlgorithmsPayload process(
                is_dst_locally_owned,
                is_dst_remote_potentially_relevant,
                communicator,
                owning_ranks_of_ghosts,
                false);

            Utilities::MPI::ConsensusAlgorithms::Selector<
              std::vector<
                std::pair<types::global_cell_index, types::global_cell_index>>,
              std::vector<unsigned int>>
              consensus_algorithm;
            consensus_algorithm.run(process, communicator);
          }

          for (unsigned i = 0;
               i < is_dst_remote_potentially_relevant.n_elements();
               ++i)
            if (owning_ranks_of_ghosts[i] != numbers::invalid_unsigned_int)
              is_dst_remote.add_index(
                is_dst_remote_potentially_relevant.nth_index_in_set(i));
        }

      // determine owner of remote cells
      std::vector<unsigned int> is_dst_remote_owners(
        is_dst_remote.n_elements());

      Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
        process(is_dst_locally_owned,
                is_dst_remote,
                communicator,
                is_dst_remote_owners,
                true);

      Utilities::MPI::ConsensusAlgorithms::Selector<
        std::vector<
          std::pair<types::global_cell_index, types::global_cell_index>>,
        std::vector<unsigned int>>
        consensus_algorithm;
      consensus_algorithm.run(process, communicator);

      this->is_dst_locally_owned = is_dst_locally_owned;
      this->is_dst_remote        = is_dst_remote;
      this->is_src_locally_owned = is_src_locally_owned;

      this->is_extended_locally_owned =
        mg_level_fine == numbers::invalid_unsigned_int ?
          dof_handler_fine.locally_owned_dofs() :
          dof_handler_fine.locally_owned_mg_dofs(mg_level_fine);

      const auto targets_with_indexset = process.get_requesters();
      std::vector<types::global_dof_index> ghost_indices;

#ifndef DEAL_II_WITH_MPI
      Assert(targets_with_indexset.empty(), ExcInternalError());
#else

      std::map<unsigned int, std::vector<types::global_dof_index>>
                               indices_to_be_sent;
      std::vector<MPI_Request> requests;
      requests.reserve(targets_with_indexset.size());
      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(communicator);

      {
        std::vector<types::global_dof_index> indices;

        for (const auto &i : targets_with_indexset)
          {
            if (i.first == my_rank)
              continue;

            indices_to_be_sent[i.first] = {};
            std::vector<types::global_dof_index> &buffer =
              indices_to_be_sent[i.first];

            for (auto cell_id : i.second)
              {
                typename DoFHandler<dim>::cell_iterator cell(
                  *dof_handler_fine.get_triangulation().create_cell_iterator(
                    cell_id_translator.to_cell_id(cell_id)),
                  &dof_handler_fine);

                indices.resize(cell->get_fe().n_dofs_per_cell());

                if (mg_level_fine == numbers::invalid_unsigned_int)
                  cell->get_dof_indices(indices);
                else
                  cell->get_mg_dof_indices(indices);

                buffer.push_back(cell->active_fe_index());
                buffer.insert(buffer.end(), indices.begin(), indices.end());
              }

            requests.emplace_back();

            const auto ierr_1 = MPI_Isend(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              i.first,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &requests.back());
            AssertThrowMPI(ierr_1);
          }
      }

      // process local cells
      {
        std::vector<types::global_dof_index> indices;

        for (const auto id : is_dst_locally_owned)
          {
            const auto cell_id = cell_id_translator.to_cell_id(id);

            typename DoFHandler<dim>::cell_iterator cell_(
              *dof_handler_fine.get_triangulation().create_cell_iterator(
                cell_id),
              &dof_handler_fine);

            indices.resize(cell_->get_fe().n_dofs_per_cell());

            if (mg_level_fine == numbers::invalid_unsigned_int)
              cell_->get_dof_indices(indices);
            else
              cell_->get_mg_dof_indices(indices);

            for (const auto i : indices)
              if (!is_extended_locally_owned.is_element(i))
                ghost_indices.push_back(i);
          }
      }

      {
        std::map<unsigned int, std::vector<types::global_dof_index>>
          rank_to_ids;
        for (unsigned int i = 0; i < is_dst_remote_owners.size(); ++i)
          rank_to_ids[is_dst_remote_owners[i]].push_back(
            is_dst_remote.nth_index_in_set(i));

        for (const auto &pair : rank_to_ids)
          {
            // above we skip messages sent to myself, so also skip MPI_Probe
            if (pair.first == my_rank)
              continue;

            MPI_Status status;
            const int  ierr_1 = MPI_Probe(
              MPI_ANY_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &status);
            AssertThrowMPI(ierr_1);

            std::vector<types::global_dof_index> buffer;

            int       message_length;
            const int ierr_2 = MPI_Get_count(
              &status,
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              &message_length);
            AssertThrowMPI(ierr_2);

            buffer.resize(message_length);

            const int ierr_3 = MPI_Recv(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              status.MPI_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr_3);

            for (unsigned int i = 0; i < buffer.size();)
              {
                const unsigned int dofs_per_cell =
                  dof_handler_fine.get_fe(buffer[i]).n_dofs_per_cell();
                ++i;
                for (unsigned int j = 0; j < dofs_per_cell; ++j, ++i)
                  {
                    AssertIndexRange(i, buffer.size());
                    if (!is_extended_locally_owned.is_element(buffer[i]))
                      ghost_indices.push_back(buffer[i]);
                  }
              }

            const unsigned int rank = status.MPI_SOURCE;

            const auto ids = rank_to_ids[rank];

            std::vector<types::global_dof_index> indices;

            for (unsigned int i = 0, k = 0; i < ids.size(); ++i)
              {
                const unsigned int active_fe_index = buffer[k++];

                indices.resize(
                  dof_handler_fine.get_fe(active_fe_index).n_dofs_per_cell());

                for (unsigned int j = 0; j < indices.size(); ++j, ++k)
                  indices[j] = buffer[k];
                map[ids[i]] = {active_fe_index, indices};
              }
          }

        if (!requests.empty())
          {
            const int ierr_1 = MPI_Waitall(requests.size(),
                                           requests.data(),
                                           MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr_1);
          }
      }
#endif

      std::sort(ghost_indices.begin(), ghost_indices.end());

      this->is_extended_ghosts =
        IndexSet(mg_level_fine == numbers::invalid_unsigned_int ?
                   dof_handler_fine.n_dofs() :
                   dof_handler_fine.n_dofs(mg_level_fine));
      this->is_extended_ghosts.add_indices(ghost_indices.begin(),
                                           ghost_indices.end());
      this->is_extended_ghosts.subtract_set(this->is_extended_locally_owned);
    }

    FineDoFHandlerViewCell
    get_cell(const typename DoFHandler<dim>::cell_iterator &cell) const
    {
      const auto id = this->cell_id_translator.translate(cell);

      const bool is_cell_locally_owned =
        this->is_dst_locally_owned.is_element(id);
      const bool is_cell_remotly_owned = this->is_dst_remote.is_element(id);

      const bool has_cell_any_children = [&]() {
        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
             ++i)
          {
            const auto j = this->cell_id_translator.translate(cell, i);

            if (this->is_dst_locally_owned.is_element(j))
              return true;

            if (this->is_dst_remote.is_element(j))
              return true;
          }

        return false;
      }();

      return FineDoFHandlerViewCell(
        [has_cell_any_children]() { return has_cell_any_children; },
        [cell, is_cell_locally_owned, is_cell_remotly_owned, id, this](
          auto &dof_indices) {
          if (is_cell_locally_owned)
            {
              typename DoFHandler<dim>::cell_iterator cell_fine(
                *dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell->id()),
                &dof_handler_fine);
              if (mg_level_fine == numbers::invalid_unsigned_int)
                cell_fine->get_dof_indices(dof_indices);
              else
                cell_fine->get_mg_dof_indices(dof_indices);
            }
          else if (is_cell_remotly_owned)
            {
              dof_indices = map.at(id).second;
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
            }
        },
        [cell, is_cell_locally_owned, is_cell_remotly_owned, id, this]()
          -> unsigned int {
          if (is_cell_locally_owned)
            {
              return (typename DoFHandler<dim>::cell_iterator(
                        *dof_handler_fine.get_triangulation()
                           .create_cell_iterator(cell->id()),
                        &dof_handler_fine))
                ->active_fe_index();
            }
          else if (is_cell_remotly_owned)
            {
              return map.at(id).first;
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
              return 0;
            }
        });
    }

    FineDoFHandlerViewCell
    get_cell(const typename DoFHandler<dim>::cell_iterator &cell,
             const unsigned int                             c) const
    {
      const auto id = this->cell_id_translator.translate(cell, c);

      const bool is_cell_locally_owned =
        this->is_dst_locally_owned.is_element(id);
      const bool is_cell_remotely_owned = this->is_dst_remote.is_element(id);

      return FineDoFHandlerViewCell(
        [cell]() {
          AssertThrow(false, ExcNotImplemented()); // currently we do not need
                                                   // children of children

          return false;
        },
        [cell, is_cell_locally_owned, is_cell_remotely_owned, c, id, this](
          auto &dof_indices) {
          if (is_cell_locally_owned)
            {
              const auto cell_fine =
                (typename DoFHandler<dim>::cell_iterator(
                   *dof_handler_fine.get_triangulation().create_cell_iterator(
                     cell->id()),
                   &dof_handler_fine))
                  ->child(c);

              if (mg_level_fine == numbers::invalid_unsigned_int)
                cell_fine->get_dof_indices(dof_indices);
              else
                cell_fine->get_mg_dof_indices(dof_indices);
            }
          else if (is_cell_remotely_owned)
            {
              dof_indices = map.at(id).second;
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
            }
        },
        []() -> unsigned int {
          AssertThrow(false, ExcNotImplemented()); // currently we do not need
                                                   // active_fe_index() for
                                                   // children

          return 0;
        });
    }

    const IndexSet &
    locally_owned_dofs() const
    {
      return is_extended_locally_owned;
    }

    const IndexSet &
    locally_relevant_dofs() const
    {
      return is_extended_ghosts;
    }

  private:
    const DoFHandler<dim> &dof_handler_fine;
    const DoFHandler<dim> &dof_handler_coarse;
    const unsigned int     mg_level_fine;
    const MPI_Comm         communicator;

  protected:
    const CellIDTranslator<dim> cell_id_translator;

  private:
    IndexSet is_dst_locally_owned;
    IndexSet is_dst_remote;
    IndexSet is_src_locally_owned;


    IndexSet is_extended_locally_owned;
    IndexSet is_extended_ghosts;

    std::map<unsigned int,
             std::pair<unsigned int, std::vector<types::global_dof_index>>>
      map;
  };

  template <int dim>
  class GlobalCoarseningFineDoFHandlerView : public FineDoFHandlerView<dim>
  {
  public:
    GlobalCoarseningFineDoFHandlerView(const DoFHandler<dim> &dof_handler_dst,
                                       const DoFHandler<dim> &dof_handler_src,
                                       const unsigned int     mg_level_fine,
                                       const unsigned int     mg_level_coarse)
      : FineDoFHandlerView<dim>(dof_handler_dst, dof_handler_src, mg_level_fine)
    {
      Assert((mg_level_fine == numbers::invalid_unsigned_int &&
              mg_level_coarse == numbers::invalid_unsigned_int) ||
               (mg_level_coarse + 1 == mg_level_fine),
             ExcNotImplemented());

      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      const auto fine_operation = [&](const auto &cell) {
        is_dst_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
      };

      const auto coarse_operation = [&](const auto &cell) {
        is_src_locally_owned.add_index(
          this->cell_id_translator.translate(cell));

        // in the case of global coarsening identity transfer is possible
        if (mg_level_coarse == numbers::invalid_unsigned_int)
          is_dst_remote.add_index(this->cell_id_translator.translate(cell));

        if (cell->level() + 1u == tria_dst.n_global_levels())
          return;

        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
             ++i)
          is_dst_remote.add_index(this->cell_id_translator.translate(cell, i));
      };

      loop_over_active_or_level_cells(tria_dst, mg_level_fine, fine_operation);
      loop_over_active_or_level_cells(tria_src,
                                      mg_level_coarse,
                                      coarse_operation);

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   true);
    }
  };

  template <int dim>
  class PermutationFineDoFHandlerView : public internal::FineDoFHandlerView<dim>
  {
  public:
    PermutationFineDoFHandlerView(const DoFHandler<dim> &dof_handler_dst,
                                  const DoFHandler<dim> &dof_handler_src,
                                  const unsigned int     mg_level_fine,
                                  const unsigned int     mg_level_coarse)
      : internal::FineDoFHandlerView<dim>(dof_handler_dst,
                                          dof_handler_src,
                                          mg_level_fine)
    {
      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      const auto fine_operation = [&](const auto &cell) {
        is_dst_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
      };

      const auto coarse_operation = [&](const auto &cell) {
        is_src_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
        is_dst_remote.add_index(this->cell_id_translator.translate(cell));
      };

      loop_over_active_or_level_cells(tria_dst, mg_level_fine, fine_operation);
      loop_over_active_or_level_cells(tria_src,
                                      mg_level_coarse,
                                      coarse_operation);

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   false);
    }
  };

  class MGTwoLevelTransferImplementation
  {
    template <int dim, typename Number>
    static std::shared_ptr<const Utilities::MPI::Partitioner>
    create_coarse_partitioner(
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_coarse)
    {
      IndexSet locally_relevant_dofs =
        (mg_level_coarse == numbers::invalid_unsigned_int) ?
          DoFTools::extract_locally_active_dofs(dof_handler_coarse) :
          DoFTools::extract_locally_active_level_dofs(dof_handler_coarse,
                                                      mg_level_coarse);

      std::vector<types::global_dof_index> locally_relevant_dofs_temp;

      for (const auto i : locally_relevant_dofs)
        {
          if (locally_relevant_dofs.is_element(i) == false)
            locally_relevant_dofs_temp.emplace_back(i);

          const auto constraints = constraints_coarse.get_constraint_entries(i);

          if (constraints)
            for (const auto &p : *constraints)
              if (locally_relevant_dofs.is_element(p.first) == false)
                locally_relevant_dofs_temp.emplace_back(p.first);
        }

      std::sort(locally_relevant_dofs_temp.begin(),
                locally_relevant_dofs_temp.end());
      locally_relevant_dofs.add_indices(locally_relevant_dofs_temp.begin(),
                                        locally_relevant_dofs_temp.end());

      return std::make_shared<Utilities::MPI::Partitioner>(
        mg_level_coarse == numbers::invalid_unsigned_int ?
          dof_handler_coarse.locally_owned_dofs() :
          dof_handler_coarse.locally_owned_mg_dofs(mg_level_coarse),
        locally_relevant_dofs,
        dof_handler_coarse.get_communicator());
    }



    template <int dim, typename Number>
    static void
    compute_weights(
      const DoFHandler<dim> &                  dof_handler_fine,
      const unsigned int                       mg_level_fine,
      const dealii::AffineConstraints<Number> &constraints_fine,
      const MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &                                         transfer,
      LinearAlgebra::distributed::Vector<Number> &touch_count)
    {
      touch_count.reinit(transfer.partitioner_fine);

      for (const auto i : transfer.constraint_info_fine.dof_indices)
        touch_count.local_element(i) += 1;
      touch_count.compress(VectorOperation::add);

      IndexSet locally_relevant_dofs;
      if (mg_level_fine == numbers::invalid_unsigned_int)
        DoFTools::extract_locally_relevant_dofs(dof_handler_fine,
                                                locally_relevant_dofs);
      else
        DoFTools::extract_locally_relevant_level_dofs(dof_handler_fine,
                                                      mg_level_fine,
                                                      locally_relevant_dofs);

      const auto partitioner_fine_ =
        std::make_shared<Utilities::MPI::Partitioner>(
          mg_level_fine == numbers::invalid_unsigned_int ?
            dof_handler_fine.locally_owned_dofs() :
            dof_handler_fine.locally_owned_mg_dofs(mg_level_fine),
          locally_relevant_dofs,
          dof_handler_fine.get_communicator());
      transfer.vec_fine.reinit(transfer.partitioner_fine);

      LinearAlgebra::distributed::Vector<Number> touch_count_;
      touch_count_.reinit(partitioner_fine_);

      touch_count_.copy_locally_owned_data_from(touch_count);

      for (unsigned int i = 0; i < touch_count_.locally_owned_size(); ++i)
        touch_count_.local_element(i) =
          constraints_fine.is_constrained(
            touch_count_.get_partitioner()->local_to_global(i)) ?
            Number(0.) :
            Number(1.) / touch_count_.local_element(i);

      touch_count_.update_ghost_values();

      // copy weights to other indexset
      touch_count.reinit(transfer.partitioner_fine);
      touch_count.copy_locally_owned_data_from(touch_count_);
      touch_count.update_ghost_values();
    }



    /**
     * Try to compress weights to Utilities::pow(3, dim) doubles.
     * Weights of a cell can be compressed if all components have the
     * same weights.
     */
    template <int dim, typename Number>
    static void
    compress_weights(
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      unsigned int n_cells = 0;
      for (const auto &scheme : transfer.schemes)
        n_cells += scheme.n_coarse_cells;

      std::vector<std::array<Number, Utilities::pow(3, dim)>> masks(n_cells);

      const Number *weight_ptr = transfer.weights.data();
      std::array<Number, Utilities::pow(3, dim)> *mask_ptr = masks.data();

      // (try to) compress weights for each cell
      for (const auto &scheme : transfer.schemes)
        for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
             ++cell, weight_ptr += scheme.n_dofs_per_cell_fine, ++mask_ptr)
          if (!compute_weights_fe_q_dofs_by_entity<dim, -1, Number>(
                weight_ptr,
                transfer.n_components,
                scheme.degree_fine + 1,
                mask_ptr->data()))
            return;

      // vectorize weights
      AlignedVector<VectorizedArray<Number>> masks_vectorized;

      const unsigned int n_lanes = VectorizedArray<Number>::size();

      unsigned int c = 0;

      for (const auto &scheme : transfer.schemes)
        {
          for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
               cell += n_lanes)
            {
              const unsigned int n_lanes_filled =
                (cell + n_lanes > scheme.n_coarse_cells) ?
                  (scheme.n_coarse_cells - cell) :
                  n_lanes;

              std::array<VectorizedArray<Number>, Utilities::pow(3, dim)> mask;
              mask.fill(VectorizedArray<Number>());

              for (unsigned int v = 0; v < n_lanes_filled; ++v, ++c)
                {
                  for (unsigned int i = 0;
                       i < Utilities::pow<unsigned int>(3, dim);
                       ++i)
                    mask[i][v] = masks[c][i];
                }

              masks_vectorized.insert_back(mask.begin(), mask.end());
            }
        }

      // copy result
      transfer.weights_compressed = masks_vectorized;
    }



  public:
    template <int dim, typename Number>
    static void
    reinit_geometric_transfer(
      const DoFHandler<dim> &                  dof_handler_fine,
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_fine,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_fine,
      const unsigned int                       mg_level_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      Assert((mg_level_fine == numbers::invalid_unsigned_int &&
              mg_level_coarse == numbers::invalid_unsigned_int) ||
               (mg_level_coarse + 1 == mg_level_fine),
             ExcNotImplemented());

      // if (mg_level_fine != numbers::invalid_unsigned_int)
      //   AssertIndexRange(mg_level_fine,
      //                    MGTools::max_level_for_coarse_mesh(
      //                      dof_handler_fine.get_triangulation()) +
      //                      1);
      //
      // if (mg_level_coarse != numbers::invalid_unsigned_int)
      //   AssertIndexRange(mg_level_coarse,
      //                    MGTools::max_level_for_coarse_mesh(
      //                      dof_handler_coarse.get_triangulation()) +
      //                      1);

      const GlobalCoarseningFineDoFHandlerView<dim> view(dof_handler_fine,
                                                         dof_handler_coarse,
                                                         mg_level_fine,
                                                         mg_level_coarse);

      // gather ranges for active FE indices on both fine and coarse dofhandlers
      std::array<unsigned int, 2> min_active_fe_indices = {
        {std::numeric_limits<unsigned int>::max(),
         std::numeric_limits<unsigned int>::max()}};
      std::array<unsigned int, 2> max_active_fe_indices = {{0, 0}};

      loop_over_active_or_level_cells(
        dof_handler_fine, mg_level_fine, [&](const auto &cell) {
          min_active_fe_indices[0] =
            std::min<unsigned int>(min_active_fe_indices[0],
                                   cell->active_fe_index());
          max_active_fe_indices[0] =
            std::max<unsigned int>(max_active_fe_indices[0],
                                   cell->active_fe_index());
        });

      loop_over_active_or_level_cells(
        dof_handler_coarse, mg_level_coarse, [&](const auto &cell) {
          min_active_fe_indices[1] =
            std::min<unsigned int>(min_active_fe_indices[1],
                                   cell->active_fe_index());
          max_active_fe_indices[1] =
            std::max<unsigned int>(max_active_fe_indices[1],
                                   cell->active_fe_index());
        });

      const auto comm = dof_handler_fine.get_communicator();

      Assert(comm == dof_handler_coarse.get_communicator(),
             ExcNotImplemented());

      ArrayView<unsigned int> temp_min(min_active_fe_indices);
      ArrayView<unsigned int> temp_max(max_active_fe_indices);
      Utilities::MPI::min(temp_min, comm, temp_min);
      Utilities::MPI::max(temp_max, comm, temp_max);

      // make sure that hp is used neither on the coarse nor on the fine
      // dofhandler
      AssertDimension(min_active_fe_indices[0], max_active_fe_indices[0]);
      AssertDimension(min_active_fe_indices[1], max_active_fe_indices[1]);

      // set up two mg-schemes
      //   (0) no refinement -> identity
      //   (1) h-refinement
      transfer.schemes.resize(2);

      const unsigned int fe_index_fine   = min_active_fe_indices[0];
      const unsigned int fe_index_coarse = min_active_fe_indices[1];

      const auto &fe_fine   = dof_handler_fine.get_fe(fe_index_fine);
      const auto &fe_coarse = dof_handler_coarse.get_fe(fe_index_coarse);

      // extract number of components
      AssertDimension(fe_fine.n_components(), fe_coarse.n_components());

      transfer.n_components = fe_fine.n_components();

      const auto reference_cell = dof_handler_fine.get_fe(0).reference_cell();

      // create partitioners and vectors for internal purposes
      {
        // ... for fine mesh
        {
          transfer.partitioner_fine =
            std::make_shared<Utilities::MPI::Partitioner>(
              view.locally_owned_dofs(),
              view.locally_relevant_dofs(),
              dof_handler_fine.get_communicator());
          transfer.vec_fine.reinit(transfer.partitioner_fine);
        }

        // ... coarse mesh (needed since user vector might be const)
        {
          transfer.partitioner_coarse =
            create_coarse_partitioner(dof_handler_coarse,
                                      constraints_coarse,
                                      mg_level_coarse);
          transfer.vec_coarse.reinit(transfer.partitioner_coarse);
        }
      }

      // helper function: to process the fine level cells; function @p fu_non_refined is
      // performed on cells that are not refined and @fu_refined is performed on
      // children of cells that are refined
      const auto process_cells = [&](const auto &fu_non_refined,
                                     const auto &fu_refined) {
        loop_over_active_or_level_cells(
          dof_handler_coarse, mg_level_coarse, [&](const auto &cell_coarse) {
            if (mg_level_coarse == numbers::invalid_unsigned_int)
              {
                // get a reference to the equivalent cell on the fine
                // triangulation
                const auto cell_coarse_on_fine_mesh =
                  view.get_cell(cell_coarse);

                // check if cell has children
                if (cell_coarse_on_fine_mesh.has_children())
                  // ... cell has children -> process children
                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       c++)
                    fu_refined(cell_coarse, view.get_cell(cell_coarse, c), c);
                else // ... cell has no children -> process cell
                  fu_non_refined(cell_coarse, cell_coarse_on_fine_mesh);
              }
            else
              {
                // check if cell has children
                if (cell_coarse->has_children())
                  // ... cell has children -> process children
                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       c++)
                    fu_refined(cell_coarse, view.get_cell(cell_coarse, c), c);
              }
          });
      };

      // check if FE is the same
      AssertDimension(fe_coarse.n_dofs_per_cell(), fe_fine.n_dofs_per_cell());


      const bool is_feq =
        fe_fine.n_base_elements() == 1 &&
        (dynamic_cast<const FE_Q<dim> *>(&fe_fine.base_element(0)) != nullptr);

      // number of dofs on coarse and fine cells
      transfer.schemes[0].n_dofs_per_cell_coarse =
        transfer.schemes[0].n_dofs_per_cell_fine =
          transfer.schemes[1].n_dofs_per_cell_coarse =
            fe_coarse.n_dofs_per_cell();
      transfer.schemes[1].n_dofs_per_cell_fine =
        is_feq ? (fe_fine.n_components() *
                  Utilities::pow(2 * fe_fine.degree + 1, dim)) :
                 (fe_coarse.n_dofs_per_cell() *
                  GeometryInfo<dim>::max_children_per_cell);

      // degree of FE on coarse and fine cell
      transfer.schemes[0].degree_coarse   = transfer.schemes[0].degree_fine =
        transfer.schemes[1].degree_coarse = fe_coarse.degree;
      transfer.schemes[1].degree_fine =
        is_feq ? (fe_coarse.degree * 2) : (fe_coarse.degree * 2 + 1);

      // continuous or discontinuous
      transfer.fine_element_is_continuous = fe_fine.n_dofs_per_vertex() > 0;

      // count coarse cells for each scheme (0, 1)
      {
        transfer.schemes[0].n_coarse_cells = 0; // reset
        transfer.schemes[1].n_coarse_cells = 0;

        // count by looping over all coarse cells
        process_cells(
          [&](const auto &, const auto &) {
            transfer.schemes[0].n_coarse_cells++;
          },
          [&](const auto &, const auto &, const auto c) {
            if (c == 0)
              transfer.schemes[1].n_coarse_cells++;
          });
      }


      const auto cell_local_children_indices =
        (reference_cell == ReferenceCells::get_hypercube<dim>()) ?
          get_child_offsets<dim>(transfer.schemes[0].n_dofs_per_cell_coarse,
                                 is_feq ? fe_fine.degree : (fe_fine.degree + 1),
                                 fe_fine.degree) :
          get_child_offsets_general<dim>(
            transfer.schemes[0].n_dofs_per_cell_coarse);

      std::vector<unsigned int> n_dof_indices_fine(transfer.schemes.size() + 1);
      std::vector<unsigned int> n_dof_indices_coarse(transfer.schemes.size() +
                                                     1);

      for (unsigned int i = 0; i < transfer.schemes.size(); ++i)
        {
          n_dof_indices_fine[i + 1] = transfer.schemes[i].n_dofs_per_cell_fine *
                                      transfer.schemes[i].n_coarse_cells;
          n_dof_indices_coarse[i + 1] =
            transfer.schemes[i].n_dofs_per_cell_coarse *
            transfer.schemes[i].n_coarse_cells;
        }

      for (unsigned int i = 0; i < transfer.schemes.size(); ++i)
        {
          n_dof_indices_fine[i + 1] += n_dof_indices_fine[i];
          n_dof_indices_coarse[i + 1] += n_dof_indices_coarse[i];
        }

      // indices
      {
        std::vector<types::global_dof_index> local_dof_indices(
          transfer.schemes[0].n_dofs_per_cell_coarse);

        // ---------------------- lexicographic_numbering ----------------------
        std::vector<unsigned int> lexicographic_numbering_fine;
        std::vector<unsigned int> lexicographic_numbering_coarse;
        if (reference_cell == ReferenceCells::get_hypercube<dim>())
          {
            const Quadrature<1> dummy_quadrature(
              std::vector<Point<1>>(1, Point<1>()));
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering_fine = shape_info.lexicographic_numbering;
            shape_info.reinit(dummy_quadrature, fe_coarse, 0);
            lexicographic_numbering_coarse = shape_info.lexicographic_numbering;
          }
        else
          {
            const auto dummy_quadrature =
              reference_cell.template get_gauss_type_quadrature<dim>(1);
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering_fine = shape_info.lexicographic_numbering;
            shape_info.reinit(dummy_quadrature, fe_coarse, 0);
            lexicographic_numbering_coarse = shape_info.lexicographic_numbering;
          }

        // ------------------------------ indices ------------------------------
        std::vector<types::global_dof_index> level_dof_indices_fine_0(
          transfer.schemes[0].n_dofs_per_cell_fine);
        std::vector<types::global_dof_index> level_dof_indices_fine_1(
          transfer.schemes[1].n_dofs_per_cell_fine);

        unsigned int cell_no_0 = 0;
        unsigned int cell_no_1 = transfer.schemes[0].n_coarse_cells;

        transfer.constraint_info_coarse.reinit(
          dof_handler_coarse,
          transfer.schemes[0].n_coarse_cells +
            transfer.schemes[1].n_coarse_cells,
          constraints_coarse.n_constraints() > 0 &&
            use_fast_hanging_node_algorithm(dof_handler_coarse,
                                            mg_level_coarse));

        transfer.constraint_info_fine.reinit(
          transfer.schemes[0].n_coarse_cells +
          transfer.schemes[1].n_coarse_cells);

        process_cells(
          [&](const auto &cell_coarse, const auto &cell_fine) {
            // parent
            {
              transfer.constraint_info_coarse.read_dof_indices(
                cell_no_0,
                mg_level_coarse,
                cell_coarse,
                constraints_coarse,
                transfer.partitioner_coarse);
            }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[0].n_dofs_per_cell_coarse;
                   i++)
                level_dof_indices_fine_0[i] =
                  local_dof_indices[lexicographic_numbering_fine[i]];

              transfer.constraint_info_fine.read_dof_indices(
                cell_no_0, level_dof_indices_fine_0, transfer.partitioner_fine);
            }

            // move pointers
            {
              cell_no_0++;
            }
          },
          [&](const auto &cell_coarse, const auto &cell_fine, const auto c) {
            // parent (only once at the beginning)
            if (c == 0)
              {
                transfer.constraint_info_coarse.read_dof_indices(
                  cell_no_1,
                  mg_level_coarse,
                  cell_coarse,
                  constraints_coarse,
                  transfer.partitioner_coarse);

                level_dof_indices_fine_1.assign(level_dof_indices_fine_1.size(),
                                                numbers::invalid_dof_index);
              }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[1].n_dofs_per_cell_coarse;
                   ++i)
                {
                  const auto index =
                    local_dof_indices[lexicographic_numbering_fine[i]];

                  Assert(level_dof_indices_fine_1
                               [cell_local_children_indices[c][i]] ==
                             numbers::invalid_dof_index ||
                           level_dof_indices_fine_1
                               [cell_local_children_indices[c][i]] == index,
                         ExcInternalError());

                  level_dof_indices_fine_1[cell_local_children_indices[c][i]] =
                    index;
                }
            }

            // move pointers (only once at the end)
            if (c + 1 == GeometryInfo<dim>::max_children_per_cell)
              {
                transfer.constraint_info_fine.read_dof_indices(
                  cell_no_1,
                  level_dof_indices_fine_1,
                  transfer.partitioner_fine);

                cell_no_1++;
              }
          });

        transfer.constraint_info_coarse.finalize();
        transfer.constraint_info_fine.finalize();
      }


      // ------------- prolongation matrix (0) -> identity matrix --------------

      // nothing to do since for identity prolongation matrices a short-cut
      // code path is used during prolongation/restriction

      // ----------------------- prolongation matrix (1) -----------------------
      {
        AssertDimension(fe_fine.n_base_elements(), 1);
        if (reference_cell == ReferenceCells::get_hypercube<dim>())
          {
            const auto fe = create_1D_fe(fe_fine.base_element(0));

            std::vector<unsigned int> renumbering(fe->n_dofs_per_cell());
            {
              AssertIndexRange(fe->n_dofs_per_vertex(), 2);
              renumbering[0] = 0;
              for (unsigned int i = 0; i < fe->dofs_per_line; ++i)
                renumbering[i + fe->n_dofs_per_vertex()] =
                  GeometryInfo<1>::vertices_per_cell * fe->n_dofs_per_vertex() +
                  i;
              if (fe->n_dofs_per_vertex() > 0)
                renumbering[fe->n_dofs_per_cell() - fe->n_dofs_per_vertex()] =
                  fe->n_dofs_per_vertex();
            }

            // TODO: data structures are saved in form of DG data structures
            // here
            const unsigned int shift =
              is_feq ? (fe->n_dofs_per_cell() - fe->n_dofs_per_vertex()) :
                       (fe->n_dofs_per_cell());
            const unsigned int n_child_dofs_1d =
              is_feq ? (fe->n_dofs_per_cell() * 2 - fe->n_dofs_per_vertex()) :
                       (fe->n_dofs_per_cell() * 2);

            {
              transfer.schemes[1].prolongation_matrix_1d.resize(
                fe->n_dofs_per_cell() * n_child_dofs_1d);

              for (unsigned int c = 0;
                   c < GeometryInfo<1>::max_children_per_cell;
                   ++c)
                for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
                  for (unsigned int j = 0; j < fe->n_dofs_per_cell(); ++j)
                    transfer.schemes[1]
                      .prolongation_matrix_1d[i * n_child_dofs_1d + j +
                                              c * shift] =
                      fe->get_prolongation_matrix(c)(renumbering[j],
                                                     renumbering[i]);
            }
            {
              transfer.schemes[1].restriction_matrix_1d.resize(
                fe->n_dofs_per_cell() * n_child_dofs_1d);

              for (unsigned int c = 0;
                   c < GeometryInfo<1>::max_children_per_cell;
                   ++c)
                {
                  const auto matrix = get_restriction_matrix(*fe, c);
                  for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
                    for (unsigned int j = 0; j < fe->n_dofs_per_cell(); ++j)
                      transfer.schemes[1]
                        .restriction_matrix_1d[i * n_child_dofs_1d + j +
                                               c * shift] +=
                        matrix(renumbering[i], renumbering[j]);
                }
            }
          }
        else
          {
            const auto &       fe              = fe_fine.base_element(0);
            const unsigned int n_dofs_per_cell = fe.n_dofs_per_cell();

            {
              transfer.schemes[1].prolongation_matrix.resize(
                n_dofs_per_cell * n_dofs_per_cell *
                GeometryInfo<dim>::max_children_per_cell);

              for (unsigned int c = 0;
                   c < GeometryInfo<dim>::max_children_per_cell;
                   ++c)
                for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                    transfer.schemes[1].prolongation_matrix
                      [i * n_dofs_per_cell *
                         GeometryInfo<dim>::max_children_per_cell +
                       j + c * n_dofs_per_cell] =
                      fe.get_prolongation_matrix(c)(j, i);
            }
            {
              transfer.schemes[1].restriction_matrix.resize(
                n_dofs_per_cell * n_dofs_per_cell *
                GeometryInfo<dim>::max_children_per_cell);

              for (unsigned int c = 0;
                   c < GeometryInfo<dim>::max_children_per_cell;
                   ++c)
                {
                  const auto matrix = get_restriction_matrix(fe, c);
                  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                      transfer.schemes[1].restriction_matrix
                        [i * n_dofs_per_cell *
                           GeometryInfo<dim>::max_children_per_cell +
                         j + c * n_dofs_per_cell] += matrix(i, j);
                }
            }
          }
      }


      // ------------------------------- weights -------------------------------
      if (transfer.fine_element_is_continuous)
        {
          // compute weights globally
          LinearAlgebra::distributed::Vector<Number> weight_vector;
          compute_weights(dof_handler_fine,
                          mg_level_fine,
                          constraints_fine,
                          transfer,
                          weight_vector);

          // ... and store them cell-wise a DG format
          transfer.weights.resize(n_dof_indices_fine.back());

          Number *weights_0 = transfer.weights.data() + n_dof_indices_fine[0];
          Number *weights_1 = transfer.weights.data() + n_dof_indices_fine[1];
          unsigned int *dof_indices_fine_0 =
            transfer.constraint_info_fine.dof_indices.data() +
            n_dof_indices_fine[0];
          unsigned int *dof_indices_fine_1 =
            transfer.constraint_info_fine.dof_indices.data() +
            n_dof_indices_fine[1];

          process_cells(
            [&](const auto &, const auto &) {
              for (unsigned int i = 0;
                   i < transfer.schemes[0].n_dofs_per_cell_fine;
                   ++i)
                weights_0[i] =
                  weight_vector.local_element(dof_indices_fine_0[i]);

              dof_indices_fine_0 += transfer.schemes[0].n_dofs_per_cell_fine;
              weights_0 += transfer.schemes[0].n_dofs_per_cell_fine;
            },
            [&](const auto &, const auto &, const auto c) {
              for (unsigned int i = 0;
                   i < transfer.schemes[1].n_dofs_per_cell_coarse;
                   ++i)
                weights_1[cell_local_children_indices[c][i]] =
                  weight_vector.local_element(
                    dof_indices_fine_1[cell_local_children_indices[c][i]]);

              // move pointers (only once at the end)
              if (c + 1 == GeometryInfo<dim>::max_children_per_cell)
                {
                  dof_indices_fine_1 +=
                    transfer.schemes[1].n_dofs_per_cell_fine;
                  weights_1 += transfer.schemes[1].n_dofs_per_cell_fine;
                }
            });

          if (is_feq)
            compress_weights(transfer);
        }
    }



    template <int dim, typename Number>
    static void
    reinit_polynomial_transfer(
      const DoFHandler<dim> &                  dof_handler_fine,
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_fine,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_fine,
      const unsigned int                       mg_level_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      Assert(
        mg_level_fine == numbers::invalid_unsigned_int ||
          mg_level_fine <= MGTools::max_level_for_coarse_mesh(
                             dof_handler_fine.get_triangulation()),
        ExcMessage(
          "Polynomial transfer is only allowed on the active level "
          "(numbers::invalid_unsigned_int) or on refinement levels without "
          "hanging nodes."));
      Assert(
        mg_level_coarse == numbers::invalid_unsigned_int ||
          mg_level_coarse <= MGTools::max_level_for_coarse_mesh(
                               dof_handler_coarse.get_triangulation()),
        ExcMessage(
          "Polynomial transfer is only allowed on the active level "
          "(numbers::invalid_unsigned_int) or on refinement levels without "
          "hanging nodes."));

      const PermutationFineDoFHandlerView<dim> view(dof_handler_fine,
                                                    dof_handler_coarse,
                                                    mg_level_fine,
                                                    mg_level_coarse);

      // TODO: adjust assert
      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_active_cells(),
        dof_handler_coarse.get_triangulation().n_global_active_cells());

      // extract number of components
      AssertDimension(dof_handler_fine.get_fe_collection().n_components(),
                      dof_handler_coarse.get_fe_collection().n_components());

      transfer.n_components =
        dof_handler_fine.get_fe_collection().n_components();

      transfer.fine_element_is_continuous =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) {
                      return fe.n_dofs_per_cell() == 0 ||
                             fe.n_dofs_per_vertex() > 0;
                    });

#ifdef DEBUG
      const bool fine_element_is_discontinuous =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) {
                      return fe.n_dofs_per_cell() == 0 ||
                             fe.n_dofs_per_vertex() == 0;
                    });

      Assert(transfer.fine_element_is_continuous !=
               fine_element_is_discontinuous,
             ExcNotImplemented());
#endif

      const bool is_feq =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) {
                      return fe.n_base_elements() == 1 &&
                             (dynamic_cast<const FE_Q<dim> *>(
                                &fe.base_element(0)) != nullptr);
                    });

      const auto process_cells = [&](const auto &fu) {
        loop_over_active_or_level_cells(
          dof_handler_coarse, mg_level_coarse, [&](const auto &cell_coarse) {
            const auto cell_coarse_on_fine_mesh = view.get_cell(cell_coarse);
            fu(cell_coarse, &cell_coarse_on_fine_mesh);
          });
      };

      std::map<std::pair<unsigned int, unsigned int>, unsigned int>
        fe_index_pairs;

      process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
        fe_index_pairs.emplace(
          std::pair<unsigned int, unsigned int>(cell_coarse->active_fe_index(),
                                                cell_fine->active_fe_index()),
          0);
      });

      unsigned int counter = 0;
      for (auto &f : fe_index_pairs)
        f.second = counter++;

      transfer.schemes.resize(fe_index_pairs.size());

      // extract number of coarse cells
      {
        for (auto &scheme : transfer.schemes)
          scheme.n_coarse_cells = 0;
        process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
          transfer
            .schemes[fe_index_pairs[std::pair<unsigned int, unsigned int>(
              cell_coarse->active_fe_index(), cell_fine->active_fe_index())]]
            .n_coarse_cells++;
        });
      }

      for (const auto &fe_index_pair : fe_index_pairs)
        {
          transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first)
              .n_dofs_per_cell();
          transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second)
              .n_dofs_per_cell();

          transfer.schemes[fe_index_pair.second].degree_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first).degree;
          transfer.schemes[fe_index_pair.second].degree_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second).degree;
        }

      {
        transfer.partitioner_fine =
          std::make_shared<Utilities::MPI::Partitioner>(
            view.locally_owned_dofs(),
            view.locally_relevant_dofs(),
            dof_handler_fine.get_communicator());
        transfer.vec_fine.reinit(transfer.partitioner_fine);
      }
      {
        transfer.partitioner_coarse =
          create_coarse_partitioner(dof_handler_coarse,
                                    constraints_coarse,
                                    mg_level_coarse);
        transfer.vec_coarse.reinit(transfer.partitioner_coarse);
      }

      std::vector<unsigned int> n_dof_indices_fine(fe_index_pairs.size() + 1);
      std::vector<unsigned int> n_dof_indices_coarse(fe_index_pairs.size() + 1);
      std::vector<unsigned int> cell_no(fe_index_pairs.size() + 1, 0);

      {
        std::vector<std::vector<unsigned int>> lexicographic_numbering_fine(
          fe_index_pairs.size());
        std::vector<std::vector<unsigned int>> lexicographic_numbering_coarse(
          fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_coarse(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_coarse_lex(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_fine(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_fine_lex(fe_index_pairs.size());

        for (const auto &fe_index_pair : fe_index_pairs)
          {
            local_dof_indices_coarse[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse);
            local_dof_indices_coarse_lex[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse);
            local_dof_indices_fine[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine);
            local_dof_indices_fine_lex[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine);

            n_dof_indices_fine[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine *
              transfer.schemes[fe_index_pair.second].n_coarse_cells;
            n_dof_indices_coarse[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse *
              transfer.schemes[fe_index_pair.second].n_coarse_cells;
            cell_no[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_coarse_cells;

            const auto reference_cell =
              dof_handler_fine.get_fe(fe_index_pair.first.second)
                .reference_cell();

            Assert(reference_cell ==
                     dof_handler_coarse.get_fe(fe_index_pair.first.first)
                       .reference_cell(),
                   ExcNotImplemented());

            // ------------------- lexicographic_numbering  --------------------
            if (reference_cell == ReferenceCells::get_hypercube<dim>())
              {
                const Quadrature<1> dummy_quadrature(
                  std::vector<Point<1>>(1, Point<1>()));
                internal::MatrixFreeFunctions::ShapeInfo<
                  VectorizedArray<Number>>
                  shape_info;
                shape_info.reinit(dummy_quadrature,
                                  dof_handler_fine.get_fe(
                                    fe_index_pair.first.second),
                                  0);
                lexicographic_numbering_fine[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;

                shape_info.reinit(dummy_quadrature,
                                  dof_handler_coarse.get_fe(
                                    fe_index_pair.first.first),
                                  0);
                lexicographic_numbering_coarse[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;
              }
            else
              {
                const auto dummy_quadrature =
                  reference_cell.template get_gauss_type_quadrature<dim>(1);

                internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
                shape_info.reinit(dummy_quadrature,
                                  dof_handler_fine.get_fe(
                                    fe_index_pair.first.second),
                                  0);
                lexicographic_numbering_fine[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;

                shape_info.reinit(dummy_quadrature,
                                  dof_handler_coarse.get_fe(
                                    fe_index_pair.first.first),
                                  0);
                lexicographic_numbering_coarse[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;
              }
          }

        for (unsigned int i = 0; i < fe_index_pairs.size(); ++i)
          {
            n_dof_indices_fine[i + 1] += n_dof_indices_fine[i];
            n_dof_indices_coarse[i + 1] += n_dof_indices_coarse[i];
            cell_no[i + 1] += cell_no[i];
          }

        // ------------------------------ indices  -----------------------------

        transfer.constraint_info_coarse.reinit(
          dof_handler_coarse,
          cell_no.back(),
          constraints_coarse.n_constraints() > 0 &&
            use_fast_hanging_node_algorithm(dof_handler_coarse,
                                            mg_level_coarse));

        transfer.constraint_info_fine.reinit(cell_no.back());

        process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
          const auto fe_pair_no =
            fe_index_pairs[std::pair<unsigned int, unsigned int>(
              cell_coarse->active_fe_index(), cell_fine->active_fe_index())];

          // parent
          {
            transfer.constraint_info_coarse.read_dof_indices(
              cell_no[fe_pair_no],
              mg_level_coarse,
              cell_coarse,
              constraints_coarse,
              transfer.partitioner_coarse);
          }

          // child
          {
            cell_fine->get_dof_indices(local_dof_indices_fine[fe_pair_no]);

            for (unsigned int i = 0;
                 i < transfer.schemes[fe_pair_no].n_dofs_per_cell_fine;
                 ++i)
              local_dof_indices_fine_lex[fe_pair_no][i] = local_dof_indices_fine
                [fe_pair_no][lexicographic_numbering_fine[fe_pair_no][i]];

            transfer.constraint_info_fine.read_dof_indices(
              cell_no[fe_pair_no],
              local_dof_indices_fine_lex[fe_pair_no],
              transfer.partitioner_fine);
          }

          // move pointers
          {
            cell_no[fe_pair_no]++;
          }
        });

        transfer.constraint_info_coarse.finalize();
        transfer.constraint_info_fine.finalize();
      }

      // ------------------------- prolongation matrix -------------------------
      for (auto const &fe_index_pair_ : fe_index_pairs)
        {
          const auto &fe_index_pair = fe_index_pair_.first;
          const auto &fe_index_no   = fe_index_pair_.second;

          AssertDimension(
            dof_handler_fine.get_fe(fe_index_pair.second).n_base_elements(), 1);
          AssertDimension(
            dof_handler_coarse.get_fe(fe_index_pair.first).n_base_elements(),
            1);

          const auto reference_cell =
            dof_handler_fine.get_fe(fe_index_pair_.first.second)
              .reference_cell();

          Assert(reference_cell ==
                   dof_handler_coarse.get_fe(fe_index_pair_.first.first)
                     .reference_cell(),
                 ExcNotImplemented());

          if (reference_cell == ReferenceCells::get_hypercube<dim>() &&
              (dof_handler_coarse.get_fe(fe_index_pair.first) !=
               dof_handler_fine.get_fe(fe_index_pair.second)) &&
              (dof_handler_coarse.get_fe(fe_index_pair.first)
                   .n_dofs_per_cell() != 0 &&
               dof_handler_fine.get_fe(fe_index_pair.second)
                   .n_dofs_per_cell() != 0))
            {
              const auto fe_fine = create_1D_fe(
                dof_handler_fine.get_fe(fe_index_pair.second).base_element(0));

              std::vector<unsigned int> renumbering_fine(
                fe_fine->n_dofs_per_cell());
              {
                AssertIndexRange(fe_fine->n_dofs_per_vertex(), 2);
                renumbering_fine[0] = 0;
                for (unsigned int i = 0; i < fe_fine->dofs_per_line; ++i)
                  renumbering_fine[i + fe_fine->n_dofs_per_vertex()] =
                    GeometryInfo<1>::vertices_per_cell *
                      fe_fine->n_dofs_per_vertex() +
                    i;
                if (fe_fine->n_dofs_per_vertex() > 0)
                  renumbering_fine[fe_fine->n_dofs_per_cell() -
                                   fe_fine->n_dofs_per_vertex()] =
                    fe_fine->n_dofs_per_vertex();
              }

              const auto fe_coarse = create_1D_fe(
                dof_handler_coarse.get_fe(fe_index_pair.first).base_element(0));

              std::vector<unsigned int> renumbering_coarse(
                fe_coarse->n_dofs_per_cell());
              {
                AssertIndexRange(fe_coarse->n_dofs_per_vertex(), 2);
                renumbering_coarse[0] = 0;
                for (unsigned int i = 0; i < fe_coarse->dofs_per_line; ++i)
                  renumbering_coarse[i + fe_coarse->n_dofs_per_vertex()] =
                    GeometryInfo<1>::vertices_per_cell *
                      fe_coarse->n_dofs_per_vertex() +
                    i;
                if (fe_coarse->n_dofs_per_vertex() > 0)
                  renumbering_coarse[fe_coarse->n_dofs_per_cell() -
                                     fe_coarse->n_dofs_per_vertex()] =
                    fe_coarse->n_dofs_per_vertex();
              }

              {
                FullMatrix<double> matrix(fe_fine->n_dofs_per_cell(),
                                          fe_coarse->n_dofs_per_cell());
                FETools::get_projection_matrix(*fe_coarse, *fe_fine, matrix);
                transfer.schemes[fe_index_no].prolongation_matrix_1d.resize(
                  fe_fine->n_dofs_per_cell() * fe_coarse->n_dofs_per_cell());

                for (unsigned int i = 0, k = 0;
                     i < fe_coarse->n_dofs_per_cell();
                     ++i)
                  for (unsigned int j = 0; j < fe_fine->n_dofs_per_cell();
                       ++j, ++k)
                    transfer.schemes[fe_index_no].prolongation_matrix_1d[k] =
                      matrix(renumbering_fine[j], renumbering_coarse[i]);
              }

              {
                FullMatrix<double> matrix(fe_coarse->n_dofs_per_cell(),
                                          fe_fine->n_dofs_per_cell());
                FETools::get_projection_matrix(*fe_fine, *fe_coarse, matrix);
                transfer.schemes[fe_index_no].restriction_matrix_1d.resize(
                  fe_fine->n_dofs_per_cell() * fe_coarse->n_dofs_per_cell());

                for (unsigned int i = 0, k = 0;
                     i < fe_coarse->n_dofs_per_cell();
                     ++i)
                  for (unsigned int j = 0; j < fe_fine->n_dofs_per_cell();
                       ++j, ++k)
                    transfer.schemes[fe_index_no].restriction_matrix_1d[k] =
                      matrix(renumbering_coarse[i], renumbering_fine[j]);
              }
            }
          else if (reference_cell != ReferenceCells::get_hypercube<dim>() &&
                   (dof_handler_coarse.get_fe(fe_index_pair.first) !=
                    dof_handler_fine.get_fe(fe_index_pair.second)) &&
                   (dof_handler_coarse.get_fe(fe_index_pair.first)
                        .n_dofs_per_cell() != 0 &&
                    dof_handler_fine.get_fe(fe_index_pair.second)
                        .n_dofs_per_cell() != 0))
            {
              const auto &fe_fine =
                dof_handler_fine.get_fe(fe_index_pair.second).base_element(0);

              const auto &fe_coarse =
                dof_handler_coarse.get_fe(fe_index_pair.first).base_element(0);

              {
                FullMatrix<double> matrix(fe_fine.n_dofs_per_cell(),
                                          fe_coarse.n_dofs_per_cell());
                FETools::get_projection_matrix(fe_coarse, fe_fine, matrix);
                transfer.schemes[fe_index_no].prolongation_matrix.resize(
                  fe_fine.n_dofs_per_cell() * fe_coarse.n_dofs_per_cell());

                for (unsigned int i = 0, k = 0; i < fe_coarse.n_dofs_per_cell();
                     ++i)
                  for (unsigned int j = 0; j < fe_fine.n_dofs_per_cell();
                       ++j, ++k)
                    transfer.schemes[fe_index_no].prolongation_matrix[k] =
                      matrix(j, i);
              }

              {
                FullMatrix<double> matrix(fe_coarse.n_dofs_per_cell(),
                                          fe_fine.n_dofs_per_cell());
                FETools::get_projection_matrix(fe_fine, fe_coarse, matrix);
                transfer.schemes[fe_index_no].restriction_matrix.resize(
                  fe_fine.n_dofs_per_cell() * fe_coarse.n_dofs_per_cell());

                for (unsigned int i = 0, k = 0; i < fe_coarse.n_dofs_per_cell();
                     ++i)
                  for (unsigned int j = 0; j < fe_fine.n_dofs_per_cell();
                       ++j, ++k)
                    transfer.schemes[fe_index_no].restriction_matrix[k] =
                      matrix(i, j);
              }
            }
        }

      // ------------------------------- weights -------------------------------
      if (transfer.fine_element_is_continuous)
        {
          // compute weights globally
          LinearAlgebra::distributed::Vector<Number> weight_vector;
          compute_weights(dof_handler_fine,
                          mg_level_fine,
                          constraints_fine,
                          transfer,
                          weight_vector);

          // ... and store them cell-wise a DG format
          transfer.weights.resize(n_dof_indices_fine.back());

          std::vector<unsigned int *> level_dof_indices_fine_(
            fe_index_pairs.size());
          std::vector<Number *> weights_(fe_index_pairs.size());

          for (unsigned int i = 0; i < fe_index_pairs.size(); ++i)
            {
              level_dof_indices_fine_[i] =
                transfer.constraint_info_fine.dof_indices.data() +
                n_dof_indices_fine[i];
              weights_[i] = transfer.weights.data() + n_dof_indices_fine[i];
            }

          process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
            const auto fe_pair_no =
              fe_index_pairs[std::pair<unsigned int, unsigned int>(
                cell_coarse->active_fe_index(), cell_fine->active_fe_index())];

            for (unsigned int i = 0;
                 i < transfer.schemes[fe_pair_no].n_dofs_per_cell_fine;
                 i++)
              weights_[fe_pair_no][i] = weight_vector.local_element(
                level_dof_indices_fine_[fe_pair_no][i]);

            level_dof_indices_fine_[fe_pair_no] +=
              transfer.schemes[fe_pair_no].n_dofs_per_cell_fine;
            weights_[fe_pair_no] +=
              transfer.schemes[fe_pair_no].n_dofs_per_cell_fine;
          });

          if (is_feq)
            compress_weights(transfer);
        }
    }
  };



  template <typename Number>
  struct SimpleVectorDataExchange
  {
    SimpleVectorDataExchange(
      const std::shared_ptr<const Utilities::MPI::Partitioner>
        &                    embedded_partitioner,
      AlignedVector<Number> &buffer)
      : embedded_partitioner(embedded_partitioner)
      , buffer(buffer)
    {}

    template <typename VectorType>
    void
    update_ghost_values(const VectorType &vec) const
    {
      update_ghost_values_start(vec);
      update_ghost_values_finish(vec);
    }

    template <typename VectorType>
    void
    update_ghost_values_start(const VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      buffer.resize_fast(embedded_partitioner->n_import_indices());

      embedded_partitioner
        ->template export_to_ghosted_array_start<Number, MemorySpace::Host>(
          0,
          dealii::ArrayView<const Number>(
            vec.begin(), embedded_partitioner->locally_owned_size()),
          dealii::ArrayView<Number>(buffer.begin(), buffer.size()),
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    update_ghost_values_finish(const VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      embedded_partitioner
        ->template export_to_ghosted_array_finish<Number, MemorySpace::Host>(
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);

      vec.set_ghost_state(true);
#endif
    }

    template <typename VectorType>
    void
    compress(VectorType &vec) const
    {
      compress_start(vec);
      compress_finish(vec);
    }

    template <typename VectorType>
    void
    compress_start(VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      buffer.resize_fast(embedded_partitioner->n_import_indices());

      embedded_partitioner
        ->template import_from_ghosted_array_start<Number, MemorySpace::Host>(
          VectorOperation::add,
          0,
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          dealii::ArrayView<Number>(buffer.begin(), buffer.size()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    compress_finish(VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      embedded_partitioner
        ->template import_from_ghosted_array_finish<Number, MemorySpace::Host>(
          VectorOperation::add,
          dealii::ArrayView<const Number>(buffer.begin(), buffer.size()),
          dealii::ArrayView<Number>(vec.begin(),
                                    embedded_partitioner->locally_owned_size()),
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    zero_out_ghost_values(const VectorType &vec) const
    {
      const auto &vector_partitioner = vec.get_partitioner();

      ArrayView<Number> ghost_array(
        const_cast<LinearAlgebra::distributed::Vector<Number> &>(vec).begin() +
          vector_partitioner->locally_owned_size(),
        vector_partitioner->n_ghost_indices());

      for (const auto &my_ghosts :
           embedded_partitioner->ghost_indices_within_larger_ghost_set())
        for (unsigned int j = my_ghosts.first; j < my_ghosts.second; ++j)
          ghost_array[j] = 0.;

      vec.set_ghost_state(false);
    }

  private:
    const std::shared_ptr<const Utilities::MPI::Partitioner>
                                     embedded_partitioner;
    dealii::AlignedVector<Number> &  buffer;
    mutable std::vector<MPI_Request> requests;
  };

} // namespace internal



namespace MGTransferGlobalCoarseningTools
{
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &fine_triangulation_in)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

    coarse_grid_triangulations.back().reset(&fine_triangulation_in, [](auto *) {
      // empty deleter, since fine_triangulation_in is an external field
      // and its destructor is called somewhere else
    });

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    Assert(
      (dynamic_cast<
         const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
         &fine_triangulation_in) == nullptr),
      ExcMessage(
        "Triangulations of type parallel::fullydistributed::Triangulation are "
        "not supported by this function!"));

    const auto create_new_empty_triangulation =
      [&]() -> std::shared_ptr<Triangulation<dim, spacedim>> {
#ifdef DEAL_II_WITH_P4EST
      if (const auto fine_triangulation = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(
            &fine_triangulation_in))
        return std::make_shared<
          parallel::distributed::Triangulation<dim, spacedim>>(
          fine_triangulation->get_communicator(),
          fine_triangulation->get_mesh_smoothing());
      else
#endif
#ifdef DEAL_II_WITH_MPI
        if (const auto fine_triangulation = dynamic_cast<
              const parallel::shared::Triangulation<dim, spacedim> *>(
              &fine_triangulation_in))
        return std::make_shared<parallel::shared::Triangulation<dim, spacedim>>(
          fine_triangulation->get_communicator(),
          fine_triangulation->get_mesh_smoothing(),
          fine_triangulation->with_artificial_cells());
      else
#endif
        return std::make_shared<Triangulation<dim, spacedim>>(
          fine_triangulation_in.get_mesh_smoothing());
    };

    const unsigned int max_level = fine_triangulation_in.n_global_levels() - 1;

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // copy triangulation
        auto new_tria = create_new_empty_triangulation();

        new_tria->copy_triangulation(*coarse_grid_triangulations[l]);

        // coarsen mesh
        new_tria->coarsen_global();

        // save mesh
        coarse_grid_triangulations[l - 1] = new_tria;
      }

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim> &                        fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool keep_fine_triangulation,
    const bool repartition_fine_triangulation)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

#ifndef DEAL_II_WITH_P4EST
    Assert(false, ExcNotImplemented());
    (void)policy;
    (void)keep_fine_triangulation;
    (void)repartition_fine_triangulation;
#else
    const auto fine_triangulation =
      dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
        &fine_triangulation_in);

    Assert(fine_triangulation, ExcNotImplemented());

    const auto comm = fine_triangulation->get_communicator();

    parallel::distributed::Triangulation<dim, spacedim> temp_triangulation(
      comm, fine_triangulation->get_mesh_smoothing());

    if (keep_fine_triangulation == true &&
        repartition_fine_triangulation == false)
      {
        coarse_grid_triangulations.back().reset(&fine_triangulation_in,
                                                [](auto *) {
                                                  // empty deleter, since
                                                  // fine_triangulation_in is an
                                                  // external field and its
                                                  // destructor is called
                                                  // somewhere else
                                                });
      }
    else
      {
        // create triangulation description
        const auto construction_data =
          repartition_fine_triangulation ?
            TriangulationDescription::Utilities::
              create_description_from_triangulation(
                *fine_triangulation, policy.partition(*fine_triangulation)) :
            TriangulationDescription::Utilities::
              create_description_from_triangulation(*fine_triangulation, comm);

        // create new triangulation
        const auto new_fine_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            new_fine_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        new_fine_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations.back() = new_fine_triangulation;
      }

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    if (keep_fine_triangulation == true)
      temp_triangulation.copy_triangulation(*fine_triangulation);

    auto *temp_triangulation_ptr =
      keep_fine_triangulation ? &temp_triangulation : fine_triangulation;

    const unsigned int max_level = fine_triangulation->n_global_levels() - 1;

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // coarsen mesh
        temp_triangulation_ptr->coarsen_global();

        // create triangulation description
        const auto construction_data = TriangulationDescription::Utilities::
          create_description_from_triangulation(
            *temp_triangulation_ptr, policy.partition(*temp_triangulation_ptr));

        // create new triangulation
        const auto level_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            level_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        level_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations[l - 1] = level_triangulation;
      }
#endif

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &                  fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation)
  {
    // remove const and convert it to flag
    return create_geometric_coarsening_sequence(
      const_cast<Triangulation<dim, spacedim> &>(fine_triangulation_in),
      policy,
      true,
      repartition_fine_triangulation);
  }

} // namespace MGTransferGlobalCoarseningTools



template <typename Number>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::
  prolongate_and_add(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const bool  use_dst_inplace = this->vec_fine.size() == 0;
  auto *const vec_fine_ptr    = use_dst_inplace ? &dst : &this->vec_fine;
  Assert(vec_fine_ptr->get_partitioner().get() == this->partitioner_fine.get(),
         ExcInternalError());

  const bool        use_src_inplace = this->vec_coarse.size() == 0;
  const auto *const vec_coarse_ptr = use_src_inplace ? &src : &this->vec_coarse;
  Assert(vec_coarse_ptr->get_partitioner().get() ==
           this->partitioner_coarse.get(),
         ExcInternalError());

  if (use_src_inplace == false)
    this->vec_coarse.copy_locally_owned_data_from(src);

  this->update_ghost_values(*vec_coarse_ptr);

  if (use_dst_inplace == false)
    *vec_fine_ptr = Number(0.);

  this->prolongate_and_add_internal(*vec_fine_ptr, *vec_coarse_ptr);

  if (fine_element_is_continuous || use_dst_inplace == false)
    this->compress(*vec_fine_ptr, VectorOperation::add);

  if (use_dst_inplace == false)
    dst += this->vec_fine;

  if (use_src_inplace)
    this->zero_out_ghost_values(*vec_coarse_ptr);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int n_lanes = VectorizedArrayType::size();

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  const Number *             weights            = nullptr;
  const VectorizedArrayType *weights_compressed = nullptr;

  if (this->fine_element_is_continuous)
    {
      weights            = this->weights.data();
      weights_compressed = this->weights_compressed.data();
    }

  unsigned int cell_counter = 0;

  for (const auto &scheme : schemes)
    {
      if (scheme.n_coarse_cells == 0)
        continue;

      const bool needs_interpolation =
        (scheme.prolongation_matrix.empty() &&
         scheme.prolongation_matrix_1d.empty()) == false;

      evaluation_data_fine.clear();
      evaluation_data_coarse.clear();

      const unsigned int max_n_dofs_per_cell =
        std::max(scheme.n_dofs_per_cell_fine, scheme.n_dofs_per_cell_coarse);
      evaluation_data_fine.resize(max_n_dofs_per_cell);
      evaluation_data_coarse.resize(max_n_dofs_per_cell);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.n_dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.n_dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from src vector (similar to FEEvaluation::read_dof_values())
          internal::VectorReader<Number, VectorizedArrayType> reader;
          constraint_info_coarse.read_write_operation(
            reader,
            src,
            evaluation_data_coarse.data(),
            cell_counter,
            n_lanes_filled,
            scheme.n_dofs_per_cell_coarse,
            true);
          constraint_info_coarse.apply_hanging_node_constraints(
            cell_counter, n_lanes_filled, false, evaluation_data_coarse);

          // ---------------------------- coarse -------------------------------
          if (needs_interpolation)
            for (int c = n_components - 1; c >= 0; --c)
              {
                CellProlongator<dim, VectorizedArrayType> cell_prolongator(
                  scheme.prolongation_matrix,
                  scheme.prolongation_matrix_1d,
                  evaluation_data_coarse.begin() + c * n_scalar_dofs_coarse,
                  evaluation_data_fine.begin() + c * n_scalar_dofs_fine);

                if (scheme.prolongation_matrix_1d.size() > 0)
                  cell_transfer.run(cell_prolongator);
                else
                  cell_prolongator.run_full(n_scalar_dofs_fine,
                                            n_scalar_dofs_coarse);
              }
          else
            evaluation_data_fine = evaluation_data_coarse; // TODO
          // ------------------------------ fine -------------------------------

          // weight
          if (this->fine_element_is_continuous &&
              this->weights_compressed.size() > 0)
            {
              internal::
                weight_fe_q_dofs_by_entity<dim, -1, VectorizedArrayType>(
                  weights_compressed,
                  n_components,
                  scheme.degree_fine + 1,
                  evaluation_data_fine.begin());
              weights_compressed += Utilities::pow(3, dim);
            }
          else if (this->fine_element_is_continuous)
            {
              for (unsigned int v = 0; v < n_lanes_filled; ++v)
                {
                  for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine; ++i)
                    evaluation_data_fine[i][v] *= weights[i];
                  weights += scheme.n_dofs_per_cell_fine;
                }
            }

          // add into dst vector
          internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
            writer;
          constraint_info_fine.read_write_operation(writer,
                                                    dst,
                                                    evaluation_data_fine.data(),
                                                    cell_counter,
                                                    n_lanes_filled,
                                                    scheme.n_dofs_per_cell_fine,
                                                    false);

          cell_counter += n_lanes_filled;
        }
    }
}



template <typename Number>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::
  restrict_and_add(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const bool        use_src_inplace = this->vec_fine.size() == 0;
  const auto *const vec_fine_ptr    = use_src_inplace ? &src : &this->vec_fine;
  Assert(vec_fine_ptr->get_partitioner().get() == this->partitioner_fine.get(),
         ExcInternalError());

  const bool  use_dst_inplace = this->vec_coarse.size() == 0;
  auto *const vec_coarse_ptr  = use_dst_inplace ? &dst : &this->vec_coarse;
  Assert(vec_coarse_ptr->get_partitioner().get() ==
           this->partitioner_coarse.get(),
         ExcInternalError());

  if (use_src_inplace == false)
    this->vec_fine.copy_locally_owned_data_from(src);

  if (fine_element_is_continuous || use_src_inplace == false)
    this->update_ghost_values(*vec_fine_ptr);

  if (use_dst_inplace == false)
    *vec_coarse_ptr = 0.0;

  this->zero_out_ghost_values(
    *vec_coarse_ptr); // since we might add into the
                      // ghost values and call compress

  this->restrict_and_add_internal(*vec_coarse_ptr, *vec_fine_ptr);

  // clean up related to update_ghost_values()
  if (fine_element_is_continuous == false && use_src_inplace == false)
    this->zero_out_ghost_values(*vec_fine_ptr); // internal vector (DG)
  else if (fine_element_is_continuous && use_src_inplace == false)
    vec_fine_ptr->set_ghost_state(false); // internal vector (CG)
  else if (fine_element_is_continuous)
    this->zero_out_ghost_values(*vec_fine_ptr); // external vector

  this->compress(*vec_coarse_ptr, VectorOperation::add);

  if (use_dst_inplace == false)
    dst += this->vec_coarse;
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int n_lanes = VectorizedArrayType::size();

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  const Number *             weights            = nullptr;
  const VectorizedArrayType *weights_compressed = nullptr;

  if (this->fine_element_is_continuous)
    {
      weights            = this->weights.data();
      weights_compressed = this->weights_compressed.data();
    }

  unsigned int cell_counter = 0;

  for (const auto &scheme : schemes)
    {
      if (scheme.n_coarse_cells == 0)
        continue;

      const bool needs_interpolation =
        (scheme.prolongation_matrix.empty() &&
         scheme.prolongation_matrix_1d.empty()) == false;

      evaluation_data_fine.clear();
      evaluation_data_coarse.clear();

      const unsigned int max_n_dofs_per_cell =
        std::max(scheme.n_dofs_per_cell_fine, scheme.n_dofs_per_cell_coarse);
      evaluation_data_fine.resize(max_n_dofs_per_cell);
      evaluation_data_coarse.resize(max_n_dofs_per_cell);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.n_dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.n_dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from source vector
          internal::VectorReader<Number, VectorizedArrayType> reader;
          constraint_info_fine.read_write_operation(reader,
                                                    src,
                                                    evaluation_data_fine.data(),
                                                    cell_counter,
                                                    n_lanes_filled,
                                                    scheme.n_dofs_per_cell_fine,
                                                    false);

          // weight
          if (this->fine_element_is_continuous &&
              this->weights_compressed.size() > 0)
            {
              internal::
                weight_fe_q_dofs_by_entity<dim, -1, VectorizedArrayType>(
                  weights_compressed,
                  n_components,
                  scheme.degree_fine + 1,
                  evaluation_data_fine.begin());
              weights_compressed += Utilities::pow(3, dim);
            }
          else if (this->fine_element_is_continuous)
            {
              for (unsigned int v = 0; v < n_lanes_filled; ++v)
                {
                  for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine; ++i)
                    evaluation_data_fine[i][v] *= weights[i];
                  weights += scheme.n_dofs_per_cell_fine;
                }
            }

          // ------------------------------ fine -------------------------------
          if (needs_interpolation)
            for (int c = n_components - 1; c >= 0; --c)
              {
                CellRestrictor<dim, VectorizedArrayType> cell_restrictor(
                  scheme.prolongation_matrix,
                  scheme.prolongation_matrix_1d,
                  evaluation_data_fine.begin() + c * n_scalar_dofs_fine,
                  evaluation_data_coarse.begin() + c * n_scalar_dofs_coarse);

                if (scheme.prolongation_matrix_1d.size() > 0)
                  cell_transfer.run(cell_restrictor);
                else
                  cell_restrictor.run_full(n_scalar_dofs_fine,
                                           n_scalar_dofs_coarse);
              }
          else
            evaluation_data_coarse = evaluation_data_fine; // TODO
          // ----------------------------- coarse ------------------------------

          // write into dst vector (similar to
          // FEEvaluation::distribute_global_to_local())
          internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
            writer;
          constraint_info_coarse.apply_hanging_node_constraints(
            cell_counter, n_lanes_filled, true, evaluation_data_coarse);
          constraint_info_coarse.read_write_operation(
            writer,
            dst,
            evaluation_data_coarse.data(),
            cell_counter,
            n_lanes_filled,
            scheme.n_dofs_per_cell_coarse,
            true);

          cell_counter += n_lanes_filled;
        }
    }
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  interpolate(LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int n_lanes = VectorizedArrayType::size();

  const bool        use_src_inplace = this->vec_fine.size() == 0;
  const auto *const vec_fine_ptr    = use_src_inplace ? &src : &this->vec_fine;
  Assert(vec_fine_ptr->get_partitioner().get() == this->partitioner_fine.get(),
         ExcInternalError());

  const bool  use_dst_inplace = this->vec_coarse.size() == 0;
  auto *const vec_coarse_ptr  = use_dst_inplace ? &dst : &this->vec_coarse;
  Assert(vec_coarse_ptr->get_partitioner().get() ==
           this->partitioner_coarse.get(),
         ExcInternalError());

  if (use_src_inplace == false)
    this->vec_fine.copy_locally_owned_data_from(src);

  if (this->fine_element_is_continuous || use_src_inplace == false)
    this->update_ghost_values(*vec_fine_ptr);

  *vec_coarse_ptr = 0.0;

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  unsigned int cell_counter = 0;

  for (const auto &scheme : schemes)
    {
      if (scheme.n_coarse_cells == 0)
        continue;

      if (scheme.n_dofs_per_cell_fine == 0 ||
          scheme.n_dofs_per_cell_coarse == 0)
        {
          cell_counter += scheme.n_coarse_cells;

          continue;
        }

      const bool needs_interpolation =
        (scheme.prolongation_matrix.empty() &&
         scheme.prolongation_matrix_1d.empty()) == false;

      // general case -> local restriction is needed
      evaluation_data_fine.resize(scheme.n_dofs_per_cell_fine);
      evaluation_data_coarse.resize(scheme.n_dofs_per_cell_fine);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.n_dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.n_dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from source vector
          internal::VectorReader<Number, VectorizedArrayType> reader;
          constraint_info_fine.read_write_operation(reader,
                                                    *vec_fine_ptr,
                                                    evaluation_data_fine.data(),
                                                    cell_counter,
                                                    n_lanes_filled,
                                                    scheme.n_dofs_per_cell_fine,
                                                    false);

          // ------------------------------ fine -------------------------------
          if (needs_interpolation)
            for (int c = n_components - 1; c >= 0; --c)
              {
                CellRestrictor<dim, VectorizedArrayType> cell_restrictor(
                  scheme.restriction_matrix,
                  scheme.restriction_matrix_1d,
                  evaluation_data_fine.begin() + c * n_scalar_dofs_fine,
                  evaluation_data_coarse.begin() + c * n_scalar_dofs_coarse);

                if (scheme.restriction_matrix_1d.size() > 0)
                  cell_transfer.run(cell_restrictor);
                else
                  cell_restrictor.run_full(n_scalar_dofs_fine,
                                           n_scalar_dofs_coarse);
              }
          else
            evaluation_data_coarse = evaluation_data_fine; // TODO
          // ----------------------------- coarse ------------------------------

          // write into dst vector (similar to
          // FEEvaluation::set_dof_values_plain())
          internal::VectorSetter<Number, VectorizedArrayType> writer;
          constraint_info_coarse.read_write_operation(
            writer,
            *vec_coarse_ptr,
            evaluation_data_coarse.data(),
            cell_counter,
            n_lanes_filled,
            scheme.n_dofs_per_cell_coarse,
            false);

          cell_counter += n_lanes_filled;
        }
    }

  // clean up related to update_ghost_values()
  if (use_src_inplace == false)
    vec_fine_ptr->set_ghost_state(false); // internal vector
  else if (this->fine_element_is_continuous)
    this->zero_out_ghost_values(*vec_fine_ptr); // external vector

  if (use_dst_inplace == false)
    dst.copy_locally_owned_data_from(this->vec_coarse);
}


namespace internal
{
  namespace
  {
    bool
    is_partitioner_contained(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::shared_ptr<const Utilities::MPI::Partitioner>
        &external_partitioner)
    {
      // no external partitioner has been given
      if (external_partitioner.get() == nullptr)
        return false;

      // check if locally owned ranges are the same
      if (external_partitioner->size() != partitioner->size())
        return false;

      if (external_partitioner->locally_owned_range() !=
          partitioner->locally_owned_range())
        return false;

      const int ghosts_locally_contained =
        ((external_partitioner->ghost_indices() &
          partitioner->ghost_indices()) == partitioner->ghost_indices()) ?
          1 :
          0;

      // check if ghost values are contained in external partititioner
      return Utilities::MPI::min(ghosts_locally_contained,
                                 partitioner->get_mpi_communicator()) == 1;
    }

    std::shared_ptr<Utilities::MPI::Partitioner>
    create_embedded_partitioner(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::shared_ptr<const Utilities::MPI::Partitioner>
        &larger_partitioner)
    {
      auto embedded_partitioner = std::make_shared<Utilities::MPI::Partitioner>(
        larger_partitioner->locally_owned_range(),
        larger_partitioner->get_mpi_communicator());

      embedded_partitioner->set_ghost_indices(
        partitioner->ghost_indices(), larger_partitioner->ghost_indices());

      return embedded_partitioner;
    }
  } // namespace
} // namespace internal

template <typename Number>
template <int dim, std::size_t width>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::
  internal_enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_fine,
    internal::MatrixFreeFunctions::ConstraintInfo<
      dim,
      VectorizedArray<Number, width>> &constraint_info_coarse,
    std::vector<unsigned int> &        dof_indices_fine)
{
  if (this->partitioner_coarse->is_globally_compatible(
        *external_partitioner_coarse))
    {
      this->vec_coarse.reinit(0);
      this->partitioner_coarse = external_partitioner_coarse;
    }
  else if (internal::is_partitioner_contained(this->partitioner_coarse,
                                              external_partitioner_coarse))
    {
      this->vec_coarse.reinit(0);

      for (auto &i : constraint_info_coarse.dof_indices)
        i = external_partitioner_coarse->global_to_local(
          this->partitioner_coarse->local_to_global(i));

      for (auto &i : constraint_info_coarse.plain_dof_indices)
        i = external_partitioner_coarse->global_to_local(
          this->partitioner_coarse->local_to_global(i));

      this->partitioner_coarse_embedded =
        internal::create_embedded_partitioner(this->partitioner_coarse,
                                              external_partitioner_coarse);

      this->partitioner_coarse = external_partitioner_coarse;
    }

  if (this->partitioner_fine->is_globally_compatible(
        *external_partitioner_fine))
    {
      this->vec_fine.reinit(0);
      this->partitioner_fine = external_partitioner_fine;
    }
  else if (internal::is_partitioner_contained(this->partitioner_fine,
                                              external_partitioner_fine))
    {
      this->vec_fine.reinit(0);

      for (auto &i : dof_indices_fine)
        i = external_partitioner_fine->global_to_local(
          this->partitioner_fine->local_to_global(i));

      this->partitioner_fine_embedded =
        internal::create_embedded_partitioner(this->partitioner_fine,
                                              external_partitioner_fine);

      this->partitioner_fine = external_partitioner_fine;
    }
}

template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_fine)
{
  this->internal_enable_inplace_operations_if_possible(
    external_partitioner_coarse,
    external_partitioner_fine,
    constraint_info_coarse,
    constraint_info_fine.dof_indices);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  reinit_geometric_transfer(const DoFHandler<dim> &          dof_handler_fine,
                            const DoFHandler<dim> &          dof_handler_coarse,
                            const AffineConstraints<Number> &constraints_fine,
                            const AffineConstraints<Number> &constraints_coarse,
                            const unsigned int               mg_level_fine,
                            const unsigned int               mg_level_coarse)
{
  internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraints_fine,
    constraints_coarse,
    mg_level_fine,
    mg_level_coarse,
    *this);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  reinit_polynomial_transfer(
    const DoFHandler<dim> &          dof_handler_fine,
    const DoFHandler<dim> &          dof_handler_coarse,
    const AffineConstraints<Number> &constraints_fine,
    const AffineConstraints<Number> &constraints_coarse,
    const unsigned int               mg_level_fine,
    const unsigned int               mg_level_coarse)
{
  internal::MGTwoLevelTransferImplementation::reinit_polynomial_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraints_fine,
    constraints_coarse,
    mg_level_fine,
    mg_level_coarse,
    *this);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::reinit(
  const DoFHandler<dim> &          dof_handler_fine,
  const DoFHandler<dim> &          dof_handler_coarse,
  const AffineConstraints<Number> &constraints_fine,
  const AffineConstraints<Number> &constraints_coarse,
  const unsigned int               mg_level_fine,
  const unsigned int               mg_level_coarse)
{
  // determine if polynomial transfer can be performed via the following two
  // criteria:
  // 1) multigrid levels can be only used with polynomial transfer
  bool do_polynomial_transfer =
    (mg_level_fine != numbers::invalid_unsigned_int) ||
    (mg_level_coarse != numbers::invalid_unsigned_int);

  // 2) the meshes are identical
  if (do_polynomial_transfer == false)
    {
      const internal::CellIDTranslator<dim> cell_id_translator(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_fine.get_triangulation().n_global_levels());

      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_coarse.get_triangulation().n_global_coarse_cells());
      AssertIndexRange(dof_handler_coarse.get_triangulation().n_global_levels(),
                       dof_handler_fine.get_triangulation().n_global_levels() +
                         1);

      IndexSet is_locally_owned_fine(cell_id_translator.size());
      IndexSet is_locally_owned_coarse(cell_id_translator.size());

      for (const auto &cell : dof_handler_fine.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        is_locally_owned_fine.add_index(cell_id_translator.translate(cell));

      for (const auto &cell : dof_handler_coarse.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        is_locally_owned_coarse.add_index(cell_id_translator.translate(cell));

      const MPI_Comm communicator = dof_handler_fine.get_communicator();

      std::vector<unsigned int> owning_ranks(
        is_locally_owned_coarse.n_elements());

      Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
        process(is_locally_owned_fine,
                is_locally_owned_coarse,
                communicator,
                owning_ranks,
                false);

      Utilities::MPI::ConsensusAlgorithms::Selector<
        std::vector<
          std::pair<types::global_cell_index, types::global_cell_index>>,
        std::vector<unsigned int>>
        consensus_algorithm;
      consensus_algorithm.run(process, communicator);

      bool all_cells_found = true;

      for (unsigned i = 0; i < is_locally_owned_coarse.n_elements(); ++i)
        all_cells_found &= (owning_ranks[i] != numbers::invalid_unsigned_int);

      do_polynomial_transfer =
        Utilities::MPI::min(static_cast<unsigned int>(all_cells_found),
                            communicator) == 1;
    }

  if (do_polynomial_transfer)
    internal::MGTwoLevelTransferImplementation::reinit_polynomial_transfer(
      dof_handler_fine,
      dof_handler_coarse,
      constraints_fine,
      constraints_coarse,
      mg_level_fine,
      mg_level_coarse,
      *this);
  else
    internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
      dof_handler_fine,
      dof_handler_coarse,
      constraints_fine,
      constraints_coarse,
      mg_level_fine,
      mg_level_coarse,
      *this);
}



template <int dim, typename Number>
bool
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  fast_polynomial_transfer_supported(const unsigned int fe_degree_fine,
                                     const unsigned int fe_degree_coarse)
{
  CellTransferFactory cell_transfer(fe_degree_fine, fe_degree_coarse);
  CellProlongatorTest cell_transfer_test;

  return cell_transfer.run(cell_transfer_test);
}



template <int dim, typename Number>
std::size_t
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  memory_consumption() const
{
  std::size_t size = 0;

  for (const auto &scheme : schemes)
    {
      size += scheme.prolongation_matrix.memory_consumption();
      size += scheme.prolongation_matrix_1d.memory_consumption();
      size += scheme.restriction_matrix.memory_consumption();
      size += scheme.restriction_matrix_1d.memory_consumption();
    }

  size += this->partitioner_fine->memory_consumption();
  size += this->partitioner_coarse->memory_consumption();
  size += this->vec_fine.memory_consumption();
  size += this->vec_coarse.memory_consumption();
  size += MemoryConsumption::memory_consumption(weights);
  size += constraint_info_coarse.memory_consumption();
  size += constraint_info_fine.memory_consumption();

  return size;
}



template <typename Number>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::
  update_ghost_values(
    const LinearAlgebra::distributed::Vector<Number> &vec) const
{
  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .update_ghost_values(vec);
  else if ((vec.get_partitioner().get() == this->partitioner_fine.get()) &&
           (this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .update_ghost_values(vec);
  else
    vec.update_ghost_values();
}



template <typename Number>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::compress(
  LinearAlgebra::distributed::Vector<Number> &vec,
  const VectorOperation::values               op) const
{
  Assert(op == VectorOperation::add, ExcNotImplemented());

  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .compress(vec);
  else if ((vec.get_partitioner().get() == this->partitioner_fine.get()) &&
           (this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .compress(vec);
  else
    vec.compress(op);
}



template <typename Number>
void
MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>::
  zero_out_ghost_values(
    const LinearAlgebra::distributed::Vector<Number> &vec) const
{
  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .zero_out_ghost_values(vec);
  else if ((vec.get_partitioner().get() == (this->partitioner_fine.get()) &&
            this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .zero_out_ghost_values(vec);
  else
    vec.zero_out_ghost_values();
}



template <int dim, typename Number>
MGTransferMF<dim, Number>::MGTransferMF(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::initialize_constraints(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->mg_constrained_dofs = &mg_constrained_dofs;
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::intitialize_internal_transfer(
  const DoFHandler<dim> &                      dof_handler,
  const SmartPointer<const MGConstrainedDoFs> &mg_constrained_dofs)
{
  const unsigned int min_level = 0;
  const unsigned int max_level =
    dof_handler.get_triangulation().n_global_levels() - 1;

  MGLevelObject<AffineConstraints<typename VectorType::value_type>> constraints(
    min_level, max_level);

  if (mg_constrained_dofs)
    for (unsigned int l = min_level; l <= max_level; ++l)
      mg_constrained_dofs->merge_constraints(
        constraints[l],
        l,
        /*add_boundary_indices*/ true,
        /*add_refinement_edge_indices*/ false,
        /*add_level_constraints*/ true,
        /*add_user_constraints*/ true);

  this->internal_transfer.resize(min_level, max_level);

  for (unsigned int l = min_level; l < max_level; ++l)
    internal_transfer[l + 1].reinit_geometric_transfer(
      dof_handler, dof_handler, constraints[l + 1], constraints[l], l + 1, l);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::build(
  const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    &external_partitioners)
{
  this->external_partitioners = external_partitioners;

  if (this->external_partitioners.size() > 0)
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();

      AssertDimension(this->external_partitioners.size(), transfer.n_levels());

      for (unsigned int l = min_level + 1; l <= max_level; ++l)
        transfer[l]->enable_inplace_operations_if_possible(
          this->external_partitioners[l - 1 - min_level],
          this->external_partitioners[l - min_level]);
    }
  else
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();

      for (unsigned int l = min_level + 1; l <= max_level; ++l)
        {
          if (l == min_level + 1)
            this->external_partitioners.push_back(
              transfer[l]->partitioner_coarse);

          this->external_partitioners.push_back(transfer[l]->partitioner_fine);
        }
    }

  this->perform_plain_copy            = true;
  this->perform_renumbered_plain_copy = false;
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::build(
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  if (initialize_dof_vector)
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();
      const unsigned int n_levels  = transfer.n_levels();

      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        external_partitioners(n_levels);

      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          LinearAlgebra::distributed::Vector<typename VectorType::value_type>
            vector;
          initialize_dof_vector(l, vector);
          external_partitioners[l - min_level] = vector.get_partitioner();
        }

      this->build(external_partitioners);
    }
  else
    {
      this->build();
    }
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::build(
  const DoFHandler<dim> &dof_handler,
  const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    &external_partitioners)
{
  this->intitialize_internal_transfer(dof_handler, this->mg_constrained_dofs);
  this->intitialize_transfer_references(internal_transfer);
  this->build(external_partitioners);
  this->fill_and_communicate_copy_indices(dof_handler);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::build(
  const DoFHandler<dim> &dof_handler,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  this->intitialize_internal_transfer(dof_handler, this->mg_constrained_dofs);
  this->intitialize_transfer_references(internal_transfer);
  this->build(initialize_dof_vector);
  this->fill_and_communicate_copy_indices(dof_handler);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::prolongate(const unsigned int to_level,
                                      VectorType &       dst,
                                      const VectorType & src) const
{
  dst = 0;
  prolongate_and_add(to_level, dst, src);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::prolongate_and_add(const unsigned int to_level,
                                              VectorType &       dst,
                                              const VectorType & src) const
{
  this->transfer[to_level]->prolongate_and_add(dst, src);
}



template <int dim, typename Number>
void
MGTransferMF<dim, Number>::restrict_and_add(const unsigned int from_level,
                                            VectorType &       dst,
                                            const VectorType & src) const
{
  this->transfer[from_level]->restrict_and_add(dst, src);
}



template <int dim, typename Number>
std::size_t
MGTransferMF<dim, Number>::memory_consumption() const
{
  std::size_t size = 0;

  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  for (unsigned int l = min_level + 1; l <= max_level; ++l)
    size += this->transfer[l]->memory_consumption();

  return size;
}



template <int dim, typename Number>
inline unsigned int
MGTransferMF<dim, Number>::min_level() const
{
  return transfer.min_level();
}



template <int dim, typename Number>
inline unsigned int
MGTransferMF<dim, Number>::max_level() const
{
  return transfer.max_level();
}


template <int dim, typename Number>
inline void
MGTransferMF<dim, Number>::clear()
{
  MGLevelGlobalTransfer<VectorType>::clear();

  internal_transfer.clear();
  transfer.clear();
  external_partitioners.clear();
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const MGTransferMF<dim, Number> &transfer_operator)
  : MGTransferBlockMatrixFreeBase<dim, Number, MGTransferMF<dim, Number>>(true)
{
  this->transfer_operators = {&transfer_operator};
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const MGConstrainedDoFs &mg_constrained_dofs)
  : MGTransferBlockMatrixFreeBase<dim, Number, MGTransferMF<dim, Number>>(true)
{
  initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const std::vector<MGConstrainedDoFs> &mg_constrained_dofs)
  : MGTransferBlockMatrixFreeBase<dim, Number, MGTransferMF<dim, Number>>(false)
{
  initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::initialize_constraints(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->transfer_operators_internal.clear();
  this->transfer_operators.clear();

  Assert(this->same_for_all,
         ExcMessage("This object was initialized with support for usage with "
                    "one DoFHandler for each block, but this method assumes "
                    "that the same DoFHandler is used for all the blocks!"));

  this->transfer_operators_internal.emplace_back(mg_constrained_dofs);
  this->transfer_operators = {&transfer_operators_internal.back()};
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::initialize_constraints(
  const std::vector<MGConstrainedDoFs> &mg_constrained_dofs)
{
  this->transfer_operators_internal.clear();
  this->transfer_operators.clear();

  Assert(!this->same_for_all,
         ExcMessage("This object was initialized with support for using "
                    "the same DoFHandler for all the blocks, but this "
                    "method assumes that there is a separate DoFHandler "
                    "for each block!"));

  for (const auto &dofs : mg_constrained_dofs)
    this->transfer_operators_internal.emplace_back(dofs);

  for (const auto &transfer : this->transfer_operators_internal)
    this->transfer_operators.emplace_back(&transfer);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::build(const DoFHandler<dim> &dof_handler)
{
  AssertDimension(transfer_operators.size(), 1);
  this->transfer_operators_internal[0].build(dof_handler);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::build(
  const std::vector<const DoFHandler<dim> *> &dof_handler)
{
  AssertDimension(transfer_operators.size(), dof_handler.size());
  AssertDimension(transfer_operators_internal.size(), dof_handler.size());

  for (unsigned int i = 0; i < dof_handler.size(); ++i)
    this->transfer_operators_internal[i].build(*dof_handler[i]);
}



template <int dim, typename Number>
const MGTransferMF<dim, Number> &
MGTransferBlockMF<dim, Number>::get_matrix_free_transfer(
  const unsigned int b) const
{
  AssertIndexRange(b, transfer_operators.size());
  return *transfer_operators[b];
}

namespace internal
{
  namespace
  {
    template <int dim, typename Number>
    std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>>
    fill_mapping_info(const Utilities::MPI::RemotePointEvaluation<dim> &rpe)
    {
      const auto &cell_data = rpe.get_cell_data();

      std::vector<typename Triangulation<dim>::active_cell_iterator>
                                           cell_iterators;
      std::vector<std::vector<Point<dim>>> unit_points_vector;

      for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
        {
          typename Triangulation<dim>::active_cell_iterator cell(
            &rpe.get_triangulation(),
            cell_data.cells[i].first,
            cell_data.cells[i].second);

          const ArrayView<const Point<dim>> unit_points(
            cell_data.reference_point_values.data() +
              cell_data.reference_point_ptrs[i],
            cell_data.reference_point_ptrs[i + 1] -
              cell_data.reference_point_ptrs[i]);

          cell_iterators.emplace_back(cell);
          unit_points_vector.emplace_back(unit_points.begin(),
                                          unit_points.end());
        }

      auto mapping_info =
        std::make_shared<NonMatching::MappingInfo<dim, dim, Number>>(
          rpe.get_mapping(), update_values);
      mapping_info->reinit_cells(cell_iterators, unit_points_vector);

      return mapping_info;
    }

    /**
     * This function provides information which DoF index is associated with
     * a support point.
     *
     * @param[in] dof_handler DoFHandler with @c FE_DGQ or @c FE_Q elements
     * providing DoF indices which are collected at support points.
     * @param[in] dof_handler_support_points DoFHandler with one component used
     * to determine support point indices (the underlying finite element is @c
     * FE_Q or @c FE_DGQ in case of polynomial degree 0).
     * @param[in] constraint AffineConstrains associated with @p dof_handler.
     *   Only unconstrained DoFs are considered
     * @return a tuple containing 0) local support point indices,
     *   1) pointers to global DoF indices, and 2) global DoF indices.
     */
    template <int dim, int spacedim, typename Number>
    std::tuple<std::vector<unsigned int>,
               std::vector<unsigned int>,
               std::vector<types::global_dof_index>>
    support_point_indices_to_dof_indices(
      const DoFHandler<dim, spacedim> &        dof_handler,
      const DoFHandler<dim, spacedim> &        dof_handler_support_points,
      const dealii::AffineConstraints<Number> &constraint)
    {
      // in case a FE_DGQ space of order 0 is provided, DoFs indices are always
      // uniquely assigned to support points (they are always defined in the
      // center of the element) and are never shared at vertices or faces.
      Assert((dynamic_cast<const FE_DGQ<dim, spacedim> *>(
                &dof_handler.get_fe()) != nullptr) ||
               (dynamic_cast<const FE_Q<dim, spacedim> *>(
                  &dof_handler.get_fe()) != nullptr),
             ExcMessage("Function expects FE_DGQ of FE_Q in dof_handler."));

      Assert(
        (dynamic_cast<const FE_Q<dim, spacedim> *>(
           &dof_handler_support_points.get_fe()) != nullptr) ||
          ((dynamic_cast<const FE_DGQ<dim, spacedim> *>(
              &dof_handler_support_points.get_fe()) != nullptr) &&
           dof_handler_support_points.get_fe().degree == 0),
        ExcMessage(
          "Function expects FE_DGQ&&degree==0 or FE_Q in dof_handler_support_points."));

      Assert(
        dof_handler_support_points.get_fe().n_components() == 1,
        ExcMessage(
          "dof_handler_support_points needs element with exactly one component."));

      Assert(&dof_handler.get_triangulation() ==
               &dof_handler_support_points.get_triangulation(),
             ExcMessage("DoFHandlers need the same underlying triangulation."));

      Assert(dof_handler.get_fe().degree ==
               dof_handler_support_points.get_fe().degree,
             ExcMessage("DoFHandlers need the same degree."));

      Assert(
        dof_handler.get_fe().n_dofs_per_cell() ==
          dof_handler_support_points.get_fe().n_dofs_per_cell(),
        ExcMessage(
          "For now multiple components are not considered and therefore DoFHandlers need the same degree."));

      const auto degree        = dof_handler.get_fe().degree;
      const auto dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
      const auto support_points_per_cell =
        dof_handler_support_points.get_fe().n_dofs_per_cell();

      std::vector<std::pair<unsigned int, types::global_dof_index>>
        support_point_dofs;
      support_point_dofs.reserve(dof_handler.n_locally_owned_dofs());

      // fill support_point_dofs
      {
        // Support points have a hierarchic numbering, L2 DoFs have
        // lexicographic numbering. Therefore, we need to convert the DoF
        // indices if DoFHander is L2 conforming and has degree > 0.
        const bool needs_conversion =
          dof_handler.get_fe().conforming_space ==
            FiniteElementData<dim>::Conformity::L2 &&
          (dof_handler.get_fe().degree > 0);
        std::vector<unsigned int> lexicographic_to_hierarchic;
        if (needs_conversion)
          lexicographic_to_hierarchic =
            FETools::lexicographic_to_hierarchic_numbering<dim>(degree);

        const Utilities::MPI::Partitioner partitioner_support_points(
          dof_handler_support_points.locally_owned_dofs(),
          dof_handler_support_points.get_communicator());

        const Utilities::MPI::Partitioner partitioner_dof(
          dof_handler.locally_owned_dofs(),
          DoFTools::extract_locally_relevant_dofs(dof_handler),
          dof_handler.get_communicator());

        std::vector<bool> dof_processed(partitioner_dof.locally_owned_size() +
                                          partitioner_dof.n_ghost_indices(),
                                        false);


        std::vector<types::global_dof_index> support_point_indices(
          support_points_per_cell);
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          {
            if (cell->is_locally_owned() || cell->is_ghost())
              {
                const auto cell_support_point =
                  cell->as_dof_handler_iterator(dof_handler_support_points);

                cell_support_point->get_dof_indices(support_point_indices);
                cell->get_dof_indices(dof_indices);

                // collect unconstrained DoFs for support point (assuming
                // support_points_per_cell==dofs_per_cell, i.e. n_components==1)
                for (unsigned int i = 0; i < support_point_indices.size(); ++i)
                  if (partitioner_support_points.in_local_range(
                        support_point_indices[i]))
                    {
                      const auto global_dof_idx =
                        needs_conversion ?
                          dof_indices[lexicographic_to_hierarchic[i]] :
                          dof_indices[i];

                      const auto local_dof_idx =
                        partitioner_dof.global_to_local(global_dof_idx);

                      AssertIndexRange(local_dof_idx, dof_processed.size());

                      if (dof_processed[local_dof_idx] == false)
                        {
                          if (!constraint.is_constrained(global_dof_idx))
                            support_point_dofs.emplace_back(
                              partitioner_support_points.global_to_local(
                                support_point_indices[i]),
                              global_dof_idx);

                          dof_processed[local_dof_idx] = true;
                        }
                    }
              }
          }
      }

      // sort for support points (stable sort needed for multiple components)
      std::stable_sort(support_point_dofs.begin(),
                       support_point_dofs.end(),
                       [](const auto &a, const auto &b) {
                         return a.first < b.first;
                       });

      // convert to CRS format
      std::vector<types::global_dof_index> dof_indices;
      dof_indices.reserve(support_point_dofs.size());
      std::vector<unsigned int> dof_ptrs;
      dof_ptrs.reserve(dof_handler_support_points.n_locally_owned_dofs() + 1);
      dof_ptrs.push_back(0);
      std::vector<unsigned int> support_point_indices;
      support_point_indices.reserve(
        dof_handler_support_points.n_locally_owned_dofs());

      auto it = support_point_dofs.begin();
      while (it != support_point_dofs.end())
        {
          const unsigned int index = std::get<0>(*it);
          while (it != support_point_dofs.end() && it->first == index)
            {
              dof_indices.push_back(it->second);
              ++it;
            }
          support_point_indices.push_back(index);
          dof_ptrs.push_back(dof_indices.size());
        }

      return std::make_tuple(std::move(support_point_indices),
                             std::move(dof_ptrs),
                             std::move(dof_indices));
    }


    /**
     * Create DoFHandler with unique support points.
     */
    template <int dim, int spacedim>
    std::shared_ptr<const DoFHandler<dim, spacedim>>
    create_support_point_dof_handler(
      const DoFHandler<dim, spacedim> &dof_handler)
    {
      const auto &fe           = dof_handler.get_fe();
      const auto &tria         = dof_handler.get_triangulation();
      const auto  degree       = fe.degree;
      const auto  n_components = fe.n_components();

      if (n_components == 1 &&
          (fe.conforming_space == FiniteElementData<dim>::Conformity::H1 ||
           degree == 0))
        {
          // in case a DG space of order 0 is provided, DoFs indices are always
          // uniquely assigned to support points (they are always defined in the
          // center of the element) and are never shared at vertices or faces.
          return std::shared_ptr<const DoFHandler<dim, spacedim>>(&dof_handler,
                                                                  [](auto *) {
                                                                  });
        }
      else
        {
          // Create dummy dof handler for support point numbering.
          // Unique support points are generally numbered according to FE_Q with
          // n_components==1. If degree==0 we use FE_DGQ which ensures a unique
          // support point numbering since the support point is located in the
          // center of the cell.
          auto dof_handler_support_points =
            std::make_shared<DoFHandler<dim, spacedim>>(tria);

          if (degree == 0)
            dof_handler_support_points->distribute_dofs(
              FE_DGQ<dim, spacedim>(degree));
          else
            dof_handler_support_points->distribute_dofs(
              FE_Q<dim, spacedim>(degree));

          return dof_handler_support_points;
        }
    }

    // Loop over cells and collect unique set of points
    template <int dim, typename Number>
    std::tuple<std::vector<Point<dim>>,
               std::vector<unsigned int>,
               std::vector<types::global_dof_index>>
    collect_unconstrained_unique_support_points(
      const DoFHandler<dim> &                  dof_handler,
      const Mapping<dim> &                     mapping,
      const dealii::AffineConstraints<Number> &constraint)
    {
      AssertThrow(dof_handler.get_fe().has_support_points(),
                  ExcNotImplemented());

      // create DoFHandler for support points
      const auto dof_handler_support_points =
        create_support_point_dof_handler(dof_handler);

      // compute mapping: index of locally owned support points to (global) DoF
      // indices
      const auto support_point_dofs_crs =
        support_point_indices_to_dof_indices(dof_handler,
                                             *dof_handler_support_points,
                                             constraint);

      const std::vector<unsigned int> &local_support_point_indices =
        std::get<0>(support_point_dofs_crs);

      // compute locally owned support points
      std::vector<Point<dim>> points;
      points.resize(local_support_point_indices.size());

      const auto locally_owned_support_point =
        dof_handler_support_points->locally_owned_dofs();
      std::vector<unsigned int> indices_state(
        locally_owned_support_point.n_elements(),
        numbers::invalid_unsigned_int);

      AssertIndexRange(local_support_point_indices.size(),
                       indices_state.size() + 1);

      for (unsigned int i = 0; i < local_support_point_indices.size(); ++i)
        indices_state[local_support_point_indices[i]] = i;

      const auto &  fe_support_point = dof_handler_support_points->get_fe();
      FEValues<dim> fe_values(mapping,
                              fe_support_point,
                              Quadrature<dim>(
                                fe_support_point.get_unit_support_points()),
                              update_quadrature_points);

      std::vector<types::global_dof_index> dof_indices(
        fe_support_point.n_dofs_per_cell());

      for (const auto &cell :
           dof_handler_support_points->active_cell_iterators() |
             IteratorFilters::LocallyOwnedCell())
        {
          fe_values.reinit(cell);
          cell->get_dof_indices(dof_indices);

          for (const unsigned int q : fe_values.quadrature_point_indices())
            if (locally_owned_support_point.is_element(dof_indices[q]))
              {
                const auto index =
                  locally_owned_support_point.index_within_set(dof_indices[q]);

                if (indices_state[index] != numbers::invalid_unsigned_int)
                  {
                    points[indices_state[index]] =
                      fe_values.quadrature_point(q);
                    indices_state[index] = numbers::invalid_unsigned_int;
                  }
              }
        }

      return std::make_tuple(
        std::move(points),
        std::move(std::get<1>(support_point_dofs_crs)),  // global_dofs_ptrs
        std::move(std::get<2>(support_point_dofs_crs))); // global_dofs_indices
    }

  } // namespace
} // namespace internal


template <int dim, typename Number>
void
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  reinit(const DoFHandler<dim> &          dof_handler_fine,
         const DoFHandler<dim> &          dof_handler_coarse,
         const Mapping<dim> &             mapping_fine,
         const Mapping<dim> &             mapping_coarse,
         const AffineConstraints<Number> &constraint_fine,
         const AffineConstraints<Number> &constraint_coarse)
{
  AssertThrow(dof_handler_coarse.get_fe().has_support_points(),
              ExcNotImplemented());
  Assert(dof_handler_coarse.get_fe().n_components() > 0 &&
           dof_handler_fine.get_fe().n_components() > 0,
         ExcNotImplemented());

  this->fine_element_is_continuous =
    dof_handler_fine.get_fe().n_dofs_per_vertex() > 0;

  // collect points, ptrs, and global indices
  const auto points_ptrs_indices =
    internal::collect_unconstrained_unique_support_points(dof_handler_fine,
                                                          mapping_fine,
                                                          constraint_fine);
  const auto &global_dof_indices = std::get<2>(points_ptrs_indices);

  // create partitioners and internal vectors
  {
    IndexSet locally_active_dofs =
      DoFTools::extract_locally_active_dofs(dof_handler_coarse);
    this->partitioner_coarse.reset(
      new Utilities::MPI::Partitioner(dof_handler_coarse.locally_owned_dofs(),
                                      locally_active_dofs,
                                      dof_handler_coarse.get_communicator()));

    this->vec_coarse.reinit(this->partitioner_coarse);
  }
  {
    // in case a DG space of order 0 is provided, DoFs indices are never defined
    // on element faces or vertices and therefore, the partitioner is fine
    IndexSet locally_relevant_dofs(dof_handler_fine.n_dofs());
    if (!this->fine_element_is_continuous &&
        dof_handler_fine.get_fe().degree != 0)
      locally_relevant_dofs.add_indices(global_dof_indices.begin(),
                                        global_dof_indices.end());

    this->partitioner_fine.reset(
      new Utilities::MPI::Partitioner(dof_handler_fine.locally_owned_dofs(),
                                      locally_relevant_dofs,
                                      dof_handler_fine.get_communicator()));

    this->vec_fine.reinit(this->partitioner_fine);
  }

  const auto &points = std::get<0>(points_ptrs_indices);

  // using level_dof_indices_fine_ptrs always works but in case of CG or DG
  // with degree==0 and n_components==1 support points to dof mapping is unique
  // and we don't need it.
  if (dof_handler_fine.get_fe().n_components() == 1 &&
      (this->fine_element_is_continuous ||
       dof_handler_fine.get_fe().degree == 0))
    this->level_dof_indices_fine_ptrs.clear();
  else
    this->level_dof_indices_fine_ptrs = std::get<1>(points_ptrs_indices);

  // fill level_dof_indices_fine with local indices
  this->level_dof_indices_fine.resize(global_dof_indices.size());
  for (unsigned int i = 0; i < global_dof_indices.size(); ++i)
    this->level_dof_indices_fine[i] =
      this->partitioner_fine->global_to_local(global_dof_indices[i]);

  // hand points over to RPE
  rpe.reinit(points, dof_handler_coarse.get_triangulation(), mapping_coarse);

  // set up MappingInfo for easier data access
  mapping_info = internal::fill_mapping_info<dim, Number>(rpe);

  // set up constraints
  const auto &cell_data = rpe.get_cell_data();

  constraint_info.reinit(dof_handler_coarse,
                         cell_data.cells.size(),
                         false /*TODO*/);

  for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
    {
      typename DoFHandler<dim>::active_cell_iterator cell(
        &rpe.get_triangulation(),
        cell_data.cells[i].first,
        cell_data.cells[i].second,
        &dof_handler_coarse);

      constraint_info.read_dof_indices(i,
                                       numbers::invalid_unsigned_int,
                                       cell,
                                       constraint_coarse,
                                       this->partitioner_coarse);
    }

  constraint_info.finalize();

  const auto &fe_base = dof_handler_coarse.get_fe().base_element(0);

  if (const auto fe = dynamic_cast<const FE_Q<dim> *>(&fe_base))
    fe_coarse = std::make_unique<FE_DGQ<dim>>(fe->get_degree());
  else if (const auto fe = dynamic_cast<const FE_DGQ<dim> *>(&fe_base))
    fe_coarse = fe->clone();
  else
    AssertThrow(false, ExcMessage(dof_handler_coarse.get_fe().get_name()));
}



template <int dim, typename Number>
void
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  std::vector<Number> evaluation_point_results;
  std::vector<Number> buffer;

  const auto evaluation_function = [&](auto &values, const auto &cell_data) {
    std::vector<Number> solution_values;

    FEPointEvaluation<1, dim, dim, Number> evaluator(*mapping_info, *fe_coarse);

    for (unsigned int cell = 0; cell < cell_data.cells.size(); ++cell)
      {
        solution_values.resize(fe_coarse->n_dofs_per_cell());

        // gather and resolve constraints
        internal::VectorReader<Number, VectorizedArrayType> reader;
        constraint_info.read_write_operation(
          reader,
          src,
          reinterpret_cast<VectorizedArrayType *>(solution_values.data()),
          cell,
          1,
          solution_values.size(),
          true);

        // evaluate and scatter
        evaluator.reinit(cell);

        evaluator.evaluate(solution_values, dealii::EvaluationFlags::values);

        for (const auto q : evaluator.quadrature_point_indices())
          values[q + cell_data.reference_point_ptrs[cell]] =
            evaluator.get_value(q);
      }
  };

  rpe.template evaluate_and_process<Number>(evaluation_point_results,
                                            buffer,
                                            evaluation_function);

  // Weight operator in case some points are owned by multiple cells.
  if (rpe.is_map_unique() == false)
    {
      const auto evaluation_point_results_temp = evaluation_point_results;
      evaluation_point_results.clear();
      evaluation_point_results.reserve(rpe.get_point_ptrs().size() - 1);

      const auto &ptr = rpe.get_point_ptrs();

      for (unsigned int i = 0; i < ptr.size() - 1; ++i)
        {
          const auto n_entries = ptr[i + 1] - ptr[i];

          Number result = 0.0;

          if (n_entries > 0)
            {
              for (unsigned int j = 0; j < n_entries; ++j)
                result += evaluation_point_results_temp[ptr[i] + j];
              result /= n_entries;
            }
          evaluation_point_results.push_back(result);
        }
    }

  for (unsigned int j = 0; j < evaluation_point_results.size(); ++j)
    {
      if (level_dof_indices_fine_ptrs.empty())
        {
          dst.local_element(this->level_dof_indices_fine[j]) +=
            evaluation_point_results[j];
        }
      else
        {
          for (unsigned int i = this->level_dof_indices_fine_ptrs[j];
               i < this->level_dof_indices_fine_ptrs[j + 1];
               ++i)
            dst.local_element(this->level_dof_indices_fine[i]) +=
              evaluation_point_results[j];
        }
    }
}



template <int dim, typename Number>
void
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  std::vector<Number> evaluation_point_results;
  std::vector<Number> buffer;

  evaluation_point_results.resize(rpe.get_point_ptrs().size() - 1);

  for (unsigned int j = 0; j < evaluation_point_results.size(); ++j)
    {
      if (level_dof_indices_fine_ptrs.empty())
        {
          evaluation_point_results[j] =
            src.local_element(this->level_dof_indices_fine[j]);
        }
      else
        {
          evaluation_point_results[j] = 0.0;
          for (unsigned int i = this->level_dof_indices_fine_ptrs[j];
               i < this->level_dof_indices_fine_ptrs[j + 1];
               ++i)
            evaluation_point_results[j] +=
              src.local_element(this->level_dof_indices_fine[i]);
        }
    }

  // Weight operator in case some points are owned by multiple cells.
  if (rpe.is_map_unique() == false)
    {
      const auto &ptr = rpe.get_point_ptrs();

      for (unsigned int i = 0; i < ptr.size() - 1; ++i)
        {
          const auto n_entries = ptr[i + 1] - ptr[i];
          if (n_entries == 0)
            continue;

          evaluation_point_results[i] /= n_entries;
        }
    }

  const auto evaluation_function = [&](const auto &values,
                                       const auto &cell_data) {
    std::vector<Number>                    solution_values;
    FEPointEvaluation<1, dim, dim, Number> evaluator(*mapping_info, *fe_coarse);

    for (unsigned int cell = 0; cell < cell_data.cells.size(); ++cell)
      {
        solution_values.resize(fe_coarse->n_dofs_per_cell());

        // gather and integrate
        evaluator.reinit(cell);

        for (const auto q : evaluator.quadrature_point_indices())
          evaluator.submit_value(
            values[q + cell_data.reference_point_ptrs[cell]], q);

        evaluator.test_and_sum(solution_values, EvaluationFlags::values);

        // resolve constraints and scatter
        internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
          writer;
        constraint_info.read_write_operation(
          writer,
          dst,
          reinterpret_cast<VectorizedArrayType *>(solution_values.data()),
          cell,
          1,
          solution_values.size(),
          true);
      }
  };

  rpe.template process_and_evaluate<Number>(evaluation_point_results,
                                            buffer,
                                            evaluation_function);
}



template <int dim, typename Number>
void
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  interpolate(LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const
{
  AssertThrow(false, ExcNotImplemented());
  (void)dst;
  (void)src;
}



template <int dim, typename Number>
void
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_fine)
{
  this->internal_enable_inplace_operations_if_possible(
    external_partitioner_coarse,
    external_partitioner_fine,
    constraint_info,
    this->level_dof_indices_fine);
}



template <int dim, typename Number>
std::size_t
MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>::
  memory_consumption() const
{
  std::size_t size = 0;

  size += this->partitioner_coarse->memory_consumption();
  size += this->vec_coarse.memory_consumption();
  size += MemoryConsumption::memory_consumption(this->level_dof_indices_fine);
  // TODO: add consumption for rpe, mapping_info and constraint_info.

  return size;
}


DEAL_II_NAMESPACE_CLOSE

#endif
