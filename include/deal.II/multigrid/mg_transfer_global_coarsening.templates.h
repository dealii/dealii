// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/matrix_free/evaluation_kernels.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

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
      if ((degree_fine == (2 * deg + 1)) && (degree_coarse == deg))
        fu.template run<2 * deg + 1, deg>(); // h-MG
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg / 2, 1u)))
        fu.template run<deg, std::max(deg / 2u, 1u)>(); // p-MG: bisection
      else if ((degree_fine == deg) && (degree_coarse == deg))
        fu.template run<deg, deg>(); // identity (nothing to do)
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg - 1, 1u)))
        fu.template run<deg, std::max(deg - 1u, 1u)>(); // p-MG: --
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
        degree_fine + 1,
        Number,
        Number>::do_forward(1,
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
        0,
        Number,
        Number>::do_forward(1,
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
        degree_fine == 0 ? -1 : (degree_fine + 1),
        Number,
        Number>::do_backward(1,
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
        0,
        Number,
        Number>::do_backward(1,
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
    get_child_offsets(const unsigned int dofs_per_cell_coarse,
                      const unsigned int fe_shift_1d,
                      const unsigned int fe_degree)
    {
      std::vector<std::vector<unsigned int>> cell_local_chilren_indices(
        GeometryInfo<dim>::max_children_per_cell,
        std::vector<unsigned int>(dofs_per_cell_coarse));
      for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
           c++)
        get_child_offset<dim>(c,
                              fe_shift_1d,
                              fe_degree,
                              cell_local_chilren_indices[c]);
      return cell_local_chilren_indices;
    }

    template <int dim>
    std::vector<std::vector<unsigned int>>
    get_child_offsets_general(const unsigned int dofs_per_cell_coarse)
    {
      std::vector<std::vector<unsigned int>> cell_local_chilren_indices(
        GeometryInfo<dim>::max_children_per_cell,
        std::vector<unsigned int>(dofs_per_cell_coarse));
      for (unsigned int c = 0, k = 0;
           c < GeometryInfo<dim>::max_children_per_cell;
           c++)
        for (unsigned int d = 0; d < dofs_per_cell_coarse; ++d, ++k)
          cell_local_chilren_indices[c][d] = k;
      return cell_local_chilren_indices;
    }

    template <int dim>
    types::coarse_cell_id
    convert_cell_id_binary_type_to_level_coarse_cell_id(
      const typename CellId::binary_type &binary_representation)
    {
      // exploiting the structure of CellId::binary_type
      // see also the documentation of CellId

      // actual coarse-grid id
      const unsigned int coarse_cell_id  = binary_representation[0];
      const unsigned int n_child_indices = binary_representation[1] >> 2;

      const unsigned int children_per_value =
        sizeof(CellId::binary_type::value_type) * 8 / dim;
      unsigned int child_level  = 0;
      unsigned int binary_entry = 2;

      // path to the get to the cell
      std::vector<unsigned int> cell_indices;
      while (child_level < n_child_indices)
        {
          Assert(binary_entry < binary_representation.size(),
                 ExcInternalError());

          for (unsigned int j = 0; j < children_per_value; ++j)
            {
              unsigned int cell_index =
                (((binary_representation[binary_entry] >> (j * dim))) &
                 (GeometryInfo<dim>::max_children_per_cell - 1));
              cell_indices.push_back(cell_index);
              ++child_level;
              if (child_level == n_child_indices)
                break;
            }
          ++binary_entry;
        }

      // compute new coarse-grid id: c_{i+1} = c_{i}*2^dim + q;
      types::coarse_cell_id level_coarse_cell_id = coarse_cell_id;
      for (auto i : cell_indices)
        level_coarse_cell_id =
          level_coarse_cell_id * GeometryInfo<dim>::max_children_per_cell + i;

      return level_coarse_cell_id;
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
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (fe.restriction_is_additive(i) == true)
                continue;

              bool do_zero = false;
              for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                if (matrix_other(i, j) != 0.)
                  do_zero = true;

              if (do_zero)
                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                  matrix(i, j) = 0.0;
            }
        }
      return matrix;
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
  class CellIDTranslator
  {
  public:
    CellIDTranslator(const types::global_cell_index n_coarse_cells,
                     const types::global_cell_index n_global_levels)
      : n_coarse_cells(n_coarse_cells)
      , n_global_levels(n_global_levels)
    {
      tree_sizes.push_back(0);
      for (unsigned int i = 0; i < n_global_levels; ++i)
        tree_sizes.push_back(
          tree_sizes.back() +
          Utilities::pow(GeometryInfo<dim>::max_children_per_cell, i) *
            n_coarse_cells);
    }

    types::global_cell_index
    size() const
    {
      return n_coarse_cells *
             (Utilities::pow(GeometryInfo<dim>::max_children_per_cell,
                             n_global_levels) -
              1);
    }

    template <typename T>
    types::global_cell_index
    translate(const T &cell) const
    {
      types::global_cell_index id = 0;

      id += convert_cell_id_binary_type_to_level_coarse_cell_id<dim>(
        cell->id().template to_binary<dim>());

      id += tree_sizes[cell->level()];

      return id;
    }

    template <typename T>
    types::global_cell_index
    translate(const T &cell, const types::global_cell_index i) const
    {
      return (translate(cell) - tree_sizes[cell->level()]) *
               GeometryInfo<dim>::max_children_per_cell +
             i + tree_sizes[cell->level() + 1];
    }

    CellId
    to_cell_id(const types::global_cell_index id) const
    {
      std::vector<std::uint8_t> child_indices;

      types::global_cell_index id_temp = id;

      types::global_cell_index level = 0;

      for (; level < n_global_levels; ++level)
        if (id < tree_sizes[level])
          break;
      level -= 1;

      id_temp -= tree_sizes[level];

      for (types::global_cell_index l = 0; l < level; ++l)
        {
          child_indices.push_back(id_temp %
                                  GeometryInfo<dim>::max_children_per_cell);
          id_temp /= GeometryInfo<dim>::max_children_per_cell;
        }

      std::reverse(child_indices.begin(), child_indices.end());

      return {id_temp, child_indices}; // TODO
    }

  private:
    const types::global_cell_index        n_coarse_cells;
    const types::global_cell_index        n_global_levels;
    std::vector<types::global_cell_index> tree_sizes;
  };

  template <int dim>
  class FineDoFHandlerView
  {
  public:
    FineDoFHandlerView(const DoFHandler<dim> &mesh_fine,
                       const DoFHandler<dim> &mesh_coarse,
                       const unsigned int     mg_level_fine)
      : mesh_fine(mesh_fine)
      , mesh_coarse(mesh_coarse)
      , mg_level_fine(mg_level_fine)
      , communicator(
          mesh_fine.get_communicator() /*TODO: fix for different comms*/)
      , cell_id_translator(n_coarse_cells(mesh_fine),
                           n_global_levels(mesh_fine))
    {
      AssertDimension(n_coarse_cells(mesh_fine), n_coarse_cells(mesh_coarse));
      AssertIndexRange(n_global_levels(mesh_coarse),
                       n_global_levels(mesh_fine) + 1);
    }

    void
    reinit(const IndexSet &is_dst_locally_owned,
           const IndexSet &is_dst_remote_input,
           const IndexSet &is_src_locally_owned,
           const bool      check_if_elements_in_is_dst_remote_exist = false)
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)is_dst_locally_owned;
      (void)is_dst_remote_input;
      (void)is_src_locally_owned;
      (void)check_if_elements_in_is_dst_remote_exist;
#else
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
              std::pair<types::global_cell_index, types::global_cell_index>,
              unsigned int>
              consensus_algorithm(process, communicator);
            consensus_algorithm.run();
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
        std::pair<types::global_cell_index, types::global_cell_index>,
        unsigned int>
        consensus_algorithm(process, communicator);
      consensus_algorithm.run();

      this->is_dst_locally_owned = is_dst_locally_owned;
      this->is_dst_remote        = is_dst_remote;
      this->is_src_locally_owned = is_src_locally_owned;

#  if false
        std::cout << "IS_SRC_LOCALLY_OWNED" << std::endl;
        for (auto i : is_src_locally_owned)
          std::cout << cell_id_translator.to_cell_id(i) << std::endl;
        std::cout << std::endl << std::endl << std::endl;


        std::cout << "IS_DST_LOCALLY_OWNED" << std::endl;
        for (auto i : is_dst_locally_owned)
          std::cout << cell_id_translator.to_cell_id(i) << std::endl;
        std::cout << std::endl << std::endl << std::endl;
#  endif

      const auto &dof_handler_dst = mesh_fine; // TODO: remove
      const auto &tria_dst        = mesh_fine.get_triangulation(); // TODO

      const auto targets_with_indexset = process.get_requesters();

      std::map<unsigned int, std::vector<types::global_dof_index>>
                               indices_to_be_sent;
      std::vector<MPI_Request> requests;
      requests.reserve(targets_with_indexset.size());

      {
        std::vector<types::global_dof_index> indices;

        for (const auto &i : targets_with_indexset)
          {
            indices_to_be_sent[i.first] = {};
            std::vector<types::global_dof_index> &buffer =
              indices_to_be_sent[i.first];

            for (auto cell_id : i.second)
              {
                typename DoFHandler<dim>::cell_iterator cell(
                  *tria_dst.create_cell_iterator(
                    cell_id_translator.to_cell_id(cell_id)),
                  &dof_handler_dst);

                indices.resize(cell->get_fe().n_dofs_per_cell());

                if (mg_level_fine == numbers::invalid_unsigned_int)
                  cell->get_dof_indices(indices);
                else
                  cell->get_mg_dof_indices(indices);

                buffer.push_back(cell->active_fe_index());
                buffer.insert(buffer.end(), indices.begin(), indices.end());
              }

            requests.resize(requests.size() + 1);

            const auto ierr_1 = MPI_Isend(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::internal::mpi_type_id(buffer.data()),
              i.first,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &requests.back());
            AssertThrowMPI(ierr_1);
          }
      }

      std::vector<types::global_dof_index> ghost_indices;

      // process local cells
      {
        std::vector<types::global_dof_index> indices;

        for (const auto id : is_dst_locally_owned)
          {
            const auto cell_id = cell_id_translator.to_cell_id(id);

            typename DoFHandler<dim>::cell_iterator cell_(
              *tria_dst.create_cell_iterator(cell_id), &dof_handler_dst);

            indices.resize(cell_->get_fe().n_dofs_per_cell());

            if (mg_level_fine == numbers::invalid_unsigned_int)
              cell_->get_dof_indices(indices);
            else
              cell_->get_mg_dof_indices(indices);

            ghost_indices.insert(ghost_indices.end(),
                                 indices.begin(),
                                 indices.end());
          }
      }

      {
        std::map<unsigned int, std::vector<types::global_dof_index>>
          rank_to_ids;

        std::set<unsigned int> ranks;

        for (auto i : is_dst_remote_owners)
          ranks.insert(i);

        for (auto i : ranks)
          rank_to_ids[i] = {};

        for (unsigned int i = 0; i < is_dst_remote_owners.size(); ++i)
          rank_to_ids[is_dst_remote_owners[i]].push_back(
            is_dst_remote.nth_index_in_set(i));


        for (unsigned int i = 0; i < ranks.size(); ++i)
          {
            MPI_Status status;
            const int  ierr_1 = MPI_Probe(
              MPI_ANY_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &status);
            AssertThrowMPI(ierr_1);

            std::vector<types::global_dof_index> buffer;

            int       message_length;
            const int ierr_2 =
              MPI_Get_count(&status,
                            Utilities::MPI::internal::mpi_type_id(
                              buffer.data()),
                            &message_length);
            AssertThrowMPI(ierr_2);

            buffer.resize(message_length);

            const int ierr_3 = MPI_Recv(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::internal::mpi_type_id(buffer.data()),
              status.MPI_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr_3);

            for (unsigned int i = 0; i < buffer.size();
                 i += dof_handler_dst.get_fe(buffer[i]).n_dofs_per_cell() + 1)
              ghost_indices.insert(
                ghost_indices.end(),
                buffer.begin() + i + 1,
                buffer.begin() + i + 1 +
                  dof_handler_dst.get_fe(buffer[i]).n_dofs_per_cell());

            const unsigned int rank = status.MPI_SOURCE;

            const auto ids = rank_to_ids[rank];

            std::vector<types::global_dof_index> indices;

            for (unsigned int i = 0, k = 0; i < ids.size(); ++i)
              {
                const unsigned int active_fe_index = buffer[k++];

                indices.resize(
                  dof_handler_dst.get_fe(active_fe_index).n_dofs_per_cell());

                for (unsigned int j = 0; j < indices.size(); ++j, ++k)
                  indices[j] = buffer[k];
                map[ids[i]] = {active_fe_index, indices};
              }
          }

        const int ierr_1 =
          MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr_1);
      }

      std::sort(ghost_indices.begin(), ghost_indices.end());
      ghost_indices.erase(std::unique(ghost_indices.begin(),
                                      ghost_indices.end()),
                          ghost_indices.end());

      this->is_extended_locally_owned = dof_handler_dst.locally_owned_dofs();

      this->is_extended_ghosts = IndexSet(dof_handler_dst.n_dofs());
      this->is_extended_ghosts.add_indices(ghost_indices.begin(),
                                           ghost_indices.end());
      this->is_extended_ghosts.subtract_set(this->is_extended_locally_owned);

#endif
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
                *mesh_fine.get_triangulation().create_cell_iterator(cell->id()),
                &mesh_fine);
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
                        *mesh_fine.get_triangulation().create_cell_iterator(
                          cell->id()),
                        &mesh_fine))
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
      const bool is_cell_remotly_owned = this->is_dst_remote.is_element(id);

      return FineDoFHandlerViewCell(
        [cell]() {
          AssertThrow(false, ExcNotImplemented()); // currently we do not need
                                                   // children of children

          return false;
        },
        [cell, is_cell_locally_owned, is_cell_remotly_owned, c, id, this](
          auto &dof_indices) {
          if (is_cell_locally_owned)
            {
              const auto cell_fine =
                (typename DoFHandler<dim>::cell_iterator(
                   *mesh_fine.get_triangulation().create_cell_iterator(
                     cell->id()),
                   &mesh_fine))
                  ->child(c);

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
    const DoFHandler<dim> &mesh_fine;
    const DoFHandler<dim> &mesh_coarse;
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

    static unsigned int
    n_coarse_cells(const DoFHandler<dim> &mesh)
    {
      types::coarse_cell_id n_coarse_cells = 0;

      for (const auto &cell : mesh.get_triangulation().active_cell_iterators())
        if (!cell->is_artificial())
          n_coarse_cells =
            std::max(n_coarse_cells, cell->id().get_coarse_cell_id());

      return Utilities::MPI::max(n_coarse_cells, mesh.get_communicator()) + 1;
    }

    static unsigned int
    n_global_levels(const DoFHandler<dim> &mesh)
    {
      return mesh.get_triangulation().n_global_levels();
    }
  };

  template <int dim>
  class GlobalCoarseningFineDoFHandlerView : public FineDoFHandlerView<dim>
  {
  public:
    GlobalCoarseningFineDoFHandlerView(const DoFHandler<dim> &dof_handler_dst,
                                       const DoFHandler<dim> &dof_handler_src)
      : FineDoFHandlerView<
          dim>(dof_handler_dst, dof_handler_src, numbers::invalid_unsigned_int /*global coarsening only possible on active levels*/)
    {
      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      for (const auto &cell : tria_dst.active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          is_dst_locally_owned.add_index(
            this->cell_id_translator.translate(cell));


      for (const auto &cell : tria_src.active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          {
            is_src_locally_owned.add_index(
              this->cell_id_translator.translate(cell));
            is_dst_remote.add_index(this->cell_id_translator.translate(cell));

            if (cell->level() + 1u == tria_dst.n_global_levels())
              continue;

            for (unsigned int i = 0;
                 i < GeometryInfo<dim>::max_children_per_cell;
                 ++i)
              is_dst_remote.add_index(
                this->cell_id_translator.translate(cell, i));
          }

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

      const auto loop = [&](const auto &tria,
                            const auto  level,
                            const auto &op) {
        if (level == numbers::invalid_unsigned_int)
          {
            for (const auto &cell : tria.active_cell_iterators())
              if (!cell->is_artificial() && cell->is_locally_owned())
                op(cell);
          }
        else
          {
            for (const auto &cell : tria.cell_iterators_on_level(level))
              if (cell->level_subdomain_id() == tria.locally_owned_subdomain())
                op(cell);
          }
      };

      loop(tria_dst, mg_level_fine, fine_operation);
      loop(tria_src, mg_level_coarse, coarse_operation);

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   false);
    }
  };

  class MGTwoLevelTransferImplementation
  {
    template <int dim, typename Number>
    static void
    precompute_constraints(
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraint_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      std::vector<types::global_dof_index> dependencies;

      const auto &locally_owned_dofs = dof_handler_coarse.locally_owned_dofs();

      for (const auto i : locally_owned_dofs)
        {
          if (constraint_coarse.is_constrained(i) == false)
            continue;

          const auto constraints = constraint_coarse.get_constraint_entries(i);

          if (constraints)
            for (const auto &p : *constraints)
              if (locally_owned_dofs.is_element(p.first) == false)
                dependencies.push_back(p.first);
        }

      std::sort(dependencies.begin(), dependencies.end());
      dependencies.erase(std::unique(dependencies.begin(), dependencies.end()),
                         dependencies.end());

      IndexSet locally_relevant_dofs(locally_owned_dofs.size());
      locally_relevant_dofs.add_indices(dependencies.begin(),
                                        dependencies.end());

      const auto partitioner = std::make_shared<Utilities::MPI::Partitioner>(
        locally_owned_dofs,
        locally_relevant_dofs,
        dof_handler_coarse.get_communicator());

      transfer.vec_coarse_constraints.reinit(partitioner);

      transfer.constraint_coarse_distribute_indices.clear();
      transfer.constraint_coarse_distribute_values.clear();
      transfer.constraint_coarse_distribute_ptr = {0};

      for (const auto i : partitioner->locally_owned_range())
        {
          Assert(constraint_coarse.is_inhomogeneously_constrained(i) == false,
                 ExcNotImplemented());

          if (constraint_coarse.is_constrained(i))
            {
              const auto constraints =
                constraint_coarse.get_constraint_entries(i);

              if (constraints)
                for (const auto &p : *constraints)
                  {
                    transfer.constraint_coarse_distribute_indices.emplace_back(
                      partitioner->global_to_local(p.first));
                    transfer.constraint_coarse_distribute_values.emplace_back(
                      p.second);
                  }
            }

          transfer.constraint_coarse_distribute_ptr.push_back(
            transfer.constraint_coarse_distribute_indices.size());
        }
    }



    template <int dim, typename Number>
    static void
    precompute_restriction_constraints(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const dealii::AffineConstraints<Number> &constraint_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      transfer.distribute_local_to_global_indices.clear();
      transfer.distribute_local_to_global_values.clear();
      transfer.distribute_local_to_global_ptr = {0};

      const auto fu = [&](const auto &index_set) {
        for (const auto i : index_set)
          {
            Assert(constraint_coarse.is_inhomogeneously_constrained(i) == false,
                   ExcNotImplemented());

            if (constraint_coarse.is_constrained(i))
              {
                const auto constraints =
                  constraint_coarse.get_constraint_entries(i);

                if (constraints)
                  for (const auto &p : *constraints)
                    {
                      transfer.distribute_local_to_global_indices.emplace_back(
                        partitioner->global_to_local(p.first));
                      transfer.distribute_local_to_global_values.emplace_back(
                        p.second);
                    }

                // add a dummy entry for homogeneous constraints
                if (transfer.distribute_local_to_global_indices.size() ==
                    transfer.distribute_local_to_global_ptr.back())
                  {
                    transfer.distribute_local_to_global_indices.emplace_back(
                      numbers::invalid_unsigned_int);
                    transfer.distribute_local_to_global_values.emplace_back(
                      0.0);
                  }
              }

            transfer.distribute_local_to_global_ptr.push_back(
              transfer.distribute_local_to_global_indices.size());
          }
      };

      fu(partitioner->locally_owned_range());
      fu(partitioner->ghost_indices());
    }



    template <int dim, typename Number>
    static void
    compute_weights(
      const DoFHandler<dim> &                  dof_handler_fine,
      const unsigned int                       mg_level_fine,
      const dealii::AffineConstraints<Number> &constraint_fine,
      const MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &                                         transfer,
      LinearAlgebra::distributed::Vector<Number> &touch_count)
    {
      LinearAlgebra::distributed::Vector<Number> touch_count_;
      touch_count.reinit(transfer.partitioner_fine);

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

      touch_count_.reinit(partitioner_fine_);

      std::vector<types::global_dof_index> local_dof_indices;

      if (mg_level_fine == numbers::invalid_unsigned_int)
        {
          for (const auto &cell : dof_handler_fine.active_cell_iterators())
            {
              if (cell->is_locally_owned() == false)
                continue;

              local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());

              cell->get_dof_indices(local_dof_indices);

              for (auto i : local_dof_indices)
                if (constraint_fine.is_constrained(i) == false)
                  touch_count_[i] += 1;
            }
        }
      else
        {
          for (const auto &cell :
               dof_handler_fine.mg_cell_iterators_on_level(mg_level_fine))
            {
              if (cell->level_subdomain_id() !=
                  dof_handler_fine.get_triangulation()
                    .locally_owned_subdomain())
                continue;

              local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());

              cell->get_mg_dof_indices(local_dof_indices);

              for (auto i : local_dof_indices)
                if (constraint_fine.is_constrained(i) == false)
                  touch_count_[i] += 1;
            }
        }

      touch_count_.compress(VectorOperation::add);

      for (unsigned int i = 0; i < touch_count_.local_size(); ++i)
        touch_count_.local_element(i) =
          constraint_fine.is_constrained(
            touch_count_.get_partitioner()->local_to_global(i)) ?
            Number(0.) :
            Number(1.) / touch_count_.local_element(i);

      touch_count_.update_ghost_values();

      // copy weights to other indexset
      touch_count.copy_locally_owned_data_from(touch_count_);
      touch_count.update_ghost_values();
    }



  public:
    template <int dim, typename Number>
    static void
    reinit_geometric_transfer(
      const DoFHandler<dim> &                  dof_handler_fine,
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraint_fine,
      const dealii::AffineConstraints<Number> &constraint_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      const GlobalCoarseningFineDoFHandlerView<dim> view(dof_handler_fine,
                                                         dof_handler_coarse);

      // copy constrain matrix; TODO: why only for the coarse level?
      precompute_constraints(dof_handler_coarse, constraint_coarse, transfer);

      // gather ranges for active FE indices on both fine and coarse dofhandlers
      std::array<unsigned int, 2> min_active_fe_indices = {
        {std::numeric_limits<unsigned int>::max(),
         std::numeric_limits<unsigned int>::max()}};
      std::array<unsigned int, 2> max_active_fe_indices = {{0, 0}};

      for (const auto &cell : dof_handler_fine.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            min_active_fe_indices[0] =
              std::min(min_active_fe_indices[0], cell->active_fe_index());
            max_active_fe_indices[0] =
              std::max(max_active_fe_indices[0], cell->active_fe_index());
          }

      for (const auto &cell : dof_handler_coarse.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            min_active_fe_indices[1] =
              std::min(min_active_fe_indices[1], cell->active_fe_index());
            max_active_fe_indices[1] =
              std::max(max_active_fe_indices[1], cell->active_fe_index());
          }

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

      const auto &fe_fine = dof_handler_fine.get_fe(min_active_fe_indices[0]);
      const auto &fe_coarse =
        dof_handler_coarse.get_fe(min_active_fe_indices[1]);

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
          IndexSet locally_relevant_dofs;
          DoFTools::extract_locally_relevant_dofs(dof_handler_coarse,
                                                  locally_relevant_dofs);

          transfer.partitioner_coarse =
            std::make_shared<Utilities::MPI::Partitioner>(
              dof_handler_coarse.locally_owned_dofs(),
              locally_relevant_dofs,
              dof_handler_coarse.get_communicator());
          transfer.vec_coarse.reinit(transfer.partitioner_coarse);
          precompute_restriction_constraints(transfer.partitioner_coarse,
                                             constraint_coarse,
                                             transfer);
        }
      }

      // helper function: to process the fine level cells; function @p fu_non_refined is
      // performed on cells that are not refined and @fu_refined is performed on
      // children of cells that are refined
      const auto process_cells =
        [&view, &dof_handler_coarse](const auto &fu_non_refined,
                                     const auto &fu_refined) {
          // loop over all active locally-owned cells
          for (const auto &cell_coarse :
               dof_handler_coarse.active_cell_iterators())
            {
              if (cell_coarse->is_artificial() == true ||
                  cell_coarse->is_locally_owned() == false)
                continue;

              // get a reference to the equivalent cell on the fine
              // triangulation
              const auto cell_coarse_on_fine_mesh = view.get_cell(cell_coarse);

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
        };

      // set up two mg-schemes
      //   (0) no refinement -> identity
      //   (1) h-refinement
      transfer.schemes.resize(2);

      // check if FE is the same; TODO: better way?
      AssertDimension(fe_coarse.dofs_per_cell, fe_fine.dofs_per_cell);

      // number of dofs on coarse and fine cells
      transfer.schemes[0].dofs_per_cell_coarse =
        transfer.schemes[0].dofs_per_cell_fine =
          transfer.schemes[1].dofs_per_cell_coarse = fe_coarse.dofs_per_cell;
      transfer.schemes[1].dofs_per_cell_fine =
        fe_coarse.dofs_per_cell * GeometryInfo<dim>::max_children_per_cell;

      // degree of FE on coarse and fine cell
      transfer.schemes[0].degree_coarse   = transfer.schemes[0].degree_fine =
        transfer.schemes[1].degree_coarse = fe_coarse.degree;
      transfer.schemes[1].degree_fine     = fe_coarse.degree * 2 + 1;

      // continuous or discontinuous
      transfer.fine_element_is_continuous = fe_fine.dofs_per_vertex > 0;

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


      const auto cell_local_chilren_indices =
        (reference_cell == ReferenceCells::get_hypercube<dim>()) ?
          get_child_offsets<dim>(transfer.schemes[0].dofs_per_cell_coarse,
                                 fe_fine.degree + 1,
                                 fe_fine.degree) :
          get_child_offsets_general<dim>(
            transfer.schemes[0].dofs_per_cell_coarse);


      // indices
      {
        transfer.schemes[0].level_dof_indices_fine.resize(
          transfer.schemes[0].dofs_per_cell_fine *
          transfer.schemes[0].n_coarse_cells);
        transfer.schemes[0].level_dof_indices_coarse.resize(
          transfer.schemes[0].dofs_per_cell_coarse *
          transfer.schemes[0].n_coarse_cells);

        transfer.schemes[1].level_dof_indices_fine.resize(
          transfer.schemes[1].dofs_per_cell_fine *
          transfer.schemes[1].n_coarse_cells);
        transfer.schemes[1].level_dof_indices_coarse.resize(
          transfer.schemes[1].dofs_per_cell_coarse *
          transfer.schemes[1].n_coarse_cells);

        std::vector<types::global_dof_index> local_dof_indices(
          transfer.schemes[0].dofs_per_cell_coarse);

        // ---------------------- lexicographic_numbering ----------------------
        std::vector<unsigned int> lexicographic_numbering;
        if (reference_cell == ReferenceCells::get_hypercube<dim>())
          {
            const Quadrature<1> dummy_quadrature(
              std::vector<Point<1>>(1, Point<1>()));
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering = shape_info.lexicographic_numbering;
          }
        else
          {
            const auto dummy_quadrature =
              reference_cell.template get_gauss_type_quadrature<dim>(1);
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering = shape_info.lexicographic_numbering;
          }

        // ------------------------------ indices ------------------------------
        unsigned int *level_dof_indices_coarse_0 =
          &transfer.schemes[0].level_dof_indices_coarse[0];
        unsigned int *level_dof_indices_fine_0 =
          &transfer.schemes[0].level_dof_indices_fine[0];

        unsigned int *level_dof_indices_coarse_1 =
          &transfer.schemes[1].level_dof_indices_coarse[0];
        unsigned int *level_dof_indices_fine_1 =
          &transfer.schemes[1].level_dof_indices_fine[0];

        process_cells(
          [&](const auto &cell_coarse, const auto &cell_fine) {
            // parent
            {
              cell_coarse->get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[0].dofs_per_cell_coarse;
                   i++)
                level_dof_indices_coarse_0[i] =
                  transfer.partitioner_coarse->global_to_local(
                    local_dof_indices[lexicographic_numbering[i]]);
            }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[0].dofs_per_cell_coarse;
                   i++)
                level_dof_indices_fine_0[i] =
                  transfer.partitioner_fine->global_to_local(
                    local_dof_indices[lexicographic_numbering[i]]);
            }

            // move pointers
            {
              level_dof_indices_coarse_0 +=
                transfer.schemes[0].dofs_per_cell_coarse;
              level_dof_indices_fine_0 +=
                transfer.schemes[0].dofs_per_cell_fine;
            }
          },
          [&](const auto &cell_coarse, const auto &cell_fine, const auto c) {
            // parent (only once at the beginning)
            if (c == 0)
              {
                cell_coarse->get_dof_indices(local_dof_indices);
                for (unsigned int i = 0;
                     i < transfer.schemes[1].dofs_per_cell_coarse;
                     i++)
                  level_dof_indices_coarse_1[i] =
                    transfer.partitioner_coarse->global_to_local(
                      local_dof_indices[lexicographic_numbering[i]]);
              }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[1].dofs_per_cell_coarse;
                   i++)
                level_dof_indices_fine_1[cell_local_chilren_indices[c][i]] =
                  transfer.partitioner_fine->global_to_local(
                    local_dof_indices[lexicographic_numbering[i]]);
            }

            // move pointers (only once at the end)
            if (c + 1 == GeometryInfo<dim>::max_children_per_cell)
              {
                level_dof_indices_coarse_1 +=
                  transfer.schemes[1].dofs_per_cell_coarse;
                level_dof_indices_fine_1 +=
                  transfer.schemes[1].dofs_per_cell_fine;
              }
          });
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

            std::vector<unsigned int> renumbering(fe->dofs_per_cell);
            {
              AssertIndexRange(fe->dofs_per_vertex, 2);
              renumbering[0] = 0;
              for (unsigned int i = 0; i < fe->dofs_per_line; ++i)
                renumbering[i + fe->dofs_per_vertex] =
                  GeometryInfo<1>::vertices_per_cell * fe->dofs_per_vertex + i;
              if (fe->dofs_per_vertex > 0)
                renumbering[fe->dofs_per_cell - fe->dofs_per_vertex] =
                  fe->dofs_per_vertex;
            }

            // TODO: data structures are saved in form of DG data structures
            // here
            const unsigned int shift           = fe->dofs_per_cell;
            const unsigned int n_child_dofs_1d = fe->dofs_per_cell * 2;

            {
              transfer.schemes[1].prolongation_matrix_1d.resize(
                fe->dofs_per_cell * n_child_dofs_1d);

              for (unsigned int c = 0;
                   c < GeometryInfo<1>::max_children_per_cell;
                   ++c)
                for (unsigned int i = 0; i < fe->dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
                    transfer.schemes[1]
                      .prolongation_matrix_1d[i * n_child_dofs_1d + j +
                                              c * shift] =
                      fe->get_prolongation_matrix(c)(renumbering[j],
                                                     renumbering[i]);
            }
            {
              transfer.schemes[1].restriction_matrix_1d.resize(
                fe->dofs_per_cell * n_child_dofs_1d);

              for (unsigned int c = 0;
                   c < GeometryInfo<1>::max_children_per_cell;
                   ++c)
                {
                  const auto matrix = get_restriction_matrix(*fe, c);
                  for (unsigned int i = 0; i < fe->dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
                      transfer.schemes[1]
                        .restriction_matrix_1d[i * n_child_dofs_1d + j +
                                               c * shift] =
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
                         j + c * n_dofs_per_cell] = matrix(i, j);
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
                          numbers::invalid_unsigned_int /*active level*/,
                          constraint_fine,
                          transfer,
                          weight_vector);

          // ... and store them cell-wise a DG format
          transfer.schemes[0].weights.resize(
            transfer.schemes[0].n_coarse_cells *
            transfer.schemes[0].dofs_per_cell_fine);
          transfer.schemes[1].weights.resize(
            transfer.schemes[1].n_coarse_cells *
            transfer.schemes[1].dofs_per_cell_fine);

          Number *      weights_0 = &transfer.schemes[0].weights[0];
          Number *      weights_1 = &transfer.schemes[1].weights[0];
          unsigned int *dof_indices_fine_0 =
            &transfer.schemes[0].level_dof_indices_fine[0];
          unsigned int *dof_indices_fine_1 =
            &transfer.schemes[1].level_dof_indices_fine[0];

          process_cells(
            [&](const auto &, const auto &) {
              for (unsigned int i = 0;
                   i < transfer.schemes[0].dofs_per_cell_fine;
                   i++)
                weights_0[i] =
                  weight_vector.local_element(dof_indices_fine_0[i]);

              dof_indices_fine_0 += transfer.schemes[0].dofs_per_cell_fine;
              weights_0 += transfer.schemes[0].dofs_per_cell_fine;
            },
            [&](const auto &, const auto &, const auto c) {
              for (unsigned int i = 0;
                   i < transfer.schemes[1].dofs_per_cell_coarse;
                   i++)
                weights_1[cell_local_chilren_indices[c][i]] =
                  weight_vector.local_element(
                    dof_indices_fine_1[cell_local_chilren_indices[c][i]]);

              // move pointers (only once at the end)
              if (c + 1 == GeometryInfo<dim>::max_children_per_cell)
                {
                  dof_indices_fine_1 += transfer.schemes[1].dofs_per_cell_fine;
                  weights_1 += transfer.schemes[1].dofs_per_cell_fine;
                }
            });
        }
    }



    template <int dim, typename Number>
    static void
    reinit_polynomial_transfer(
      const DoFHandler<dim> &                  dof_handler_fine,
      const DoFHandler<dim> &                  dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraint_fine,
      const dealii::AffineConstraints<Number> &constraint_coarse,
      const unsigned int                       mg_level_fine,
      const unsigned int                       mg_level_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      Assert(
        mg_level_fine == 0 || mg_level_fine == numbers::invalid_unsigned_int,
        ExcMessage(
          "Polynomial transfer is only allowed on the active level "
          "(numbers::invalid_unsigned_int) or the coarse-grid level (0)."));
      Assert(
        mg_level_coarse == 0 ||
          mg_level_coarse == numbers::invalid_unsigned_int,
        ExcMessage(
          "Polynomial transfer is only allowed on the active level "
          "(numbers::invalid_unsigned_int) or the coarse-grid level (0)."));

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
        dof_handler_fine.get_fe(0).dofs_per_vertex > 0;

      for (unsigned int i = 0; i < dof_handler_fine.get_fe_collection().size();
           ++i)
        Assert(transfer.fine_element_is_continuous ==
                 (dof_handler_fine.get_fe(i).dofs_per_vertex > 0),
               ExcInternalError());

      const auto process_cells = [&](const auto &fu) {
        if (mg_level_coarse == numbers::invalid_unsigned_int)
          {
            for (const auto &cell_coarse :
                 dof_handler_coarse.active_cell_iterators())
              {
                if (!cell_coarse->is_locally_owned())
                  continue;

                const auto cell_coarse_on_fine_mesh =
                  view.get_cell(cell_coarse);

                fu(cell_coarse, &cell_coarse_on_fine_mesh);
              }
          }
        else
          {
            for (const auto &cell_coarse :
                 dof_handler_coarse.mg_cell_iterators_on_level(mg_level_coarse))
              {
                if (cell_coarse->level_subdomain_id() !=
                    dof_handler_coarse.get_triangulation()
                      .locally_owned_subdomain())
                  continue;

                const auto cell_coarse_on_fine_mesh =
                  view.get_cell(cell_coarse);

                fu(cell_coarse, &cell_coarse_on_fine_mesh);
              }
          }
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
          transfer.schemes[fe_index_pair.second].dofs_per_cell_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first).dofs_per_cell;
          transfer.schemes[fe_index_pair.second].dofs_per_cell_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second).dofs_per_cell;

          transfer.schemes[fe_index_pair.second].degree_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first).degree;
          transfer.schemes[fe_index_pair.second].degree_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second).degree;
        }

      precompute_constraints(dof_handler_coarse, constraint_coarse, transfer);

      const auto comm = dof_handler_coarse.get_communicator();
      {
        transfer.partitioner_fine =
          std::make_shared<Utilities::MPI::Partitioner>(
            view.locally_owned_dofs(),
            view.locally_relevant_dofs(),
            dof_handler_fine.get_communicator());
        transfer.vec_fine.reinit(transfer.partitioner_fine);
      }
      {
        IndexSet locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler_coarse,
                                                locally_relevant_dofs);

        transfer.partitioner_coarse =
          std::make_shared<Utilities::MPI::Partitioner>(
            dof_handler_coarse.locally_owned_dofs(),
            locally_relevant_dofs,
            comm);
        transfer.vec_coarse.reinit(transfer.partitioner_coarse);
        precompute_restriction_constraints(transfer.partitioner_coarse,
                                           constraint_coarse,
                                           transfer);
      }

      {
        std::vector<std::vector<unsigned int>> lexicographic_numbering_fine(
          fe_index_pairs.size());
        std::vector<std::vector<unsigned int>> lexicographic_numbering_coarse(
          fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_coarse(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_fine(fe_index_pairs.size());

        for (const auto &fe_index_pair : fe_index_pairs)
          {
            local_dof_indices_coarse[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].dofs_per_cell_coarse);
            local_dof_indices_fine[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].dofs_per_cell_fine);

            transfer.schemes[fe_index_pair.second]
              .level_dof_indices_fine.resize(
                transfer.schemes[fe_index_pair.second].dofs_per_cell_fine *
                transfer.schemes[fe_index_pair.second].n_coarse_cells);
            transfer.schemes[fe_index_pair.second]
              .level_dof_indices_coarse.resize(
                transfer.schemes[fe_index_pair.second].dofs_per_cell_coarse *
                transfer.schemes[fe_index_pair.second].n_coarse_cells);

            const auto reference_cell =
              dof_handler_fine.get_fe(fe_index_pair.first.second)
                .reference_cell();

            Assert(reference_cell ==
                     dof_handler_coarse.get_fe(fe_index_pair.first.second)
                       .reference_cell(),
                   ExcNotImplemented());

            // ------------------- lexicographic_numbering  --------------------
            if (reference_cell == ReferenceCells::get_hypercube<dim>())
              {
                const Quadrature<1> dummy_quadrature(
                  std::vector<Point<1>>(1, Point<1>()));
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

        // ------------------------------ indices  -----------------------------
        std::vector<unsigned int *> level_dof_indices_coarse_(
          fe_index_pairs.size());
        std::vector<unsigned int *> level_dof_indices_fine_(
          fe_index_pairs.size());

        for (unsigned int i = 0; i < fe_index_pairs.size(); i++)
          {
            level_dof_indices_coarse_[i] =
              transfer.schemes[i].level_dof_indices_coarse.data();
            level_dof_indices_fine_[i] =
              transfer.schemes[i].level_dof_indices_fine.data();
          }

        bool     fine_indices_touch_remote_dofs = false;
        IndexSet locally_owned_dofs = dof_handler_fine.locally_owned_dofs();

        process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
          const auto fe_pair_no =
            fe_index_pairs[std::pair<unsigned int, unsigned int>(
              cell_coarse->active_fe_index(), cell_fine->active_fe_index())];

          if (mg_level_coarse == numbers::invalid_unsigned_int)
            cell_coarse->get_dof_indices(local_dof_indices_coarse[fe_pair_no]);
          else
            cell_coarse->get_mg_dof_indices(
              local_dof_indices_coarse[fe_pair_no]);

          for (unsigned int i = 0;
               i < transfer.schemes[fe_pair_no].dofs_per_cell_coarse;
               i++)
            level_dof_indices_coarse_[fe_pair_no][i] =
              transfer.partitioner_coarse->global_to_local(
                local_dof_indices_coarse
                  [fe_pair_no][lexicographic_numbering_coarse[fe_pair_no][i]]);

          cell_fine->get_dof_indices(local_dof_indices_fine[fe_pair_no]);
          for (unsigned int i = 0;
               i < transfer.schemes[fe_pair_no].dofs_per_cell_fine;
               i++)
            {
              level_dof_indices_fine_[fe_pair_no][i] =
                transfer.partitioner_fine->global_to_local(
                  local_dof_indices_fine
                    [fe_pair_no][lexicographic_numbering_fine[fe_pair_no][i]]);
              if (!locally_owned_dofs.is_element(
                    local_dof_indices_fine
                      [fe_pair_no]
                      [lexicographic_numbering_fine[fe_pair_no][i]]))
                fine_indices_touch_remote_dofs = true;
            }

          level_dof_indices_coarse_[fe_pair_no] +=
            transfer.schemes[fe_pair_no].dofs_per_cell_coarse;
          level_dof_indices_fine_[fe_pair_no] +=
            transfer.schemes[fe_pair_no].dofs_per_cell_fine;
        });

        // if all access goes to the locally owned dofs on all ranks, we do
        // not need the vec_fine vector
        if (Utilities::MPI::max(static_cast<unsigned int>(
                                  fine_indices_touch_remote_dofs),
                                comm) == 0)
          transfer.vec_fine.reinit(0);
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
                   dof_handler_coarse.get_fe(fe_index_pair_.first.second)
                     .reference_cell(),
                 ExcNotImplemented());

          if (reference_cell == ReferenceCells::get_hypercube<dim>())
            {
              const auto fe_fine = create_1D_fe(
                dof_handler_fine.get_fe(fe_index_pair.second).base_element(0));

              std::vector<unsigned int> renumbering_fine(
                fe_fine->dofs_per_cell);
              {
                AssertIndexRange(fe_fine->dofs_per_vertex, 2);
                renumbering_fine[0] = 0;
                for (unsigned int i = 0; i < fe_fine->dofs_per_line; ++i)
                  renumbering_fine[i + fe_fine->dofs_per_vertex] =
                    GeometryInfo<1>::vertices_per_cell *
                      fe_fine->dofs_per_vertex +
                    i;
                if (fe_fine->dofs_per_vertex > 0)
                  renumbering_fine[fe_fine->dofs_per_cell -
                                   fe_fine->dofs_per_vertex] =
                    fe_fine->dofs_per_vertex;
              }

              const auto fe_coarse = create_1D_fe(
                dof_handler_coarse.get_fe(fe_index_pair.first).base_element(0));

              std::vector<unsigned int> renumbering_coarse(
                fe_coarse->dofs_per_cell);
              {
                AssertIndexRange(fe_coarse->dofs_per_vertex, 2);
                renumbering_coarse[0] = 0;
                for (unsigned int i = 0; i < fe_coarse->dofs_per_line; ++i)
                  renumbering_coarse[i + fe_coarse->dofs_per_vertex] =
                    GeometryInfo<1>::vertices_per_cell *
                      fe_coarse->dofs_per_vertex +
                    i;
                if (fe_coarse->dofs_per_vertex > 0)
                  renumbering_coarse[fe_coarse->dofs_per_cell -
                                     fe_coarse->dofs_per_vertex] =
                    fe_coarse->dofs_per_vertex;
              }

              {
                FullMatrix<double> matrix(fe_fine->dofs_per_cell,
                                          fe_coarse->dofs_per_cell);
                FETools::get_projection_matrix(*fe_coarse, *fe_fine, matrix);
                transfer.schemes[fe_index_no].prolongation_matrix_1d.resize(
                  fe_fine->dofs_per_cell * fe_coarse->dofs_per_cell);

                for (unsigned int i = 0, k = 0; i < fe_coarse->dofs_per_cell;
                     ++i)
                  for (unsigned int j = 0; j < fe_fine->dofs_per_cell; ++j, ++k)
                    transfer.schemes[fe_index_no].prolongation_matrix_1d[k] =
                      matrix(renumbering_fine[j], renumbering_coarse[i]);
              }

              {
                FullMatrix<double> matrix(fe_coarse->dofs_per_cell,
                                          fe_fine->dofs_per_cell);
                FETools::get_projection_matrix(*fe_fine, *fe_coarse, matrix);
                transfer.schemes[fe_index_no].restriction_matrix_1d.resize(
                  fe_fine->dofs_per_cell * fe_coarse->dofs_per_cell);

                for (unsigned int i = 0, k = 0; i < fe_coarse->dofs_per_cell;
                     ++i)
                  for (unsigned int j = 0; j < fe_fine->dofs_per_cell; ++j, ++k)
                    transfer.schemes[fe_index_no].restriction_matrix_1d[k] =
                      matrix(renumbering_coarse[i], renumbering_fine[j]);
              }
            }
          else
            {
              const auto &fe_fine =
                dof_handler_fine.get_fe(fe_index_pair.second).base_element(0);

              const auto &fe_coarse =
                dof_handler_coarse.get_fe(fe_index_pair.first).base_element(0);

              {
                FullMatrix<double> matrix(fe_fine.dofs_per_cell,
                                          fe_coarse.dofs_per_cell);
                FETools::get_projection_matrix(fe_coarse, fe_fine, matrix);
                transfer.schemes[fe_index_no].prolongation_matrix.resize(
                  fe_fine.dofs_per_cell * fe_coarse.dofs_per_cell);

                for (unsigned int i = 0, k = 0; i < fe_coarse.dofs_per_cell;
                     ++i)
                  for (unsigned int j = 0; j < fe_fine.dofs_per_cell; ++j, ++k)
                    transfer.schemes[fe_index_no].prolongation_matrix[k] =
                      matrix(j, i);
              }

              {
                FullMatrix<double> matrix(fe_coarse.dofs_per_cell,
                                          fe_fine.dofs_per_cell);
                FETools::get_projection_matrix(fe_fine, fe_coarse, matrix);
                transfer.schemes[fe_index_no].restriction_matrix.resize(
                  fe_fine.dofs_per_cell * fe_coarse.dofs_per_cell);

                for (unsigned int i = 0, k = 0; i < fe_coarse.dofs_per_cell;
                     ++i)
                  for (unsigned int j = 0; j < fe_fine.dofs_per_cell; ++j, ++k)
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
                          constraint_fine,
                          transfer,
                          weight_vector);

          // ... and store them cell-wise a DG format
          for (auto &scheme : transfer.schemes)
            scheme.weights.resize(scheme.n_coarse_cells *
                                  scheme.dofs_per_cell_fine);

          std::vector<unsigned int *> level_dof_indices_fine_(
            fe_index_pairs.size());
          std::vector<Number *> weights_(fe_index_pairs.size());

          for (unsigned int i = 0; i < fe_index_pairs.size(); i++)
            {
              level_dof_indices_fine_[i] =
                &transfer.schemes[i].level_dof_indices_fine[0];
              weights_[i] = &transfer.schemes[i].weights[0];
            }

          process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
            const auto fe_pair_no =
              fe_index_pairs[std::pair<unsigned int, unsigned int>(
                cell_coarse->active_fe_index(), cell_fine->active_fe_index())];

            for (unsigned int i = 0;
                 i < transfer.schemes[fe_pair_no].dofs_per_cell_fine;
                 i++)
              weights_[fe_pair_no][i] = weight_vector.local_element(
                level_dof_indices_fine_[fe_pair_no][i]);

            level_dof_indices_fine_[fe_pair_no] +=
              transfer.schemes[fe_pair_no].dofs_per_cell_fine;
            weights_[fe_pair_no] +=
              transfer.schemes[fe_pair_no].dofs_per_cell_fine;
          });
        }
    }
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



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::prolongate(
  LinearAlgebra::distributed::Vector<Number> &      dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  using VectorizedArrayType = VectorizedArray<Number>;

  const unsigned int n_lanes = VectorizedArrayType::size();

  const bool use_dst_inplace = this->vec_fine.size() == 0;

  const auto vec_fine_ptr = use_dst_inplace ? &dst : &this->vec_fine;

  // the following code is equivalent to:
  // this->constraint_coarse.distribute(this->vec_coarse);
  {
    this->vec_coarse_constraints.copy_locally_owned_data_from(src);
    this->vec_coarse_constraints.update_ghost_values();

    for (unsigned int i = 0; i < constraint_coarse_distribute_ptr.size() - 1;
         ++i)
      if (constraint_coarse_distribute_ptr[i + 1] ==
          constraint_coarse_distribute_ptr[i])
        // not constrained entries
        vec_coarse.local_element(i) = vec_coarse_constraints.local_element(i);
      else
        {
          // constrained entries (ignoring homogeneous entries)
          Number val = 0.0;

          for (unsigned int j = constraint_coarse_distribute_ptr[i];
               j < constraint_coarse_distribute_ptr[i + 1];
               ++j)
            val += vec_coarse_constraints.local_element(
                     constraint_coarse_distribute_indices[j]) *
                   constraint_coarse_distribute_values[j];

          vec_coarse.local_element(i) = val;
        }
  }

  this->vec_coarse
    .update_ghost_values(); // note: make sure that ghost values are set

  if (fine_element_is_continuous)
    *vec_fine_ptr = Number(0.);

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  for (const auto &scheme : schemes)
    {
      // identity -> take short cut and work directly on global vectors
      if (scheme.degree_fine == scheme.degree_coarse)
        {
          const unsigned int *indices_fine =
            scheme.level_dof_indices_fine.data();
          const unsigned int *indices_coarse =
            scheme.level_dof_indices_coarse.data();
          const Number *weights = nullptr;

          if (fine_element_is_continuous)
            weights = scheme.weights.data();

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
            {
              if (fine_element_is_continuous)
                for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                  vec_fine_ptr->local_element(indices_fine[i]) +=
                    this->vec_coarse.local_element(indices_coarse[i]) *
                    weights[i];
              else
                for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                  vec_fine_ptr->local_element(indices_fine[i]) =
                    this->vec_coarse.local_element(indices_coarse[i]);

              indices_fine += scheme.dofs_per_cell_fine;
              indices_coarse += scheme.dofs_per_cell_coarse;

              if (fine_element_is_continuous)
                weights += scheme.dofs_per_cell_fine;
            }

          continue;
        }

      // general case -> local restriction is needed
      evaluation_data_fine.resize(scheme.dofs_per_cell_fine);
      evaluation_data_coarse.resize(scheme.dofs_per_cell_fine);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from source vector
          {
            const unsigned int *indices =
              &scheme
                 .level_dof_indices_coarse[cell * scheme.dofs_per_cell_coarse];
            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                for (unsigned int i = 0; i < scheme.dofs_per_cell_coarse; ++i)
                  evaluation_data_coarse[i][v] =
                    this->vec_coarse.local_element(indices[i]);
                indices += scheme.dofs_per_cell_coarse;
              }
          }

          // ---------------------------- coarse -----------------------------
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
          // ------------------------------ fine -----------------------------

          // weight and write into dst vector
          {
            const unsigned int *indices =
              &scheme.level_dof_indices_fine[cell * scheme.dofs_per_cell_fine];
            const Number *weights = nullptr;

            if (fine_element_is_continuous)
              weights = &scheme.weights[cell * scheme.dofs_per_cell_fine];

            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                if (fine_element_is_continuous)
                  for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                    vec_fine_ptr->local_element(indices[i]) +=
                      evaluation_data_fine[i][v] * weights[i];
                else
                  for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                    vec_fine_ptr->local_element(indices[i]) =
                      evaluation_data_fine[i][v];

                indices += scheme.dofs_per_cell_fine;

                if (fine_element_is_continuous)
                  weights += scheme.dofs_per_cell_fine;
              }
          }
        }
    }

  this->vec_coarse.zero_out_ghost_values(); // clear ghost values; else compress
                                            // in do_restrict_add does not work

  if (fine_element_is_continuous || use_dst_inplace == false)
    vec_fine_ptr->compress(VectorOperation::add);

  if (use_dst_inplace == false)
    dst.copy_locally_owned_data_from(this->vec_fine);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  restrict_and_add(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const
{
  using VectorizedArrayType = VectorizedArray<Number>;

  const unsigned int n_lanes = VectorizedArrayType::size();

  const bool use_src_inplace = this->vec_fine.size() == 0;
  const auto vec_fine_ptr    = use_src_inplace ? &src : &this->vec_fine;

  if (use_src_inplace == false)
    this->vec_fine.copy_locally_owned_data_from(src);

  if (fine_element_is_continuous || use_src_inplace == false)
    vec_fine_ptr->update_ghost_values();

  this->vec_coarse.copy_locally_owned_data_from(dst);

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  // a helper function similar to AffineConstraints::distribute_local_to_global
  // but working with local indices
  const auto distribute_local_to_global =
    [&](const auto &index, const auto &value, auto &global_vector) {
      if (distribute_local_to_global_ptr[index + 1] ==
          distribute_local_to_global_ptr[index])
        global_vector.local_element(index) += value;
      else if (((distribute_local_to_global_ptr[index + 1] -
                 distribute_local_to_global_ptr[index]) == 1 &&
                distribute_local_to_global_indices
                    [distribute_local_to_global_ptr[index]] ==
                  numbers::invalid_unsigned_int) == false)
        for (unsigned int j = distribute_local_to_global_ptr[index];
             j < distribute_local_to_global_ptr[index + 1];
             ++j)
          global_vector.local_element(distribute_local_to_global_indices[j]) +=
            value * distribute_local_to_global_values[j];
    };

  for (const auto &scheme : schemes)
    {
      // identity -> take short cut and work directly on global vectors
      if (scheme.degree_fine == scheme.degree_coarse)
        {
          const unsigned int *indices_fine =
            scheme.level_dof_indices_fine.data();
          const unsigned int *indices_coarse =
            scheme.level_dof_indices_coarse.data();
          const Number *weights = nullptr;

          if (fine_element_is_continuous)
            weights = scheme.weights.data();

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
            {
              if (fine_element_is_continuous)
                for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                  distribute_local_to_global(
                    indices_coarse[i],
                    vec_fine_ptr->local_element(indices_fine[i]) * weights[i],
                    this->vec_coarse);
              else
                for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                  distribute_local_to_global(indices_coarse[i],
                                             vec_fine_ptr->local_element(
                                               indices_fine[i]),
                                             this->vec_coarse);

              indices_fine += scheme.dofs_per_cell_fine;
              indices_coarse += scheme.dofs_per_cell_coarse;

              if (fine_element_is_continuous)
                weights += scheme.dofs_per_cell_fine;
            }

          continue;
        }

      // general case -> local restriction is needed
      evaluation_data_fine.resize(scheme.dofs_per_cell_fine);
      evaluation_data_coarse.resize(scheme.dofs_per_cell_fine);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from source vector and weight
          {
            const unsigned int *indices =
              &scheme.level_dof_indices_fine[cell * scheme.dofs_per_cell_fine];
            const Number *weights = nullptr;

            if (fine_element_is_continuous)
              weights = &scheme.weights[cell * scheme.dofs_per_cell_fine];

            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                if (fine_element_is_continuous)
                  for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                    evaluation_data_fine[i][v] =
                      vec_fine_ptr->local_element(indices[i]) * weights[i];
                else
                  for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                    evaluation_data_fine[i][v] =
                      vec_fine_ptr->local_element(indices[i]);

                indices += scheme.dofs_per_cell_fine;

                if (fine_element_is_continuous)
                  weights += scheme.dofs_per_cell_fine;
              }
          }

          // ------------------------------ fine -----------------------------
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
          // ----------------------------- coarse ----------------------------

          // write into dst vector
          {
            const unsigned int *indices =
              &scheme
                 .level_dof_indices_coarse[cell * scheme.dofs_per_cell_coarse];
            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                for (unsigned int i = 0; i < scheme.dofs_per_cell_coarse; ++i)
                  distribute_local_to_global(indices[i],
                                             evaluation_data_coarse[i][v],
                                             this->vec_coarse);
                indices += scheme.dofs_per_cell_coarse;
              }
          }
        }
    }

  this->vec_coarse.compress(VectorOperation::add);

  dst.copy_locally_owned_data_from(this->vec_coarse);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  interpolate(LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const
{
  using VectorizedArrayType = VectorizedArray<Number>;

  const unsigned int n_lanes = VectorizedArrayType::size();

  const bool use_src_inplace = this->vec_fine.size() == 0;
  const auto vec_fine_ptr    = use_src_inplace ? &src : &this->vec_fine;

  if (use_src_inplace == false)
    this->vec_fine.copy_locally_owned_data_from(src);

  if (fine_element_is_continuous || use_src_inplace == false)
    vec_fine_ptr->update_ghost_values();

  this->vec_coarse = 0.0;

  AlignedVector<VectorizedArrayType> evaluation_data_fine;
  AlignedVector<VectorizedArrayType> evaluation_data_coarse;

  for (const auto &scheme : schemes)
    {
      // identity -> take short cut and work directly on global vectors
      if (scheme.degree_fine == scheme.degree_coarse)
        {
          const unsigned int *indices_fine =
            scheme.level_dof_indices_fine.data();
          const unsigned int *indices_coarse =
            scheme.level_dof_indices_coarse.data();

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
            {
              for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                this->vec_coarse.local_element(indices_coarse[i]) =
                  vec_fine_ptr->local_element(indices_fine[i]);

              indices_fine += scheme.dofs_per_cell_fine;
              indices_coarse += scheme.dofs_per_cell_coarse;
            }

          continue;
        }

      // general case -> local restriction is needed
      evaluation_data_fine.resize(scheme.dofs_per_cell_fine);
      evaluation_data_coarse.resize(scheme.dofs_per_cell_fine);

      CellTransferFactory cell_transfer(scheme.degree_fine,
                                        scheme.degree_coarse);

      const unsigned int n_scalar_dofs_fine =
        scheme.dofs_per_cell_fine / n_components;
      const unsigned int n_scalar_dofs_coarse =
        scheme.dofs_per_cell_coarse / n_components;

      for (unsigned int cell = 0; cell < scheme.n_coarse_cells; cell += n_lanes)
        {
          const unsigned int n_lanes_filled =
            (cell + n_lanes > scheme.n_coarse_cells) ?
              (scheme.n_coarse_cells - cell) :
              n_lanes;

          // read from source vector and weight
          {
            const unsigned int *indices =
              &scheme.level_dof_indices_fine[cell * scheme.dofs_per_cell_fine];

            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                for (unsigned int i = 0; i < scheme.dofs_per_cell_fine; ++i)
                  evaluation_data_fine[i][v] =
                    vec_fine_ptr->local_element(indices[i]);

                indices += scheme.dofs_per_cell_fine;
              }
          }

          // ------------------------------ fine -----------------------------
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
          // ----------------------------- coarse ----------------------------

          // write into dst vector
          {
            const unsigned int *indices =
              &scheme
                 .level_dof_indices_coarse[cell * scheme.dofs_per_cell_coarse];
            for (unsigned int v = 0; v < n_lanes_filled; ++v)
              {
                for (unsigned int i = 0; i < scheme.dofs_per_cell_coarse; ++i)
                  this->vec_coarse.local_element(indices[i]) =
                    evaluation_data_coarse[i][v];
                indices += scheme.dofs_per_cell_coarse;
              }
          }
        }
    }

  dst.copy_locally_owned_data_from(this->vec_coarse);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  reinit_geometric_transfer(const DoFHandler<dim> &          dof_handler_fine,
                            const DoFHandler<dim> &          dof_handler_coarse,
                            const AffineConstraints<Number> &constraint_fine,
                            const AffineConstraints<Number> &constraint_coarse)
{
  internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraint_fine,
    constraint_coarse,
    *this);
}



template <int dim, typename Number>
void
MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>::
  reinit_polynomial_transfer(const DoFHandler<dim> &dof_handler_fine,
                             const DoFHandler<dim> &dof_handler_coarse,
                             const AffineConstraints<Number> &constraint_fine,
                             const AffineConstraints<Number> &constraint_coarse,
                             const unsigned int               mg_level_fine,
                             const unsigned int               mg_level_coarse)
{
  internal::MGTwoLevelTransferImplementation::reinit_polynomial_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraint_fine,
    constraint_coarse,
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

DEAL_II_NAMESPACE_CLOSE

#endif
