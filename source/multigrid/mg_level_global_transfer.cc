// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2019 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer.templates.h>
#include <deal.II/multigrid/mg_transfer_internal.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


/* ------------------ MGLevelGlobalTransfer<VectorType> ----------------- */


template <typename VectorType>
template <int dim, int spacedim>
void
MGLevelGlobalTransfer<VectorType>::fill_and_communicate_copy_indices(
  const DoFHandler<dim, spacedim> &mg_dof)
{
  internal::MGTransfer::fill_copy_indices(mg_dof,
                                          mg_constrained_dofs,
                                          copy_indices,
                                          copy_indices_global_mine,
                                          copy_indices_level_mine);

  // check if we can run a plain copy operation between the global DoFs and
  // the finest level.
  bool my_perform_plain_copy =
    (copy_indices.back().size() == mg_dof.locally_owned_dofs().n_elements()) &&
    (mg_dof.locally_owned_dofs().n_elements() ==
     mg_dof
       .locally_owned_mg_dofs(mg_dof.get_triangulation().n_global_levels() - 1)
       .n_elements());
  if (my_perform_plain_copy)
    {
      AssertDimension(copy_indices_global_mine.back().size(), 0);
      AssertDimension(copy_indices_level_mine.back().size(), 0);

      // check whether there is a renumbering of degrees of freedom on
      // either the finest level or the global dofs, which means that we
      // cannot apply a plain copy
      for (unsigned int i = 0; i < copy_indices.back().size(); ++i)
        if (copy_indices.back()[i].first != copy_indices.back()[i].second)
          {
            my_perform_plain_copy = false;
            break;
          }
    }

  // now do a global reduction over all processors to see what operation
  // they can agree upon
  if (const parallel::TriangulationBase<dim, spacedim> *ptria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &mg_dof.get_triangulation()))
    perform_plain_copy = (Utilities::MPI::min(my_perform_plain_copy ? 1 : 0,
                                              ptria->get_communicator()) == 1);
  else
    perform_plain_copy = my_perform_plain_copy;
}



template <typename VectorType>
void
MGLevelGlobalTransfer<VectorType>::clear()
{
  sizes.resize(0);
  copy_indices.clear();
  copy_indices_global_mine.clear();
  copy_indices_level_mine.clear();
  component_to_block_map.resize(0);
  mg_constrained_dofs = nullptr;
  perform_plain_copy  = false;
}



template <typename VectorType>
void
MGLevelGlobalTransfer<VectorType>::print_indices(std::ostream &os) const
{
  for (unsigned int level = 0; level < copy_indices.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices[level].size(); ++i)
        os << "copy_indices[" << level << "]\t" << copy_indices[level][i].first
           << '\t' << copy_indices[level][i].second << std::endl;
    }

  for (unsigned int level = 0; level < copy_indices_level_mine.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices_level_mine[level].size(); ++i)
        os << "copy_ifrom  [" << level << "]\t"
           << copy_indices_level_mine[level][i].first << '\t'
           << copy_indices_level_mine[level][i].second << std::endl;
    }
  for (unsigned int level = 0; level < copy_indices_global_mine.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices_global_mine[level].size(); ++i)
        os << "copy_ito    [" << level << "]\t"
           << copy_indices_global_mine[level][i].first << '\t'
           << copy_indices_global_mine[level][i].second << std::endl;
    }
}



template <typename VectorType>
std::size_t
MGLevelGlobalTransfer<VectorType>::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(sizes);
  result += MemoryConsumption::memory_consumption(copy_indices);
  result += MemoryConsumption::memory_consumption(copy_indices_global_mine);
  result += MemoryConsumption::memory_consumption(copy_indices_level_mine);

  return result;
}



/* ------------------ MGLevelGlobalTransfer<VectorType> ----------------- */

namespace
{
  template <int dim, int spacedim, typename Number>
  void
  fill_internal(
    const DoFHandler<dim, spacedim> &mg_dof,
    SmartPointer<
      const MGConstrainedDoFs,
      MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>>
                                                mg_constrained_dofs,
    const MPI_Comm                              mpi_communicator,
    const bool                                  transfer_solution_vectors,
    std::vector<Table<2, unsigned int>> &       copy_indices,
    std::vector<Table<2, unsigned int>> &       copy_indices_global_mine,
    std::vector<Table<2, unsigned int>> &       copy_indices_level_mine,
    LinearAlgebra::distributed::Vector<Number> &ghosted_global_vector,
    MGLevelObject<LinearAlgebra::distributed::Vector<Number>>
      &ghosted_level_vector)
  {
    // first go to the usual routine...
    std::vector<
      std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
      my_copy_indices;
    std::vector<
      std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
      my_copy_indices_global_mine;
    std::vector<
      std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
      my_copy_indices_level_mine;

    internal::MGTransfer::fill_copy_indices(mg_dof,
                                            mg_constrained_dofs,
                                            my_copy_indices,
                                            my_copy_indices_global_mine,
                                            my_copy_indices_level_mine,
                                            !transfer_solution_vectors);

    // get all degrees of freedom that we need read access to in copy_to_mg
    // and copy_from_mg, respectively. We fill an IndexSet once on each level
    // (for the global_mine indices accessing remote level indices) and once
    // globally (for the level_mine indices accessing remote global indices).

    // the variables index_set and level_index_set are going to define the
    // ghost indices of the respective vectors (due to construction, these are
    // precisely the indices that we need)

    IndexSet index_set(mg_dof.locally_owned_dofs().size());
    std::vector<types::global_dof_index> accessed_indices;
    ghosted_level_vector.resize(0,
                                mg_dof.get_triangulation().n_global_levels() -
                                  1);
    std::vector<IndexSet> level_index_set(
      mg_dof.get_triangulation().n_global_levels());
    for (unsigned int l = 0; l < mg_dof.get_triangulation().n_global_levels();
         ++l)
      {
        for (const auto &indices : my_copy_indices_level_mine[l])
          accessed_indices.push_back(indices.first);
        std::vector<types::global_dof_index> accessed_level_indices;
        for (const auto &indices : my_copy_indices_global_mine[l])
          accessed_level_indices.push_back(indices.second);
        std::sort(accessed_level_indices.begin(), accessed_level_indices.end());
        level_index_set[l].set_size(mg_dof.locally_owned_mg_dofs(l).size());
        level_index_set[l].add_indices(accessed_level_indices.begin(),
                                       accessed_level_indices.end());
        level_index_set[l].compress();
        ghosted_level_vector[l].reinit(mg_dof.locally_owned_mg_dofs(l),
                                       level_index_set[l],
                                       mpi_communicator);
      }
    std::sort(accessed_indices.begin(), accessed_indices.end());
    index_set.add_indices(accessed_indices.begin(), accessed_indices.end());
    index_set.compress();
    ghosted_global_vector.reinit(mg_dof.locally_owned_dofs(),
                                 index_set,
                                 mpi_communicator);

    // localize the copy indices for faster access. Since all access will be
    // through the ghosted vector in 'data', we can use this (much faster)
    // option
    copy_indices.resize(mg_dof.get_triangulation().n_global_levels());
    copy_indices_level_mine.resize(
      mg_dof.get_triangulation().n_global_levels());
    copy_indices_global_mine.resize(
      mg_dof.get_triangulation().n_global_levels());
    for (unsigned int level = 0;
         level < mg_dof.get_triangulation().n_global_levels();
         ++level)
      {
        const Utilities::MPI::Partitioner &global_partitioner =
          *ghosted_global_vector.get_partitioner();
        const Utilities::MPI::Partitioner &level_partitioner =
          *ghosted_level_vector[level].get_partitioner();

        auto translate_indices =
          [&](const std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>
                &                     global_copy_indices,
              Table<2, unsigned int> &local_copy_indices) {
            local_copy_indices.reinit(2, global_copy_indices.size());
            for (unsigned int i = 0; i < global_copy_indices.size(); ++i)
              {
                local_copy_indices(0, i) = global_partitioner.global_to_local(
                  global_copy_indices[i].first);
                local_copy_indices(1, i) = level_partitioner.global_to_local(
                  global_copy_indices[i].second);
              }
          };

        // owned-owned case
        translate_indices(my_copy_indices[level], copy_indices[level]);

        // remote-owned case
        translate_indices(my_copy_indices_level_mine[level],
                          copy_indices_level_mine[level]);

        // owned-remote case
        translate_indices(my_copy_indices_global_mine[level],
                          copy_indices_global_mine[level]);
      }
  }
} // namespace

template <typename Number>
template <int dim, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>::
  fill_and_communicate_copy_indices(const DoFHandler<dim, spacedim> &mg_dof)
{
  const parallel::TriangulationBase<dim, spacedim> *ptria =
    dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
      &mg_dof.get_triangulation());
  const MPI_Comm mpi_communicator =
    ptria != nullptr ? ptria->get_communicator() : MPI_COMM_SELF;

  fill_internal(mg_dof,
                mg_constrained_dofs,
                mpi_communicator,
                false,
                this->copy_indices,
                this->copy_indices_global_mine,
                this->copy_indices_level_mine,
                ghosted_global_vector,
                ghosted_level_vector);

  fill_internal(mg_dof,
                mg_constrained_dofs,
                mpi_communicator,
                true,
                this->solution_copy_indices,
                this->solution_copy_indices_global_mine,
                this->solution_copy_indices_level_mine,
                solution_ghosted_global_vector,
                solution_ghosted_level_vector);

  bool my_perform_renumbered_plain_copy =
    (this->copy_indices.back().n_cols() ==
     mg_dof.locally_owned_dofs().n_elements());
  bool my_perform_plain_copy = false;
  if (my_perform_renumbered_plain_copy)
    {
      my_perform_plain_copy = true;
      AssertDimension(this->copy_indices_global_mine.back().n_rows(), 0);
      AssertDimension(this->copy_indices_level_mine.back().n_rows(), 0);

      // check whether there is a renumbering of degrees of freedom on
      // either the finest level or the global dofs, which means that we
      // cannot apply a plain copy
      for (unsigned int i = 0; i < this->copy_indices.back().n_cols(); ++i)
        if (this->copy_indices.back()(0, i) != this->copy_indices.back()(1, i))
          {
            my_perform_plain_copy = false;
            break;
          }
    }

  // now do a global reduction over all processors to see what operation
  // they can agree upon
  perform_plain_copy =
    Utilities::MPI::min(static_cast<int>(my_perform_plain_copy),
                        mpi_communicator);
  perform_renumbered_plain_copy =
    Utilities::MPI::min(static_cast<int>(my_perform_renumbered_plain_copy),
                        mpi_communicator);

  // if we do a plain copy, no need to hold additional ghosted vectors
  if (perform_renumbered_plain_copy)
    {
      for (unsigned int i = 0; i < this->copy_indices.back().n_cols(); ++i)
        AssertDimension(this->copy_indices.back()(0, i), i);

      ghosted_global_vector.reinit(0);
      ghosted_level_vector.resize(0, 0);
      solution_ghosted_global_vector.reinit(0);
      solution_ghosted_level_vector.resize(0, 0);
    }
}



template <typename Number>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>::clear()
{
  sizes.resize(0);
  copy_indices.clear();
  copy_indices_global_mine.clear();
  copy_indices_level_mine.clear();
  component_to_block_map.resize(0);
  mg_constrained_dofs = nullptr;
  ghosted_global_vector.reinit(0);
  ghosted_level_vector.resize(0, 0);
  perform_plain_copy            = false;
  perform_renumbered_plain_copy = false;
}



template <typename Number>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>::
  print_indices(std::ostream &os) const
{
  for (unsigned int level = 0; level < copy_indices.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices[level].n_cols(); ++i)
        os << "copy_indices[" << level << "]\t" << copy_indices[level](0, i)
           << '\t' << copy_indices[level](1, i) << std::endl;
    }

  for (unsigned int level = 0; level < copy_indices_level_mine.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices_level_mine[level].n_cols(); ++i)
        os << "copy_ifrom  [" << level << "]\t"
           << copy_indices_level_mine[level](0, i) << '\t'
           << copy_indices_level_mine[level](1, i) << std::endl;
    }
  for (unsigned int level = 0; level < copy_indices_global_mine.size(); ++level)
    {
      for (unsigned int i = 0; i < copy_indices_global_mine[level].n_cols();
           ++i)
        os << "copy_ito    [" << level << "]\t"
           << copy_indices_global_mine[level](0, i) << '\t'
           << copy_indices_global_mine[level](1, i) << std::endl;
    }
}



template <typename Number>
std::size_t
MGLevelGlobalTransfer<
  LinearAlgebra::distributed::Vector<Number>>::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(sizes);
  result += MemoryConsumption::memory_consumption(copy_indices);
  result += MemoryConsumption::memory_consumption(copy_indices_global_mine);
  result += MemoryConsumption::memory_consumption(copy_indices_level_mine);
  result += ghosted_global_vector.memory_consumption();
  for (unsigned int i = ghosted_level_vector.min_level();
       i <= ghosted_level_vector.max_level();
       ++i)
    result += ghosted_level_vector[i].memory_consumption();

  return result;
}



// explicit instantiation
#include "mg_level_global_transfer.inst"

// create an additional instantiation currently not supported by the automatic
// template instantiation scheme
template class MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<float>>;


DEAL_II_NAMESPACE_CLOSE
