// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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


#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_accessor.h>
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
MGLevelGlobalTransfer<VectorType>::fill_and_communicate_copy_indices
(const DoFHandler<dim,spacedim> &mg_dof)
{
  internal::MGTransfer::fill_copy_indices(mg_dof, mg_constrained_dofs, copy_indices,
                                          copy_indices_global_mine, copy_indices_level_mine);

  // check if we can run a plain copy operation between the global DoFs and
  // the finest level.
  perform_plain_copy =
    (copy_indices.back().size() == mg_dof.locally_owned_dofs().n_elements())
    &&
    (mg_dof.locally_owned_dofs().n_elements() ==
     mg_dof.locally_owned_mg_dofs(mg_dof.get_triangulation().n_global_levels()-1).n_elements());
  if (perform_plain_copy)
    {
      AssertDimension(copy_indices_global_mine.back().size(), 0);
      AssertDimension(copy_indices_level_mine.back().size(), 0);

      // check whether there is a renumbering of degrees of freedom on
      // either the finest level or the global dofs, which means that we
      // cannot apply a plain copy
      for (unsigned int i=0; i<copy_indices.back().size(); ++i)
        if (copy_indices.back()[i].first != copy_indices.back()[i].second)
          {
            perform_plain_copy = false;
            break;
          }
    }
  const parallel::Triangulation<dim, spacedim> *ptria =
    dynamic_cast<const parallel::Triangulation<dim, spacedim> *>
    (&mg_dof.get_tria());
  const MPI_Comm mpi_communicator = ptria != nullptr ? ptria->get_communicator() :
                                    MPI_COMM_SELF;
  perform_plain_copy =
    Utilities::MPI::min(static_cast<int>(perform_plain_copy),
                        mpi_communicator);

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
}



template <typename VectorType>
void
MGLevelGlobalTransfer<VectorType>::print_indices (std::ostream &os) const
{
  for (unsigned int level = 0; level<copy_indices.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices[level].size(); ++i)
        os << "copy_indices[" << level
           << "]\t" << copy_indices[level][i].first << '\t' << copy_indices[level][i].second << std::endl;
    }

  for (unsigned int level = 0; level<copy_indices_level_mine.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_level_mine[level].size(); ++i)
        os << "copy_ifrom  [" << level
           << "]\t" << copy_indices_level_mine[level][i].first << '\t' << copy_indices_level_mine[level][i].second << std::endl;
    }
  for (unsigned int level = 0; level<copy_indices_global_mine.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_global_mine[level].size(); ++i)
        os << "copy_ito    [" << level
           << "]\t" << copy_indices_global_mine[level][i].first << '\t' << copy_indices_global_mine[level][i].second << std::endl;
    }
}



template <typename VectorType>
std::size_t
MGLevelGlobalTransfer<VectorType>::memory_consumption () const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(sizes);
  result += MemoryConsumption::memory_consumption(copy_indices);
  result += MemoryConsumption::memory_consumption(copy_indices_global_mine);
  result += MemoryConsumption::memory_consumption(copy_indices_level_mine);

  return result;
}



/* ------------------ MGLevelGlobalTransfer<VectorType> ----------------- */


template <typename Number>
template <int dim, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::fill_and_communicate_copy_indices
(const DoFHandler<dim,spacedim> &mg_dof)
{
  // first go to the usual routine...
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  my_copy_indices;
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  my_copy_indices_global_mine;
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  my_copy_indices_level_mine;

  internal::MGTransfer::fill_copy_indices(mg_dof, mg_constrained_dofs, my_copy_indices,
                                          my_copy_indices_global_mine, my_copy_indices_level_mine);

  // get all degrees of freedom that we need read access to in copy_to_mg
  // and copy_from_mg, respectively. We fill an IndexSet once on each level
  // (for the global_mine indices accessing remote level indices) and once
  // globally (for the level_mine indices accessing remote global indices).

  // the variables index_set and level_index_set are going to define the
  // ghost indices of the respective vectors (due to construction, these are
  // precisely the indices that we need)
  const parallel::Triangulation<dim, spacedim> *ptria =
    dynamic_cast<const parallel::Triangulation<dim, spacedim> *>
    (&mg_dof.get_tria());
  const MPI_Comm mpi_communicator = ptria != nullptr ? ptria->get_communicator() :
                                    MPI_COMM_SELF;

  IndexSet index_set(mg_dof.locally_owned_dofs().size());
  std::vector<types::global_dof_index> accessed_indices;
  ghosted_level_vector.resize(0, mg_dof.get_triangulation().n_global_levels()-1);
  std::vector<IndexSet> level_index_set(mg_dof.get_triangulation().n_global_levels());
  for (unsigned int l=0; l<mg_dof.get_triangulation().n_global_levels(); ++l)
    {
      for (unsigned int i=0; i<my_copy_indices_level_mine[l].size(); ++i)
        accessed_indices.push_back(my_copy_indices_level_mine[l][i].first);
      std::vector<types::global_dof_index> accessed_level_indices;
      for (unsigned int i=0; i<my_copy_indices_global_mine[l].size(); ++i)
        accessed_level_indices.push_back(my_copy_indices_global_mine[l][i].second);
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
  this->copy_indices.resize(mg_dof.get_triangulation().n_global_levels());
  this->copy_indices_level_mine.resize(mg_dof.get_triangulation().n_global_levels());
  this->copy_indices_global_mine.resize(mg_dof.get_triangulation().n_global_levels());
  for (unsigned int level=0; level<mg_dof.get_triangulation().n_global_levels(); ++level)
    {
      const Utilities::MPI::Partitioner &global_partitioner =
        *ghosted_global_vector.get_partitioner();
      const Utilities::MPI::Partitioner &level_partitioner =
        *ghosted_level_vector[level].get_partitioner();
      // owned-owned case: the locally owned indices are going to control
      // the local index
      this->copy_indices[level].resize(my_copy_indices[level].size());
      for (unsigned int i=0; i<my_copy_indices[level].size(); ++i)
        this->copy_indices[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(my_copy_indices[level][i].first),
           level_partitioner.global_to_local(my_copy_indices[level][i].second));

      // remote-owned case: the locally owned indices for the level and the
      // ghost dofs for the global indices set the local index
      this->copy_indices_level_mine[level].
      resize(my_copy_indices_level_mine[level].size());
      for (unsigned int i=0; i<my_copy_indices_level_mine[level].size(); ++i)
        this->copy_indices_level_mine[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(my_copy_indices_level_mine[level][i].first),
           level_partitioner.global_to_local(my_copy_indices_level_mine[level][i].second));

      // owned-remote case: the locally owned indices for the global dofs
      // and the ghost dofs for the level indices set the local index
      this->copy_indices_global_mine[level].
      resize(my_copy_indices_global_mine[level].size());
      for (unsigned int i=0; i<my_copy_indices_global_mine[level].size(); ++i)
        this->copy_indices_global_mine[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(my_copy_indices_global_mine[level][i].first),
           level_partitioner.global_to_local(my_copy_indices_global_mine[level][i].second));
    }

  perform_plain_copy = this->copy_indices.back().size()
                       == mg_dof.locally_owned_dofs().n_elements();
  if (perform_plain_copy)
    {
      AssertDimension(this->copy_indices_global_mine.back().size(), 0);
      AssertDimension(this->copy_indices_level_mine.back().size(), 0);

      // check whether there is a renumbering of degrees of freedom on
      // either the finest level or the global dofs, which means that we
      // cannot apply a plain copy
      for (unsigned int i=0; i<this->copy_indices.back().size(); ++i)
        if (this->copy_indices.back()[i].first !=
            this->copy_indices.back()[i].second)
          {
            perform_plain_copy = false;
            break;
          }
    }
  perform_plain_copy =
    Utilities::MPI::min(static_cast<int>(perform_plain_copy),
                        mpi_communicator);

  // if we do a plain copy, no need to hold additional ghosted vectors
  if (perform_plain_copy)
    {
      ghosted_global_vector.reinit(0);
      ghosted_level_vector.resize(0, 0);
    }
}



template <typename Number>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::clear()
{
  sizes.resize(0);
  copy_indices.clear();
  copy_indices_global_mine.clear();
  copy_indices_level_mine.clear();
  component_to_block_map.resize(0);
  mg_constrained_dofs = nullptr;
  ghosted_global_vector.reinit(0);
  ghosted_level_vector.resize(0, 0);
}



template <typename Number>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::print_indices (std::ostream &os) const
{
  for (unsigned int level = 0; level<copy_indices.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices[level].size(); ++i)
        os << "copy_indices[" << level
           << "]\t" << copy_indices[level][i].first << '\t' << copy_indices[level][i].second << std::endl;
    }

  for (unsigned int level = 0; level<copy_indices_level_mine.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_level_mine[level].size(); ++i)
        os << "copy_ifrom  [" << level
           << "]\t" << copy_indices_level_mine[level][i].first << '\t' << copy_indices_level_mine[level][i].second << std::endl;
    }
  for (unsigned int level = 0; level<copy_indices_global_mine.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_global_mine[level].size(); ++i)
        os << "copy_ito    [" << level
           << "]\t" << copy_indices_global_mine[level][i].first << '\t' << copy_indices_global_mine[level][i].second << std::endl;
    }
}



template <typename Number>
std::size_t
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::memory_consumption () const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(sizes);
  result += MemoryConsumption::memory_consumption(copy_indices);
  result += MemoryConsumption::memory_consumption(copy_indices_global_mine);
  result += MemoryConsumption::memory_consumption(copy_indices_level_mine);
  result += ghosted_global_vector.memory_consumption();
  for (unsigned int i=ghosted_level_vector.min_level();
       i<=ghosted_level_vector.max_level(); ++i)
    result += ghosted_level_vector[i].memory_consumption();

  return result;
}



// explicit instantiation
#include "mg_level_global_transfer.inst"

// create an additional instantiation currently not supported by the automatic
// template instantiation scheme
template class MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<float> >;


DEAL_II_NAMESPACE_CLOSE
