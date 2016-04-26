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
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer.templates.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * Internal data structure that is used in the MPI communication in fill_and_communicate_copy_indices().
   * It represents an entry in the copy_indices* map, that associates a level dof index with a global dof index.
   */
  struct DoFPair
  {
    unsigned int level;
    types::global_dof_index global_dof_index;
    types::global_dof_index level_dof_index;

    DoFPair(const unsigned int level,
            const types::global_dof_index global_dof_index,
            const types::global_dof_index level_dof_index)
      :
      level(level), global_dof_index(global_dof_index), level_dof_index(level_dof_index)
    {}

    DoFPair()
    {}
  };



  /**
   * Internal function for filling the copy indices from global to level indices
   */
  template <int dim, int spacedim>
  void fill_copy_indices(const DoFHandler<dim,spacedim> &mg_dof,
                         const MGConstrainedDoFs        *mg_constrained_dofs,
                         std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices,
                         std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_global_mine,
                         std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_level_mine)
  {
    // Now we are filling the variables copy_indices*, which are essentially
    // maps from global to mgdof for each level stored as a std::vector of
    // pairs. We need to split this map on each level depending on the ownership
    // of the global and mgdof, so that we later not access non-local elements
    // in copy_to/from_mg.
    // We keep track in the bitfield dof_touched which global dof has
    // been processed already (on the current level). This is the same as
    // the multigrid running in serial.

    // map cpu_index -> vector of data
    // that will be copied into copy_indices_level_mine
    std::vector<DoFPair> send_data_temp;

    const unsigned int n_levels = mg_dof.get_triangulation().n_global_levels();
    copy_indices.resize(n_levels);
    copy_indices_global_mine.resize(n_levels);
    copy_indices_level_mine.resize(n_levels);
    IndexSet globally_relevant;
    DoFTools::extract_locally_relevant_dofs(mg_dof, globally_relevant);

    const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> global_dof_indices (dofs_per_cell);
    std::vector<types::global_dof_index> level_dof_indices  (dofs_per_cell);

    for (unsigned int level=0; level<n_levels; ++level)
      {
        std::vector<bool> dof_touched(globally_relevant.n_elements(), false);
        copy_indices[level].clear();
        copy_indices_level_mine[level].clear();
        copy_indices_global_mine[level].clear();

        typename DoFHandler<dim,spacedim>::active_cell_iterator
        level_cell = mg_dof.begin_active(level);
        const typename DoFHandler<dim,spacedim>::active_cell_iterator
        level_end  = mg_dof.end_active(level);

        for (; level_cell!=level_end; ++level_cell)
          {
            if (mg_dof.get_triangulation().locally_owned_subdomain()!=numbers::invalid_subdomain_id
                &&  (level_cell->level_subdomain_id()==numbers::artificial_subdomain_id
                     ||  level_cell->subdomain_id()==numbers::artificial_subdomain_id)
               )
              continue;

            // get the dof numbers of this cell for the global and the level-wise
            // numbering
            level_cell->get_dof_indices (global_dof_indices);
            level_cell->get_mg_dof_indices (level_dof_indices);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                // we need to ignore if the DoF is on a refinement edge (hanging node)
                if (mg_constrained_dofs != 0
                    && mg_constrained_dofs->at_refinement_edge(level, level_dof_indices[i]))
                  continue;
                types::global_dof_index global_idx = globally_relevant.index_within_set(global_dof_indices[i]);
                //skip if we did this global dof already (on this or a coarser level)
                if (dof_touched[global_idx])
                  continue;
                bool global_mine = mg_dof.locally_owned_dofs().is_element(global_dof_indices[i]);
                bool level_mine = mg_dof.locally_owned_mg_dofs(level).is_element(level_dof_indices[i]);


                if (global_mine && level_mine)
                  {
                    copy_indices[level].push_back(
                      std::make_pair (global_dof_indices[i], level_dof_indices[i]));
                  }
                else if (global_mine)
                  {
                    copy_indices_global_mine[level].push_back(
                      std::make_pair (global_dof_indices[i], level_dof_indices[i]));

                    //send this to the owner of the level_dof:
                    send_data_temp.push_back(DoFPair(level, global_dof_indices[i], level_dof_indices[i]));
                  }
                else
                  {
                    // somebody will send those to me
                  }

                dof_touched[global_idx] = true;
              }
          }
      }

    const dealii::parallel::distributed::Triangulation<dim,spacedim> *tria =
      (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
       (&mg_dof.get_triangulation()));
    AssertThrow(send_data_temp.size()==0 || tria!=NULL, ExcMessage("parallel Multigrid only works with a distributed Triangulation!"));

#ifdef DEAL_II_WITH_MPI
    if (tria)
      {
        // TODO: Searching the owner for every single DoF becomes quite
        // inefficient. Please fix this, Timo.
        // The list of neighbors is symmetric (our neighbors have us as a neighbor),
        // so we can use it to send and to know how many messages we will get.
        std::set<unsigned int> neighbors = tria->level_ghost_owners();
        std::map<int, std::vector<DoFPair> > send_data;

        // * find owners of the level dofs and insert into send_data accordingly
        for (typename std::vector<DoFPair>::iterator dofpair=send_data_temp.begin(); dofpair != send_data_temp.end(); ++dofpair)
          {
            std::set<unsigned int>::iterator it;
            for (it = neighbors.begin(); it != neighbors.end(); ++it)
              {
                if (mg_dof.locally_owned_mg_dofs_per_processor(dofpair->level)[*it].is_element(dofpair->level_dof_index))
                  {
                    send_data[*it].push_back(*dofpair);
                    break;
                  }
              }
            // Is this level dof not owned by any of our neighbors? That
            // would certainly be a bug!
            Assert(it!=neighbors.end(), ExcMessage("could not find DoF owner."));
          }

        // * send
        std::vector<MPI_Request> requests;
        {
          for (std::set<unsigned int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
            {
              requests.push_back(MPI_Request());
              unsigned int dest = *it;
              std::vector<DoFPair> &data = send_data[dest];
              // If there is nothing to send, we still need to send a message, because
              // the receiving end will be waitng. In that case we just send
              // an empty message.
              if (data.size())
                MPI_Isend(&data[0], data.size()*sizeof(data[0]), MPI_BYTE, dest, 71, tria->get_communicator(), &*requests.rbegin());
              else
                MPI_Isend(NULL, 0, MPI_BYTE, dest, 71, tria->get_communicator(), &*requests.rbegin());
            }
        }

        // * receive
        {
          // We should get one message from each of our neighbors
          std::vector<DoFPair> receive_buffer;
          for (unsigned int counter=0; counter<neighbors.size(); ++counter)
            {
              MPI_Status status;
              int len;
              MPI_Probe(MPI_ANY_SOURCE, 71, tria->get_communicator(), &status);
              MPI_Get_count(&status, MPI_BYTE, &len);

              if (len==0)
                {
                  int err = MPI_Recv(NULL, 0, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                     tria->get_communicator(), &status);
                  AssertThrow(err==MPI_SUCCESS, ExcInternalError());
                  continue;
                }

              int count = len / sizeof(DoFPair);
              Assert(static_cast<int>(count * sizeof(DoFPair)) == len, ExcInternalError());
              receive_buffer.resize(count);

              void *ptr = &receive_buffer[0];
              int err = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                 tria->get_communicator(), &status);
              AssertThrow(err==MPI_SUCCESS, ExcInternalError());

              for (unsigned int i=0; i<receive_buffer.size(); ++i)
                {
                  copy_indices_level_mine[receive_buffer[i].level].push_back(
                    std::make_pair (receive_buffer[i].global_dof_index, receive_buffer[i].level_dof_index)
                  );
                }
            }
        }

        // * wait for all MPI_Isend to complete
        if (requests.size() > 0)
          {
            MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
            requests.clear();
          }
#ifdef DEBUG
        // Make sure in debug mode, that everybody sent/received all packages
        // on this level. If a deadlock occurs here, the list of expected
        // senders is not computed correctly.
        MPI_Barrier(tria->get_communicator());
#endif
      }
#endif

    // Sort the indices. This will produce more reliable debug output for regression tests
    // and likely won't hurt performance even in release mode.
    std::less<std::pair<types::global_dof_index, types::global_dof_index> > compare;
    for (unsigned int level=0; level<copy_indices.size(); ++level)
      std::sort(copy_indices[level].begin(), copy_indices[level].end(), compare);
    for (unsigned int level=0; level<copy_indices_level_mine.size(); ++level)
      std::sort(copy_indices_level_mine[level].begin(), copy_indices_level_mine[level].end(), compare);
    for (unsigned int level=0; level<copy_indices_global_mine.size(); ++level)
      std::sort(copy_indices_global_mine[level].begin(), copy_indices_global_mine[level].end(), compare);
  }
}



/* ------------------ MGLevelGlobalTransfer<VectorType> ----------------- */


template <typename VectorType>
template <int dim, int spacedim>
void
MGLevelGlobalTransfer<VectorType>::fill_and_communicate_copy_indices
(const DoFHandler<dim,spacedim> &mg_dof)
{
  fill_copy_indices(mg_dof, mg_constrained_dofs, copy_indices,
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
  const MPI_Comm mpi_communicator = ptria != 0 ? ptria->get_communicator() :
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
  mg_constrained_dofs = 0;
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
MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::fill_and_communicate_copy_indices
(const DoFHandler<dim,spacedim> &mg_dof)
{
  // first go to the usual routine...
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  copy_indices;
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  copy_indices_global_mine;
  std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > >
  copy_indices_level_mine;

  fill_copy_indices(mg_dof, mg_constrained_dofs, copy_indices,
                    copy_indices_global_mine, copy_indices_level_mine);

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
  const MPI_Comm mpi_communicator = ptria != 0 ? ptria->get_communicator() :
                                    MPI_COMM_SELF;

  IndexSet index_set(mg_dof.locally_owned_dofs().size());
  std::vector<types::global_dof_index> accessed_indices;
  ghosted_level_vector.resize(0, mg_dof.get_triangulation().n_global_levels()-1);
  std::vector<IndexSet> level_index_set(mg_dof.get_triangulation().n_global_levels());
  for (unsigned int l=0; l<mg_dof.get_triangulation().n_global_levels(); ++l)
    {
      for (unsigned int i=0; i<copy_indices_level_mine[l].size(); ++i)
        accessed_indices.push_back(copy_indices_level_mine[l][i].first);
      std::vector<types::global_dof_index> accessed_level_indices;
      for (unsigned int i=0; i<copy_indices_global_mine[l].size(); ++i)
        accessed_level_indices.push_back(copy_indices_global_mine[l][i].second);
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
      this->copy_indices[level].resize(copy_indices[level].size());
      for (unsigned int i=0; i<copy_indices[level].size(); ++i)
        this->copy_indices[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(copy_indices[level][i].first),
           level_partitioner.global_to_local(copy_indices[level][i].second));

      // remote-owned case: the locally owned indices for the level and the
      // ghost dofs for the global indices set the local index
      this->copy_indices_level_mine[level].
      resize(copy_indices_level_mine[level].size());
      for (unsigned int i=0; i<copy_indices_level_mine[level].size(); ++i)
        this->copy_indices_level_mine[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(copy_indices_level_mine[level][i].first),
           level_partitioner.global_to_local(copy_indices_level_mine[level][i].second));

      // owned-remote case: the locally owned indices for the global dofs
      // and the ghost dofs for the level indices set the local index
      this->copy_indices_global_mine[level].
      resize(copy_indices_global_mine[level].size());
      for (unsigned int i=0; i<copy_indices_global_mine[level].size(); ++i)
        this->copy_indices_global_mine[level][i] =
          std::pair<unsigned int,unsigned int>
          (global_partitioner.global_to_local(copy_indices_global_mine[level][i].first),
           level_partitioner.global_to_local(copy_indices_global_mine[level][i].second));
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
MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::clear()
{
  sizes.resize(0);
  copy_indices.clear();
  copy_indices_global_mine.clear();
  copy_indices_level_mine.clear();
  component_to_block_map.resize(0);
  mg_constrained_dofs = 0;
  ghosted_global_vector.reinit(0);
  ghosted_level_vector.resize(0, 0);
}



template <typename Number>
void
MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::print_indices (std::ostream &os) const
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
MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::memory_consumption () const
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
template class MGLevelGlobalTransfer<parallel::distributed::Vector<float> >;


DEAL_II_NAMESPACE_CLOSE
