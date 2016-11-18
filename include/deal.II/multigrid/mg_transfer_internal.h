// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


#ifndef dealii__mg_transfer_internal_h
#define dealii__mg_transfer_internal_h

#include <deal.II/dofs/dof_tools.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MGTransfer
  {
    /**
     * Internal data structure that is used in the MPI communication in
     * fill_and_communicate_copy_indices().  It represents an entry in the
     * copy_indices* map, that associates a level dof index with a global dof
     * index.
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
        level(level),
        global_dof_index(global_dof_index),
        level_dof_index(level_dof_index)
      {}

      DoFPair()
        :
        level (numbers::invalid_unsigned_int),
        global_dof_index (numbers::invalid_dof_index),
        level_dof_index (numbers::invalid_dof_index)
      {}
    };


    /**
     * Internal function for filling the copy indices from global to level
     * indices
     */
    template <int dim, int spacedim>
    void fill_copy_indices(const dealii::DoFHandler<dim,spacedim> &mg_dof,
                           const MGConstrainedDoFs        *mg_constrained_dofs,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_global_mine,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_level_mine)
    {
      // Now we are filling the variables copy_indices*, which are essentially
      // maps from global to mgdof for each level stored as a std::vector of
      // pairs. We need to split this map on each level depending on the
      // ownership of the global and mgdof, so that we later not access
      // non-local elements in copy_to/from_mg.
      // We keep track in the bitfield dof_touched which global dof has been
      // processed already (on the current level). This is the same as the
      // multigrid running in serial.

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

          typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
          level_cell = mg_dof.begin_active(level);
          const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
          level_end  = mg_dof.end_active(level);

          for (; level_cell!=level_end; ++level_cell)
            {
              if (mg_dof.get_triangulation().locally_owned_subdomain()!=numbers::invalid_subdomain_id
                  &&  (level_cell->level_subdomain_id()==numbers::artificial_subdomain_id
                       ||  level_cell->subdomain_id()==numbers::artificial_subdomain_id)
                 )
                continue;

              // get the dof numbers of this cell for the global and the
              // level-wise numbering
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
          // The list of neighbors is symmetric (our neighbors have us as a
          // neighbor), so we can use it to send and to know how many messages
          // we will get.
          std::set<types::subdomain_id> neighbors = tria->level_ghost_owners();
          std::map<int, std::vector<DoFPair> > send_data;

          // * find owners of the level dofs and insert into send_data accordingly
          for (typename std::vector<DoFPair>::iterator dofpair=send_data_temp.begin(); dofpair != send_data_temp.end(); ++dofpair)
            {
              std::set<types::subdomain_id>::iterator it;
              for (it = neighbors.begin(); it != neighbors.end(); ++it)
                {
                  if (mg_dof.locally_owned_mg_dofs_per_processor(dofpair->level)[*it].is_element(dofpair->level_dof_index))
                    {
                      send_data[*it].push_back(*dofpair);
                      break;
                    }
                }
              // Is this level dof not owned by any of our neighbors? That would
              // certainly be a bug!
              Assert(it!=neighbors.end(), ExcMessage("could not find DoF owner."));
            }

          // * send
          std::vector<MPI_Request> requests;
          {
            for (std::set<types::subdomain_id>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
              {
                requests.push_back(MPI_Request());
                unsigned int dest = *it;
                std::vector<DoFPair> &data = send_data[dest];
                // If there is nothing to send, we still need to send a message,
                // because the receiving end will be waitng. In that case we
                // just send an empty message.
                if (data.size())
                  {
                    const int ierr = MPI_Isend(&data[0], data.size()*sizeof(data[0]),
                                               MPI_BYTE, dest, 71, tria->get_communicator(),
                                               &*requests.rbegin());
                    AssertThrowMPI(ierr);
                  }
                else
                  {
                    const int ierr = MPI_Isend(NULL, 0, MPI_BYTE, dest, 71,
                                               tria->get_communicator(), &*requests.rbegin());
                    AssertThrowMPI(ierr);
                  }
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
                int ierr = MPI_Probe(MPI_ANY_SOURCE, 71, tria->get_communicator(), &status);
                AssertThrowMPI(ierr);
                ierr = MPI_Get_count(&status, MPI_BYTE, &len);
                AssertThrowMPI(ierr);

                if (len==0)
                  {
                    ierr = MPI_Recv(NULL, 0, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                    tria->get_communicator(), &status);
                    AssertThrowMPI(ierr);
                    continue;
                  }

                int count = len / sizeof(DoFPair);
                Assert(static_cast<int>(count * sizeof(DoFPair)) == len, ExcInternalError());
                receive_buffer.resize(count);

                void *ptr = &receive_buffer[0];
                ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                tria->get_communicator(), &status);
                AssertThrowMPI(ierr);

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
              const int ierr = MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
              requests.clear();
            }
#ifdef DEBUG
          // Make sure in debug mode, that everybody sent/received all packages
          // on this level. If a deadlock occurs here, the list of expected
          // senders is not computed correctly.
          const int ierr = MPI_Barrier(tria->get_communicator());
          AssertThrowMPI(ierr);
#endif
        }
#endif

      // Sort the indices. This will produce more reliable debug output for
      // regression tests and likely won't hurt performance even in release
      // mode.
      std::less<std::pair<types::global_dof_index, types::global_dof_index> > compare;
      for (unsigned int level=0; level<copy_indices.size(); ++level)
        std::sort(copy_indices[level].begin(), copy_indices[level].end(), compare);
      for (unsigned int level=0; level<copy_indices_level_mine.size(); ++level)
        std::sort(copy_indices_level_mine[level].begin(), copy_indices_level_mine[level].end(), compare);
      for (unsigned int level=0; level<copy_indices_global_mine.size(); ++level)
        std::sort(copy_indices_global_mine[level].begin(), copy_indices_global_mine[level].end(), compare);
    }

  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
