// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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


#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector_memory.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  template <int dim, int spacedim>
  TriangulationBase<dim, spacedim>::TriangulationBase(
    MPI_Comm mpi_communicator,
    const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
               smooth_grid,
    const bool check_for_distorted_cells)
    : dealii::Triangulation<dim, spacedim>(smooth_grid,
                                           check_for_distorted_cells)
    , mpi_communicator(mpi_communicator)
    , my_subdomain(Utilities::MPI::this_mpi_process(this->mpi_communicator))
    , n_subdomains(Utilities::MPI::n_mpi_processes(this->mpi_communicator))
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false,
           ExcMessage("You compiled deal.II without MPI support, for "
                      "which parallel::TriangulationBase is not available."));
#endif
  }



  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::copy_triangulation(
    const dealii::Triangulation<dim, spacedim> &other_tria)
  {
#ifndef DEAL_II_WITH_MPI
    (void)other_tria;
    Assert(false,
           ExcMessage("You compiled deal.II without MPI support, for "
                      "which parallel::TriangulationBase is not available."));
#else
    dealii::Triangulation<dim, spacedim>::copy_triangulation(other_tria);

    if (const dealii::parallel::TriangulationBase<dim, spacedim> *other_tria_x =
          dynamic_cast<const dealii::parallel::TriangulationBase<dim, spacedim>
                         *>(&other_tria))
      {
        mpi_communicator = other_tria_x->get_communicator();

        // release unused vector memory because we will have very different
        // vectors now
        ::dealii::internal::GrowingVectorMemoryImplementation::
          release_all_unused_memory();
      }
#endif
  }



  template <int dim, int spacedim>
  std::size_t
  TriangulationBase<dim, spacedim>::memory_consumption() const
  {
    std::size_t mem =
      this->dealii::Triangulation<dim, spacedim>::memory_consumption() +
      MemoryConsumption::memory_consumption(mpi_communicator) +
      MemoryConsumption::memory_consumption(my_subdomain) +
      MemoryConsumption::memory_consumption(
        number_cache.n_global_active_cells) +
      MemoryConsumption::memory_consumption(number_cache.n_global_levels);
    return mem;
  }

  template <int dim, int spacedim>
  TriangulationBase<dim, spacedim>::~TriangulationBase()
  {
    // release unused vector memory because the vector layout is going to look
    // very different now
    ::dealii::internal::GrowingVectorMemoryImplementation::
      release_all_unused_memory();
  }

  template <int dim, int spacedim>
  TriangulationBase<dim, spacedim>::NumberCache::NumberCache()
    : n_locally_owned_active_cells(0)
    , n_global_levels(0)
  {}

  template <int dim, int spacedim>
  unsigned int
  TriangulationBase<dim, spacedim>::n_locally_owned_active_cells() const
  {
    return number_cache.n_locally_owned_active_cells;
  }

  template <int dim, int spacedim>
  unsigned int
  TriangulationBase<dim, spacedim>::n_global_levels() const
  {
    return number_cache.n_global_levels;
  }

  template <int dim, int spacedim>
  types::global_dof_index
  TriangulationBase<dim, spacedim>::n_global_active_cells() const
  {
    return number_cache.n_global_active_cells;
  }

  template <int dim, int spacedim>
  std::vector<unsigned int>
  TriangulationBase<dim, spacedim>::
    compute_n_locally_owned_active_cells_per_processor() const
  {
    ;
#ifdef DEAL_II_WITH_MPI
    std::vector<unsigned int> n_locally_owned_active_cells_per_processor(
      Utilities::MPI::n_mpi_processes(this->mpi_communicator), 0);

    if (this->n_levels() > 0)
      {
        const int ierr =
          MPI_Allgather(&number_cache.n_locally_owned_active_cells,
                        1,
                        MPI_UNSIGNED,
                        n_locally_owned_active_cells_per_processor.data(),
                        1,
                        MPI_UNSIGNED,
                        this->mpi_communicator);
        AssertThrowMPI(ierr);
      }

    return n_locally_owned_active_cells_per_processor;
#else
    return {number_cache.n_locally_owned_active_cells};
#endif
  }

  template <int dim, int spacedim>
  MPI_Comm
  TriangulationBase<dim, spacedim>::get_communicator() const
  {
    return mpi_communicator;
  }

#ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::update_number_cache()
  {
    number_cache.ghost_owners.clear();
    number_cache.level_ghost_owners.clear();
    number_cache.n_locally_owned_active_cells = 0;

    if (this->n_levels() == 0)
      {
        // Skip communication done below if we do not have any cells
        // (meaning the Triangulation is empty on all processors). This will
        // happen when called from the destructor of Triangulation, which
        // can get called during exception handling causing a hang in this
        // function.
        number_cache.n_global_active_cells = 0;
        number_cache.n_global_levels       = 0;
        return;
      }


    {
      // find ghost owners
      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
             this->begin_active();
           cell != this->end();
           ++cell)
        if (cell->is_ghost())
          number_cache.ghost_owners.insert(cell->subdomain_id());

      Assert(number_cache.ghost_owners.size() <
               Utilities::MPI::n_mpi_processes(mpi_communicator),
             ExcInternalError());
    }

    if (this->n_levels() > 0)
      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
             this->begin_active();
           cell != this->end();
           ++cell)
        if (cell->subdomain_id() == my_subdomain)
          ++number_cache.n_locally_owned_active_cells;

    number_cache.n_global_active_cells =
      Utilities::MPI::sum(number_cache.n_locally_owned_active_cells,
                          this->mpi_communicator);
    number_cache.n_global_levels =
      Utilities::MPI::max(this->n_levels(), this->mpi_communicator);
  }



  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::fill_level_ghost_owners()
  {
    number_cache.level_ghost_owners.clear();

    // if there is nothing to do, then do nothing
    if (this->n_levels() == 0)
      return;

    // find level ghost owners
    for (typename Triangulation<dim, spacedim>::cell_iterator cell =
           this->begin();
         cell != this->end();
         ++cell)
      if (cell->level_subdomain_id() != numbers::artificial_subdomain_id &&
          cell->level_subdomain_id() != this->locally_owned_subdomain())
        this->number_cache.level_ghost_owners.insert(
          cell->level_subdomain_id());

#  ifdef DEBUG
    // Check that level_ghost_owners is symmetric by sending a message to
    // everyone
    {
      int ierr = MPI_Barrier(this->mpi_communicator);
      AssertThrowMPI(ierr);

      // important: preallocate to avoid (re)allocation:
      std::vector<MPI_Request> requests(
        this->number_cache.level_ghost_owners.size());
      unsigned int dummy       = 0;
      unsigned int req_counter = 0;

      for (std::set<types::subdomain_id>::iterator it =
             this->number_cache.level_ghost_owners.begin();
           it != this->number_cache.level_ghost_owners.end();
           ++it, ++req_counter)
        {
          ierr = MPI_Isend(&dummy,
                           1,
                           MPI_UNSIGNED,
                           *it,
                           9001,
                           this->mpi_communicator,
                           &requests[req_counter]);
          AssertThrowMPI(ierr);
        }

      for (std::set<types::subdomain_id>::iterator it =
             this->number_cache.level_ghost_owners.begin();
           it != this->number_cache.level_ghost_owners.end();
           ++it)
        {
          unsigned int dummy;
          ierr = MPI_Recv(&dummy,
                          1,
                          MPI_UNSIGNED,
                          *it,
                          9001,
                          this->mpi_communicator,
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
        }

      if (requests.size() > 0)
        {
          ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }

      ierr = MPI_Barrier(this->mpi_communicator);
      AssertThrowMPI(ierr);
    }
#  endif

    Assert(this->number_cache.level_ghost_owners.size() <
             Utilities::MPI::n_mpi_processes(this->mpi_communicator),
           ExcInternalError());
  }

#else

  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::update_number_cache()
  {
    Assert(false, ExcNotImplemented());
  }

  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::fill_level_ghost_owners()
  {
    Assert(false, ExcNotImplemented());
  }

#endif

  template <int dim, int spacedim>
  types::subdomain_id
  TriangulationBase<dim, spacedim>::locally_owned_subdomain() const
  {
    return my_subdomain;
  }



  template <int dim, int spacedim>
  const std::set<types::subdomain_id> &
  TriangulationBase<dim, spacedim>::ghost_owners() const
  {
    return number_cache.ghost_owners;
  }



  template <int dim, int spacedim>
  const std::set<types::subdomain_id> &
  TriangulationBase<dim, spacedim>::level_ghost_owners() const
  {
    return number_cache.level_ghost_owners;
  }



  template <int dim, int spacedim>
  std::map<unsigned int, std::set<dealii::types::subdomain_id>>
  TriangulationBase<dim, spacedim>::compute_vertices_with_ghost_neighbors()
    const
  {
    // 1) collect for each vertex on periodic faces all vertices it coincides
    //    with
    std::map<unsigned int, std::vector<unsigned int>> coinciding_vertex_groups;
    std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;

    {
      static const int lookup_table_2d[2][2] =
        //           flip:
        {
          {0, 1}, // false
          {1, 0}  // true
        };

      static const int lookup_table_3d[2][2][2][4] =
        //                   orientation flip  rotation
        {{{
            {0, 2, 1, 3}, // false       false false
            {2, 3, 0, 1}  // false       false true
          },
          {
            {3, 1, 2, 0}, // false       true  false
            {1, 0, 3, 2}  // false       true  true
          }},
         {{
            {0, 1, 2, 3}, // true        false false
            {1, 3, 0, 2}  // true        false true
          },
          {
            {3, 2, 1, 0}, // true        true  false
            {2, 0, 3, 1}  // true        true  true
          }}};

      // loop over all periodic face pairs
      for (const auto &pair : this->get_periodic_face_map())
        {
          const auto face_a = pair.first.first->face(pair.first.second);
          const auto face_b =
            pair.second.first.first->face(pair.second.first.second);
          const auto mask = pair.second.second;

          // loop over all vertices on face
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face;
               ++i)
            {
              const bool face_orientation = mask[0];
              const bool face_flip        = mask[1];
              const bool face_rotation    = mask[2];

              // find the right local vertex index for the second face
              unsigned int j = 0;
              switch (dim)
                {
                  case 1:
                    j = i;
                    break;
                  case 2:
                    j = lookup_table_2d[face_flip][i];
                    break;
                  case 3:
                    j = lookup_table_3d[face_orientation][face_flip]
                                       [face_rotation][i];
                    break;
                  default:
                    AssertThrow(false, ExcNotImplemented());
                }

              // get vertex indices and store in map
              const auto   vertex_a = face_a->vertex_index(i);
              const auto   vertex_b = face_b->vertex_index(j);
              unsigned int temp     = std::min(vertex_a, vertex_b);

              auto it_a = vertex_to_coinciding_vertex_group.find(vertex_a);
              if (it_a != vertex_to_coinciding_vertex_group.end())
                temp = std::min(temp, it_a->second);

              auto it_b = vertex_to_coinciding_vertex_group.find(vertex_b);
              if (it_b != vertex_to_coinciding_vertex_group.end())
                temp = std::min(temp, it_b->second);

              if (it_a != vertex_to_coinciding_vertex_group.end())
                it_a->second = temp;
              else
                vertex_to_coinciding_vertex_group[vertex_a] = temp;

              if (it_b != vertex_to_coinciding_vertex_group.end())
                it_b->second = temp;
              else
                vertex_to_coinciding_vertex_group[vertex_b] = temp;
            }
        }

      // 1b) compress map: let vertices point to the coinciding vertex with
      //     the smallest index
      for (auto &p : vertex_to_coinciding_vertex_group)
        {
          if (p.first == p.second)
            continue;
          unsigned int temp = p.second;
          while (temp != vertex_to_coinciding_vertex_group[temp])
            temp = vertex_to_coinciding_vertex_group[temp];
          p.second = temp;
        }

      // 1c) create a map: smallest index of coinciding index -> all
      // coinciding indices
      for (auto p : vertex_to_coinciding_vertex_group)
        coinciding_vertex_groups[p.second] = {};

      for (auto p : vertex_to_coinciding_vertex_group)
        coinciding_vertex_groups[p.second].push_back(p.first);
    }

    // 2) collect vertices belonging to local cells
    std::vector<bool> vertex_of_own_cell(this->n_vertices(), false);
    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_locally_owned())
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          vertex_of_own_cell[cell->vertex_index(v)] = true;

    // 3) for each vertex belonging to a locally owned cell all ghost
    //    neighbors (including the periodic own)
    std::map<unsigned int, std::set<types::subdomain_id>> result;

    // loop over all active ghost cells
    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_ghost())
        {
          const types::subdomain_id owner = cell->subdomain_id();

          // loop over all its vertices
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              // set owner if vertex belongs to a local cell
              if (vertex_of_own_cell[cell->vertex_index(v)])
                result[cell->vertex_index(v)].insert(owner);

              // mark also nodes coinciding due to periodicity
              auto coinciding_vertex_group =
                vertex_to_coinciding_vertex_group.find(cell->vertex_index(v));
              if (coinciding_vertex_group !=
                  vertex_to_coinciding_vertex_group.end())
                for (auto coinciding_vertex :
                     coinciding_vertex_groups[coinciding_vertex_group->second])
                  if (vertex_of_own_cell[coinciding_vertex])
                    result[coinciding_vertex].insert(owner);
            }
        }

    return result;
  }



  template <int dim, int spacedim>
  DistributedTriangulationBase<dim, spacedim>::DistributedTriangulationBase(
    MPI_Comm mpi_communicator,
    const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
               smooth_grid,
    const bool check_for_distorted_cells)
    : dealii::parallel::TriangulationBase<dim, spacedim>(
        mpi_communicator,
        smooth_grid,
        check_for_distorted_cells)
  {}

} // end namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "tria_base.inst"

DEAL_II_NAMESPACE_CLOSE
