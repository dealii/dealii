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
    number_cache.n_locally_owned_active_cells.resize(n_subdomains);
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
        number_cache.n_locally_owned_active_cells) +
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
    : n_global_active_cells(0)
    , n_global_levels(0)
  {}

  template <int dim, int spacedim>
  unsigned int
  TriangulationBase<dim, spacedim>::n_locally_owned_active_cells() const
  {
    return number_cache.n_locally_owned_active_cells[my_subdomain];
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
  const std::vector<unsigned int> &
  TriangulationBase<dim, spacedim>::n_locally_owned_active_cells_per_processor()
    const
  {
    return number_cache.n_locally_owned_active_cells;
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
    Assert(number_cache.n_locally_owned_active_cells.size() ==
             Utilities::MPI::n_mpi_processes(this->mpi_communicator),
           ExcInternalError());

    std::fill(number_cache.n_locally_owned_active_cells.begin(),
              number_cache.n_locally_owned_active_cells.end(),
              0);

    number_cache.ghost_owners.clear();
    number_cache.level_ghost_owners.clear();

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
          ++number_cache.n_locally_owned_active_cells[my_subdomain];

    unsigned int send_value =
      number_cache.n_locally_owned_active_cells[my_subdomain];
    const int ierr =
      MPI_Allgather(&send_value,
                    1,
                    MPI_UNSIGNED,
                    number_cache.n_locally_owned_active_cells.data(),
                    1,
                    MPI_UNSIGNED,
                    this->mpi_communicator);
    AssertThrowMPI(ierr);

    number_cache.n_global_active_cells =
      std::accumulate(number_cache.n_locally_owned_active_cells.begin(),
                      number_cache.n_locally_owned_active_cells.end(),
                      /* ensure sum is computed with correct data type:*/
                      static_cast<types::global_dof_index>(0));
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
    std::map<unsigned int, std::set<unsigned int>>
      vertex_to_coinciding_vertices;
    {
      // 1a) collect nodes coinciding due to periodicity
      std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex;
      for (auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned() || cell->is_ghost())
          for (auto i = 0u; i < GeometryInfo<dim>::faces_per_cell; ++i)
            if (cell->has_periodic_neighbor(i) &&
                cell->periodic_neighbor(i)->active())
              {
                auto face_t = cell->face(i);
                auto face_n = cell->periodic_neighbor(i)->face(
                  cell->periodic_neighbor_face_no(i));
                for (auto j = 0u; j < GeometryInfo<dim>::vertices_per_face; ++j)
                  {
                    auto         v_t  = face_t->vertex_index(j);
                    auto         v_n  = face_n->vertex_index(j);
                    unsigned int temp = std::min(v_t, v_n);
                    {
                      auto it = vertex_to_coinciding_vertex.find(v_t);
                      if (it != vertex_to_coinciding_vertex.end())
                        temp = std::min(temp, it->second);
                    }
                    {
                      auto it = vertex_to_coinciding_vertex.find(v_n);
                      if (it != vertex_to_coinciding_vertex.end())
                        temp = std::min(temp, it->second);
                    }
                    vertex_to_coinciding_vertex[v_t] = temp;
                    vertex_to_coinciding_vertex[v_n] = temp;
                  }
              }

      // 1b) compress map: let vertices point to the coinciding vertex with
      //     the smallest index
      for (auto &p : vertex_to_coinciding_vertex)
        {
          if (p.first == p.second)
            continue;
          unsigned int temp = p.second;
          while (temp != vertex_to_coinciding_vertex[temp])
            temp = vertex_to_coinciding_vertex[temp];
          p.second = temp;
        }

#ifdef DEBUG
      // check if map is actually compressed
      for (auto p : vertex_to_coinciding_vertex)
        {
          if (p.first == p.second)
            continue;
          auto pp = vertex_to_coinciding_vertex.find(p.second);
          if (pp->first == pp->second)
            continue;
          AssertThrow(false, ExcMessage("Map has to be compressed!"));
        }
#endif

      // 1c) create a map: smallest index of coinciding index -> all
      // coinciding indices
      std::map<unsigned int, std::set<unsigned int>>
        smallest_coinciding_vertex_to_coinciding_vertices;
      for (auto p : vertex_to_coinciding_vertex)
        smallest_coinciding_vertex_to_coinciding_vertices[p.second] =
          std::set<unsigned int>();

      for (auto p : vertex_to_coinciding_vertex)
        smallest_coinciding_vertex_to_coinciding_vertices[p.second].insert(
          p.first);

      // 1d) create a map: vertex -> all coinciding indices
      for (auto &s : smallest_coinciding_vertex_to_coinciding_vertices)
        for (auto &ss : s.second)
          vertex_to_coinciding_vertices[ss] = s.second;
    }

    // 2) collect vertices belonging to local cells
    std::vector<bool> vertex_of_own_cell(this->n_vertices(), false);
    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_locally_owned())
        for (auto v = 0u; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          vertex_of_own_cell[cell->vertex_index(v)] = true;

    // 3) for for each vertex belonging to a locally owned cell all ghost
    //    neighbors (including the periodic own)
    std::map<unsigned int, std::set<dealii::types::subdomain_id>> result;

    // loop over all active ghost cells
    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_ghost())
        {
          const auto owner = cell->subdomain_id();

          // loop over all its vertices
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              // set owner if vertex belongs to a local cell
              if (vertex_of_own_cell[cell->vertex_index(v)])
                result[cell->vertex_index(v)].insert(owner);

              // mark also nodes coinciding due to periodicity
              auto coinciding_vertices =
                vertex_to_coinciding_vertices.find(cell->vertex_index(v));
              if (coinciding_vertices != vertex_to_coinciding_vertices.end())
                for (auto coinciding_vertex : coinciding_vertices->second)
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
