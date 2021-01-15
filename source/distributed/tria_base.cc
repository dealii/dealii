// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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
#include <deal.II/base/mpi.templates.h>
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
    const MPI_Comm &mpi_communicator,
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
    Assert(false, ExcNeedsMPI());
#endif
  }



  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::copy_triangulation(
    const dealii::Triangulation<dim, spacedim> &other_tria)
  {
#ifndef DEAL_II_WITH_MPI
    (void)other_tria;
    Assert(false, ExcNeedsMPI());
#else
    dealii::Triangulation<dim, spacedim>::copy_triangulation(other_tria);

    if (const dealii::parallel::TriangulationBase<dim, spacedim> *other_tria_x =
          dynamic_cast<const dealii::parallel::TriangulationBase<dim, spacedim>
                         *>(&other_tria))
      {
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
  types::global_cell_index
  TriangulationBase<dim, spacedim>::n_global_active_cells() const
  {
    return number_cache.n_global_active_cells;
  }

  template <int dim, int spacedim>
  const MPI_Comm &
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
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_ghost())
          number_cache.ghost_owners.insert(cell->subdomain_id());

      Assert(number_cache.ghost_owners.size() <
               Utilities::MPI::n_mpi_processes(mpi_communicator),
             ExcInternalError());
    }

    if (this->n_levels() > 0)
      for (const auto &cell : this->active_cell_iterators())
        if (cell->subdomain_id() == my_subdomain)
          ++number_cache.n_locally_owned_active_cells;

    // Potentially cast to a 64 bit type before accumulating to avoid overflow:
    number_cache.n_global_active_cells =
      Utilities::MPI::sum(static_cast<types::global_cell_index>(
                            number_cache.n_locally_owned_active_cells),
                          this->mpi_communicator);

    number_cache.n_global_levels =
      Utilities::MPI::max(this->n_levels(), this->mpi_communicator);

    // Store MPI ranks of level ghost owners of this processor on all levels.
    if (this->is_multilevel_hierarchy_constructed() == true)
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

          const int mpi_tag = Utilities::MPI::internal::Tags::
            triangulation_base_fill_level_ghost_owners;

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
                               mpi_tag,
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
                              mpi_tag,
                              this->mpi_communicator,
                              MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr);
            }

          if (requests.size() > 0)
            {
              ierr = MPI_Waitall(requests.size(),
                                 requests.data(),
                                 MPI_STATUSES_IGNORE);
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

    // reset global cell ids
    this->reset_global_cell_indices();
  }

#else

  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::update_number_cache()
  {
    Assert(false, ExcNotImplemented());
  }

#endif

  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::update_reference_cell_types()
  {
    // run algorithm for locally-owned cells
    dealii::Triangulation<dim, spacedim>::update_reference_cell_types();

    // translate ReferenceCell::Type to unsigned int (needed by
    // Utilities::MPI::compute_set_union)
    std::vector<unsigned int> reference_cell_types_ui;

    for (const auto &i : this->reference_cell_types)
      reference_cell_types_ui.push_back(static_cast<unsigned int>(i));

    // create union
    reference_cell_types_ui =
      Utilities::MPI::compute_set_union(reference_cell_types_ui,
                                        this->mpi_communicator);

    // transform back and store result
    this->reference_cell_types.clear();
    for (const auto &i : reference_cell_types_ui)
      this->reference_cell_types.emplace_back(
        static_cast<ReferenceCell::Type::CellKinds>(i));
  }


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
    return GridTools::compute_vertices_with_ghost_neighbors(*this);
  }



  template <int dim, int spacedim>
  std::vector<types::boundary_id>
  TriangulationBase<dim, spacedim>::get_boundary_ids() const
  {
    return Utilities::MPI::compute_set_union(
      dealii::Triangulation<dim, spacedim>::get_boundary_ids(),
      this->mpi_communicator);
  }



  template <int dim, int spacedim>
  std::vector<types::manifold_id>
  TriangulationBase<dim, spacedim>::get_manifold_ids() const
  {
    return Utilities::MPI::compute_set_union(
      dealii::Triangulation<dim, spacedim>::get_manifold_ids(),
      this->mpi_communicator);
  }



  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::reset_global_cell_indices()
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false, ExcNeedsMPI());
#else

    // currently only implemented for distributed triangulations
    if (dynamic_cast<const parallel::DistributedTriangulationBase<dim, spacedim>
                       *>(this) == nullptr)
      return;

    // 1) determine number of active locally-owned cells
    const types::global_cell_index n_locally_owned_cells =
      this->n_locally_owned_active_cells();

    // 2) determine the offset of each process
    types::global_cell_index cell_index = 0;

    MPI_Exscan(&n_locally_owned_cells,
               &cell_index,
               1,
               Utilities::MPI::internal::mpi_type_id(&n_locally_owned_cells),
               MPI_SUM,
               this->mpi_communicator);

    // 3) give global indices to locally-owned cells and mark all other cells as
    //    invalid
    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_global_active_cell_index(cell_index++);
      else
        cell->set_global_active_cell_index(numbers::invalid_dof_index);

    // 4) determine the global indices of ghost cells
    GridTools::exchange_cell_data_to_ghosts<types::global_cell_index>(
      *this,
      [](const auto &cell) { return cell->global_active_cell_index(); },
      [](const auto &cell, const auto &id) {
        cell->set_global_active_cell_index(id);
      });

    // 5) set up new partitioner
    IndexSet is_local(this->n_global_active_cells());
    IndexSet is_ghost(this->n_global_active_cells());

    for (const auto &cell : this->active_cell_iterators())
      if (!cell->is_artificial())
        {
          const auto index = cell->global_active_cell_index();

          if (index == numbers::invalid_dof_index)
            continue;

          if (cell->is_locally_owned())
            is_local.add_index(index);
          else
            is_ghost.add_index(index);
        }

    number_cache.active_cell_index_partitioner =
      Utilities::MPI::Partitioner(is_local, is_ghost, this->mpi_communicator);

    // 6) proceed with multigrid levels if requested
    if (this->is_multilevel_hierarchy_constructed() == true)
      {
        // 1) determine number of locally-owned cells on levels
        std::vector<types::global_cell_index> n_locally_owned_cells(
          this->n_global_levels(), 0);

        for (auto cell : this->cell_iterators())
          if (cell->level_subdomain_id() == this->locally_owned_subdomain())
            n_locally_owned_cells[cell->level()]++;

        // 2) determine the offset of each process
        std::vector<types::global_cell_index> cell_index(
          this->n_global_levels(), 0);

        MPI_Exscan(n_locally_owned_cells.data(),
                   cell_index.data(),
                   this->n_global_levels(),
                   Utilities::MPI::internal::mpi_type_id(
                     n_locally_owned_cells.data()),
                   MPI_SUM,
                   this->mpi_communicator);

        // 3) determine global number of "active" cells on each level
        std::vector<types::global_cell_index> n_cells_level(
          this->n_global_levels(), 0);

        for (unsigned int l = 0; l < this->n_global_levels(); ++l)
          n_cells_level[l] = n_locally_owned_cells[l] + cell_index[l];

        MPI_Bcast(n_cells_level.data(),
                  this->n_global_levels(),
                  Utilities::MPI::internal::mpi_type_id(n_cells_level.data()),
                  this->n_subdomains - 1,
                  this->mpi_communicator);

        // 4) give global indices to locally-owned cells on level and mark
        //    all other cells as invalid
        for (auto cell : this->cell_iterators())
          if (cell->level_subdomain_id() == this->locally_owned_subdomain())
            cell->set_global_level_cell_index(cell_index[cell->level()]++);
          else
            cell->set_global_level_cell_index(numbers::invalid_dof_index);

        // 5) update the numbers of ghost level cells
        GridTools::exchange_cell_data_to_level_ghosts<
          types::global_cell_index,
          dealii::Triangulation<dim, spacedim>>(
          *this,
          [](const auto &cell) { return cell->global_level_cell_index(); },
          [](const auto &cell, const auto &id) {
            return cell->set_global_level_cell_index(id);
          });

        number_cache.level_cell_index_partitioners.resize(
          this->n_global_levels());

        // 6) set up cell partitioners for each level
        for (unsigned int l = 0; l < this->n_global_levels(); ++l)
          {
            IndexSet is_local(n_cells_level[l]);
            IndexSet is_ghost(n_cells_level[l]);

            for (const auto &cell : this->cell_iterators_on_level(l))
              if (cell->level_subdomain_id() !=
                  dealii::numbers::artificial_subdomain_id)
                {
                  const auto index = cell->global_level_cell_index();

                  if (index == numbers::invalid_dof_index)
                    continue;

                  if (cell->level_subdomain_id() ==
                      this->locally_owned_subdomain())
                    is_local.add_index(index);
                  else
                    is_ghost.add_index(index);
                }

            number_cache.level_cell_index_partitioners[l] =
              Utilities::MPI::Partitioner(is_local,
                                          is_ghost,
                                          this->mpi_communicator);
          }
      }

#endif
  }



  template <int dim, int spacedim>
  const Utilities::MPI::Partitioner &
  TriangulationBase<dim, spacedim>::global_active_cell_index_partitioner() const
  {
    return number_cache.active_cell_index_partitioner;
  }

  template <int dim, int spacedim>
  const Utilities::MPI::Partitioner &
  TriangulationBase<dim, spacedim>::global_level_cell_index_partitioner(
    const unsigned int level) const
  {
    Assert(this->is_multilevel_hierarchy_constructed(), ExcNotImplemented());
    AssertIndexRange(level, this->n_global_levels());

    return number_cache.level_cell_index_partitioners[level];
  }

  template <int dim, int spacedim>
  DistributedTriangulationBase<dim, spacedim>::DistributedTriangulationBase(
    const MPI_Comm &mpi_communicator,
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
