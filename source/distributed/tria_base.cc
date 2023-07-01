// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II authors
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
#include <deal.II/base/mpi_large_count.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector_memory.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  TriangulationBase<dim, spacedim>::TriangulationBase(
    const MPI_Comm mpi_communicator,
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
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::copy_triangulation(
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
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  std::size_t TriangulationBase<dim, spacedim>::memory_consumption() const
  {
    std::size_t mem =
      this->dealii::Triangulation<dim, spacedim>::memory_consumption() +
      MemoryConsumption::memory_consumption(this->mpi_communicator) +
      MemoryConsumption::memory_consumption(my_subdomain) +
      MemoryConsumption::memory_consumption(
        number_cache.n_global_active_cells) +
      MemoryConsumption::memory_consumption(number_cache.n_global_levels);
    return mem;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  TriangulationBase<dim, spacedim>::~TriangulationBase()
  {
    // release unused vector memory because the vector layout is going to look
    // very different now
    ::dealii::internal::GrowingVectorMemoryImplementation::
      release_all_unused_memory();
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  TriangulationBase<dim, spacedim>::NumberCache::NumberCache()
    : n_locally_owned_active_cells(0)
    , number_of_global_coarse_cells(0)
    , n_global_levels(0)
  {}



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  unsigned int TriangulationBase<dim, spacedim>::n_locally_owned_active_cells()
    const
  {
    return number_cache.n_locally_owned_active_cells;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  unsigned int TriangulationBase<dim, spacedim>::n_global_levels() const
  {
    return number_cache.n_global_levels;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  types::global_cell_index
    TriangulationBase<dim, spacedim>::n_global_active_cells() const
  {
    return number_cache.n_global_active_cells;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  MPI_Comm TriangulationBase<dim, spacedim>::get_communicator() const
  {
    return mpi_communicator;
  }



#ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::update_number_cache()
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
               Utilities::MPI::n_mpi_processes(this->mpi_communicator),
             ExcInternalError());
    }

    if (this->n_levels() > 0)
      number_cache.n_locally_owned_active_cells = std::count_if(
        this->begin_active(),
        typename Triangulation<dim, spacedim>::active_cell_iterator(
          this->end()),
        [](const auto &i) { return i.is_locally_owned(); });
    else
      number_cache.n_locally_owned_active_cells = 0;

    // Potentially cast to a 64 bit type before accumulating to avoid
    // overflow:
    number_cache.n_global_active_cells =
      Utilities::MPI::sum(static_cast<types::global_cell_index>(
                            number_cache.n_locally_owned_active_cells),
                          this->mpi_communicator);

    number_cache.n_global_levels =
      Utilities::MPI::max(this->n_levels(), this->mpi_communicator);

    // Store MPI ranks of level ghost owners of this processor on all
    // levels.
    if (this->is_multilevel_hierarchy_constructed() == true)
      {
        number_cache.level_ghost_owners.clear();

        // if there is nothing to do, then do nothing
        if (this->n_levels() == 0)
          return;

        // find level ghost owners
        for (const auto &cell : this->cell_iterators())
          if (cell->level_subdomain_id() != numbers::artificial_subdomain_id &&
              cell->level_subdomain_id() != this->locally_owned_subdomain())
            this->number_cache.level_ghost_owners.insert(
              cell->level_subdomain_id());

#  ifdef DEBUG
        // Check that level_ghost_owners is symmetric by sending a message
        // to everyone
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

          for (const auto &it : this->number_cache.level_ghost_owners)
            {
              ierr = MPI_Isend(&dummy,
                               1,
                               MPI_UNSIGNED,
                               it,
                               mpi_tag,
                               this->mpi_communicator,
                               &requests[req_counter]);
              AssertThrowMPI(ierr);
              ++req_counter;
            }

          for (const auto &it : this->number_cache.level_ghost_owners)
            {
              unsigned int dummy;
              ierr = MPI_Recv(&dummy,
                              1,
                              MPI_UNSIGNED,
                              it,
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

    this->number_cache.number_of_global_coarse_cells = this->n_cells(0);

    // reset global cell ids
    this->reset_global_cell_indices();
  }

#else

  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::update_number_cache()
  {
    Assert(false, ExcNeedsMPI());
  }

#endif

  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::update_reference_cells()
  {
    // run algorithm for locally-owned cells
    dealii::Triangulation<dim, spacedim>::update_reference_cells();

    // translate ReferenceCell to unsigned int (needed by
    // Utilities::MPI::compute_set_union)
    std::vector<unsigned int> reference_cells_ui;

    reference_cells_ui.reserve(this->reference_cells.size());
    for (const auto &i : this->reference_cells)
      reference_cells_ui.push_back(static_cast<unsigned int>(i));

    // create union
    reference_cells_ui =
      Utilities::MPI::compute_set_union(reference_cells_ui,
                                        this->mpi_communicator);

    // transform back and store result
    this->reference_cells.clear();
    for (const auto &i : reference_cells_ui)
      this->reference_cells.emplace_back(
        dealii::internal::make_reference_cell_from_int(i));
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  types::subdomain_id
    TriangulationBase<dim, spacedim>::locally_owned_subdomain() const
  {
    return my_subdomain;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  const std::set<types::subdomain_id>
    &TriangulationBase<dim, spacedim>::ghost_owners() const
  {
    return number_cache.ghost_owners;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  const std::set<types::subdomain_id>
    &TriangulationBase<dim, spacedim>::level_ghost_owners() const
  {
    return number_cache.level_ghost_owners;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  std::vector<types::boundary_id> TriangulationBase<dim, spacedim>::
    get_boundary_ids() const
  {
    return Utilities::MPI::compute_set_union(
      dealii::Triangulation<dim, spacedim>::get_boundary_ids(),
      this->mpi_communicator);
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  std::vector<types::manifold_id> TriangulationBase<dim, spacedim>::
    get_manifold_ids() const
  {
    return Utilities::MPI::compute_set_union(
      dealii::Triangulation<dim, spacedim>::get_manifold_ids(),
      this->mpi_communicator);
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::reset_global_cell_indices()
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false, ExcNeedsMPI());
#else
    if (const auto pst =
          dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            this))
      if (pst->with_artificial_cells() == false)
        {
          // Specialization for parallel::shared::Triangulation without
          // artificial cells. The code below only works if a halo of a single
          // ghost cells is needed.

          std::vector<unsigned int> cell_counter(n_subdomains + 1);

          // count number of cells of each process
          for (const auto &cell : this->active_cell_iterators())
            cell_counter[cell->subdomain_id() + 1]++;

          // take prefix sum to obtain offset of each process
          for (unsigned int i = 0; i < n_subdomains; ++i)
            cell_counter[i + 1] += cell_counter[i];

          AssertDimension(cell_counter.back(), this->n_active_cells());

          // create partitioners
          IndexSet is_local(this->n_active_cells());
          is_local.add_range(cell_counter[my_subdomain],
                             cell_counter[my_subdomain + 1]);
          number_cache.active_cell_index_partitioner =
            std::make_shared<const Utilities::MPI::Partitioner>(
              is_local,
              complete_index_set(this->n_active_cells()),
              this->mpi_communicator);

          // set global active cell indices and increment process-local counters
          for (const auto &cell : this->active_cell_iterators())
            cell->set_global_active_cell_index(
              cell_counter[cell->subdomain_id()]++);

          Assert(this->is_multilevel_hierarchy_constructed() == false,
                 ExcNotImplemented());

          return;
        }

    // 1) determine number of active locally-owned cells
    const types::global_cell_index n_locally_owned_cells =
      this->n_locally_owned_active_cells();

    // 2) determine the offset of each process
    types::global_cell_index cell_index = 0;

    const int ierr = MPI_Exscan(
      &n_locally_owned_cells,
      &cell_index,
      1,
      Utilities::MPI::mpi_type_id_for_type<decltype(n_locally_owned_cells)>,
      MPI_SUM,
      this->mpi_communicator);
    AssertThrowMPI(ierr);

    // 3) give global indices to locally-owned cells and mark all other cells as
    //    invalid
    std::pair<types::global_cell_index, types::global_cell_index> my_range;
    my_range.first = cell_index;

    for (const auto &cell : this->active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_global_active_cell_index(cell_index++);
      else
        cell->set_global_active_cell_index(numbers::invalid_dof_index);

    my_range.second = cell_index;

    // 4) determine the global indices of ghost cells
    std::vector<types::global_dof_index> is_ghost_vector;
    GridTools::exchange_cell_data_to_ghosts<types::global_cell_index>(
      static_cast<dealii::Triangulation<dim, spacedim> &>(*this),
      [](const auto &cell) { return cell->global_active_cell_index(); },
      [&is_ghost_vector](const auto &cell, const auto &id) {
        cell->set_global_active_cell_index(id);
        is_ghost_vector.push_back(id);
      });

    // 5) set up new partitioner
    IndexSet is_local(this->n_global_active_cells());
    is_local.add_range(my_range.first, my_range.second);

    std::sort(is_ghost_vector.begin(), is_ghost_vector.end());
    IndexSet is_ghost(this->n_global_active_cells());
    is_ghost.add_indices(is_ghost_vector.begin(), is_ghost_vector.end());

    number_cache.active_cell_index_partitioner =
      std::make_shared<const Utilities::MPI::Partitioner>(
        is_local, is_ghost, this->mpi_communicator);

    // 6) proceed with multigrid levels if requested
    if (this->is_multilevel_hierarchy_constructed() == true)
      {
        // 1) determine number of locally-owned cells on levels
        std::vector<types::global_cell_index> n_cells_level(
          this->n_global_levels(), 0);

        for (auto cell : this->cell_iterators())
          if (cell->level_subdomain_id() == this->locally_owned_subdomain())
            n_cells_level[cell->level()]++;

        // 2) determine the offset of each process
        std::vector<types::global_cell_index> cell_index(
          this->n_global_levels(), 0);

        int ierr = MPI_Exscan(
          n_cells_level.data(),
          cell_index.data(),
          this->n_global_levels(),
          Utilities::MPI::mpi_type_id_for_type<decltype(*n_cells_level.data())>,
          MPI_SUM,
          this->mpi_communicator);
        AssertThrowMPI(ierr);

        // 3) determine global number of "active" cells on each level
        Utilities::MPI::sum(n_cells_level,
                            this->mpi_communicator,
                            n_cells_level);

        // 4) give global indices to locally-owned cells on level and mark
        //    all other cells as invalid
        std::vector<
          std::pair<types::global_cell_index, types::global_cell_index>>
          my_ranges(this->n_global_levels());
        for (unsigned int l = 0; l < this->n_global_levels(); ++l)
          my_ranges[l].first = cell_index[l];

        for (auto cell : this->cell_iterators())
          if (cell->level_subdomain_id() == this->locally_owned_subdomain())
            cell->set_global_level_cell_index(cell_index[cell->level()]++);
          else
            cell->set_global_level_cell_index(numbers::invalid_dof_index);

        for (unsigned int l = 0; l < this->n_global_levels(); ++l)
          my_ranges[l].second = cell_index[l];

        // 5) update the numbers of ghost level cells
        std::vector<std::vector<types::global_dof_index>> is_ghost_vectors(
          this->n_global_levels());
        GridTools::exchange_cell_data_to_level_ghosts<
          types::global_cell_index,
          dealii::Triangulation<dim, spacedim>>(
          *this,
          [](const auto &cell) { return cell->global_level_cell_index(); },
          [&is_ghost_vectors](const auto &cell, const auto &id) {
            cell->set_global_level_cell_index(id);
            is_ghost_vectors[cell->level()].push_back(id);
          });

        number_cache.level_cell_index_partitioners.resize(
          this->n_global_levels());

        // 6) set up cell partitioners for each level
        for (unsigned int l = 0; l < this->n_global_levels(); ++l)
          {
            IndexSet is_local(n_cells_level[l]);
            is_local.add_range(my_ranges[l].first, my_ranges[l].second);

            IndexSet is_ghost(n_cells_level[l]);
            std::sort(is_ghost_vectors[l].begin(), is_ghost_vectors[l].end());
            is_ghost.add_indices(is_ghost_vectors[l].begin(),
                                 is_ghost_vectors[l].end());

            number_cache.level_cell_index_partitioners[l] =
              std::make_shared<const Utilities::MPI::Partitioner>(
                is_local, is_ghost, this->mpi_communicator);
          }
      }

#endif
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void TriangulationBase<dim, spacedim>::communicate_locally_moved_vertices(
    const std::vector<bool> &vertex_locally_moved)
  {
    AssertDimension(vertex_locally_moved.size(), this->n_vertices());
#ifdef DEBUG
    {
      const std::vector<bool> locally_owned_vertices =
        dealii::GridTools::get_locally_owned_vertices(*this);
      for (unsigned int i = 0; i < locally_owned_vertices.size(); ++i)
        Assert((vertex_locally_moved[i] == false) ||
                 (locally_owned_vertices[i] == true),
               ExcMessage("The vertex_locally_moved argument must not "
                          "contain vertices that are not locally owned"));
    }
#endif

    Point<spacedim> invalid_point;
    for (unsigned int d = 0; d < spacedim; ++d)
      invalid_point[d] = std::numeric_limits<double>::quiet_NaN();

    const auto pack = [&](const auto &cell) {
      std::vector<Point<spacedim>> vertices(cell->n_vertices());

      for (const auto v : cell->vertex_indices())
        if (vertex_locally_moved[cell->vertex_index(v)])
          vertices[v] = cell->vertex(v);
        else
          vertices[v] = invalid_point;

      return vertices;
    };

    const auto unpack = [&](const auto &cell, const auto &vertices) {
      for (const auto v : cell->vertex_indices())
        if (numbers::is_nan(vertices[v][0]) == false)
          cell->vertex(v) = vertices[v];
    };

    if (this->is_multilevel_hierarchy_constructed())
      GridTools::exchange_cell_data_to_level_ghosts<
        std::vector<Point<spacedim>>>(
        static_cast<dealii::Triangulation<dim, spacedim> &>(*this),
        pack,
        unpack);
    else
      GridTools::exchange_cell_data_to_ghosts<std::vector<Point<spacedim>>>(
        static_cast<dealii::Triangulation<dim, spacedim> &>(*this),
        pack,
        unpack);
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  std::weak_ptr<const Utilities::MPI::Partitioner> TriangulationBase<
    dim,
    spacedim>::global_active_cell_index_partitioner() const
  {
    return number_cache.active_cell_index_partitioner;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  std::weak_ptr<const Utilities::MPI::Partitioner> TriangulationBase<dim,
                                                                     spacedim>::
    global_level_cell_index_partitioner(const unsigned int level) const
  {
    Assert(this->is_multilevel_hierarchy_constructed(), ExcNotImplemented());
    AssertIndexRange(level, this->n_global_levels());

    return number_cache.level_cell_index_partitioners[level];
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  types::coarse_cell_id
    TriangulationBase<dim, spacedim>::n_global_coarse_cells() const
  {
    return number_cache.number_of_global_coarse_cells;
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  DistributedTriangulationBase<dim, spacedim>::DistributedTriangulationBase(
    const MPI_Comm mpi_communicator,
    const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
               smooth_grid,
    const bool check_for_distorted_cells)
    : dealii::parallel::TriangulationBase<dim, spacedim>(
        mpi_communicator,
        smooth_grid,
        check_for_distorted_cells)
  {}



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  void DistributedTriangulationBase<dim, spacedim>::clear()
  {
    dealii::Triangulation<dim, spacedim>::clear();
  }



  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  bool DistributedTriangulationBase<dim, spacedim>::has_hanging_nodes() const
  {
    if (this->n_global_levels() <= 1)
      return false; // can not have hanging nodes without refined cells

    // if there are any active cells with level less than n_global_levels()-1,
    // then there is obviously also one with level n_global_levels()-1, and
    // consequently there must be a hanging node somewhere.
    //
    // The problem is that we cannot just ask for the first active cell, but
    // instead need to filter over locally owned cells.
    const bool have_coarser_cell =
      std::any_of(this->begin_active(this->n_global_levels() - 2),
                  this->end_active(this->n_global_levels() - 2),
                  [](const CellAccessor<dim, spacedim> &cell) {
                    return cell.is_locally_owned();
                  });

    // return true if at least one process has a coarser cell
    return Utilities::MPI::max(have_coarser_cell ? 1 : 0,
                               this->mpi_communicator) != 0;
  }


} // end namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "tria_base.inst"

DEAL_II_NAMESPACE_CLOSE
