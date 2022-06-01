// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2022 by the deal.II authors
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

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector_memory.h>

#include <algorithm>
#include <cstdint>
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
      MemoryConsumption::memory_consumption(this->mpi_communicator) +
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
    , number_of_global_coarse_cells(0)
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
  void
  TriangulationBase<dim, spacedim>::update_number_cache()
  {
    Assert(false, ExcNeedsMPI());
  }

#endif

  template <int dim, int spacedim>
  void
  TriangulationBase<dim, spacedim>::update_reference_cells()
  {
    // run algorithm for locally-owned cells
    dealii::Triangulation<dim, spacedim>::update_reference_cells();

    // translate ReferenceCell to unsigned int (needed by
    // Utilities::MPI::compute_set_union)
    std::vector<unsigned int> reference_cells_ui;

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
        dealii::internal::ReferenceCell::make_reference_cell_from_int(i));
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
      *this,
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
  void
  TriangulationBase<dim, spacedim>::communicate_locally_moved_vertices(
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
        if (std::isnan(vertices[v][0]) == false)
          cell->vertex(v) = vertices[v];
    };

    if (this->is_multilevel_hierarchy_constructed())
      GridTools::exchange_cell_data_to_level_ghosts<
        std::vector<Point<spacedim>>>(*this, pack, unpack);
    else
      GridTools::exchange_cell_data_to_ghosts<std::vector<Point<spacedim>>>(
        *this, pack, unpack);
  }



  template <int dim, int spacedim>
  const std::weak_ptr<const Utilities::MPI::Partitioner>
  TriangulationBase<dim, spacedim>::global_active_cell_index_partitioner() const
  {
    return number_cache.active_cell_index_partitioner;
  }



  template <int dim, int spacedim>
  const std::weak_ptr<const Utilities::MPI::Partitioner>
  TriangulationBase<dim, spacedim>::global_level_cell_index_partitioner(
    const unsigned int level) const
  {
    Assert(this->is_multilevel_hierarchy_constructed(), ExcNotImplemented());
    AssertIndexRange(level, this->n_global_levels());

    return number_cache.level_cell_index_partitioners[level];
  }



  template <int dim, int spacedim>
  types::coarse_cell_id
  TriangulationBase<dim, spacedim>::n_global_coarse_cells() const
  {
    return number_cache.number_of_global_coarse_cells;
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
    , cell_attached_data({0, 0, {}, {}})
    , data_transfer(mpi_communicator)
  {}



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::clear()
  {
    cell_attached_data = {0, 0, {}, {}};
    data_transfer.clear();

    dealii::Triangulation<dim, spacedim>::clear();
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::save_attached_data(
    const unsigned int global_first_cell,
    const unsigned int global_num_cells,
    const std::string &filename) const
  {
    // cast away constness
    auto tria = const_cast<
      dealii::parallel::DistributedTriangulationBase<dim, spacedim> *>(this);

    if (this->cell_attached_data.n_attached_data_sets > 0)
      {
        // pack attached data first
        tria->data_transfer.pack_data(
          tria->local_cell_relations,
          tria->cell_attached_data.pack_callbacks_fixed,
          tria->cell_attached_data.pack_callbacks_variable);

        // then store buffers in file
        tria->data_transfer.save(global_first_cell, global_num_cells, filename);

        // and release the memory afterwards
        tria->data_transfer.clear();
      }

    // clear all of the callback data, as explained in the documentation of
    // register_data_attach()
    {
      tria->cell_attached_data.n_attached_data_sets = 0;
      tria->cell_attached_data.pack_callbacks_fixed.clear();
      tria->cell_attached_data.pack_callbacks_variable.clear();
    }
  }



  template <int dim, int spacedim>
  bool
  DistributedTriangulationBase<dim, spacedim>::has_hanging_nodes() const
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



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::load_attached_data(
    const unsigned int global_first_cell,
    const unsigned int global_num_cells,
    const unsigned int local_num_cells,
    const std::string &filename,
    const unsigned int n_attached_deserialize_fixed,
    const unsigned int n_attached_deserialize_variable)
  {
    // load saved data, if any was stored
    if (this->cell_attached_data.n_attached_deserialize > 0)
      {
        this->data_transfer.load(global_first_cell,
                                 global_num_cells,
                                 local_num_cells,
                                 filename,
                                 n_attached_deserialize_fixed,
                                 n_attached_deserialize_variable);

        this->data_transfer.unpack_cell_status(this->local_cell_relations);

        // the CellStatus of all stored cells should always be CELL_PERSIST.
        for (const auto &cell_rel : this->local_cell_relations)
          {
            (void)cell_rel;
            Assert(
              (cell_rel.second == // cell_status
               parallel::DistributedTriangulationBase<dim,
                                                      spacedim>::CELL_PERSIST),
              ExcInternalError());
          }
      }
  }



  template <int dim, int spacedim>
  unsigned int
  DistributedTriangulationBase<dim, spacedim>::register_data_attach(
    const std::function<std::vector<char>(const cell_iterator &,
                                          const CellStatus)> &pack_callback,
    const bool returns_variable_size_data)
  {
    unsigned int handle = numbers::invalid_unsigned_int;

    // Add new callback function to the corresponding register.
    // Encode handles according to returns_variable_size_data.
    if (returns_variable_size_data)
      {
        handle = 2 * cell_attached_data.pack_callbacks_variable.size();
        cell_attached_data.pack_callbacks_variable.push_back(pack_callback);
      }
    else
      {
        handle = 2 * cell_attached_data.pack_callbacks_fixed.size() + 1;
        cell_attached_data.pack_callbacks_fixed.push_back(pack_callback);
      }

    // Increase overall counter.
    ++cell_attached_data.n_attached_data_sets;

    return handle;
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::notify_ready_to_unpack(
    const unsigned int handle,
    const std::function<
      void(const cell_iterator &,
           const CellStatus,
           const boost::iterator_range<std::vector<char>::const_iterator> &)>
      &unpack_callback)
  {
    // perform unpacking
    data_transfer.unpack_data(local_cell_relations, handle, unpack_callback);

    // decrease counters
    --cell_attached_data.n_attached_data_sets;
    if (cell_attached_data.n_attached_deserialize > 0)
      --cell_attached_data.n_attached_deserialize;

    // important: only remove data if we are not in the deserialization
    // process. There, each SolutionTransfer registers and unpacks before
    // the next one does this, so n_attached_data_sets is only 1 here.  This
    // would destroy the saved data before the second SolutionTransfer can
    // get it. This created a bug that is documented in
    // tests/mpi/p4est_save_03 with more than one SolutionTransfer.

    if (cell_attached_data.n_attached_data_sets == 0 &&
        cell_attached_data.n_attached_deserialize == 0)
      {
        // everybody got their data, time for cleanup!
        cell_attached_data.pack_callbacks_fixed.clear();
        cell_attached_data.pack_callbacks_variable.clear();
        data_transfer.clear();

        // reset all cell_status entries after coarsening/refinement
        for (auto &cell_rel : local_cell_relations)
          cell_rel.second =
            parallel::DistributedTriangulationBase<dim, spacedim>::CELL_PERSIST;
      }
  }



  /* ------------------ class DataTransfer<dim,spacedim> ----------------- */


  template <int dim, int spacedim>
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::DataTransfer(
    const MPI_Comm &mpi_communicator)
    : variable_size_data_stored(false)
    , mpi_communicator(mpi_communicator)
  {}



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::pack_data(
    const std::vector<cell_relation_t> &cell_relations,
    const std::vector<typename CellAttachedData::pack_callback_t>
      &pack_callbacks_fixed,
    const std::vector<typename CellAttachedData::pack_callback_t>
      &pack_callbacks_variable)
  {
    Assert(src_data_fixed.size() == 0,
           ExcMessage("Previously packed data has not been released yet!"));
    Assert(src_sizes_variable.size() == 0, ExcInternalError());

    const unsigned int n_callbacks_fixed    = pack_callbacks_fixed.size();
    const unsigned int n_callbacks_variable = pack_callbacks_variable.size();

    // Store information that we packed variable size data in
    // a member variable for later.
    variable_size_data_stored = (n_callbacks_variable > 0);

    // If variable transfer is scheduled, we will store the data size that
    // each variable size callback function writes in this auxiliary
    // container. The information will be stored by each cell in this vector
    // temporarily.
    std::vector<unsigned int> cell_sizes_variable_cumulative(
      n_callbacks_variable);

    // Prepare the buffer structure, in which each callback function will
    // store its data for each active cell.
    // The outmost shell in this container construct corresponds to the
    // data packed per cell. The next layer resembles the data that
    // each callback function packs on the corresponding cell. These
    // buffers are chains of chars stored in an std::vector<char>.
    // A visualisation of the data structure:
    /* clang-format off */
      // |             cell_1                | |             cell_2                | ...
      // ||  callback_1  ||  callback_2  |...| ||  callback_1  ||  callback_2  |...| ...
      // |||char|char|...|||char|char|...|...| |||char|char|...|||char|char|...|...| ...
    /* clang-format on */
    std::vector<std::vector<std::vector<char>>> packed_fixed_size_data(
      cell_relations.size());
    std::vector<std::vector<std::vector<char>>> packed_variable_size_data(
      variable_size_data_stored ? cell_relations.size() : 0);

    //
    // --------- Pack data for fixed and variable size transfer ---------
    //
    // Iterate over all cells, call all callback functions on each cell,
    // and store their data in the corresponding buffer scope.
    {
      auto cell_rel_it           = cell_relations.cbegin();
      auto data_cell_fixed_it    = packed_fixed_size_data.begin();
      auto data_cell_variable_it = packed_variable_size_data.begin();
      for (; cell_rel_it != cell_relations.cend(); ++cell_rel_it)
        {
          const auto &dealii_cell = cell_rel_it->first;
          const auto &cell_status = cell_rel_it->second;

          // Assertions about the tree structure.
          switch (cell_status)
            {
              case parallel::DistributedTriangulationBase<dim, spacedim>::
                CELL_PERSIST:
              case parallel::DistributedTriangulationBase<dim, spacedim>::
                CELL_REFINE:
                // double check the condition that we will only ever attach
                // data to active cells when we get here
                Assert(dealii_cell->is_active(), ExcInternalError());
                break;

              case parallel::DistributedTriangulationBase<dim, spacedim>::
                CELL_COARSEN:
                // double check the condition that we will only ever attach
                // data to cells with children when we get here. however, we
                // can only tolerate one level of coarsening at a time, so
                // check that the children are all active
                Assert(dealii_cell->is_active() == false, ExcInternalError());
                for (unsigned int c = 0;
                     c < GeometryInfo<dim>::max_children_per_cell;
                     ++c)
                  Assert(dealii_cell->child(c)->is_active(),
                         ExcInternalError());
                break;

              case parallel::DistributedTriangulationBase<dim, spacedim>::
                CELL_INVALID:
                // do nothing on invalid cells
                break;

              default:
                Assert(false, ExcInternalError());
                break;
            }

          // Reserve memory corresponding to the number of callback
          // functions that will be called.
          // If variable size transfer is scheduled, we need to leave
          // room for an array that holds information about how many
          // bytes each of the variable size callback functions will
          // write.
          // On cells flagged with CELL_INVALID, only its CellStatus
          // will be stored.
          const unsigned int n_fixed_size_data_sets_on_cell =
            1 +
            ((cell_status ==
              parallel::DistributedTriangulationBase<dim,
                                                     spacedim>::CELL_INVALID) ?
               0 :
               ((variable_size_data_stored ? 1 : 0) + n_callbacks_fixed));
          data_cell_fixed_it->resize(n_fixed_size_data_sets_on_cell);

          // We continue with packing all data on this specific cell.
          auto data_fixed_it = data_cell_fixed_it->begin();

          // First, we pack the CellStatus information.
          // to get consistent data sizes on each cell for the fixed size
          // transfer, we won't allow compression
          *data_fixed_it =
            Utilities::pack(cell_status, /*allow_compression=*/false);
          ++data_fixed_it;

          // Proceed with all registered callback functions.
          // Skip cells with the CELL_INVALID flag.
          if (cell_status !=
              parallel::DistributedTriangulationBase<dim,
                                                     spacedim>::CELL_INVALID)
            {
              // Pack fixed size data.
              for (auto callback_it = pack_callbacks_fixed.cbegin();
                   callback_it != pack_callbacks_fixed.cend();
                   ++callback_it, ++data_fixed_it)
                {
                  *data_fixed_it = (*callback_it)(dealii_cell, cell_status);
                }

              // Pack variable size data.
              // If we store variable size data, we need to transfer
              // the sizes of each corresponding callback function
              // via fixed size transfer as well.
              if (variable_size_data_stored)
                {
                  const unsigned int n_variable_size_data_sets_on_cell =
                    ((cell_status ==
                      parallel::DistributedTriangulationBase<dim, spacedim>::
                        CELL_INVALID) ?
                       0 :
                       n_callbacks_variable);
                  data_cell_variable_it->resize(
                    n_variable_size_data_sets_on_cell);

                  auto callback_it      = pack_callbacks_variable.cbegin();
                  auto data_variable_it = data_cell_variable_it->begin();
                  auto sizes_variable_it =
                    cell_sizes_variable_cumulative.begin();
                  for (; callback_it != pack_callbacks_variable.cend();
                       ++callback_it, ++data_variable_it, ++sizes_variable_it)
                    {
                      *data_variable_it =
                        (*callback_it)(dealii_cell, cell_status);

                      // Store data sizes for each callback function first.
                      // Make it cumulative below.
                      *sizes_variable_it = data_variable_it->size();
                    }

                  // Turn size vector into its cumulative representation.
                  std::partial_sum(cell_sizes_variable_cumulative.begin(),
                                   cell_sizes_variable_cumulative.end(),
                                   cell_sizes_variable_cumulative.begin());

                  // Serialize cumulative variable size vector
                  // value-by-value. This way we can circumvent the overhead
                  // of storing the container object as a whole, since we
                  // know its size by the number of registered callback
                  // functions.
                  data_fixed_it->resize(n_callbacks_variable *
                                        sizeof(unsigned int));
                  for (unsigned int i = 0; i < n_callbacks_variable; ++i)
                    std::memcpy(&(data_fixed_it->at(i * sizeof(unsigned int))),
                                &(cell_sizes_variable_cumulative.at(i)),
                                sizeof(unsigned int));

                  ++data_fixed_it;
                }

              // Double check that we packed everything we wanted
              // in the fixed size buffers.
              Assert(data_fixed_it == data_cell_fixed_it->end(),
                     ExcInternalError());
            }

          ++data_cell_fixed_it;

          // Increment the variable size data iterator
          // only if we actually pack this kind of data
          // to avoid getting out of bounds.
          if (variable_size_data_stored)
            ++data_cell_variable_it;
        } // loop over cell_relations
    }

    //
    // ----------- Gather data sizes for fixed size transfer ------------
    //
    // Generate a vector which stores the sizes of each callback function,
    // including the packed CellStatus transfer.
    // Find the very first cell that we wrote to with all callback
    // functions (i.e. a cell that was not flagged with CELL_INVALID)
    // and store the sizes of each buffer.
    //
    // To deal with the case that at least one of the processors does not
    // own any cell at all, we will exchange the information about the data
    // sizes among them later. The code in between is still well-defined,
    // since the following loops will be skipped.
    std::vector<unsigned int> local_sizes_fixed(
      1 + n_callbacks_fixed + (variable_size_data_stored ? 1 : 0));
    for (const auto &data_cell : packed_fixed_size_data)
      {
        if (data_cell.size() == local_sizes_fixed.size())
          {
            auto sizes_fixed_it = local_sizes_fixed.begin();
            auto data_fixed_it  = data_cell.cbegin();
            for (; data_fixed_it != data_cell.cend();
                 ++data_fixed_it, ++sizes_fixed_it)
              {
                *sizes_fixed_it = data_fixed_it->size();
              }

            break;
          }
      }

    // Check if all cells have valid sizes.
    for (auto data_cell_fixed_it = packed_fixed_size_data.cbegin();
         data_cell_fixed_it != packed_fixed_size_data.cend();
         ++data_cell_fixed_it)
      {
        Assert((data_cell_fixed_it->size() == 1) ||
                 (data_cell_fixed_it->size() == local_sizes_fixed.size()),
               ExcInternalError());
      }

    // Share information about the packed data sizes
    // of all callback functions across all processors, in case one
    // of them does not own any cells at all.
    std::vector<unsigned int> global_sizes_fixed(local_sizes_fixed.size());
    Utilities::MPI::max(local_sizes_fixed,
                        this->mpi_communicator,
                        global_sizes_fixed);

    // Construct cumulative sizes, since this is the only information
    // we need from now on.
    sizes_fixed_cumulative.resize(global_sizes_fixed.size());
    std::partial_sum(global_sizes_fixed.begin(),
                     global_sizes_fixed.end(),
                     sizes_fixed_cumulative.begin());

    //
    // ---------- Gather data sizes for variable size transfer ----------
    //
    if (variable_size_data_stored)
      {
        src_sizes_variable.reserve(packed_variable_size_data.size());
        for (const auto &data_cell : packed_variable_size_data)
          {
            int variable_data_size_on_cell = 0;

            for (const auto &data : data_cell)
              variable_data_size_on_cell += data.size();

            src_sizes_variable.push_back(variable_data_size_on_cell);
          }
      }

    //
    // ------------------------ Build buffers ---------------------------
    //
    const unsigned int expected_size_fixed =
      cell_relations.size() * sizes_fixed_cumulative.back();
    const unsigned int expected_size_variable =
      std::accumulate(src_sizes_variable.begin(),
                      src_sizes_variable.end(),
                      std::vector<int>::size_type(0));

    // Move every piece of packed fixed size data into the consecutive
    // buffer.
    src_data_fixed.reserve(expected_size_fixed);
    for (const auto &data_cell_fixed : packed_fixed_size_data)
      {
        // Move every fraction of packed data into the buffer
        // reserved for this particular cell.
        for (const auto &data_fixed : data_cell_fixed)
          std::move(data_fixed.begin(),
                    data_fixed.end(),
                    std::back_inserter(src_data_fixed));

        // If we only packed the CellStatus information
        // (i.e. encountered a cell flagged CELL_INVALID),
        // fill the remaining space with invalid entries.
        // We can skip this if there is nothing else to pack.
        if ((data_cell_fixed.size() == 1) &&
            (sizes_fixed_cumulative.size() > 1))
          {
            const std::size_t bytes_skipped =
              sizes_fixed_cumulative.back() - sizes_fixed_cumulative.front();

            src_data_fixed.insert(src_data_fixed.end(),
                                  bytes_skipped,
                                  static_cast<char>(-1)); // invalid_char
          }
      }

    // Move every piece of packed variable size data into the consecutive
    // buffer.
    if (variable_size_data_stored)
      {
        src_data_variable.reserve(expected_size_variable);
        for (const auto &data_cell : packed_variable_size_data)
          {
            // Move every fraction of packed data into the buffer
            // reserved for this particular cell.
            for (const auto &data : data_cell)
              std::move(data.begin(),
                        data.end(),
                        std::back_inserter(src_data_variable));
          }
      }

    // Double check that we packed everything correctly.
    Assert(src_data_fixed.size() == expected_size_fixed, ExcInternalError());
    Assert(src_data_variable.size() == expected_size_variable,
           ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::unpack_cell_status(
    std::vector<cell_relation_t> &cell_relations) const
  {
    Assert(sizes_fixed_cumulative.size() > 0,
           ExcMessage("No data has been packed!"));
    if (cell_relations.size() > 0)
      {
        Assert(dest_data_fixed.size() > 0,
               ExcMessage("No data has been received!"));
      }

    // Size of CellStatus object that will be unpacked on each cell.
    const unsigned int size = sizes_fixed_cumulative.front();

    // Iterate over all cells and overwrite the CellStatus
    // information from the transferred data.
    // Proceed buffer iterator position to next cell after
    // each iteration.
    auto cell_rel_it   = cell_relations.begin();
    auto dest_fixed_it = dest_data_fixed.cbegin();
    for (; cell_rel_it != cell_relations.end();
         ++cell_rel_it, dest_fixed_it += sizes_fixed_cumulative.back())
      {
        cell_rel_it->second = // cell_status
          Utilities::unpack<typename parallel::DistributedTriangulationBase<
            dim,
            spacedim>::CellStatus>(dest_fixed_it,
                                   dest_fixed_it + size,
                                   /*allow_compression=*/false);
      }
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::unpack_data(
    const std::vector<cell_relation_t> &cell_relations,
    const unsigned int                  handle,
    const std::function<
      void(const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
           const typename dealii::Triangulation<dim, spacedim>::CellStatus &,
           const boost::iterator_range<std::vector<char>::const_iterator> &)>
      &unpack_callback) const
  {
    // We decode the handle returned by register_data_attach() back into
    // a format we can use. All even handles belong to those callback
    // functions which write/read variable size data, all odd handles
    // interact with fixed size buffers.
    const bool         callback_variable_transfer = (handle % 2 == 0);
    const unsigned int callback_index             = handle / 2;

    // Cells will always receive fixed size data (i.e., CellStatus
    // information), but not necessarily variable size data (e.g., with a
    // ParticleHandler a cell might not contain any particle at all).
    // Thus it is sufficient to check if fixed size data has been received.
    Assert(sizes_fixed_cumulative.size() > 0,
           ExcMessage("No data has been packed!"));
    if (cell_relations.size() > 0)
      {
        Assert(dest_data_fixed.size() > 0,
               ExcMessage("No data has been received!"));
      }

    std::vector<char>::const_iterator dest_data_it;
    std::vector<char>::const_iterator dest_sizes_cell_it;

    // Depending on whether our callback function unpacks fixed or
    // variable size data, we have to pursue different approaches
    // to localize the correct fraction of the buffer from which
    // we are allowed to read.
    unsigned int offset         = numbers::invalid_unsigned_int;
    unsigned int size           = numbers::invalid_unsigned_int;
    unsigned int data_increment = numbers::invalid_unsigned_int;

    if (callback_variable_transfer)
      {
        // For the variable size data, we need to extract the
        // data size from the fixed size buffer on each cell.
        //
        // We packed this information last, so the last packed
        // object in the fixed size buffer corresponds to the
        // variable data sizes.
        //
        // The last entry of sizes_fixed_cumulative corresponds
        // to the size of all fixed size data packed on the cell.
        // To get the offset for the last packed object, we need
        // to get the next-to-last entry.
        const unsigned int offset_variable_data_sizes =
          sizes_fixed_cumulative[sizes_fixed_cumulative.size() - 2];

        // This iterator points to the data size that the
        // callback_function packed for each specific cell.
        // Adjust buffer iterator to the offset of the callback
        // function so that we only have to advance its position
        // to the next cell after each iteration.
        dest_sizes_cell_it = dest_data_fixed.cbegin() +
                             offset_variable_data_sizes +
                             callback_index * sizeof(unsigned int);

        // Let the data iterator point to the correct buffer.
        dest_data_it = dest_data_variable.cbegin();
      }
    else
      {
        // For the fixed size data, we can get the information about
        // the buffer location on each cell directly from the
        // sizes_fixed_cumulative vector.
        offset         = sizes_fixed_cumulative[callback_index];
        size           = sizes_fixed_cumulative[callback_index + 1] - offset;
        data_increment = sizes_fixed_cumulative.back();

        // Let the data iterator point to the correct buffer.
        // Adjust buffer iterator to the offset of the callback
        // function so that we only have to advance its position
        // to the next cell after each iteration.
        if (cell_relations.begin() != cell_relations.end())
          dest_data_it = dest_data_fixed.cbegin() + offset;
      }

    // Iterate over all cells and unpack the transferred data.
    auto cell_rel_it   = cell_relations.begin();
    auto dest_sizes_it = dest_sizes_variable.cbegin();
    for (; cell_rel_it != cell_relations.end(); ++cell_rel_it)
      {
        const auto &dealii_cell = cell_rel_it->first;
        const auto &cell_status = cell_rel_it->second;

        if (callback_variable_transfer)
          {
            // Update the increment according to the whole data size
            // of the current cell.
            data_increment = *dest_sizes_it;

            if (cell_status !=
                parallel::DistributedTriangulationBase<dim,
                                                       spacedim>::CELL_INVALID)
              {
                // Extract the corresponding values for offset and size from
                // the cumulative sizes array stored in the fixed size
                // buffer.
                if (callback_index == 0)
                  offset = 0;
                else
                  std::memcpy(&offset,
                              &(*(dest_sizes_cell_it - sizeof(unsigned int))),
                              sizeof(unsigned int));

                std::memcpy(&size,
                            &(*dest_sizes_cell_it),
                            sizeof(unsigned int));

                size -= offset;

                // Move the data iterator to the corresponding position
                // of the callback function and adjust the increment
                // accordingly.
                dest_data_it += offset;
                data_increment -= offset;
              }

            // Advance data size iterators to the next cell, avoid iterating
            // past the end of dest_sizes_cell_it
            if (cell_rel_it != cell_relations.end() - 1)
              dest_sizes_cell_it += sizes_fixed_cumulative.back();
            ++dest_sizes_it;
          }

        switch (cell_status)
          {
            case parallel::DistributedTriangulationBase<dim,
                                                        spacedim>::CELL_PERSIST:
            case parallel::DistributedTriangulationBase<dim,
                                                        spacedim>::CELL_COARSEN:
              unpack_callback(dealii_cell,
                              cell_status,
                              boost::make_iterator_range(dest_data_it,
                                                         dest_data_it + size));
              break;

            case parallel::DistributedTriangulationBase<dim,
                                                        spacedim>::CELL_REFINE:
              unpack_callback(dealii_cell->parent(),
                              cell_status,
                              boost::make_iterator_range(dest_data_it,
                                                         dest_data_it + size));
              break;

            case parallel::DistributedTriangulationBase<dim,
                                                        spacedim>::CELL_INVALID:
              // Skip this cell.
              break;

            default:
              Assert(false, ExcInternalError());
              break;
          }

        if (cell_rel_it != cell_relations.end() - 1)
          dest_data_it += data_increment;
      }
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::save(
    const unsigned int global_first_cell,
    const unsigned int global_num_cells,
    const std::string &filename) const
  {
#ifdef DEAL_II_WITH_MPI
    // Large fractions of this function have been copied from
    // DataOutInterface::write_vtu_in_parallel.
    // TODO: Write general MPIIO interface.

    Assert(sizes_fixed_cumulative.size() > 0,
           ExcMessage("No data has been packed!"));

    const int myrank = Utilities::MPI::this_mpi_process(mpi_communicator);

    const unsigned int bytes_per_cell = sizes_fixed_cumulative.back();

    //
    // ---------- Fixed size data ----------
    //
    {
      const std::string fname_fixed = std::string(filename) + "_fixed.data";

      MPI_Info info;
      int      ierr = MPI_Info_create(&info);
      AssertThrowMPI(ierr);

      MPI_File fh;
      ierr = MPI_File_open(mpi_communicator,
                           fname_fixed.c_str(),
                           MPI_MODE_CREATE | MPI_MODE_WRONLY,
                           info,
                           &fh);
      AssertThrowMPI(ierr);

      ierr = MPI_File_set_size(fh, 0); // delete the file contents
      AssertThrowMPI(ierr);
      // this barrier is necessary, because otherwise others might already
      // write while one core is still setting the size to zero.
      ierr = MPI_Barrier(mpi_communicator);
      AssertThrowMPI(ierr);
      ierr = MPI_Info_free(&info);
      AssertThrowMPI(ierr);
      // ------------------

      // Write cumulative sizes to file.
      // Since each processor owns the same information about the data
      // sizes, it is sufficient to let only the first processor perform
      // this task.
      if (myrank == 0)
        {
          ierr = Utilities::MPI::LargeCount::File_write_at_c(
            fh,
            0,
            sizes_fixed_cumulative.data(),
            sizes_fixed_cumulative.size(),
            MPI_UNSIGNED,
            MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
        }

      // Write packed data to file simultaneously.
      const MPI_Offset size_header =
        sizes_fixed_cumulative.size() * sizeof(unsigned int);

      // Make sure we do the following computation in 64bit integers to be
      // able to handle 4GB+ files:
      const MPI_Offset my_global_file_position =
        size_header +
        static_cast<MPI_Offset>(global_first_cell) * bytes_per_cell;

      ierr =
        Utilities::MPI::LargeCount::File_write_at_c(fh,
                                                    my_global_file_position,
                                                    src_data_fixed.data(),
                                                    src_data_fixed.size(),
                                                    MPI_BYTE,
                                                    MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);

      ierr = MPI_File_close(&fh);
      AssertThrowMPI(ierr);
    }



    //
    // ---------- Variable size data ----------
    //
    if (variable_size_data_stored)
      {
        const std::string fname_variable =
          std::string(filename) + "_variable.data";

        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        MPI_File fh;
        ierr = MPI_File_open(mpi_communicator,
                             fname_variable.c_str(),
                             MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             info,
                             &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_File_set_size(fh, 0); // delete the file contents
        AssertThrowMPI(ierr);
        // this barrier is necessary, because otherwise others might already
        // write while one core is still setting the size to zero.
        ierr = MPI_Barrier(mpi_communicator);
        AssertThrowMPI(ierr);
        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);

        // Write sizes of each cell into file simultaneously.
        {
          const MPI_Offset my_global_file_position =
            static_cast<MPI_Offset>(global_first_cell) * sizeof(unsigned int);

          // It is very unlikely that a single process has more than
          // 2 billion cells, but we might as well check.
          AssertThrow(src_sizes_variable.size() <
                        static_cast<std::size_t>(
                          std::numeric_limits<int>::max()),
                      ExcNotImplemented());

          ierr = Utilities::MPI::LargeCount::File_write_at_c(
            fh,
            my_global_file_position,
            src_sizes_variable.data(),
            src_sizes_variable.size(),
            MPI_INT,
            MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
        }

        // Gather size of data in bytes we want to store from this
        // processor and compute the prefix sum. We do this in 64 bit
        // to avoid overflow for files larger than 4GB:
        const std::uint64_t size_on_proc = src_data_variable.size();
        std::uint64_t       prefix_sum   = 0;
        ierr                             = MPI_Exscan(&size_on_proc,
                          &prefix_sum,
                          1,
                          MPI_UINT64_T,
                          MPI_SUM,
                          mpi_communicator);
        AssertThrowMPI(ierr);

        const MPI_Offset my_global_file_position =
          static_cast<MPI_Offset>(global_num_cells) * sizeof(unsigned int) +
          prefix_sum;

        // Write data consecutively into file.
        ierr =
          Utilities::MPI::LargeCount::File_write_at_c(fh,
                                                      my_global_file_position,
                                                      src_data_variable.data(),
                                                      src_data_variable.size(),
                                                      MPI_BYTE,
                                                      MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);


        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);
      }
#else
    (void)global_first_cell;
    (void)global_num_cells;
    (void)filename;

    AssertThrow(false, ExcNeedsMPI());
#endif
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::load(
    const unsigned int global_first_cell,
    const unsigned int global_num_cells,
    const unsigned int local_num_cells,
    const std::string &filename,
    const unsigned int n_attached_deserialize_fixed,
    const unsigned int n_attached_deserialize_variable)
  {
#ifdef DEAL_II_WITH_MPI
    // Large fractions of this function have been copied from
    // DataOutInterface::write_vtu_in_parallel.
    // TODO: Write general MPIIO interface.

    Assert(dest_data_fixed.size() == 0,
           ExcMessage("Previously loaded data has not been released yet!"));

    variable_size_data_stored = (n_attached_deserialize_variable > 0);

    //
    // ---------- Fixed size data ----------
    //
    {
      const std::string fname_fixed = std::string(filename) + "_fixed.data";

      MPI_Info info;
      int      ierr = MPI_Info_create(&info);
      AssertThrowMPI(ierr);

      MPI_File fh;
      ierr = MPI_File_open(
        mpi_communicator, fname_fixed.c_str(), MPI_MODE_RDONLY, info, &fh);
      AssertThrowMPI(ierr);

      ierr = MPI_Info_free(&info);
      AssertThrowMPI(ierr);

      // Read cumulative sizes from file.
      // Since all processors need the same information about the data
      // sizes, let each of them retrieve it by reading from the same
      // location in the file.
      sizes_fixed_cumulative.resize(1 + n_attached_deserialize_fixed +
                                    (variable_size_data_stored ? 1 : 0));
      ierr = Utilities::MPI::LargeCount::File_read_at_c(
        fh,
        0,
        sizes_fixed_cumulative.data(),
        sizes_fixed_cumulative.size(),
        MPI_UNSIGNED,
        MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);

      // Allocate sufficient memory.
      const unsigned int bytes_per_cell = sizes_fixed_cumulative.back();
      dest_data_fixed.resize(static_cast<size_t>(local_num_cells) *
                             bytes_per_cell);

      // Read packed data from file simultaneously.
      const MPI_Offset size_header =
        sizes_fixed_cumulative.size() * sizeof(unsigned int);

      // Make sure we do the following computation in 64bit integers to be
      // able to handle 4GB+ files:
      const MPI_Offset my_global_file_position =
        size_header +
        static_cast<MPI_Offset>(global_first_cell) * bytes_per_cell;

      ierr = Utilities::MPI::LargeCount::File_read_at_c(fh,
                                                        my_global_file_position,
                                                        dest_data_fixed.data(),
                                                        dest_data_fixed.size(),
                                                        MPI_BYTE,
                                                        MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);


      ierr = MPI_File_close(&fh);
      AssertThrowMPI(ierr);
    }

    //
    // ---------- Variable size data ----------
    //
    if (variable_size_data_stored)
      {
        const std::string fname_variable =
          std::string(filename) + "_variable.data";

        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        MPI_File fh;
        ierr = MPI_File_open(
          mpi_communicator, fname_variable.c_str(), MPI_MODE_RDONLY, info, &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);

        // Read sizes of all locally owned cells.
        dest_sizes_variable.resize(local_num_cells);

        const MPI_Offset my_global_file_position_sizes =
          static_cast<MPI_Offset>(global_first_cell) * sizeof(unsigned int);

        ierr = Utilities::MPI::LargeCount::File_read_at_c(
          fh,
          my_global_file_position_sizes,
          dest_sizes_variable.data(),
          dest_sizes_variable.size(),
          MPI_INT,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);


        // Compute my data size in bytes and compute prefix sum. We do this
        // in 64 bit to avoid overflow for files larger than 4 GB:
        const std::uint64_t size_on_proc =
          std::accumulate(dest_sizes_variable.begin(),
                          dest_sizes_variable.end(),
                          0ULL);

        std::uint64_t prefix_sum = 0;
        ierr                     = MPI_Exscan(&size_on_proc,
                          &prefix_sum,
                          1,
                          MPI_UINT64_T,
                          MPI_SUM,
                          mpi_communicator);
        AssertThrowMPI(ierr);

        const MPI_Offset my_global_file_position =
          static_cast<MPI_Offset>(global_num_cells) * sizeof(unsigned int) +
          prefix_sum;

        dest_data_variable.resize(size_on_proc);

        ierr =
          Utilities::MPI::LargeCount::File_read_at_c(fh,
                                                     my_global_file_position,
                                                     dest_data_variable.data(),
                                                     dest_data_variable.size(),
                                                     MPI_BYTE,
                                                     MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);
      }
#else
    (void)global_first_cell;
    (void)global_num_cells;
    (void)local_num_cells;
    (void)filename;
    (void)n_attached_deserialize_fixed;
    (void)n_attached_deserialize_variable;

    AssertThrow(false, ExcNeedsMPI());
#endif
  }



  template <int dim, int spacedim>
  void
  DistributedTriangulationBase<dim, spacedim>::DataTransfer::clear()
  {
    variable_size_data_stored = false;

    // free information about data sizes
    sizes_fixed_cumulative.clear();
    sizes_fixed_cumulative.shrink_to_fit();

    // free fixed size transfer data
    src_data_fixed.clear();
    src_data_fixed.shrink_to_fit();

    dest_data_fixed.clear();
    dest_data_fixed.shrink_to_fit();

    // free variable size transfer data
    src_sizes_variable.clear();
    src_sizes_variable.shrink_to_fit();

    src_data_variable.clear();
    src_data_variable.shrink_to_fit();

    dest_sizes_variable.clear();
    dest_sizes_variable.shrink_to_fit();

    dest_data_variable.clear();
    dest_data_variable.shrink_to_fit();
  }

} // end namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "tria_base.inst"

DEAL_II_NAMESPACE_CLOSE
