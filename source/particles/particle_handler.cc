// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/signaling_nan.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/particle_handler.h>

#include <limits>
#include <memory>
#include <utility>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  namespace
  {
    template <int dim, int spacedim>
    std::vector<char>
    pack_particles(std::vector<ParticleIterator<dim, spacedim>> &particles)
    {
      std::vector<char> buffer;

      if (particles.empty())
        return buffer;

      buffer.resize(particles.size() *
                    particles.front()->serialized_size_in_bytes());
      void *current_data = buffer.data();

      for (const auto &particle : particles)
        {
          current_data = particle->write_particle_data_to_memory(current_data);
        }

      return buffer;
    }
  } // namespace



  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::ParticleHandler()
    : triangulation()
    , mapping()
    , property_pool(std::make_unique<PropertyPool<dim, spacedim>>(0))
    , global_number_of_particles(0)
    , number_of_locally_owned_particles(0)
    , global_max_particles_per_cell(0)
    , next_free_particle_index(0)
    , size_callback()
    , store_callback()
    , load_callback()
    , tria_attached_data_index(numbers::invalid_unsigned_int)
    , tria_listeners()
  {
    reset_particle_container(particles);
  }



  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::ParticleHandler(
    const Triangulation<dim, spacedim> &triangulation,
    const Mapping<dim, spacedim>       &mapping,
    const unsigned int                  n_properties)
    : triangulation(&triangulation, typeid(*this).name())
    , mapping(&mapping, typeid(*this).name())
    , property_pool(std::make_unique<PropertyPool<dim, spacedim>>(n_properties))
    , cells_to_particle_cache(triangulation.n_active_cells(), particles.end())
    , global_number_of_particles(0)
    , number_of_locally_owned_particles(0)
    , global_max_particles_per_cell(0)
    , next_free_particle_index(0)
    , size_callback()
    , store_callback()
    , load_callback()
    , tria_attached_data_index(numbers::invalid_unsigned_int)
    , triangulation_cache(
        std::make_unique<GridTools::Cache<dim, spacedim>>(triangulation,
                                                          mapping))
    , tria_listeners()
  {
    reset_particle_container(particles);
    connect_to_triangulation_signals();
  }



  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::~ParticleHandler()
  {
    clear_particles();

    for (const auto &connection : tria_listeners)
      connection.disconnect();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::initialize(
    const Triangulation<dim, spacedim> &new_triangulation,
    const Mapping<dim, spacedim>       &new_mapping,
    const unsigned int                  n_properties)
  {
    clear();

    triangulation = &new_triangulation;
    mapping       = &new_mapping;

    reset_particle_container(particles);

    // Create the memory pool that will store all particle properties
    property_pool = std::make_unique<PropertyPool<dim, spacedim>>(n_properties);

    // Create the grid cache to cache the information about the triangulation
    // that is used to locate the particles into subdomains and cells
    triangulation_cache =
      std::make_unique<GridTools::Cache<dim, spacedim>>(new_triangulation,
                                                        new_mapping);

    cells_to_particle_cache.resize(triangulation->n_active_cells(),
                                   particles.end());

    connect_to_triangulation_signals();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::copy_from(
    const ParticleHandler<dim, spacedim> &particle_handler)
  {
    const unsigned int n_properties =
      particle_handler.property_pool->n_properties_per_slot();
    initialize(*particle_handler.triangulation,
               *particle_handler.mapping,
               n_properties);

    property_pool = std::make_unique<PropertyPool<dim, spacedim>>(
      *(particle_handler.property_pool));

    // copy static members
    global_number_of_particles = particle_handler.global_number_of_particles;
    number_of_locally_owned_particles =
      particle_handler.number_of_locally_owned_particles;

    global_max_particles_per_cell =
      particle_handler.global_max_particles_per_cell;
    next_free_particle_index = particle_handler.next_free_particle_index;

    // Manually copy over the particles because we do not want to touch the
    // anchor iterators set by initialize()
    particles.insert(particle_container_owned_end(),
                     particle_handler.particle_container_owned_begin(),
                     particle_handler.particle_container_owned_end());
    particles.insert(particle_container_ghost_end(),
                     particle_handler.particle_container_ghost_begin(),
                     particle_handler.particle_container_ghost_end());

    for (auto it = particles.begin(); it != particles.end(); ++it)
      if (!it->particles.empty())
        cells_to_particle_cache[it->cell->active_cell_index()] = it;

    ghost_particles_cache.ghost_particles_by_domain =
      particle_handler.ghost_particles_cache.ghost_particles_by_domain;
    tria_attached_data_index = particle_handler.tria_attached_data_index;
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::clear()
  {
    clear_particles();
    global_number_of_particles        = 0;
    number_of_locally_owned_particles = 0;
    next_free_particle_index          = 0;
    global_max_particles_per_cell     = 0;
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::clear_particles()
  {
    for (auto &particles_in_cell : particles)
      for (auto &particle : particles_in_cell.particles)
        if (particle != PropertyPool<dim, spacedim>::invalid_handle)
          property_pool->deregister_particle(particle);

    cells_to_particle_cache.clear();
    reset_particle_container(particles);
    if (triangulation != nullptr)
      cells_to_particle_cache.resize(triangulation->n_active_cells(),
                                     particles.end());

    // the particle properties have already been deleted by their destructor,
    // but the memory is still allocated. Return the memory as well.
    property_pool->clear();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::reserve(const std::size_t n_particles)
  {
    property_pool->reserve(n_particles);
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::reset_particle_container(
    particle_container &given_particles)
  {
    // Make sure to set a valid past-the-end iterator also in case we have no
    // triangulation
    const typename Triangulation<dim, spacedim>::cell_iterator
      past_the_end_iterator =
        triangulation != nullptr ?
          triangulation->end() :
          typename Triangulation<dim, spacedim>::cell_iterator(nullptr, -1, -1);

    given_particles.clear();
    for (unsigned int i = 0; i < 3; ++i)
      given_particles.emplace_back(
        std::vector<typename PropertyPool<dim, spacedim>::Handle>(),
        past_the_end_iterator);

    // Set the end of owned particles to the middle of the three elements
    const_cast<typename particle_container::iterator &>(owned_particles_end) =
      ++given_particles.begin();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::update_cached_numbers()
  {
    // first sort the owned particles by the active cell index
    bool sort_is_necessary = false;
    {
      auto previous = particle_container_owned_begin();
      for (auto next = previous; next != particle_container_owned_end(); ++next)
        {
          if (previous->cell.state() == IteratorState::valid &&
              next->cell.state() == IteratorState::valid &&
              previous->cell > next->cell)
            {
              sort_is_necessary = true;
              break;
            }
          previous = next;
        }
    }
    if (sort_is_necessary)
      {
        // we could have tried to call std::list::sort with a custom
        // comparator to get things sorted, but things are complicated by the
        // three anchor entries that we do not want to move, and we would
        // hence pay an O(N log(N)) algorithm with a large constant in front
        // of it (on the order of 20+ instructions with many
        // difficult-to-predict branches). Therefore, we simply copy the list
        // into a new one (keeping alive the possibly large vectors with
        // particles on cells) into a new container.
        particle_container sorted_particles;

        // note that this call updates owned_particles_end, so that
        // particle_container_owned_end() below already points to the
        // new container
        reset_particle_container(sorted_particles);

        // iterate over cells and insert the entries in the new order
        for (const auto &cell : triangulation->active_cell_iterators())
          if (!cell->is_artificial())
            if (cells_to_particle_cache[cell->active_cell_index()] !=
                particles.end())
              {
                // before we move the sorted_particles into particles
                // particle_container_ghost_end() still points to the
                // old particles container. Therefore this condition looks
                // quirky.
                typename particle_container::iterator insert_position =
                  cell->is_locally_owned() ? particle_container_owned_end() :
                                             --sorted_particles.end();
                typename particle_container::iterator new_entry =
                  sorted_particles.insert(
                    insert_position, typename particle_container::value_type());
                new_entry->cell = cell;
                new_entry->particles =
                  std::move(cells_to_particle_cache[cell->active_cell_index()]
                              ->particles);
              }
        particles = std::move(sorted_particles);

        // refresh cells_to_particle_cache
        cells_to_particle_cache.clear();
        cells_to_particle_cache.resize(triangulation->n_active_cells(),
                                       particles.end());
        for (auto it = particles.begin(); it != particles.end(); ++it)
          if (!it->particles.empty())
            cells_to_particle_cache[it->cell->active_cell_index()] = it;
      }

    // Ensure that we did not accidentally modify the anchor entries with
    // special purpose.
    Assert(particles.front().cell.state() == IteratorState::past_the_end &&
             particles.front().particles.empty() &&
             particles.back().cell.state() == IteratorState::past_the_end &&
             particles.back().particles.empty() &&
             owned_particles_end->cell.state() == IteratorState::past_the_end &&
             owned_particles_end->particles.empty(),
           ExcInternalError());

    if constexpr (running_in_debug_mode())
      {
        // check that no cache element hits the three anchor states in the list
        // of particles
        for (const auto &it : cells_to_particle_cache)
          Assert(it != particles.begin() && it != owned_particles_end &&
                   it != --(particles.end()),
                 ExcInternalError());

        // check that we only have locally owned particles in the first region
        // of cells; note that we skip the very first anchor element
        for (auto it = particle_container_owned_begin();
             it != particle_container_owned_end();
             ++it)
          Assert(it->cell->is_locally_owned(), ExcInternalError());

        // check that the cache is consistent with the iterators
        std::vector<typename particle_container::iterator> verify_cache(
          triangulation->n_active_cells(), particles.end());
        for (auto it = particles.begin(); it != particles.end(); ++it)
          if (!it->particles.empty())
            verify_cache[it->cell->active_cell_index()] = it;

        for (unsigned int i = 0; i < verify_cache.size(); ++i)
          Assert(verify_cache[i] == cells_to_particle_cache[i],
                 ExcInternalError());
      }

    // now compute local result with the function above and then compute the
    // collective results
    number_of_locally_owned_particles = 0;

    types::particle_index result[2] = {};
    for (const auto &particles_in_cell : particles)
      {
        const types::particle_index n_particles_in_cell =
          particles_in_cell.particles.size();

        // local_max_particles_per_cell
        result[0] = std::max(result[0], n_particles_in_cell);

        // number of locally owned particles
        if (n_particles_in_cell > 0 &&
            particles_in_cell.cell->is_locally_owned())
          number_of_locally_owned_particles += n_particles_in_cell;

        // local_max_particle_index
        for (const auto &particle : particles_in_cell.particles)
          result[1] = std::max(result[1], property_pool->get_id(particle));
      }

    global_number_of_particles =
      dealii::Utilities::MPI::sum(number_of_locally_owned_particles,
                                  triangulation->get_mpi_communicator());

    if (global_number_of_particles == 0)
      {
        next_free_particle_index      = 0;
        global_max_particles_per_cell = 0;
      }
    else
      {
        Utilities::MPI::max(result,
                            triangulation->get_mpi_communicator(),
                            result);

        next_free_particle_index      = result[1] + 1;
        global_max_particles_per_cell = result[0];
      }
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::n_particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
    const
  {
    if (cells_to_particle_cache.empty())
      return 0;

    if (cell->is_artificial() == false)
      {
        return cells_to_particle_cache[cell->active_cell_index()] !=
                   particles.end() ?
                 cells_to_particle_cache[cell->active_cell_index()]
                   ->particles.size() :
                 0;
      }
    else
      AssertThrow(false,
                  ExcMessage("You can't ask for the particles on an artificial "
                             "cell since we don't know what exists on these "
                             "kinds of cells."));

    return numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator_range
  ParticleHandler<dim, spacedim>::particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
    const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))
      ->particles_in_cell(cell);
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator_range
  ParticleHandler<dim, spacedim>::particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
  {
    const unsigned int active_cell_index = cell->active_cell_index();

    if (cell->is_artificial() == false)
      {
        if (cells_to_particle_cache[active_cell_index] == particles.end())
          {
            return boost::make_iterator_range(
              particle_iterator(particles.begin(), *property_pool, 0),
              particle_iterator(particles.begin(), *property_pool, 0));
          }
        else
          {
            const typename particle_container::iterator
              particles_in_current_cell =
                cells_to_particle_cache[active_cell_index];
            typename particle_container::iterator particles_in_next_cell =
              particles_in_current_cell;
            ++particles_in_next_cell;
            return boost::make_iterator_range(
              particle_iterator(particles_in_current_cell, *property_pool, 0),
              particle_iterator(particles_in_next_cell, *property_pool, 0));
          }
      }
    else
      AssertThrow(false,
                  ExcMessage("You can't ask for the particles on an artificial "
                             "cell since we don't know what exists on these "
                             "kinds of cells."));

    return {};
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::remove_particle(
    const ParticleHandler<dim, spacedim>::particle_iterator &particle)
  {
    auto &particles_on_cell = particle->particles_in_cell->particles;

    // if the particle has an invalid handle (e.g. because it has
    // been duplicated before calling this function) do not try
    // to deallocate its memory again
    auto handle = particle->get_handle();
    if (handle != PropertyPool<dim, spacedim>::invalid_handle)
      property_pool->deregister_particle(handle);

    // need to reduce the cached number before deleting, because the iterator
    // may be invalid after removing the particle even if only
    // accessing the cell
    const auto cell       = particle->get_surrounding_cell();
    const bool owned_cell = cell->is_locally_owned();
    if (owned_cell)
      --number_of_locally_owned_particles;

    if (particles_on_cell.size() > 1)
      {
        particles_on_cell[particle->particle_index_within_cell] =
          std::move(particles_on_cell.back());
        particles_on_cell.resize(particles_on_cell.size() - 1);
      }
    else
      {
        particles.erase(particle->particles_in_cell);
        cells_to_particle_cache[cell->active_cell_index()] = particles.end();
      }
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::remove_particles(
    const std::vector<ParticleHandler<dim, spacedim>::particle_iterator>
      &particles_to_remove)
  {
    // We need to remove particles backwards on each cell to keep the particle
    // iterators alive as we keep removing particles on the same cell. To
    // ensure that this is safe, we either check if we already have sorted
    // iterators or if we need to manually sort
    const auto check_greater = [](const particle_iterator &a,
                                  const particle_iterator &b) {
      return a->particles_in_cell->cell > b->particles_in_cell->cell ||
             (a->particles_in_cell->cell == b->particles_in_cell->cell &&
              a->particle_index_within_cell > b->particle_index_within_cell);
    };

    bool particles_are_sorted = true;
    auto previous             = particles_to_remove.begin();
    for (auto next = previous; next != particles_to_remove.end(); ++next)
      {
        if (check_greater(*previous, *next))
          {
            particles_are_sorted = false;
            break;
          }
        previous = next;
      }
    if (particles_are_sorted)
      {
        // pass along backwards in array
        for (auto it = particles_to_remove.rbegin();
             it != particles_to_remove.rend();
             ++it)
          remove_particle(*it);
      }
    else
      {
        std::vector<ParticleHandler<dim, spacedim>::particle_iterator>
          sorted_particles(particles_to_remove);
        std::sort(sorted_particles.begin(),
                  sorted_particles.end(),
                  check_greater);

        for (const auto &particle : sorted_particles)
          remove_particle(particle);
      }

    update_cached_numbers();
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::insert_particle(
    const Particle<dim, spacedim>                                     &particle,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
  {
    return insert_particle(particle.get_location(),
                           particle.get_reference_location(),
                           particle.get_id(),
                           cell,
                           particle.get_properties());
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::insert_particle(
    const typename PropertyPool<dim, spacedim>::Handle          handle,
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    const unsigned int active_cell_index = cell->active_cell_index();
    typename particle_container::iterator &cache =
      cells_to_particle_cache[active_cell_index];
    if (cache == particles.end())
      {
        const typename particle_container::iterator insert_position =
          cell->is_locally_owned() ? particle_container_owned_end() :
                                     particle_container_ghost_end();
        cache = particles.emplace(
          insert_position,
          std::vector<typename PropertyPool<dim, spacedim>::Handle>{handle},
          cell);
      }
    else
      {
        cache->particles.push_back(handle);
        Assert(cache->cell == cell, ExcInternalError());
      }
    return particle_iterator(cache,
                             *property_pool,
                             cache->particles.size() - 1);
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::insert_particle(
    const void                                                       *&data,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
  {
    Assert(triangulation != nullptr, ExcInternalError());
    Assert(cells_to_particle_cache.size() == triangulation->n_active_cells(),
           ExcInternalError());
    Assert(cell->is_locally_owned(),
           ExcMessage("You tried to insert particles into a cell that is not "
                      "locally owned. This is not supported."));

    particle_iterator particle_it =
      insert_particle(property_pool->register_particle(), cell);

    data = particle_it->read_particle_data_from_memory(data);

    ++number_of_locally_owned_particles;

    return particle_it;
  }



  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::insert_particle(
    const Point<spacedim>      &position,
    const Point<dim>           &reference_position,
    const types::particle_index particle_index,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const ArrayView<const double> &properties)
  {
    Assert(triangulation != nullptr, ExcInternalError());
    Assert(cells_to_particle_cache.size() == triangulation->n_active_cells(),
           ExcInternalError());
    Assert(cell.state() == IteratorState::valid, ExcInternalError());
    Assert(cell->is_locally_owned(),
           ExcMessage("You tried to insert particles into a cell that is not "
                      "locally owned. This is not supported."));

    particle_iterator particle_it =
      insert_particle(property_pool->register_particle(), cell);

    particle_it->set_location(position);
    particle_it->set_reference_location(reference_position);
    particle_it->set_id(particle_index);

    if (properties.size() != 0)
      particle_it->set_properties(properties);

    ++number_of_locally_owned_particles;

    return particle_it;
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::insert_particles(
    const std::multimap<
      typename Triangulation<dim, spacedim>::active_cell_iterator,
      Particle<dim, spacedim>> &new_particles)
  {
    reserve(n_locally_owned_particles() + new_particles.size());
    for (const auto &cell_and_particle : new_particles)
      insert_particle(cell_and_particle.second, cell_and_particle.first);

    update_cached_numbers();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::insert_particles(
    const std::vector<Point<spacedim>> &positions)
  {
    Assert(triangulation != nullptr, ExcInternalError());

    update_cached_numbers();
    reserve(n_locally_owned_particles() + positions.size());

    // Determine the starting particle index of this process, which
    // is the highest currently existing particle index plus the sum
    // of the number of newly generated particles of all
    // processes with a lower rank if in a parallel computation.
    const types::particle_index local_next_particle_index =
      get_next_free_particle_index();
    types::particle_index local_start_index = 0;

#ifdef DEAL_II_WITH_MPI
    if (const auto parallel_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &*triangulation))
      {
        types::particle_index particles_to_add_locally = positions.size();
        const int             ierr =
          MPI_Scan(&particles_to_add_locally,
                   &local_start_index,
                   1,
                   Utilities::MPI::mpi_type_id_for_type<types::particle_index>,
                   MPI_SUM,
                   parallel_triangulation->get_mpi_communicator());
        AssertThrowMPI(ierr);
        local_start_index -= particles_to_add_locally;
      }
#endif

    local_start_index += local_next_particle_index;

    auto point_locations =
      GridTools::compute_point_locations_try_all(*triangulation_cache,
                                                 positions);

    auto &cells           = std::get<0>(point_locations);
    auto &local_positions = std::get<1>(point_locations);
    auto &index_map       = std::get<2>(point_locations);
    auto &missing_points  = std::get<3>(point_locations);
    // If a point was not found, throwing an error, as the old
    // implementation of compute_point_locations would have done
    AssertThrow(missing_points.empty(),
                VectorTools::ExcPointNotAvailableHere());

    (void)missing_points;

    for (unsigned int i = 0; i < cells.size(); ++i)
      for (unsigned int p = 0; p < local_positions[i].size(); ++p)
        insert_particle(positions[index_map[i][p]],
                        local_positions[i][p],
                        local_start_index + index_map[i][p],
                        cells[i]);

    update_cached_numbers();
  }



  template <int dim, int spacedim>
  std::map<unsigned int, IndexSet>
  ParticleHandler<dim, spacedim>::insert_global_particles(
    const std::vector<Point<spacedim>> &positions,
    const std::vector<std::vector<BoundingBox<spacedim>>>
                                             &global_bounding_boxes,
    const std::vector<std::vector<double>>   &properties,
    const std::vector<types::particle_index> &ids)
  {
    if (!properties.empty())
      {
        AssertDimension(properties.size(), positions.size());
        if constexpr (running_in_debug_mode())
          {
            for (const auto &p : properties)
              AssertDimension(p.size(), n_properties_per_particle());
          }
      }

    if (!ids.empty())
      AssertDimension(ids.size(), positions.size());

    const auto comm = triangulation->get_mpi_communicator();

    const auto n_mpi_processes = Utilities::MPI::n_mpi_processes(comm);

    // Compute the global number of properties
    const auto n_global_properties =
      Utilities::MPI::sum(properties.size(), comm);

    // Gather the number of points per processor
    const auto n_particles_per_proc =
      Utilities::MPI::all_gather(comm, positions.size());

    // Calculate all starting points locally
    std::vector<unsigned int> particle_start_indices(n_mpi_processes);

    unsigned int particle_start_index = get_next_free_particle_index();
    for (unsigned int process = 0; process < particle_start_indices.size();
         ++process)
      {
        particle_start_indices[process] = particle_start_index;
        particle_start_index += n_particles_per_proc[process];
      }

    // Get all local information
    const auto cells_positions_and_index_maps =
      GridTools::distributed_compute_point_locations(*triangulation_cache,
                                                     positions,
                                                     global_bounding_boxes);

    // Unpack the information into several vectors:
    // All cells that contain at least one particle
    const auto &local_cells_containing_particles =
      std::get<0>(cells_positions_and_index_maps);

    // The reference position of every particle in the local part of the
    // triangulation.
    const auto &local_reference_positions =
      std::get<1>(cells_positions_and_index_maps);
    // The original index in the positions vector for each particle in the
    // local part of the triangulation
    const auto &original_indices_of_local_particles =
      std::get<2>(cells_positions_and_index_maps);
    // The real spatial position of every particle in the local part of the
    // triangulation.
    const auto &local_positions = std::get<3>(cells_positions_and_index_maps);
    // The MPI process that inserted each particle
    const auto &calling_process_indices =
      std::get<4>(cells_positions_and_index_maps);

    // Create the map of cpu to indices, indicating who sent us what particle
    std::map<unsigned int, std::vector<unsigned int>>
      original_process_to_local_particle_indices_tmp;
    for (unsigned int i_cell = 0;
         i_cell < local_cells_containing_particles.size();
         ++i_cell)
      {
        for (unsigned int i_particle = 0;
             i_particle < local_positions[i_cell].size();
             ++i_particle)
          {
            const unsigned int local_id_on_calling_process =
              original_indices_of_local_particles[i_cell][i_particle];
            const unsigned int calling_process =
              calling_process_indices[i_cell][i_particle];

            original_process_to_local_particle_indices_tmp[calling_process]
              .push_back(local_id_on_calling_process);
          }
      }
    std::map<unsigned int, IndexSet> original_process_to_local_particle_indices;
    for (auto &process_and_particle_indices :
         original_process_to_local_particle_indices_tmp)
      {
        const unsigned int calling_process = process_and_particle_indices.first;
        original_process_to_local_particle_indices.insert(
          {calling_process, IndexSet(n_particles_per_proc[calling_process])});
        std::sort(process_and_particle_indices.second.begin(),
                  process_and_particle_indices.second.end());
        original_process_to_local_particle_indices[calling_process].add_indices(
          process_and_particle_indices.second.begin(),
          process_and_particle_indices.second.end());
        original_process_to_local_particle_indices[calling_process].compress();
      }

    // A map from mpi process to properties, ordered as in the IndexSet.
    // Notice that this ordering may be different from the ordering in the
    // vectors above, since no local ordering is guaranteed by the
    // distribute_compute_point_locations() call.
    // This is only filled if n_global_properties is > 0
    std::map<unsigned int, std::vector<std::vector<double>>>
      locally_owned_properties_from_other_processes;

    // A map from mpi process to ids, ordered as in the IndexSet.
    // Notice that this ordering may be different from the ordering in the
    // vectors above, since no local ordering is guaranteed by the
    // distribute_compute_point_locations() call.
    // This is only filled if ids.size() is > 0
    std::map<unsigned int, std::vector<types::particle_index>>
      locally_owned_ids_from_other_processes;

    if (n_global_properties > 0 || !ids.empty())
      {
        // Gather whom I sent my own particles to, to decide whom to send
        // the particle properties or the ids
        auto send_to_cpu = Utilities::MPI::some_to_some(
          comm, original_process_to_local_particle_indices);

        // Prepare the vector of properties to send
        if (n_global_properties > 0)
          {
            std::map<unsigned int, std::vector<std::vector<double>>>
              non_locally_owned_properties;

            for (const auto &it : send_to_cpu)
              {
                std::vector<std::vector<double>> properties_to_send(
                  it.second.n_elements(),
                  std::vector<double>(n_properties_per_particle()));
                unsigned int index = 0;
                for (const auto el : it.second)
                  properties_to_send[index++] = properties[el];
                non_locally_owned_properties.insert(
                  {it.first, properties_to_send});
              }

            // Send the non locally owned properties to each mpi process
            // that needs them
            locally_owned_properties_from_other_processes =
              Utilities::MPI::some_to_some(comm, non_locally_owned_properties);

            AssertDimension(
              locally_owned_properties_from_other_processes.size(),
              original_process_to_local_particle_indices.size());
          }

        if (!ids.empty())
          {
            std::map<unsigned int, std::vector<types::particle_index>>
              non_locally_owned_ids;
            for (const auto &it : send_to_cpu)
              {
                std::vector<types::particle_index> ids_to_send(
                  it.second.n_elements());
                unsigned int index = 0;
                for (const auto el : it.second)
                  ids_to_send[index++] = ids[el];
                non_locally_owned_ids.insert({it.first, ids_to_send});
              }

            // Send the non locally owned ids to each mpi process
            // that needs them
            locally_owned_ids_from_other_processes =
              Utilities::MPI::some_to_some(comm, non_locally_owned_ids);

            AssertDimension(locally_owned_ids_from_other_processes.size(),
                            original_process_to_local_particle_indices.size());
          }
      }

    // Now fill up the actual particles
    for (unsigned int i_cell = 0;
         i_cell < local_cells_containing_particles.size();
         ++i_cell)
      {
        for (unsigned int i_particle = 0;
             i_particle < local_positions[i_cell].size();
             ++i_particle)
          {
            const unsigned int local_id_on_calling_process =
              original_indices_of_local_particles[i_cell][i_particle];

            const unsigned int calling_process =
              calling_process_indices[i_cell][i_particle];

            const unsigned int index_within_set =
              original_process_to_local_particle_indices[calling_process]
                .index_within_set(local_id_on_calling_process);

            const unsigned int particle_id =
              ids.empty() ?
                local_id_on_calling_process +
                  particle_start_indices[calling_process] :
                locally_owned_ids_from_other_processes[calling_process]
                                                      [index_within_set];

            auto particle_it =
              insert_particle(local_positions[i_cell][i_particle],
                              local_reference_positions[i_cell][i_particle],
                              particle_id,
                              local_cells_containing_particles[i_cell]);

            if (n_global_properties > 0)
              {
                particle_it->set_properties(
                  locally_owned_properties_from_other_processes
                    [calling_process][index_within_set]);
              }
          }
      }

    update_cached_numbers();

    return original_process_to_local_particle_indices;
  }



  template <int dim, int spacedim>
  std::map<unsigned int, IndexSet>
  ParticleHandler<dim, spacedim>::insert_global_particles(
    const std::vector<Particle<dim, spacedim>> &particles,
    const std::vector<std::vector<BoundingBox<spacedim>>>
      &global_bounding_boxes)
  {
    // Store the positions in a vector of points, the ids in a vector of ids,
    // and the properties, if any, in a vector of vector of properties.
    std::vector<Point<spacedim>>       positions;
    std::vector<std::vector<double>>   properties;
    std::vector<types::particle_index> ids;
    positions.resize(particles.size());
    ids.resize(particles.size());
    if (n_properties_per_particle() > 0)
      properties.resize(particles.size(),
                        std::vector<double>(n_properties_per_particle()));

    unsigned int i = 0;
    for (const auto &p : particles)
      {
        positions[i] = p.get_location();
        ids[i]       = p.get_id();
        if (p.has_properties())
          properties[i] = {p.get_properties().begin(),
                           p.get_properties().end()};
        ++i;
      }

    return insert_global_particles(positions,
                                   global_bounding_boxes,
                                   properties,
                                   ids);
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::n_global_particles() const
  {
    return global_number_of_particles;
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::n_global_max_particles_per_cell() const
  {
    return global_max_particles_per_cell;
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::n_locally_owned_particles() const
  {
    return number_of_locally_owned_particles;
  }



  template <int dim, int spacedim>
  unsigned int
  ParticleHandler<dim, spacedim>::n_properties_per_particle() const
  {
    return property_pool->n_properties_per_slot();
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::get_next_free_particle_index() const
  {
    return next_free_particle_index;
  }



  template <int dim, int spacedim>
  IndexSet
  ParticleHandler<dim, spacedim>::locally_owned_particle_ids() const
  {
    IndexSet                           set(get_next_free_particle_index());
    std::vector<types::particle_index> indices;
    indices.reserve(n_locally_owned_particles());
    for (const auto &p : *this)
      indices.push_back(p.get_id());
    set.add_indices(indices.begin(), indices.end());
    set.compress();
    return set;
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::get_max_local_particle_index() const
  {
    return property_pool->n_slots();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::get_particle_positions(
    std::vector<Point<spacedim>> &positions,
    const bool                    add_to_output_vector)
  {
    // There should be one point per particle to gather
    AssertDimension(positions.size(), n_locally_owned_particles());

    unsigned int i = 0;
    for (auto it = begin(); it != end(); ++it, ++i)
      {
        if (add_to_output_vector)
          positions[i] = positions[i] + it->get_location();
        else
          positions[i] = it->get_location();
      }
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::set_particle_positions(
    const std::vector<Point<spacedim>> &new_positions,
    const bool                          displace_particles)
  {
    // There should be one point per particle to fix the new position
    AssertDimension(new_positions.size(), n_locally_owned_particles());

    unsigned int i = 0;
    for (auto it = begin(); it != end(); ++it, ++i)
      {
        Point<spacedim> &location = it->get_location();
        if (displace_particles)
          location += new_positions[i];
        else
          location = new_positions[i];
      }
    sort_particles_into_subdomains_and_cells();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::set_particle_positions(
    const Function<spacedim> &function,
    const bool                displace_particles)
  {
    // The function should have sufficient components to displace the
    // particles
    AssertDimension(function.n_components, spacedim);

    Vector<double> new_position(spacedim);
    for (auto &particle : *this)
      {
        Point<spacedim> &particle_location = particle.get_location();
        function.vector_value(particle_location, new_position);
        if (displace_particles)
          for (unsigned int d = 0; d < spacedim; ++d)
            particle_location[d] += new_position[d];
        else
          for (unsigned int d = 0; d < spacedim; ++d)
            particle_location[d] = new_position[d];
      }
    sort_particles_into_subdomains_and_cells();
  }



  template <int dim, int spacedim>
  PropertyPool<dim, spacedim> &
  ParticleHandler<dim, spacedim>::get_property_pool() const
  {
    return *property_pool;
  }



  namespace
  {
    /**
     * This function is used as comparison argument to std::sort to sort the
     * vector of tensors @p center_directions by its scalar product with the
     * @p particle_direction tensor. The sorted indices allow to
     * loop over @p center_directions with increasing angle between
     * @p particle_direction and @p center_directions. This function assumes
     * that @p particle_direction and @p center_directions are normalized
     * to length one before calling this function.
     */
    template <int dim>
    bool
    compare_particle_association(
      const unsigned int                 a,
      const unsigned int                 b,
      const Tensor<1, dim>              &particle_direction,
      const std::vector<Tensor<1, dim>> &center_directions)
    {
      const double scalar_product_a = center_directions[a] * particle_direction;
      const double scalar_product_b = center_directions[b] * particle_direction;

      // The function is supposed to return if a is before b. We are looking
      // for the alignment of particle direction and center direction,
      // therefore return if the scalar product of a is larger.
      return (scalar_product_a > scalar_product_b);
    }
  } // namespace



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::sort_particles_into_subdomains_and_cells()
  {
    Assert(triangulation != nullptr, ExcInternalError());
    Assert(cells_to_particle_cache.size() == triangulation->n_active_cells(),
           ExcInternalError());

    // TODO: The current algorithm only works for particles that are in
    // the local domain or in ghost cells, because it only knows the
    // subdomain_id of ghost cells, but not of artificial cells. This
    // can be extended using the distributed version of compute point
    // locations.
    // TODO: Extend this function to allow keeping particles on other
    // processes around (with an invalid cell).

    std::vector<particle_iterator> particles_out_of_cell;

    // Reserve some space for particles that need sorting to avoid frequent
    // re-allocation. Guess 25% of particles need sorting. Balance memory
    // overhead and performance.
    particles_out_of_cell.reserve(n_locally_owned_particles() / 4);

    // Now update the reference locations of the moved particles
    std::vector<Point<spacedim>> real_locations;
    std::vector<Point<dim>>      reference_locations;
    real_locations.reserve(global_max_particles_per_cell);
    reference_locations.reserve(global_max_particles_per_cell);

    for (const auto &cell : triangulation->active_cell_iterators())
      {
        // Particles can be inserted into arbitrary cells, e.g. if their cell is
        // not known. However, for artificial cells we can not evaluate
        // the reference position of particles. Do not sort particles that are
        // not locally owned, because they will be sorted by the process that
        // owns them.
        if (cell->is_locally_owned() == false)
          {
            continue;
          }

        const unsigned int n_pic = n_particles_in_cell(cell);
        auto               pic   = particles_in_cell(cell);

        real_locations.clear();
        for (const auto &particle : pic)
          real_locations.push_back(particle.get_location());

        reference_locations.resize(n_pic);
        mapping->transform_points_real_to_unit_cell(cell,
                                                    real_locations,
                                                    reference_locations);

        auto particle = pic.begin();
        for (const auto &p_unit : reference_locations)
          {
            if (numbers::is_finite(p_unit[0]) &&
                cell->reference_cell().contains_point(p_unit,
                                                      tolerance_inside_cell))
              particle->set_reference_location(p_unit);
            else
              particles_out_of_cell.push_back(particle);

            ++particle;
          }
      }

    // There are three reasons why a particle is not in its old cell:
    // It moved to another cell, to another subdomain or it left the mesh.
    // Particles that moved to another cell are updated and moved inside the
    // particles vector, particles that moved to another domain are
    // collected in the moved_particles_domain vector. Particles that left
    // the mesh completely are ignored and removed.
    std::map<types::subdomain_id, std::vector<particle_iterator>>
      moved_particles;
    std::map<
      types::subdomain_id,
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      moved_cells;

    // We do not know exactly how many particles are lost, exchanged between
    // domains, or remain on this process. Therefore we pre-allocate
    // approximate sizes for these vectors. If more space is needed an
    // automatic and relatively fast (compared to other parts of this
    // algorithm) re-allocation will happen.
    std::set<types::subdomain_id> ghost_owners;
    if (const auto parallel_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &*triangulation))
      ghost_owners = parallel_triangulation->ghost_owners();

    // Reserve some space for particles that need communication to avoid
    // frequent re-allocation. Guess 25% of particles out of their old cell need
    // communication. Balance memory overhead and performance.
    for (const auto &ghost_owner : ghost_owners)
      moved_particles[ghost_owner].reserve(particles_out_of_cell.size() / 4);
    for (const auto &ghost_owner : ghost_owners)
      moved_cells[ghost_owner].reserve(particles_out_of_cell.size() / 4);

    {
      // Create a map from vertices to adjacent cells using grid cache
      const std::vector<
        std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
        &vertex_to_cells = triangulation_cache->get_vertex_to_cell_map();

      // Create a corresponding map of vectors from vertex to cell center
      // using grid cache
      const std::vector<std::vector<Tensor<1, spacedim>>>
        &vertex_to_cell_centers =
          triangulation_cache->get_vertex_to_cell_centers_directions();

      std::vector<unsigned int> search_order;

      // Reuse these vectors below, but only with a single element.
      // Avoid resizing for every particle.
      reference_locations.resize(1, numbers::signaling_nan<Point<dim>>());
      real_locations.resize(1, numbers::signaling_nan<Point<spacedim>>());

      // Find the cells that the particles moved to.
      for (auto &out_particle : particles_out_of_cell)
        {
          // make a copy of the current cell, since we will modify the
          // variable current_cell below, but we need the original in
          // the case the particle is not found
          auto current_cell = out_particle->get_surrounding_cell();

          real_locations[0] = out_particle->get_location();

          // Record if the new cell was found
          bool found_cell = false;

          // Check if the particle is in one of the old cell's neighbors
          // that are adjacent to the closest vertex
          const unsigned int closest_vertex =
            GridTools::find_closest_vertex_of_cell<dim, spacedim>(
              current_cell, out_particle->get_location(), *mapping);
          const unsigned int closest_vertex_index =
            current_cell->vertex_index(closest_vertex);

          const auto &candidate_cells = vertex_to_cells[closest_vertex_index];
          const unsigned int n_candidate_cells = candidate_cells.size();

          // The order of searching through the candidate cells matters for
          // performance reasons. Start with a simple order.
          search_order.resize(n_candidate_cells);
          for (unsigned int i = 0; i < n_candidate_cells; ++i)
            search_order[i] = i;

          // If the particle is not on a vertex, we can do better by
          // sorting the candidate cells by alignment with
          // the vertex_to_particle direction.
          Tensor<1, spacedim> vertex_to_particle =
            out_particle->get_location() - current_cell->vertex(closest_vertex);

          // Only do this if the particle is not on a vertex, otherwise we
          // cannot normalize
          if (vertex_to_particle.norm_square() >
              1e4 * std::numeric_limits<double>::epsilon() *
                std::numeric_limits<double>::epsilon() *
                vertex_to_cell_centers[closest_vertex_index][0].norm_square())
            {
              vertex_to_particle /= vertex_to_particle.norm();
              const auto &vertex_to_cells =
                vertex_to_cell_centers[closest_vertex_index];

              std::sort(search_order.begin(),
                        search_order.end(),
                        [&vertex_to_particle,
                         &vertex_to_cells](const unsigned int a,
                                           const unsigned int b) {
                          return compare_particle_association(
                            a, b, vertex_to_particle, vertex_to_cells);
                        });
            }

          // Search all of the candidate cells according to the determined
          // order. Most likely we will find the particle in them.
          for (unsigned int i = 0; i < n_candidate_cells; ++i)
            {
              typename std::set<
                typename Triangulation<dim, spacedim>::active_cell_iterator>::
                const_iterator candidate_cell = candidate_cells.begin();

              std::advance(candidate_cell, search_order[i]);
              mapping->transform_points_real_to_unit_cell(*candidate_cell,
                                                          real_locations,
                                                          reference_locations);

              if ((*candidate_cell)
                    ->reference_cell()
                    .contains_point(reference_locations[0],
                                    tolerance_inside_cell))
                {
                  current_cell = *candidate_cell;
                  found_cell   = true;
                  break;
                }
            }

          // If we did not find a cell the particle is not in a neighbor of
          // its old cell. Look for the new cell in the whole local domain.
          // This case should be rare.
          if (!found_cell)
            {
              // For some clang-based compilers and boost versions the call to
              // RTree::query doesn't compile. We use a slower implementation as
              // workaround.
              // This is fixed in boost in
              // https://github.com/boostorg/numeric_conversion/commit/50a1eae942effb0a9b90724323ef8f2a67e7984a
#if defined(DEAL_II_WITH_BOOST_BUNDLED) ||                \
  !(defined(__clang_major__) && __clang_major__ >= 16) || \
  BOOST_VERSION >= 108100

              std::vector<std::pair<Point<spacedim>, unsigned int>>
                closest_vertex_in_domain;
              triangulation_cache->get_used_vertices_rtree().query(
                boost::geometry::index::nearest(out_particle->get_location(),
                                                1),
                std::back_inserter(closest_vertex_in_domain));

              // We should have one and only one result
              AssertDimension(closest_vertex_in_domain.size(), 1);
              const unsigned int closest_vertex_index_in_domain =
                closest_vertex_in_domain[0].second;
#else
              const unsigned int closest_vertex_index_in_domain =
                GridTools::find_closest_vertex(*mapping,
                                               *triangulation,
                                               out_particle->get_location());
#endif

              // Search all of the cells adjacent to the closest vertex of the
              // domain. Most likely we will find the particle in them.
              for (const auto &cell :
                   vertex_to_cells[closest_vertex_index_in_domain])
                {
                  mapping->transform_points_real_to_unit_cell(
                    cell, real_locations, reference_locations);

                  if (cell->reference_cell().contains_point(
                        reference_locations[0], tolerance_inside_cell))
                    {
                      current_cell = cell;
                      found_cell   = true;
                      break;
                    }
                }
            }

          if (!found_cell)
            {
              // We can find no cell for this particle. It has left the
              // domain due to an integration error or an open boundary.
              // Signal the loss and move on.
              signals.particle_lost(out_particle,
                                    out_particle->get_surrounding_cell());
              continue;
            }

          // If we are here, we found a cell and reference position for this
          // particle
          out_particle->set_reference_location(reference_locations[0]);

          // Reinsert the particle into our domain if we own its cell.
          // Mark it for MPI transfer otherwise
          if (current_cell->is_locally_owned())
            {
              typename PropertyPool<dim, spacedim>::Handle &old =
                out_particle->particles_in_cell
                  ->particles[out_particle->particle_index_within_cell];

              // Avoid deallocating the memory of this particle
              const auto old_value = old;
              old = PropertyPool<dim, spacedim>::invalid_handle;

              // Allocate particle with the old handle
              insert_particle(old_value, current_cell);
            }
          else
            {
              moved_particles[current_cell->subdomain_id()].push_back(
                out_particle);
              moved_cells[current_cell->subdomain_id()].push_back(current_cell);
            }
        }
    }

    // Exchange particles between processors if we have more than one process
#ifdef DEAL_II_WITH_MPI
    if (const auto parallel_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &*triangulation))
      {
        if (dealii::Utilities::MPI::n_mpi_processes(
              parallel_triangulation->get_mpi_communicator()) > 1)
          send_recv_particles(moved_particles, moved_cells);
      }
#endif

    // remove_particles also calls update_cached_numbers()
    remove_particles(particles_out_of_cell);

    // now make sure particle data is sorted in order of iteration
    std::vector<typename PropertyPool<dim, spacedim>::Handle> unsorted_handles;
    unsorted_handles.reserve(property_pool->n_registered_slots());

    typename PropertyPool<dim, spacedim>::Handle sorted_handle = 0;
    for (auto &particles_in_cell : particles)
      for (auto &particle : particles_in_cell.particles)
        {
          unsorted_handles.push_back(particle);
          particle = sorted_handle++;
        }

    property_pool->sort_memory_slots(unsorted_handles);

  } // namespace Particles



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::exchange_ghost_particles(
    const bool enable_cache)
  {
    // Nothing to do in serial computations
    const auto parallel_triangulation =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*triangulation);
    if (parallel_triangulation != nullptr)
      {
        if (dealii::Utilities::MPI::n_mpi_processes(
              parallel_triangulation->get_mpi_communicator()) == 1)
          return;
      }
    else
      return;

#ifndef DEAL_II_WITH_MPI
    (void)enable_cache;
#else
    // Clear ghost particles and their properties
    for (const auto &cell : triangulation->active_cell_iterators())
      if (cell->is_ghost() &&
          cells_to_particle_cache[cell->active_cell_index()] != particles.end())
        {
          Assert(cells_to_particle_cache[cell->active_cell_index()]->cell ==
                   cell,
                 ExcInternalError());
          // Clear particle properties
          for (auto &ghost_particle :
               cells_to_particle_cache[cell->active_cell_index()]->particles)
            property_pool->deregister_particle(ghost_particle);

          // Clear particles themselves
          particles.erase(cells_to_particle_cache[cell->active_cell_index()]);
          cells_to_particle_cache[cell->active_cell_index()] = particles.end();
        }

    // Clear ghost particles cache and invalidate it
    ghost_particles_cache.ghost_particles_by_domain.clear();
    ghost_particles_cache.valid = false;

    // In the case of a parallel simulation with periodic boundary conditions
    // the vertices associated with periodic boundaries are not directly
    // connected to the ghost cells but they are connected to the ghost cells
    // through their coinciding vertices. We gather this information using the
    // vertices_with_ghost_neighbors map
    const std::map<unsigned int, std::set<types::subdomain_id>>
      &vertices_with_ghost_neighbors =
        triangulation_cache->get_vertices_with_ghost_neighbors();

    const std::set<types::subdomain_id> ghost_owners =
      parallel_triangulation->ghost_owners();
    for (const auto ghost_owner : ghost_owners)
      ghost_particles_cache.ghost_particles_by_domain[ghost_owner].reserve(
        n_locally_owned_particles() / 4);

    const std::vector<std::set<unsigned int>> vertex_to_neighbor_subdomain =
      triangulation_cache->get_vertex_to_neighbor_subdomain();

    for (const auto &cell : triangulation->active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            std::set<unsigned int> cell_to_neighbor_subdomain;
            for (const unsigned int v : cell->vertex_indices())
              {
                const auto vertex_ghost_neighbors =
                  vertices_with_ghost_neighbors.find(cell->vertex_index(v));
                if (vertex_ghost_neighbors !=
                    vertices_with_ghost_neighbors.end())
                  {
                    cell_to_neighbor_subdomain.insert(
                      vertex_ghost_neighbors->second.begin(),
                      vertex_ghost_neighbors->second.end());
                  }
              }

            if (cell_to_neighbor_subdomain.size() > 0)
              {
                const particle_iterator_range particle_range =
                  particles_in_cell(cell);

                for (const auto domain : cell_to_neighbor_subdomain)
                  {
                    for (typename particle_iterator_range::iterator particle =
                           particle_range.begin();
                         particle != particle_range.end();
                         ++particle)
                      ghost_particles_cache.ghost_particles_by_domain[domain]
                        .push_back(particle);
                  }
              }
          }
      }

    send_recv_particles(
      ghost_particles_cache.ghost_particles_by_domain,
      std::map<
        types::subdomain_id,
        std::vector<
          typename Triangulation<dim, spacedim>::active_cell_iterator>>(),
      enable_cache);
#endif
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::update_ghost_particles()
  {
    // Nothing to do in serial computations
    const auto parallel_triangulation =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*triangulation);
    if (parallel_triangulation == nullptr ||
        dealii::Utilities::MPI::n_mpi_processes(
          parallel_triangulation->get_mpi_communicator()) == 1)
      {
        return;
      }


#ifdef DEAL_II_WITH_MPI
    // First clear the current ghost_particle information
    // ghost_particles.clear();
    Assert(ghost_particles_cache.valid,
           ExcMessage(
             "Ghost particles cannot be updated if they first have not been "
             "exchanged at least once with the cache enabled"));


    send_recv_particles_properties_and_location(
      ghost_particles_cache.ghost_particles_by_domain);
#endif
  }



#ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::send_recv_particles(
    const std::map<types::subdomain_id, std::vector<particle_iterator>>
      &particles_to_send,
    const std::map<
      types::subdomain_id,
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>>
              &send_cells,
    const bool build_cache)
  {
    Assert(triangulation != nullptr, ExcInternalError());
    Assert(cells_to_particle_cache.size() == triangulation->n_active_cells(),
           ExcInternalError());

    ghost_particles_cache.valid = build_cache;

    const auto parallel_triangulation =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*triangulation);
    Assert(parallel_triangulation,
           ExcMessage("This function is only implemented for "
                      "parallel::TriangulationBase objects."));

    // Determine the communication pattern
    const std::set<types::subdomain_id> ghost_owners =
      parallel_triangulation->ghost_owners();
    const std::vector<types::subdomain_id> neighbors(ghost_owners.begin(),
                                                     ghost_owners.end());
    const unsigned int                     n_neighbors = neighbors.size();

    if (send_cells.size() != 0)
      Assert(particles_to_send.size() == send_cells.size(), ExcInternalError());

    // If we do not know the subdomain this particle needs to be send to,
    // throw an error
    Assert(particles_to_send.find(numbers::artificial_subdomain_id) ==
             particles_to_send.end(),
           ExcInternalError());

    // TODO: Implement the shipping of particles to processes that are not
    // ghost owners of the local domain
    for (auto send_particles = particles_to_send.begin();
         send_particles != particles_to_send.end();
         ++send_particles)
      Assert(ghost_owners.find(send_particles->first) != ghost_owners.end(),
             ExcNotImplemented());

    std::size_t n_send_particles = 0;
    for (auto send_particles = particles_to_send.begin();
         send_particles != particles_to_send.end();
         ++send_particles)
      n_send_particles += send_particles->second.size();

    const unsigned int cellid_size = sizeof(CellId::binary_type);

    // Containers for the amount and offsets of data we will send
    // to other processors and the data itself.
    std::vector<unsigned int> n_send_data(n_neighbors, 0);
    std::vector<unsigned int> send_offsets(n_neighbors, 0);
    std::vector<char>         send_data;

    Particle<dim, spacedim> test_particle;
    test_particle.set_property_pool(*property_pool);

    const unsigned int individual_particle_data_size =
      test_particle.serialized_size_in_bytes() +
      (size_callback ? size_callback() : 0);

    const unsigned int individual_total_particle_data_size =
      individual_particle_data_size + cellid_size;

    // Only serialize things if there are particles to be send.
    // We can not return early even if no particles
    // are send, because we might receive particles from other processes
    if (n_send_particles > 0)
      {
        // Allocate space for sending particle data
        send_data.resize(n_send_particles *
                         individual_total_particle_data_size);

        void *data = static_cast<void *>(&send_data.front());

        // Serialize the data sorted by receiving process
        for (unsigned int i = 0; i < n_neighbors; ++i)
          {
            send_offsets[i] = reinterpret_cast<std::size_t>(data) -
                              reinterpret_cast<std::size_t>(&send_data.front());

            const unsigned int n_particles_to_send =
              particles_to_send.at(neighbors[i]).size();

            Assert(static_cast<std::size_t>(n_particles_to_send) *
                       individual_total_particle_data_size ==
                     static_cast<std::size_t>(
                       n_particles_to_send *
                       individual_total_particle_data_size),
                   ExcMessage("Overflow when trying to send particle "
                              "data"));

            for (unsigned int j = 0; j < n_particles_to_send; ++j)
              {
                // If no target cells are given, use the iterator
                // information
                typename Triangulation<dim, spacedim>::active_cell_iterator
                  cell;
                if (send_cells.empty())
                  cell = particles_to_send.at(neighbors[i])[j]
                           ->get_surrounding_cell();
                else
                  cell = send_cells.at(neighbors[i])[j];

                const CellId::binary_type cellid =
                  cell->id().template to_binary<dim>();
                memcpy(data, &cellid, cellid_size);
                data = static_cast<char *>(data) + cellid_size;

                data = particles_to_send.at(neighbors[i])[j]
                         ->write_particle_data_to_memory(data);
                if (store_callback)
                  data =
                    store_callback(particles_to_send.at(neighbors[i])[j], data);
              }
            n_send_data[i] = n_particles_to_send;
          }
      }

    // Containers for the data we will receive from other processors
    std::vector<unsigned int> n_recv_data(n_neighbors);
    std::vector<unsigned int> recv_offsets(n_neighbors);

    {
      const int mpi_tag = Utilities::MPI::internal::Tags::
        particle_handler_send_recv_particles_setup;

      std::vector<MPI_Request> n_requests(2 * n_neighbors);
      for (unsigned int i = 0; i < n_neighbors; ++i)
        {
          const int ierr =
            MPI_Irecv(&(n_recv_data[i]),
                      1,
                      MPI_UNSIGNED,
                      neighbors[i],
                      mpi_tag,
                      parallel_triangulation->get_mpi_communicator(),
                      &(n_requests[2 * i]));
          AssertThrowMPI(ierr);
        }
      for (unsigned int i = 0; i < n_neighbors; ++i)
        {
          const int ierr =
            MPI_Isend(&(n_send_data[i]),
                      1,
                      MPI_UNSIGNED,
                      neighbors[i],
                      mpi_tag,
                      parallel_triangulation->get_mpi_communicator(),
                      &(n_requests[2 * i + 1]));
          AssertThrowMPI(ierr);
        }
      const int ierr =
        MPI_Waitall(2 * n_neighbors, n_requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Determine how many particles and data we will receive
    unsigned int total_recv_data = 0;
    for (unsigned int neighbor_id = 0; neighbor_id < n_neighbors; ++neighbor_id)
      {
        recv_offsets[neighbor_id] = total_recv_data;
        total_recv_data +=
          n_recv_data[neighbor_id] * individual_total_particle_data_size;
      }

    // Set up the space for the received particle data
    std::vector<char> recv_data(total_recv_data);

    // Exchange the particle data between domains
    {
      std::vector<MPI_Request> requests(2 * n_neighbors);
      unsigned int             send_ops = 0;
      unsigned int             recv_ops = 0;

      const int mpi_tag = Utilities::MPI::internal::Tags::
        particle_handler_send_recv_particles_send;

      for (unsigned int i = 0; i < n_neighbors; ++i)
        if (n_recv_data[i] > 0)
          {
            const int ierr =
              MPI_Irecv(&(recv_data[recv_offsets[i]]),
                        n_recv_data[i] * individual_total_particle_data_size,
                        MPI_CHAR,
                        neighbors[i],
                        mpi_tag,
                        parallel_triangulation->get_mpi_communicator(),
                        &(requests[send_ops]));
            AssertThrowMPI(ierr);
            ++send_ops;
          }

      for (unsigned int i = 0; i < n_neighbors; ++i)
        if (n_send_data[i] > 0)
          {
            const int ierr =
              MPI_Isend(&(send_data[send_offsets[i]]),
                        n_send_data[i] * individual_total_particle_data_size,
                        MPI_CHAR,
                        neighbors[i],
                        mpi_tag,
                        parallel_triangulation->get_mpi_communicator(),
                        &(requests[send_ops + recv_ops]));
            AssertThrowMPI(ierr);
            ++recv_ops;
          }
      const int ierr =
        MPI_Waitall(send_ops + recv_ops, requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Put the received particles into the domain if they are in the
    // triangulation
    const void *recv_data_it = static_cast<const void *>(recv_data.data());

    // Store the particle iterators in the cache
    auto &ghost_particles_iterators =
      ghost_particles_cache.ghost_particles_iterators;

    if (build_cache)
      {
        ghost_particles_iterators.clear();

        auto &send_pointers_particles = ghost_particles_cache.send_pointers;
        send_pointers_particles.assign(n_neighbors + 1, 0);

        for (unsigned int i = 0; i < n_neighbors; ++i)
          send_pointers_particles[i + 1] =
            send_pointers_particles[i] +
            n_send_data[i] * individual_particle_data_size;

        auto &recv_pointers_particles = ghost_particles_cache.recv_pointers;
        recv_pointers_particles.assign(n_neighbors + 1, 0);

        for (unsigned int i = 0; i < n_neighbors; ++i)
          recv_pointers_particles[i + 1] =
            recv_pointers_particles[i] +
            n_recv_data[i] * individual_particle_data_size;

        ghost_particles_cache.neighbors = neighbors;

        ghost_particles_cache.send_data.resize(
          ghost_particles_cache.send_pointers.back());
        ghost_particles_cache.recv_data.resize(
          ghost_particles_cache.recv_pointers.back());
      }

    while (reinterpret_cast<std::size_t>(recv_data_it) -
             reinterpret_cast<std::size_t>(recv_data.data()) <
           total_recv_data)
      {
        CellId::binary_type binary_cellid;
        memcpy(&binary_cellid, recv_data_it, cellid_size);
        const CellId id(binary_cellid);
        recv_data_it = static_cast<const char *>(recv_data_it) + cellid_size;

        const typename Triangulation<dim, spacedim>::active_cell_iterator cell =
          triangulation->create_cell_iterator(id);

        insert_particle(property_pool->register_particle(), cell);
        const typename particle_container::iterator &cache =
          cells_to_particle_cache[cell->active_cell_index()];
        Assert(cache->cell == cell, ExcInternalError());

        particle_iterator particle_it(cache,
                                      *property_pool,
                                      cache->particles.size() - 1);

        recv_data_it =
          particle_it->read_particle_data_from_memory(recv_data_it);

        if (load_callback)
          recv_data_it = load_callback(particle_it, recv_data_it);

        if (build_cache) // TODO: is this safe?
          ghost_particles_iterators.push_back(particle_it);
      }

    AssertThrow(recv_data_it == recv_data.data() + recv_data.size(),
                ExcMessage(
                  "The amount of data that was read into new particles "
                  "does not match the amount of data sent around."));
  }
#endif



#ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::send_recv_particles_properties_and_location(
    const std::map<types::subdomain_id, std::vector<particle_iterator>>
      &particles_to_send)
  {
    const auto parallel_triangulation =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*triangulation);
    Assert(
      parallel_triangulation,
      ExcMessage(
        "This function is only implemented for parallel::TriangulationBase "
        "objects."));

    const auto &neighbors     = ghost_particles_cache.neighbors;
    const auto &send_pointers = ghost_particles_cache.send_pointers;
    const auto &recv_pointers = ghost_particles_cache.recv_pointers;

    std::vector<char> &send_data = ghost_particles_cache.send_data;

    // Fill data to send
    if (send_pointers.back() > 0)
      {
        void *data = static_cast<void *>(&send_data.front());

        // Serialize the data sorted by receiving process
        for (const auto i : neighbors)
          for (const auto &p : particles_to_send.at(i))
            {
              data = p->write_particle_data_to_memory(data);
              if (store_callback)
                data = store_callback(p, data);
            }
      }

    std::vector<char> &recv_data = ghost_particles_cache.recv_data;

    // Exchange the particle data between domains
    {
      std::vector<MPI_Request> requests(2 * neighbors.size());
      unsigned int             send_ops = 0;
      unsigned int             recv_ops = 0;

      const int mpi_tag = Utilities::MPI::internal::Tags::
        particle_handler_send_recv_particles_send;

      for (unsigned int i = 0; i < neighbors.size(); ++i)
        if ((recv_pointers[i + 1] - recv_pointers[i]) > 0)
          {
            const int ierr =
              MPI_Irecv(recv_data.data() + recv_pointers[i],
                        recv_pointers[i + 1] - recv_pointers[i],
                        MPI_CHAR,
                        neighbors[i],
                        mpi_tag,
                        parallel_triangulation->get_mpi_communicator(),
                        &(requests[send_ops]));
            AssertThrowMPI(ierr);
            ++send_ops;
          }

      for (unsigned int i = 0; i < neighbors.size(); ++i)
        if ((send_pointers[i + 1] - send_pointers[i]) > 0)
          {
            const int ierr =
              MPI_Isend(send_data.data() + send_pointers[i],
                        send_pointers[i + 1] - send_pointers[i],
                        MPI_CHAR,
                        neighbors[i],
                        mpi_tag,
                        parallel_triangulation->get_mpi_communicator(),
                        &(requests[send_ops + recv_ops]));
            AssertThrowMPI(ierr);
            ++recv_ops;
          }
      const int ierr =
        MPI_Waitall(send_ops + recv_ops, requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Put the received particles into the domain if they are in the
    // triangulation
    const void *recv_data_it = static_cast<const void *>(recv_data.data());

    // Gather ghost particle iterators from the cache
    auto &ghost_particles_iterators =
      ghost_particles_cache.ghost_particles_iterators;

    for (auto &recv_particle : ghost_particles_iterators)
      {
        // Update particle data using previously allocated memory space
        // for efficiency reasons
        recv_data_it =
          recv_particle->read_particle_data_from_memory(recv_data_it);

        Assert(recv_particle->particles_in_cell->cell->is_ghost(),
               ExcInternalError());

        if (load_callback)
          recv_data_it = load_callback(
            particle_iterator(recv_particle->particles_in_cell,
                              *property_pool,
                              recv_particle->particle_index_within_cell),
            recv_data_it);
      }

    AssertThrow(recv_data_it == recv_data.data() + recv_data.size(),
                ExcMessage(
                  "The amount of data that was read into new particles "
                  "does not match the amount of data sent around."));
  }
#endif

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::register_additional_store_load_functions(
    const std::function<std::size_t()>                             &size_callb,
    const std::function<void *(const particle_iterator &, void *)> &store_callb,
    const std::function<const void *(const particle_iterator &, const void *)>
      &load_callb)
  {
    size_callback  = size_callb;
    store_callback = store_callb;
    load_callback  = load_callb;
  }


  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::connect_to_triangulation_signals()
  {
    // First disconnect existing connections
    for (const auto &connection : tria_listeners)
      connection.disconnect();

    tria_listeners.clear();

    tria_listeners.push_back(triangulation->signals.create.connect([&]() {
      this->initialize(*(this->triangulation),
                       *(this->mapping),
                       this->property_pool->n_properties_per_slot());
    }));

    this->tria_listeners.push_back(
      this->triangulation->signals.clear.connect([&]() { this->clear(); }));

    // for distributed triangulations, connect to distributed signals
    if (dynamic_cast<const parallel::DistributedTriangulationBase<dim, spacedim>
                       *>(&(*triangulation)) != nullptr)
      {
        tria_listeners.push_back(
          triangulation->signals.post_distributed_refinement.connect(
            [&]() { this->post_mesh_change_action(); }));
        tria_listeners.push_back(
          triangulation->signals.post_distributed_repartition.connect(
            [&]() { this->post_mesh_change_action(); }));
        tria_listeners.push_back(
          triangulation->signals.post_distributed_load.connect(
            [&]() { this->post_mesh_change_action(); }));
      }
    else
      {
        tria_listeners.push_back(triangulation->signals.post_refinement.connect(
          [&]() { this->post_mesh_change_action(); }));
      }
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::post_mesh_change_action()
  {
    Assert(triangulation != nullptr, ExcInternalError());

    const bool distributed_triangulation =
      dynamic_cast<
        const parallel::DistributedTriangulationBase<dim, spacedim> *>(
        &(*triangulation)) != nullptr;
    (void)distributed_triangulation;

    Assert(
      distributed_triangulation || number_of_locally_owned_particles == 0,
      ExcMessage(
        "Mesh refinement in a non-distributed triangulation is not supported "
        "by the ParticleHandler class. Either insert particles after mesh "
        "creation, or use a distributed triangulation."));

    // Resize the container if it is possible without
    // transferring particles
    if (number_of_locally_owned_particles == 0)
      cells_to_particle_cache.resize(triangulation->n_active_cells(),
                                     particles.end());
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::prepare_for_coarsening_and_refinement()
  {
    register_data_attach();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::prepare_for_serialization()
  {
    register_data_attach();
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::register_data_attach()
  {
    const auto callback_function =
      [this](const typename Triangulation<dim, spacedim>::cell_iterator
                             &cell_iterator,
             const CellStatus cell_status) {
        return this->pack_callback(cell_iterator, cell_status);
      };

    tria_attached_data_index =
      const_cast<Triangulation<dim, spacedim> *>(&*triangulation)
        ->register_data_attach(callback_function,
                               /*returns_variable_size_data=*/true);
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::unpack_after_coarsening_and_refinement()
  {
    const bool serialization = false;
    notify_ready_to_unpack(serialization);
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::deserialize()
  {
    const bool serialization = true;
    notify_ready_to_unpack(serialization);
  }


  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::notify_ready_to_unpack(
    const bool serialization)
  {
    // First prepare container for insertion
    clear();

    // If we are resuming from a checkpoint, we first have to register the
    // store function again, to set the triangulation to the same state as
    // before the serialization. Only afterwards we know how to deserialize the
    // data correctly.
    if (serialization)
      register_data_attach();

    // Check if something was stored and load it
    if (tria_attached_data_index != numbers::invalid_unsigned_int)
      {
        const auto callback_function =
          [this](const typename Triangulation<dim, spacedim>::cell_iterator
                                 &cell_iterator,
                 const CellStatus cell_status,
                 const boost::iterator_range<std::vector<char>::const_iterator>
                   &range_iterator) {
            this->unpack_callback(cell_iterator, cell_status, range_iterator);
          };

        const_cast<Triangulation<dim, spacedim> *>(&*triangulation)
          ->notify_ready_to_unpack(tria_attached_data_index, callback_function);

        // Reset handle and update global numbers.
        tria_attached_data_index = numbers::invalid_unsigned_int;
        update_cached_numbers();
      }
  }



  template <int dim, int spacedim>
  std::vector<char>
  ParticleHandler<dim, spacedim>::pack_callback(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellStatus                                            status) const
  {
    std::vector<particle_iterator> stored_particles_on_cell;

    switch (status)
      {
        case CellStatus::cell_will_persist:
        case CellStatus::cell_will_be_refined:
          // If the cell persist or is refined store all particles of the
          // current cell.
          {
            const unsigned int n_particles = n_particles_in_cell(cell);
            stored_particles_on_cell.reserve(n_particles);

            for (unsigned int i = 0; i < n_particles; ++i)
              stored_particles_on_cell.push_back(particle_iterator(
                cells_to_particle_cache[cell->active_cell_index()],
                *property_pool,
                i));
          }
          break;

        case CellStatus::children_will_be_coarsened:
          // If this cell is the parent of children that will be coarsened,
          // collect the particles of all children.
          {
            for (const auto &child : cell->child_iterators())
              {
                const unsigned int n_particles = n_particles_in_cell(child);

                stored_particles_on_cell.reserve(
                  stored_particles_on_cell.size() + n_particles);

                const typename particle_container::iterator &cache =
                  cells_to_particle_cache[child->active_cell_index()];
                for (unsigned int i = 0; i < n_particles; ++i)
                  stored_particles_on_cell.push_back(
                    particle_iterator(cache, *property_pool, i));
              }
          }
          break;

        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }

    return pack_particles(stored_particles_on_cell);
  }



  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::unpack_callback(
    const typename Triangulation<dim, spacedim>::cell_iterator     &cell,
    const CellStatus                                                status,
    const boost::iterator_range<std::vector<char>::const_iterator> &data_range)
  {
    if (data_range.begin() == data_range.end())
      return;

    const auto cell_to_store_particles =
      (status != CellStatus::cell_will_be_refined) ? cell : cell->child(0);

    // deserialize particles and insert into local storage
    if (data_range.begin() != data_range.end())
      {
        const void *data = static_cast<const void *>(&(*data_range.begin()));
        const void *end  = static_cast<const void *>(
          &(*data_range.begin()) + (data_range.end() - data_range.begin()));

        while (data < end)
          {
            const void *old_data = data;
            const auto  x = insert_particle(data, cell_to_store_particles);

            // Ensure that the particle read exactly as much data as
            // it promised it needs to store its data
            const void *new_data = data;
            (void)old_data;
            (void)new_data;
            (void)x;
            AssertDimension((const char *)new_data - (const char *)old_data,
                            x->serialized_size_in_bytes());
          }

        Assert(data == end,
               ExcMessage(
                 "The particle data could not be deserialized successfully. "
                 "Check that when deserializing the particles you expect "
                 "the same number of properties that were serialized."));
      }

    auto loaded_particles_on_cell = particles_in_cell(cell_to_store_particles);

    // now update particle storage location and properties if necessary
    switch (status)
      {
        case CellStatus::cell_will_persist:
          {
            // all particles are correctly inserted
          }
          break;

        case CellStatus::children_will_be_coarsened:
          {
            // all particles are in correct cell, but their reference location
            // has changed
            for (auto &particle : loaded_particles_on_cell)
              {
                const Point<dim> p_unit =
                  mapping->transform_real_to_unit_cell(cell_to_store_particles,
                                                       particle.get_location());
                particle.set_reference_location(p_unit);
              }
          }
          break;

        case CellStatus::cell_will_be_refined:
          {
            // we need to find the correct child to store the particles and
            // their reference location has changed
            typename particle_container::iterator &cache =
              cells_to_particle_cache[cell_to_store_particles
                                        ->active_cell_index()];

            // make sure that the call above has inserted an entry
            Assert(cache != particles.end(), ExcInternalError());

            // Cannot use range-based loop, because number of particles in cell
            // is going to change
            auto particle = loaded_particles_on_cell.begin();
            for (unsigned int i = 0; i < cache->particles.size();)
              {
                bool found_new_cell = false;

                for (const auto &child : cell->child_iterators())
                  {
                    Assert(child->is_locally_owned(), ExcInternalError());

                    try
                      {
                        const Point<dim> p_unit =
                          mapping->transform_real_to_unit_cell(
                            child, particle->get_location());
                        if (cell->reference_cell().contains_point(
                              p_unit, tolerance_inside_cell))
                          {
                            found_new_cell = true;
                            particle->set_reference_location(p_unit);

                            // if the particle is not in the cell we stored it
                            // in above, its handle is in the wrong place
                            if (child != cell_to_store_particles)
                              {
                                // move handle into correct cell
                                insert_particle(cache->particles[i], child);
                                // remove handle by replacing it with last one
                                cache->particles[i] = cache->particles.back();
                                cache->particles.pop_back();
                                // no loop increment, we need to process
                                // the new i-th particle.
                              }
                            else
                              {
                                // move on to next particle
                                ++i;
                                ++particle;
                              }
                            break;
                          }
                      }
                    catch (typename Mapping<dim>::ExcTransformationFailed &)
                      {}
                  }

                if (found_new_cell == false)
                  {
                    // If we get here, we did not find the particle in any
                    // child. This case may happen for particles that are at the
                    // boundary for strongly curved cells. We apply a tolerance
                    // in the call to ReferenceCell::contains_point() to
                    // account for this, but if that is not enough, we still
                    // need to prevent an endless loop here. Delete the particle
                    // and move on.
                    signals.particle_lost(particle,
                                          particle->get_surrounding_cell());
                    if (cache->particles[i] !=
                        PropertyPool<dim, spacedim>::invalid_handle)
                      property_pool->deregister_particle(cache->particles[i]);
                    cache->particles[i] = cache->particles.back();
                    cache->particles.pop_back();
                  }
              }
            // clean up in case child 0 has no particle left
            if (cache->particles.empty())
              {
                particles.erase(cache);
                cache = particles.end();
              }
          }
          break;

        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }
  }
} // namespace Particles

#include "particles/particle_handler.inst"

DEAL_II_NAMESPACE_CLOSE
