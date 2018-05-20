// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#include <deal.II/particles/particle_handler.h>

#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <utility>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST

namespace Particles
{
  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::ParticleHandler()
    : triangulation(),
      particles(),
      ghost_particles(),
      global_number_of_particles(0),
      global_max_particles_per_cell(0),
      next_free_particle_index(0),
      property_pool(new PropertyPool(0)),
      size_callback(),
      store_callback(),
      load_callback(),
      handle(numbers::invalid_unsigned_int)
  {}

  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::ParticleHandler(
    const parallel::distributed::Triangulation<dim, spacedim>& triangulation,
    const Mapping<dim, spacedim>&                              mapping,
    const unsigned int                                         n_properties)
    : triangulation(&triangulation, typeid(*this).name()),
      mapping(&mapping, typeid(*this).name()),
      particles(),
      ghost_particles(),
      global_number_of_particles(0),
      global_max_particles_per_cell(0),
      next_free_particle_index(0),
      property_pool(new PropertyPool(n_properties)),
      size_callback(),
      store_callback(),
      load_callback(),
      handle(numbers::invalid_unsigned_int)
  {}

  template <int dim, int spacedim>
  ParticleHandler<dim, spacedim>::~ParticleHandler()
  {}

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::initialize(
    const parallel::distributed::Triangulation<dim, spacedim>& tria,
    const Mapping<dim, spacedim>&                              mapp,
    const unsigned int                                         n_properties)
  {
    triangulation = &tria;
    mapping       = &mapp;

    // Create the memory pool that will store all particle properties
    property_pool = std_cxx14::make_unique<PropertyPool>(n_properties);
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::clear()
  {
    clear_particles();
    global_number_of_particles    = 0;
    next_free_particle_index      = 0;
    global_max_particles_per_cell = 0;
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::clear_particles()
  {
    particles.clear();
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::update_cached_numbers()
  {
    types::particle_index locally_highest_index        = 0;
    unsigned int          local_max_particles_per_cell = 0;
    unsigned int          current_particles_per_cell   = 0;
    typename Triangulation<dim, spacedim>::active_cell_iterator current_cell
      = triangulation->begin_active();

    for(particle_iterator particle = begin(); particle != end(); ++particle)
      {
        locally_highest_index
          = std::max(locally_highest_index, particle->get_id());

        if(particle->get_surrounding_cell(*triangulation) != current_cell)
          {
            current_particles_per_cell = 0;
            current_cell = particle->get_surrounding_cell(*triangulation);
          }

        ++current_particles_per_cell;
        local_max_particles_per_cell
          = std::max(local_max_particles_per_cell, current_particles_per_cell);
      }

    global_number_of_particles = dealii::Utilities::MPI::sum(
      particles.size(), triangulation->get_communicator());
    next_free_particle_index
      = dealii::Utilities::MPI::max(locally_highest_index,
                                    triangulation->get_communicator())
        + 1;
    global_max_particles_per_cell = dealii::Utilities::MPI::max(
      local_max_particles_per_cell, triangulation->get_communicator());
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin() const
  {
    return (const_cast<ParticleHandler<dim, spacedim>*>(this))->begin();
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin()
  {
    return particle_iterator(particles, particles.begin());
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end() const
  {
    return (const_cast<ParticleHandler<dim, spacedim>*>(this))->end();
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end()
  {
    return particle_iterator(particles, particles.end());
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim>*>(this))->begin_ghost();
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost()
  {
    return particle_iterator(ghost_particles, ghost_particles.begin());
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim>*>(this))->end_ghost();
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost()
  {
    return particle_iterator(ghost_particles, ghost_particles.end());
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator_range
  ParticleHandler<dim, spacedim>::particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator& cell)
    const
  {
    return (const_cast<ParticleHandler<dim, spacedim>*>(this))
      ->particles_in_cell(cell);
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator_range
  ParticleHandler<dim, spacedim>::particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator& cell)
  {
    const internal::LevelInd level_index
      = std::make_pair<int, int>(cell->level(), cell->index());

    if(cell->is_ghost())
      {
        const auto particles_in_cell = ghost_particles.equal_range(level_index);
        return boost::make_iterator_range(
          particle_iterator(ghost_particles, particles_in_cell.first),
          particle_iterator(ghost_particles, particles_in_cell.second));
      }

    const auto particles_in_cell = particles.equal_range(level_index);
    return boost::make_iterator_range(
      particle_iterator(particles, particles_in_cell.first),
      particle_iterator(particles, particles_in_cell.second));
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::remove_particle(
    const ParticleHandler<dim, spacedim>::particle_iterator& particle)
  {
    particles.erase(particle->particle);
  }

  template <int dim, int spacedim>
  typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::insert_particle(
    const Particle<dim, spacedim>&                                     particle,
    const typename Triangulation<dim, spacedim>::active_cell_iterator& cell)
  {
    typename std::multimap<internal::LevelInd,
                           Particle<dim, spacedim>>::iterator it
      = particles.insert(std::make_pair(
        internal::LevelInd(cell->level(), cell->index()), particle));

    particle_iterator particle_it(particles, it);
    particle_it->set_property_pool(*property_pool);

    if(particle.has_properties())
      for(unsigned int n = 0; n < particle.get_properties().size(); ++n)
        particle_it->get_properties()[n] = particle.get_properties()[n];

    return particle_it;
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::insert_particles(
    const std::multimap<
      typename Triangulation<dim, spacedim>::active_cell_iterator,
      Particle<dim, spacedim>>& new_particles)
  {
    for(auto particle = new_particles.begin(); particle != new_particles.end();
        ++particle)
      particles.insert(
        particles.end(),
        std::make_pair(internal::LevelInd(particle->first->level(),
                                          particle->first->index()),
                       particle->second));

    update_cached_numbers();
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::insert_particles(
    const std::vector<Point<spacedim>>& positions)
  {
    update_cached_numbers();

    // Determine the starting particle index of this process, which
    // is the highest currently existing particle index plus the sum
    // of the number of newly generated particles of all
    // processes with a lower rank if in a parallel computation.
    const types::particle_index local_next_particle_index
      = get_next_free_particle_index();
    types::particle_index local_start_index = 0;

#  ifdef DEAL_II_WITH_MPI
    types::particle_index particles_to_add_locally = positions.size();
    const int             ierr = MPI_Scan(&particles_to_add_locally,
                              &local_start_index,
                              1,
                              PARTICLE_INDEX_MPI_TYPE,
                              MPI_SUM,
                              triangulation->get_communicator());
    AssertThrowMPI(ierr);
    local_start_index -= particles_to_add_locally;
#  endif

    local_start_index += local_next_particle_index;

    GridTools::Cache<dim, spacedim> cache(*triangulation, *mapping);
    auto point_locations = GridTools::compute_point_locations(cache, positions);

    auto& cells           = std::get<0>(point_locations);
    auto& local_positions = std::get<1>(point_locations);
    auto& index_map       = std::get<2>(point_locations);

    if(cells.size() == 0)
      return;

    auto hint
      = particles.find(std::make_pair(cells[0]->level(), cells[0]->index()));
    for(unsigned int i = 0; i < cells.size(); ++i)
      {
        internal::LevelInd current_cell(cells[i]->level(), cells[i]->index());
        for(unsigned int p = 0; p < local_positions[i].size(); ++p)
          {
            hint = particles.insert(
              hint,
              std::make_pair(
                current_cell,
                Particle<dim, spacedim>(positions[index_map[i][p]],
                                        local_positions[i][p],
                                        local_start_index + index_map[i][p])));
          }
      }

    update_cached_numbers();
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
    return particles.size();
  }

  template <int dim, int spacedim>
  unsigned int
  ParticleHandler<dim, spacedim>::n_properties_per_particle() const
  {
    return property_pool->n_properties_per_slot();
  }

  template <int dim, int spacedim>
  unsigned int
  ParticleHandler<dim, spacedim>::n_particles_in_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator& cell)
    const
  {
    const internal::LevelInd found_cell
      = std::make_pair<int, int>(cell->level(), cell->index());

    if(cell->is_locally_owned())
      return particles.count(found_cell);
    else if(cell->is_ghost())
      return ghost_particles.count(found_cell);
    else if(cell->is_artificial())
      AssertThrow(false, ExcInternalError());

    return 0;
  }

  template <int dim, int spacedim>
  types::particle_index
  ParticleHandler<dim, spacedim>::get_next_free_particle_index() const
  {
    return next_free_particle_index;
  }

  template <int dim, int spacedim>
  PropertyPool&
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
      const Tensor<1, dim>&              particle_direction,
      const std::vector<Tensor<1, dim>>& center_directions)
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
    // TODO: The current algorithm only works for particles that are in
    // the local domain or in ghost cells, because it only knows the
    // subdomain_id of ghost cells, but not of artificial cells. This
    // can be extended using the distributed version of compute point
    // locations.
    // TODO: Extend this function to allow keeping particles on other
    // processes around (with an invalid cell).

    std::vector<particle_iterator> particles_out_of_cell;
    particles_out_of_cell.reserve(n_locally_owned_particles());

    // Now update the reference locations of the moved particles
    for(particle_iterator it = begin(); it != end(); ++it)
      {
        const typename Triangulation<dim, spacedim>::cell_iterator cell
          = it->get_surrounding_cell(*triangulation);

        try
          {
            const Point<dim> p_unit
              = mapping->transform_real_to_unit_cell(cell, it->get_location());
            if(GeometryInfo<dim>::is_inside_unit_cell(p_unit))
              {
                it->set_reference_location(p_unit);
              }
            else
              {
                // The particle has left the cell
                particles_out_of_cell.push_back(it);
              }
          }
        catch(typename Mapping<dim>::ExcTransformationFailed&)
          {
            // The particle has left the cell
            particles_out_of_cell.push_back(it);
          }
      }

    // There are three reasons why a particle is not in its old cell:
    // It moved to another cell, to another subdomain or it left the mesh.
    // Particles that moved to another cell are updated and stored inside the
    // sorted_particles vector, particles that moved to another domain are
    // collected in the moved_particles_domain vector. Particles that left
    // the mesh completely are ignored and removed.
    std::vector<std::pair<internal::LevelInd, Particle<dim, spacedim>>>
      sorted_particles;
    std::map<types::subdomain_id, std::vector<particle_iterator>>
      moved_particles;
    std::map<
      types::subdomain_id,
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      moved_cells;

    // We do not know exactly how many particles are lost, exchanged between
    // domains, or remain on this process. Therefore we pre-allocate approximate
    // sizes for these vectors. If more space is needed an automatic and
    // relatively fast (compared to other parts of this algorithm)
    // re-allocation will happen.
    typedef typename std::vector<particle_iterator>::size_type vector_size;
    sorted_particles.reserve(
      static_cast<vector_size>(particles_out_of_cell.size() * 1.25));

    const std::set<types::subdomain_id> ghost_owners
      = triangulation->ghost_owners();

    for(auto ghost_domain_id = ghost_owners.begin();
        ghost_domain_id != ghost_owners.end();
        ++ghost_domain_id)
      moved_particles[*ghost_domain_id].reserve(
        static_cast<vector_size>(particles_out_of_cell.size() * 0.25));
    for(auto ghost_domain_id = ghost_owners.begin();
        ghost_domain_id != ghost_owners.end();
        ++ghost_domain_id)
      moved_cells[*ghost_domain_id].reserve(
        static_cast<vector_size>(particles_out_of_cell.size() * 0.25));

    {
      // Create a map from vertices to adjacent cells
      const std::vector<
        std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
        vertex_to_cells(GridTools::vertex_to_cell_map(*triangulation));

      // Create a corresponding map of vectors from vertex to cell center
      const std::vector<std::vector<Tensor<1, spacedim>>>
        vertex_to_cell_centers(GridTools::vertex_to_cell_centers_directions(
          *triangulation, vertex_to_cells));

      std::vector<unsigned int> neighbor_permutation;

      // Find the cells that the particles moved to.
      typename std::vector<particle_iterator>::iterator it
        = particles_out_of_cell.begin(),
        end_particle = particles_out_of_cell.end();

      for(; it != end_particle; ++it)
        {
          // The cell the particle is in
          Point<dim> current_reference_position;
          bool       found_cell = false;

          // Check if the particle is in one of the old cell's neighbors
          // that are adjacent to the closest vertex
          typename Triangulation<dim, spacedim>::active_cell_iterator
            current_cell
            = (*it)->get_surrounding_cell(*triangulation);

          const unsigned int closest_vertex
            = GridTools::find_closest_vertex_of_cell<dim, spacedim>(
              current_cell, (*it)->get_location());
          Tensor<1, spacedim> vertex_to_particle
            = (*it)->get_location() - current_cell->vertex(closest_vertex);
          vertex_to_particle /= vertex_to_particle.norm();

          const unsigned int closest_vertex_index
            = current_cell->vertex_index(closest_vertex);
          const unsigned int n_neighbor_cells
            = vertex_to_cells[closest_vertex_index].size();

          neighbor_permutation.resize(n_neighbor_cells);
          for(unsigned int i = 0; i < n_neighbor_cells; ++i)
            neighbor_permutation[i] = i;

          std::sort(
            neighbor_permutation.begin(),
            neighbor_permutation.end(),
            std::bind(&compare_particle_association<spacedim>,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::cref(vertex_to_particle),
                      std::cref(vertex_to_cell_centers[closest_vertex_index])));

          // Search all of the cells adjacent to the closest vertex of the previous cell
          // Most likely we will find the particle in them.
          for(unsigned int i = 0; i < n_neighbor_cells; ++i)
            {
              try
                {
                  typename std::set<typename Triangulation<dim, spacedim>::
                                      active_cell_iterator>::const_iterator cell
                    = vertex_to_cells[closest_vertex_index].begin();

                  std::advance(cell, neighbor_permutation[i]);
                  const Point<dim> p_unit
                    = mapping->transform_real_to_unit_cell(
                      *cell, (*it)->get_location());
                  if(GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                    {
                      current_cell               = *cell;
                      current_reference_position = p_unit;
                      found_cell                 = true;
                      break;
                    }
                }
              catch(typename Mapping<dim>::ExcTransformationFailed&)
                {}
            }

          if(!found_cell)
            {
              // The particle is not in a neighbor of the old cell.
              // Look for the new cell in the whole local domain.
              // This case is rare.
              try
                {
                  const std::pair<const typename Triangulation<dim, spacedim>::
                                    active_cell_iterator,
                                  Point<dim>>
                    current_cell_and_position
                    = GridTools::find_active_cell_around_point<>(
                      *mapping, *triangulation, (*it)->get_location());
                  current_cell               = current_cell_and_position.first;
                  current_reference_position = current_cell_and_position.second;
                }
              catch(GridTools::ExcPointNotFound<spacedim>&)
                {
                  // We can find no cell for this particle. It has left the
                  // domain due to an integration error or an open boundary.
                  continue;
                }
            }

          // If we are here, we found a cell and reference position for this particle
          (*it)->set_reference_location(current_reference_position);

          // Reinsert the particle into our domain if we own its cell.
          // Mark it for MPI transfer otherwise
          if(current_cell->is_locally_owned())
            {
              sorted_particles.push_back(
                std::make_pair(internal::LevelInd(current_cell->level(),
                                                  current_cell->index()),
                               (*it)->particle->second));
            }
          else
            {
              moved_particles[current_cell->subdomain_id()].push_back(*it);
              moved_cells[current_cell->subdomain_id()].push_back(current_cell);
            }
        }
    }

    // Sort the updated particles. This pre-sort speeds up inserting
    // them into particles to O(N) complexity.
    std::multimap<internal::LevelInd, Particle<dim, spacedim>>
      sorted_particles_map;

    // Exchange particles between processors if we have more than one process
#  ifdef DEAL_II_WITH_MPI
    if(dealii::Utilities::MPI::n_mpi_processes(
         triangulation->get_communicator())
       > 1)
      send_recv_particles(moved_particles, sorted_particles_map, moved_cells);
#  endif

    sorted_particles_map.insert(sorted_particles.begin(),
                                sorted_particles.end());

    for(unsigned int i = 0; i < particles_out_of_cell.size(); ++i)
      remove_particle(particles_out_of_cell[i]);

    particles.insert(sorted_particles_map.begin(), sorted_particles_map.end());
    update_cached_numbers();
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::exchange_ghost_particles()
  {
    // Nothing to do in serial computations
    if(dealii::Utilities::MPI::n_mpi_processes(
         triangulation->get_communicator())
       == 1)
      return;

#  ifdef DEAL_II_WITH_MPI
    // First clear the current ghost_particle information
    ghost_particles.clear();

    std::map<types::subdomain_id, std::vector<particle_iterator>>
      ghost_particles_by_domain;

    const std::set<types::subdomain_id> ghost_owners
      = triangulation->ghost_owners();
    for(auto ghost_domain_id = ghost_owners.begin();
        ghost_domain_id != ghost_owners.end();
        ++ghost_domain_id)
      ghost_particles_by_domain[*ghost_domain_id].reserve(
        static_cast<typename std::vector<particle_iterator>::size_type>(
          particles.size() * 0.25));

    std::vector<std::set<unsigned int>> vertex_to_neighbor_subdomain(
      triangulation->n_vertices());

    typename Triangulation<dim, spacedim>::active_cell_iterator cell
      = triangulation->begin_active(),
      endc = triangulation->end();
    for(; cell != endc; ++cell)
      {
        if(cell->is_ghost())
          for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
            vertex_to_neighbor_subdomain[cell->vertex_index(v)].insert(
              cell->subdomain_id());
      }

    cell = triangulation->begin_active();
    for(; cell != endc; ++cell)
      {
        if(!cell->is_ghost())
          {
            std::set<unsigned int> cell_to_neighbor_subdomain;
            for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                ++v)
              {
                cell_to_neighbor_subdomain.insert(
                  vertex_to_neighbor_subdomain[cell->vertex_index(v)].begin(),
                  vertex_to_neighbor_subdomain[cell->vertex_index(v)].end());
              }

            if(cell_to_neighbor_subdomain.size() > 0)
              {
                const particle_iterator_range particle_range
                  = particles_in_cell(cell);

                for(std::set<types::subdomain_id>::const_iterator domain
                    = cell_to_neighbor_subdomain.begin();
                    domain != cell_to_neighbor_subdomain.end();
                    ++domain)
                  {
                    for(typename particle_iterator_range::iterator particle
                        = particle_range.begin();
                        particle != particle_range.end();
                        ++particle)
                      ghost_particles_by_domain[*domain].push_back(particle);
                  }
              }
          }
      }

    send_recv_particles(ghost_particles_by_domain, ghost_particles);
#  endif
  }

#  ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::send_recv_particles(
    const std::map<types::subdomain_id, std::vector<particle_iterator>>&
      particles_to_send,
    std::multimap<internal::LevelInd, Particle<dim, spacedim>>&
      received_particles,
    const std::map<
      types::subdomain_id,
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>>&
      send_cells)
  {
    // Determine the communication pattern
    const std::set<types::subdomain_id> ghost_owners
      = triangulation->ghost_owners();
    const std::vector<types::subdomain_id> neighbors(ghost_owners.begin(),
                                                     ghost_owners.end());
    const unsigned int                     n_neighbors = neighbors.size();

    if(send_cells.size() != 0)
      Assert(particles_to_send.size() == send_cells.size(), ExcInternalError());

    // If we do not know the subdomain this particle needs to be send to, throw an error
    Assert(particles_to_send.find(numbers::artificial_subdomain_id)
             == particles_to_send.end(),
           ExcInternalError());

    // TODO: Implement the shipping of particles to processes that are not ghost owners of the local domain
    for(auto send_particles = particles_to_send.begin();
        send_particles != particles_to_send.end();
        ++send_particles)
      Assert(ghost_owners.find(send_particles->first) != ghost_owners.end(),
             ExcNotImplemented());

    unsigned int n_send_particles = 0;
    for(auto send_particles = particles_to_send.begin();
        send_particles != particles_to_send.end();
        ++send_particles)
      n_send_particles += send_particles->second.size();

    const unsigned int cellid_size = sizeof(CellId::binary_type);

    // Containers for the amount and offsets of data we will send
    // to other processors and the data itself.
    std::vector<unsigned int> n_send_data(n_neighbors, 0);
    std::vector<unsigned int> send_offsets(n_neighbors, 0);
    std::vector<char>         send_data;

    // Only serialize things if there are particles to be send.
    // We can not return early even if no particles
    // are send, because we might receive particles from other processes
    if(n_send_particles > 0)
      {
        // Allocate space for sending particle data
        const unsigned int particle_size
          = begin()->serialized_size_in_bytes() + cellid_size
            + (size_callback ? size_callback() : 0);
        send_data.resize(n_send_particles * particle_size);
        void* data = static_cast<void*>(&send_data.front());

        // Serialize the data sorted by receiving process
        for(unsigned int i = 0; i < n_neighbors; ++i)
          {
            send_offsets[i]
              = reinterpret_cast<std::size_t>(data)
                - reinterpret_cast<std::size_t>(&send_data.front());

            for(unsigned int j = 0;
                j < particles_to_send.at(neighbors[i]).size();
                ++j)
              {
                // If no target cells are given, use the iterator information
                typename Triangulation<dim, spacedim>::active_cell_iterator
                  cell;
                if(send_cells.size() == 0)
                  cell = particles_to_send.at(neighbors[i])[j]
                           ->get_surrounding_cell(*triangulation);
                else
                  cell = send_cells.at(neighbors[i])[j];

                const CellId::binary_type cellid
                  = cell->id().template to_binary<dim>();
                memcpy(data, &cellid, cellid_size);
                data = static_cast<char*>(data) + cellid_size;

                particles_to_send.at(neighbors[i])[j]->write_data(data);
                if(store_callback)
                  data = store_callback(particles_to_send.at(neighbors[i])[j],
                                        data);
              }
            n_send_data[i]
              = reinterpret_cast<std::size_t>(data) - send_offsets[i]
                - reinterpret_cast<std::size_t>(&send_data.front());
          }
      }

    // Containers for the data we will receive from other processors
    std::vector<unsigned int> n_recv_data(n_neighbors);
    std::vector<unsigned int> recv_offsets(n_neighbors);

    // Notify other processors how many particles we will send
    {
      std::vector<MPI_Request> n_requests(2 * n_neighbors);
      for(unsigned int i = 0; i < n_neighbors; ++i)
        {
          const int ierr = MPI_Irecv(&(n_recv_data[i]),
                                     1,
                                     MPI_INT,
                                     neighbors[i],
                                     0,
                                     triangulation->get_communicator(),
                                     &(n_requests[2 * i]));
          AssertThrowMPI(ierr);
        }
      for(unsigned int i = 0; i < n_neighbors; ++i)
        {
          const int ierr = MPI_Isend(&(n_send_data[i]),
                                     1,
                                     MPI_INT,
                                     neighbors[i],
                                     0,
                                     triangulation->get_communicator(),
                                     &(n_requests[2 * i + 1]));
          AssertThrowMPI(ierr);
        }
      const int ierr
        = MPI_Waitall(2 * n_neighbors, &n_requests[0], MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Determine how many particles and data we will receive
    unsigned int total_recv_data = 0;
    for(unsigned int neighbor_id = 0; neighbor_id < n_neighbors; ++neighbor_id)
      {
        recv_offsets[neighbor_id] = total_recv_data;
        total_recv_data += n_recv_data[neighbor_id];
      }

    // Set up the space for the received particle data
    std::vector<char> recv_data(total_recv_data);

    // Exchange the particle data between domains
    {
      std::vector<MPI_Request> requests(2 * n_neighbors);
      unsigned int             send_ops = 0;
      unsigned int             recv_ops = 0;

      for(unsigned int i = 0; i < n_neighbors; ++i)
        if(n_recv_data[i] > 0)
          {
            const int ierr = MPI_Irecv(&(recv_data[recv_offsets[i]]),
                                       n_recv_data[i],
                                       MPI_CHAR,
                                       neighbors[i],
                                       1,
                                       triangulation->get_communicator(),
                                       &(requests[send_ops]));
            AssertThrowMPI(ierr);
            send_ops++;
          }

      for(unsigned int i = 0; i < n_neighbors; ++i)
        if(n_send_data[i] > 0)
          {
            const int ierr = MPI_Isend(&(send_data[send_offsets[i]]),
                                       n_send_data[i],
                                       MPI_CHAR,
                                       neighbors[i],
                                       1,
                                       triangulation->get_communicator(),
                                       &(requests[send_ops + recv_ops]));
            AssertThrowMPI(ierr);
            recv_ops++;
          }
      const int ierr
        = MPI_Waitall(send_ops + recv_ops, &requests[0], MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Put the received particles into the domain if they are in the triangulation
    const void* recv_data_it = static_cast<const void*>(recv_data.data());

    while(reinterpret_cast<std::size_t>(recv_data_it)
            - reinterpret_cast<std::size_t>(recv_data.data())
          < total_recv_data)
      {
        CellId::binary_type binary_cellid;
        memcpy(&binary_cellid, recv_data_it, cellid_size);
        const CellId id(binary_cellid);
        recv_data_it = static_cast<const char*>(recv_data_it) + cellid_size;

        const typename Triangulation<dim, spacedim>::active_cell_iterator cell
          = id.to_cell(*triangulation);

        typename std::multimap<internal::LevelInd,
                               Particle<dim, spacedim>>::iterator recv_particle
          = received_particles.insert(std::make_pair(
            internal::LevelInd(cell->level(), cell->index()),
            Particle<dim, spacedim>(recv_data_it, property_pool.get())));

        if(load_callback)
          recv_data_it = load_callback(
            particle_iterator(received_particles, recv_particle), recv_data_it);
      }

    AssertThrow(
      recv_data_it == recv_data.data() + recv_data.size(),
      ExcMessage("The amount of data that was read into new particles "
                 "does not match the amount of data sent around."));
  }
#  endif

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::register_additional_store_load_functions(
    const std::function<std::size_t()>&                          size_callb,
    const std::function<void*(const particle_iterator&, void*)>& store_callb,
    const std::function<const void*(const particle_iterator&, const void*)>&
      load_callb)
  {
    size_callback  = size_callb;
    store_callback = store_callb;
    load_callback  = load_callb;
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::register_store_callback_function(
    const bool serialization)
  {
    parallel::distributed::Triangulation<dim, spacedim>* non_const_triangulation
      = const_cast<parallel::distributed::Triangulation<dim, spacedim>*>(
        &(*triangulation));

    // Only save and load particles if there are any, we might get here for
    // example if somebody created a ParticleHandler but generated 0 particles.
    update_cached_numbers();

    if(global_max_particles_per_cell > 0)
      {
        const std::function<void(
          const typename Triangulation<dim, spacedim>::cell_iterator&,
          const typename Triangulation<dim, spacedim>::CellStatus,
          void*)>
          callback_function
          = std::bind(&ParticleHandler<dim, spacedim>::store_particles,
                      std::cref(*this),
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3);

        // Compute the size per serialized particle. This is simple if we own
        // particles, simply ask one of them. Otherwise create a temporary particle,
        // ask it for its size and add the size of its properties.
        const std::size_t size_per_particle
          = (particles.size() > 0) ?
              begin()->serialized_size_in_bytes() :
              Particle<dim, spacedim>().serialized_size_in_bytes()
                + property_pool->n_properties_per_slot() * sizeof(double);

        // We need to transfer the number of particles for this cell and
        // the particle data itself. If we are in the process of refinement
        // (i.e. not in serialization) we need to provide 2^dim times the
        // space for the data in case a cell is coarsened and all particles
        // of the children have to be stored in the parent cell.
        const std::size_t transfer_size_per_cell
          = sizeof(unsigned int)
            + (size_per_particle * global_max_particles_per_cell)
                * (serialization ? 1 : std::pow(2, dim));

        handle = non_const_triangulation->register_data_attach(
          transfer_size_per_cell, callback_function);
      }
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::register_load_callback_function(
    const bool serialization)
  {
    // All particles have been stored, when we reach this point. Empty the
    // particle data.
    clear_particles();

    parallel::distributed::Triangulation<dim, spacedim>* non_const_triangulation
      = const_cast<parallel::distributed::Triangulation<dim, spacedim>*>(
        &(*triangulation));

    // If we are resuming from a checkpoint, we first have to register the
    // store function again, to set the triangulation in the same state as
    // before the serialization. Only by this it knows how to deserialize the
    // data correctly. Only do this if something was actually stored.
    if(serialization && (global_max_particles_per_cell > 0))
      {
        const std::function<void(
          const typename Triangulation<dim, spacedim>::cell_iterator&,
          const typename Triangulation<dim, spacedim>::CellStatus,
          void*)>
          callback_function
          = std::bind(&ParticleHandler<dim, spacedim>::store_particles,
                      std::cref(*this),
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3);

        // Compute the size per serialized particle. This is simple if we own
        // particles, simply ask one of them. Otherwise create a temporary particle,
        // ask it for its size and add the size of its properties.
        const std::size_t size_per_particle
          = (particles.size() > 0) ?
              begin()->serialized_size_in_bytes() :
              Particle<dim, spacedim>().serialized_size_in_bytes()
                + property_pool->n_properties_per_slot() * sizeof(double);

        // We need to transfer the number of particles for this cell and
        // the particle data itself and we need to provide 2^dim times the
        // space for the data in case a cell is coarsened
        const std::size_t transfer_size_per_cell
          = sizeof(unsigned int)
            + (size_per_particle * global_max_particles_per_cell);
        handle = non_const_triangulation->register_data_attach(
          transfer_size_per_cell, callback_function);
      }

    // Check if something was stored and load it
    if(handle != numbers::invalid_unsigned_int)
      {
        const std::function<void(
          const typename Triangulation<dim, spacedim>::cell_iterator&,
          const typename Triangulation<dim, spacedim>::CellStatus,
          const void*)>
          callback_function
          = std::bind(&ParticleHandler<dim, spacedim>::load_particles,
                      std::ref(*this),
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3);

        non_const_triangulation->notify_ready_to_unpack(handle,
                                                        callback_function);

        // Reset handle and update global number of particles. The number
        // can change because of discarded or newly generated particles
        handle = numbers::invalid_unsigned_int;
        update_cached_numbers();
      }
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::store_particles(
    const typename Triangulation<dim, spacedim>::cell_iterator& cell,
    const typename Triangulation<dim, spacedim>::CellStatus     status,
    void*                                                       data) const
  {
    unsigned int n_particles(0);

    // If the cell persist or is refined store all particles of the current cell.
    if(status
         == parallel::distributed::Triangulation<dim, spacedim>::CELL_PERSIST
       || status
            == parallel::distributed::Triangulation<dim, spacedim>::CELL_REFINE)
      {
        const boost::iterator_range<particle_iterator> particle_range
          = particles_in_cell(cell);
        n_particles
          = std::distance(particle_range.begin(), particle_range.end());

        unsigned int* ndata = static_cast<unsigned int*>(data);
        *ndata              = n_particles;
        data                = static_cast<void*>(ndata + 1);

        for(particle_iterator particle = particle_range.begin();
            particle != particle_range.end();
            ++particle)
          {
            particle->write_data(data);
          }
      }
    // If this cell is the parent of children that will be coarsened, collect
    // the particles of all children.
    else if(status
            == parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_COARSEN)
      {
        for(unsigned int child_index = 0;
            child_index < GeometryInfo<dim>::max_children_per_cell;
            ++child_index)
          {
            const typename Triangulation<dim, spacedim>::cell_iterator child
              = cell->child(child_index);
            n_particles += n_particles_in_cell(child);
          }

        unsigned int* ndata = static_cast<unsigned int*>(data);
        *ndata              = n_particles;

        data = static_cast<void*>(ndata + 1);

        for(unsigned int child_index = 0;
            child_index < GeometryInfo<dim>::max_children_per_cell;
            ++child_index)
          {
            const typename Triangulation<dim, spacedim>::cell_iterator child
              = cell->child(child_index);
            const boost::iterator_range<particle_iterator> particle_range
              = particles_in_cell(child);

            for(particle_iterator particle = particle_range.begin();
                particle != particle_range.end();
                ++particle)
              {
                particle->write_data(data);
              }
          }
      }
    else
      Assert(false, ExcInternalError());
  }

  template <int dim, int spacedim>
  void
  ParticleHandler<dim, spacedim>::load_particles(
    const typename Triangulation<dim, spacedim>::cell_iterator& cell,
    const typename Triangulation<dim, spacedim>::CellStatus     status,
    const void*                                                 data)
  {
    const unsigned int* n_particles_in_cell_ptr
      = static_cast<const unsigned int*>(data);
    const void* pdata
      = reinterpret_cast<const void*>(n_particles_in_cell_ptr + 1);

    if(*n_particles_in_cell_ptr == 0)
      return;

    // Load all particles from the data stream and store them in the local
    // particle map.
    if(status
       == parallel::distributed::Triangulation<dim, spacedim>::CELL_PERSIST)
      {
        typename std::multimap<internal::LevelInd,
                               Particle<dim, spacedim>>::iterator position_hint
          = particles.end();
        for(unsigned int i = 0; i < *n_particles_in_cell_ptr; ++i)
          {
            // Use std::multimap::emplace_hint to speed up insertion of
            // particles. This is a C++11 function, but not all compilers
            // that report a -std=c++11 (like gcc 4.6) implement it, so
            // require C++14 instead.
#  ifdef DEAL_II_WITH_CXX14
            position_hint = particles.emplace_hint(
              position_hint,
              std::make_pair(cell->level(), cell->index()),
              Particle<dim, spacedim>(pdata, property_pool.get()));
#  else
            position_hint = particles.insert(
              position_hint,
              std::make_pair(
                std::make_pair(cell->level(), cell->index()),
                Particle<dim, spacedim>(pdata, property_pool.get())));
#  endif
            ++position_hint;
          }
      }

    else if(status
            == parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_COARSEN)
      {
        typename std::multimap<internal::LevelInd,
                               Particle<dim, spacedim>>::iterator position_hint
          = particles.end();
        for(unsigned int i = 0; i < *n_particles_in_cell_ptr; ++i)
          {
            // Use std::multimap::emplace_hint to speed up insertion of
            // particles. This is a C++11 function, but not all compilers
            // that report a -std=c++11 (like gcc 4.6) implement it, so
            // require C++14 instead.
#  ifdef DEAL_II_WITH_CXX14
            position_hint = particles.emplace_hint(
              position_hint,
              std::make_pair(cell->level(), cell->index()),
              Particle<dim, spacedim>(pdata, property_pool.get()));
#  else
            position_hint = particles.insert(
              position_hint,
              std::make_pair(
                std::make_pair(cell->level(), cell->index()),
                Particle<dim, spacedim>(pdata, property_pool.get())));
#  endif
            const Point<dim> p_unit = mapping->transform_real_to_unit_cell(
              cell, position_hint->second.get_location());
            position_hint->second.set_reference_location(p_unit);
            ++position_hint;
          }
      }
    else if(status
            == parallel::distributed::Triangulation<dim, spacedim>::CELL_REFINE)
      {
        std::vector<typename std::multimap<internal::LevelInd,
                                           Particle<dim, spacedim>>::iterator>
          position_hints(GeometryInfo<dim>::max_children_per_cell);
        for(unsigned int child_index = 0;
            child_index < GeometryInfo<dim>::max_children_per_cell;
            ++child_index)
          {
            const typename Triangulation<dim, spacedim>::cell_iterator child
              = cell->child(child_index);
            position_hints[child_index] = particles.upper_bound(
              std::make_pair(child->level(), child->index()));
          }

        for(unsigned int i = 0; i < *n_particles_in_cell_ptr; ++i)
          {
            Particle<dim, spacedim> p(pdata, property_pool.get());

            for(unsigned int child_index = 0;
                child_index < GeometryInfo<dim>::max_children_per_cell;
                ++child_index)
              {
                const typename Triangulation<dim, spacedim>::cell_iterator child
                  = cell->child(child_index);

                try
                  {
                    const Point<dim> p_unit
                      = mapping->transform_real_to_unit_cell(child,
                                                             p.get_location());
                    if(GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                      {
                        p.set_reference_location(p_unit);
                        // Use std::multimap::emplace_hint to speed up insertion of
                        // particles. This is a C++11 function, but not all compilers
                        // that report a -std=c++11 (like gcc 4.6) implement it, so
                        // require C++14 instead.
#  ifdef DEAL_II_WITH_CXX14
                        position_hints[child_index] = particles.emplace_hint(
                          position_hints[child_index],
                          std::make_pair(child->level(), child->index()),
                          std::move(p));
#  else
                        position_hints[child_index] = particles.insert(
                          position_hints[child_index],
                          std::make_pair(
                            std::make_pair(child->level(), child->index()), p));
#  endif
                        ++position_hints[child_index];
                        break;
                      }
                  }
                catch(typename Mapping<dim>::ExcTransformationFailed&)
                  {}
              }
          }
      }
  }
} // namespace Particles

#endif

DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST
#  include "particle_handler.inst"
#endif

DEAL_II_NAMESPACE_CLOSE
