// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_particles_particle_handler_h
#define dealii_particles_particle_handler_h

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/array_view.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/fe/mapping.h>

#include <boost/serialization/map.hpp>
#include <boost/range/iterator_range.hpp>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * This class manages the storage and handling of particles. It provides
   * the data structures necessary to store particles efficiently, accessor
   * functions to iterate over particles and find particles, and algorithms
   * to distribute particles in parallel domains.
   *
   * @ingroup Particle
   */
  template <int dim, int spacedim=dim>
  class ParticleHandler: public Subscriptor
  {
  public:
    /**
     * A type that can be used to iterate over all particles in the domain.
     */
    typedef ParticleIterator<dim,spacedim> particle_iterator;

    /**
     * A type that represents a range of particles.
     */
    typedef boost::iterator_range<particle_iterator> particle_iterator_range;

    /**
     * Default constructor.
     */
    ParticleHandler();

    /**
     * Constructor that initializes the particle handler with respect to
     * a given triangulation and MPI communicator. Pointers to the
     * triangulation and the communicator are stored inside of the particle
     *
     */
    ParticleHandler(const parallel::distributed::Triangulation<dim,spacedim> &tria,
                    const Mapping<dim,spacedim> &mapping,
                    const unsigned int n_properties = 0);

    /**
     * Destructor.
     */
    ~ParticleHandler();

    /**
     * Initialize the particle handler. This function does not clear the
     * internal data structures, it just sets the connections to the
     * MPI communicator and the triangulation.
     */
    void initialize(const parallel::distributed::Triangulation<dim,spacedim> &tria,
                    const Mapping<dim,spacedim> &mapping,
                    const unsigned int n_properties = 0);

    /**
     * Clear all particle related data.
     */
    void clear();

    /**
     * Only clear particle data, but keep cache information about number
     * of particles. This is useful during reorganization of particle data
     * between processes.
     */
    void clear_particles();

    /**
     * Updates all internally cached numbers. Note that all functions that
     * modify internal data structures and act on multiple particles will
     * call this function automatically (e.g. insert_particles), while
     * functions that act on single particles will not call this function
     * (e.g. insert_particle). This is done because the update is
     * expensive compared to single operations.
     */
    void update_cached_numbers();

    /**
     * Return an iterator to the first particle.
     */
    ParticleHandler<dim,spacedim>::particle_iterator begin() const;

    /**
     * Return an iterator to the first particle.
     */
    particle_iterator begin();

    /**
     * Return an iterator past the end of the particles.
     */
    particle_iterator end() const;

    /**
     * Return an iterator past the end of the particles.
     */
    particle_iterator end();

    /**
     * Return a pair of particle iterators that mark the begin and end of
     * the particles in a particular cell. The last iterator is the first
     * particle that is no longer in the cell.
     */
    particle_iterator_range
    particles_in_cell(const typename Triangulation<dim,spacedim>::active_cell_iterator &cell);


    /**
     * Return a pair of particle iterators that mark the begin and end of
     * the particles in a particular cell. The last iterator is the first
     * particle that is no longer in the cell.
     */
    particle_iterator_range
    particles_in_cell(const typename Triangulation<dim,spacedim>::active_cell_iterator &cell) const;

    /**
     * Remove a particle pointed to by the iterator.
     */
    void
    remove_particle(const particle_iterator &particle);

    /**
     * Insert a particle into the collection of particles. Return an iterator
     * to the new position of the particle. This function involves a copy of
     * the particle and its properties. Note that this function is of $O(N \log N)$
     * complexity for $N$ particles.
     */
    particle_iterator
    insert_particle(const Particle<dim,spacedim> &particle,
                    const typename Triangulation<dim,spacedim>::active_cell_iterator &cell);

    /**
     * Insert a number of particle into the collection of particles.
     * This function involves a copy of the particles and their properties.
     * Note that this function is of O(n_existing_particles + n_particles) complexity.
     */
    void
    insert_particles(const std::multimap<types::LevelInd, Particle<dim,spacedim> > &particles);

    /**
     * This function allows to register three additional functions that are
     * called every time a particle is transferred to another process
     * (i.e. during sorting into cells, during ghost particle transfer, or
     * during serialization of all particles).
     *
     * @param size_callback A function that is called when serializing
     * particle data. The function gets no arguments and is expected to
     * return the size of the additional data that is serialized per
     * particle. Note that this currently implies the data size has to be
     * the same for every particle.
     * @param store_callback A function that is called once per particle
     * when serializing particle data. Arguments to the function are a
     * particle iterator that identifies the current particle and a void
     * pointer that points to a data block of size size_callback() in which
     * the function can store additional data. The function is expected to
     * return a void pointer pointing to a position right after its data
     * block.
     * @param load_callback A function that is called once per particle
     * when deserializing particle data. Arguments to the function are a
     * particle iterator that identifies the current particle and a void
     * pointer that points to a data block of size size_callback() in which
     * additional data was stored by the store_callback function. The
     * function is expected to return a void pointer pointing to a position
     * right after its data block.
     */
    void
    register_additional_store_load_functions(const std::function<std::size_t ()> &size_callback,
                                             const std::function<void *(const particle_iterator &,
                                                 void *)> &store_callback,
                                             const std::function<const void *(const particle_iterator &,
                                                 const void *)> &load_callback);

    /**
     * Return the total number of particles that were managed by this class
     * the last time the update_n_global_particles() function was called.
     * The actual number of particles may have changed since then if
     * particles have been added or removed.
     *
     * @return Total number of particles in simulation.
     */
    types::particle_index n_global_particles() const;

    /**
     * Return the number of particles in the local part of the
     * triangulation.
     */
    types::particle_index n_locally_owned_particles() const;

    /**
     * Return the number of particles in the local part of the
     * triangulation.
     */
    types::particle_index get_next_free_particle_index() const;

    /**
     * Return the number of properties each particle has.
     */
    unsigned int n_properties_per_particle() const;

    /**
     * Return a reference to the property pool that owns all particle
     * properties, and organizes them physically.
     */
    PropertyPool &
    get_property_pool() const;

    /**
     * Return the number of particles in the given cell.
     */
    unsigned int
    n_particles_in_cell(const typename Triangulation<dim,spacedim>::active_cell_iterator &cell) const;

    /**
     * Finds the cells containing each particle for all locally owned
     * particles. If particles moved out of the local subdomain
     * they will be sent to their new process and inserted there.
     * After this function call every particle is either on its current
     * process and in its current cell, or deleted (if it could not find
     * its new process or cell).
     *
     * TODO: Extend this to allow keeping particles on other processes
     * around (with an invalid cell).
     */
    void
    sort_particles_into_subdomains_and_cells();


    /**
     * Exchanges all particles that live in cells that are ghost cells to
     * other processes. Clears and re-populates the ghost_neighbors
     * member variable.
     */
    void
    exchange_ghost_particles();

    /**
     * Callback function that should be called before every
     * refinement and when writing checkpoints.
     * Allows registering store_particles() in the triangulation.
     */
    void
    register_store_callback_function(const bool serialization);

    /**
     * Callback function that should be called after every
     * refinement and after resuming from a checkpoint.
     * Allows registering load_particles() in the triangulation.
     */
    void
    register_load_callback_function(const bool serialization);

    /**
     * Serialize the contents of this class.
     */
    template <class Archive>
    void serialize (Archive &ar, const unsigned int version);

  private:
    /**
     * Address of the triangulation to work on.
     */
    SmartPointer<const parallel::distributed::Triangulation<dim,spacedim>,ParticleHandler<dim,spacedim> > triangulation;

    /**
     * Address of the mapping to work on.
     */
    SmartPointer<const Mapping<dim,spacedim>,ParticleHandler<dim,spacedim> > mapping;

    /**
     * Set of particles currently in the local domain, organized by
     * the level/index of the cell they are in.
     */
    std::multimap<types::LevelInd, Particle<dim,spacedim> > particles;

    /**
     * Set of particles currently in the ghost cells of the local domain,
     * organized by the level/index of the cell they are in. These
     * particles are marked read-only.
     */
    std::multimap<types::LevelInd, Particle<dim,spacedim> > ghost_particles;

    /**
     * This variable stores how many particles are stored globally. It is
     * calculated by update_n_global_particles().
     */
    types::particle_index global_number_of_particles;

    /**
     * The maximum number of particles per cell in the global domain. This
     * variable is important to store and load particle data during
     * repartition and serialization of the solution. Note that the
     * variable is only updated when it is needed, e.g. after particle
     * movement, before/after mesh refinement, before creating a
     * checkpoint and after resuming from a checkpoint.
     */
    unsigned int global_max_particles_per_cell;

    /**
     * This variable stores the next free particle index that is available
     * globally in case new particles need to be generated.
     */
    types::particle_index next_free_particle_index;

    /**
     * This object owns and organizes the memory for all particle
     * properties.
     */
    std::unique_ptr<PropertyPool> property_pool;

    /**
     * A function that can be registered by calling
     * register_additional_store_load_functions. It is called when serializing
     * particle data. The function gets no arguments and is expected to
     * return the size of the additional data that is serialized per
     * particle. Note that this currently implies the data size has to be
     * the same for every particle, but it does not have to be the same for
     * every serialization process (e.g. a serialization during particle
     * movement might include temporary data, while a serialization after
     * movement was finished does not need to transfer this data).
     */
    std::function<std::size_t ()> size_callback;

    /**
     * A function that can be registered by calling
     * register_additional_store_load_functions. It is called once per
     * particle when serializing particle data. Arguments to the function
     * are a particle iterator that identifies the current particle and a void
     * pointer that points to a data block of size size_callback() in which
     * the function can store additional data. The function is expected to
     * return a void pointer pointing to a position right after its data
     * block.
     */
    std::function<void *(const particle_iterator &,
                         void *)> store_callback;

    /**
     * A function that is called once per particle
     * when deserializing particle data. Arguments to the function are a
     * particle iterator that identifies the current particle and a void
     * pointer that points to a data block of size size_callback() from
     * which the function can load additional data. This block was filled
     * by the store_callback function during serialization. This function
     * is expected to return a void pointer pointing to a position right
     * after its data block.
     */
    std::function<const void *(const particle_iterator &,
                               const void *)> load_callback;

    /**
     * This variable is set by the register_store_callback_function()
     * function and used by the register_load_callback_function() function
     * to check where the particle data was stored.
     */
    unsigned int data_offset;

    /**
     * Transfer particles that have crossed subdomain boundaries to other
     * processors.
     * All received particles and their new cells will be appended to the
     * @p received_particles vector.
     *
     * @param [in] particles_to_send All particles that should be sent and
     * their new subdomain_ids are in this map.
     *
     * @param [in,out] received_particles Vector that stores all received
     * particles. Note that it is not required nor checked that the list
     * is empty, received particles are simply attached to the end of
     * the vector.
     *
     * @param [in] new_cells_for_particles Optional vector of cell
     * iterators with the same structure as @p particles_to_send. If this
     * parameter is given it should contain the cell iterator for every
     * particle to be send in which the particle belongs. This parameter
     * is necessary if the cell information of the particle iterator is
     * outdated (e.g. after particle movement).
     */
    void
    send_recv_particles(const std::vector<std::vector<particle_iterator> >      &particles_to_send,
                        std::multimap<types::LevelInd,Particle <dim,spacedim> > &received_particles,
                        const std::vector<std::vector<typename Triangulation<dim,spacedim>::active_cell_iterator> > &new_cells_for_particles = std::vector<std::vector<typename Triangulation<dim,spacedim>::active_cell_iterator> > ());

    /**
     * Called by listener functions from Triangulation for every cell
     * before a refinement step. All particles have to be attached to their
     * cell to be sent around to the new processes.
     */
    void
    store_particles(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                    const typename Triangulation<dim,spacedim>::CellStatus status,
                    void *data) const;

    /**
     * Called by listener functions after a refinement step. The local map
     * of particles has to be read from the triangulation user_pointer.
     */
    void
    load_particles(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                   const typename Triangulation<dim,spacedim>::CellStatus status,
                   const void *data);

    /**
     * Get a map between subdomain id and a contiguous
     * number from 0 to n_neighbors, which is interpreted as the neighbor index.
     * In other words the returned map answers the question: Given a subdomain id, which
     * neighbor of the current processor's domain owns this subdomain?
     */
    std::map<types::subdomain_id, unsigned int>
    get_subdomain_id_to_neighbor_map() const;
  };

  /* -------------------------- inline and template functions ---------------------- */

  template <int dim, int spacedim>
  template <class Archive>
  void ParticleHandler<dim,spacedim>::serialize (Archive &ar, const unsigned int)
  {
    // Note that we do not serialize the particle data itself. Instead we
    // use the serialization functionality of the triangulation class, because
    // this guarantees that data is immediately shipped to new processes if
    // the domain is distributed differently after resuming from a checkpoint.
    ar //&particles
    &global_number_of_particles
    &global_max_particles_per_cell
    &next_free_particle_index;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

