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

#ifndef dealii_particles_particle_handler_h
#define dealii_particles_particle_handler_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/function.h>
#include <deal.II/base/observer_pointer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/partitioner.h>
#include <deal.II/particles/property_pool.h>

#include <boost/range/iterator_range.hpp>
#include <boost/serialization/map.hpp>
#include <boost/signals2.hpp>


DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * This class manages the storage and handling of particles. It provides
   * the data structures necessary to store particles efficiently, accessor
   * functions to iterate over particles and find particles, and algorithms
   * to distribute particles in parallel domains. Note that the class
   * is designed in a similar way as the triangulation class. In particular,
   * we call particles in the domain of the local process local particles,
   * and particles that belong to neighbor processes and live in the ghost cells
   * around the locally owned domain "ghost particles". The class also includes
   * functionality that is similar to the DoFHandler class (it knows which
   * particles live on which cells) and the SolutionTransfer class (it knows
   * how to transfer particles between cells and subdomains).
   *
   * For examples on how to use this class to track particles, store properties
   * on particles, and let the properties on the particles influence the
   * finite-element solution see step-19, step-68, step-70, and step-83.
   */
  template <int dim, int spacedim = dim>
  class ParticleHandler : public EnableObserverPointer
  {
  public:
    /**
     * A type that can be used to iterate over all particles in the domain.
     */
    using particle_iterator = ParticleIterator<dim, spacedim>;

    /**
     * A type that represents a range of particles.
     */
    using particle_iterator_range = boost::iterator_range<particle_iterator>;

    /**
     * A type for the storage container for particles.
     */
    using particle_container =
      typename ParticleAccessor<dim, spacedim>::particle_container;

    /**
     * Default constructor.
     */
    ParticleHandler();

    /**
     * Constructor that initializes the particle handler with
     * a given triangulation and mapping. Since particles are stored with
     * respect to their surrounding cells this information is necessary to
     * correctly organize the particle collection.
     * This constructor is equivalent to calling the default constructor and
     * the initialize function.
     */
    ParticleHandler(const Triangulation<dim, spacedim> &tria,
                    const Mapping<dim, spacedim>       &mapping,
                    const unsigned int                  n_properties = 0);

    /**
     * Destructor.
     */
    virtual ~ParticleHandler();

    /**
     * Initialize the particle handler. This function clears the
     * internal data structures, and sets the triangulation and the
     * mapping to be used.
     */
    void
    initialize(const Triangulation<dim, spacedim> &tria,
               const Mapping<dim, spacedim>       &mapping,
               const unsigned int                  n_properties = 0);

    /**
     * Copy the state of particle handler @p particle_handler into the
     * current object. This will copy
     * all particles and properties and leave this object
     * as an identical copy of @p particle_handler. Existing
     * particles in this object are deleted. Be aware that this
     * does not copy functions that are connected to the signals of
     * @p particle_handler, nor does it connect the current object's member
     * functions to triangulation signals, which must be done by the caller
     * if necessary, that is if the @p particle_handler had
     * connected functions.
     *
     * This function is expensive as it has to duplicate all data
     * in @p particle_handler, and insert it into this object,
     * which may be a significant amount of data. However, it can
     * be useful to save the state of a particle
     * collection at a certain point in time and reset this
     * state later under certain conditions, for example if
     * a timestep has to be undone and repeated.
     */
    void
    copy_from(const ParticleHandler<dim, spacedim> &particle_handler);

    /**
     * Clear all particle related data but keep the handler initialized.
     */
    void
    clear();

    /**
     * Only clear particle data, but keep cache information about number
     * of particles. This is useful during reorganization of particle data
     * between processes.
     */
    void
    clear_particles();

    /**
     * This function can be used to preemptively reserve memory for particle
     * data. Calling this function before inserting particles will reduce
     * memory allocations and therefore increase the performance. Calling
     * this function is optional; if memory is not already allocated it will
     * be allocated automatically during the insertion. It is recommended to
     * use this function if you know the number of particles that will be
     * inserted, but cannot use one of the collective particle insertion
     * functions.
     *
     * @param n_particles Number of particles to reserve memory for. Note that
     * this is the total number of particles to be stored, not the number of
     * particles to be newly inserted.
     */
    void
    reserve(const std::size_t n_particles);

    /**
     * Update all internally cached numbers. Note that all functions that
     * modify internal data structures and act on multiple particles will
     * call this function automatically (e.g. insert_particles), while
     * functions that act on single particles will not call this function
     * (e.g. insert_particle). This is done because the update is
     * expensive compared to single operations.
     */
    void
    update_cached_numbers();

    /**
     * Return an iterator to the first locally owned particle.
     */
    particle_iterator
    begin() const;

    /**
     * Return an iterator to the first locally owned particle.
     */
    particle_iterator
    begin();

    /**
     * Return an iterator past the end of the particles.
     */
    particle_iterator
    end() const;

    /**
     * Return an iterator past the end of the particles.
     */
    particle_iterator
    end();

    /**
     * Return an iterator to the first ghost particle.
     */
    particle_iterator
    begin_ghost() const;

    /**
     * Return an iterator to the first ghost particle.
     */
    particle_iterator
    begin_ghost();

    /**
     * Return an iterator past the end of the ghost particles.
     */
    particle_iterator
    end_ghost() const;

    /**
     * Return an iterator past the end of the ghost particles.
     */
    particle_iterator
    end_ghost();

    /**
     * Return the number of particles that live on the given cell.
     */
    types::particle_index
    n_particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
      const;

    /**
     * Return a pair of particle iterators that mark the begin and end of
     * the particles in a particular cell. The last iterator is the first
     * particle that is no longer in the cell.
     *
     * The number of elements in the returned range equals what the
     * n_particles_in_cell() function returns.
     */
    particle_iterator_range
    particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

    /**
     * Return a pair of particle iterators that mark the begin and end of
     * the particles in a particular cell. The last iterator is the first
     * particle that is no longer in the cell.
     *
     * The number of elements in the returned range equals what the
     * n_particles_in_cell() function returns.
     */
    particle_iterator_range
    particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
      const;

    /**
     * Remove a particle pointed to by the iterator. Note that @p particle
     * and all iterators that point to other particles in the same cell
     * as @p particle will be invalidated during this call.
     */
    void
    remove_particle(const particle_iterator &particle);

    /**
     * Remove a vector of particles indicated by the particle iterators.
     * The iterators and all other particle iterators are invalidated
     * during the function call.
     */
    void
    remove_particles(const std::vector<particle_iterator> &particles);

    /**
     * Insert a particle into the collection of particles. Return an iterator
     * to the new position of the particle. This function involves a copy of
     * the particle and its properties. Note that this function is of $O(N \log
     * N)$ complexity for $N$ particles.
     */
    particle_iterator
    insert_particle(
      const Particle<dim, spacedim> &particle,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

    /**
     * Insert a particle into the collection of particles given all the
     * properties necessary for creating a particle. This function is used to
     * efficiently generate particles without the detour through a Particle
     * object.
     *
     * @param[in] position Initial position of the particle in real space.
     * @param[in] reference_position Initial position of the particle
     * in the coordinate system of the reference cell.
     * @param[in] particle_index Globally unique identifier for this particle.
     * @param[in] cell The cell in which the particle is located.
     * @param[in] properties An optional ArrayView that describes the
     * particle properties. If given this has to be of size
     * n_properties_per_particle().
     */
    particle_iterator
    insert_particle(
      const Point<spacedim>      &position,
      const Point<dim>           &reference_position,
      const types::particle_index particle_index,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
      const ArrayView<const double> &properties = {});

    /**
     * Insert a number of particles into the collection of particles.
     * This function involves a copy of the particles and their properties.
     * Note that this function is of O(n_existing_particles + n_particles)
     * complexity.
     */
    void
    insert_particles(
      const std::multimap<
        typename Triangulation<dim, spacedim>::active_cell_iterator,
        Particle<dim, spacedim>> &particles);

    /**
     * Create and insert a number of particles into the collection of particles.
     * This function takes a list of positions and creates a set of particles
     * at these positions, which are then added to the local particle
     * collection. Note that this function currently uses
     * GridTools::compute_point_locations(), which assumes all positions are
     * within the local part of the triangulation. If one of them is not in the
     * local domain this function will throw an exception.
     */
    void
    insert_particles(const std::vector<Point<spacedim>> &positions);

    /**
     * Create and insert a number of particles into the collection of particles.
     * This function takes a list of positions and creates a set of particles
     * at these positions, which are then distributed and added to the local
     * particle collection of a processor. Note that this function uses
     * GridTools::distributed_compute_point_locations(). Consequently, it can
     * require intense communications between the processors. This function
     * is used in step-70.
     *
     * This function figures out what mpi process owns the points that do not
     * fall within the locally owned part of the triangulation, it sends
     * to that process the points passed to this function on this process,
     * and receives the points that fall within the locally owned cells of
     * the triangulation from whoever received them as input.
     *
     * In order to keep track of what mpi process received what points, a map
     * from mpi process to IndexSet is returned by the function. This IndexSet
     * contains the local indices of the points that were passed to this
     * function on the calling mpi process, and that falls within the part of
     * triangulation owned by this mpi process.
     *
     * The ids of the resulting particles are assigned from the optional
     * argument @p ids. If the vector of @p ids is empty, then the ids are
     * computed automatically from the get_next_free_particle_index() onward.
     * For example, if the method get_next_free_particle_index() returns n0,
     * calling this function with two MPI processes each adding n1 and n2
     * particles will result in the n1 particles added by process zero having
     * ids equal to `[n0,n0+n1)`, and the n2 particles added by process one
     * having ids `[n0+n1, n0+n1+n2)`.
     *
     * @param[in] positions A vector of points that do not need to be on the
     * local processor, but have to be in the triangulation that is associated
     * with this ParticleHandler object.
     *
     * @param[in] global_bounding_boxes A vector of vectors of bounding boxes.
     * The bounding boxes `global_bboxes[rk]` describe which part of the mesh is
     * locally owned by the mpi process with rank `rk`. The local description
     * can be obtained from GridTools::compute_mesh_predicate_bounding_box(),
     * and the global one can be obtained by passing the local ones to
     * Utilities::MPI::all_gather().
     *
     * @param[in] properties (Optional) A vector of vector of properties
     * associated with each local point. The size of the vector should be either
     * zero (no properties will be transferred nor attached to the generated
     * particles) or it should be a vector of `positions.size()` vectors of size
     * `n_properties_per_particle()`. Notice that this function call will
     * transfer the properties from the local mpi process to the final mpi
     * process that will own each of the particles, and it may therefore be
     * communication intensive.
     *
     * @param[in] ids (Optional) A vector of ids to associate to each particle.
     * If the vector is empty, the ids are assigned as a continuous range
     * from the first available index, as documented above. If the vector is not
     * empty, then its size must match the size of the @p positions vector.
     *
     * @return A map from owner to IndexSet, that contains the local indices
     * of the points that were passed to this function on the calling mpi
     * process, and that falls within the part of triangulation owned by this
     * mpi process.
     */
    std::map<unsigned int, IndexSet>
    insert_global_particles(
      const std::vector<Point<spacedim>> &positions,
      const std::vector<std::vector<BoundingBox<spacedim>>>
                                               &global_bounding_boxes,
      const std::vector<std::vector<double>>   &properties = {},
      const std::vector<types::particle_index> &ids        = {});

    /**
     * Insert a number of particles into the collection of particles. This
     * function takes a list of particles for which we don't know the associated
     * cell iterator, and distributes them to the correct local particle
     * collection of a processor, by unpacking the locations, figuring out where
     * to send the particles by calling
     * GridTools::distributed_compute_point_locations(), and sending the
     * particles to the corresponding process.
     *
     * In order to keep track of what mpi process received what particles, a map
     * from mpi process to IndexSet is returned by the function. This IndexSet
     * contains the local indices of the particles that were passed to this
     * function on the calling mpi process, and that falls within the part of
     * the triangulation owned by this mpi process.
     *
     * @param[in] particles A vector of particles that do not need to be on the
     * local processor.
     *
     * @param[in] global_bounding_boxes A vector of vectors of bounding boxes.
     * The bounding boxes `global_bboxes[rk]` describe which part of the mesh is
     * locally owned by the mpi process with rank `rk`. The local description
     * can be obtained from GridTools::compute_mesh_predicate_bounding_box(),
     * and the global one can be obtained by passing the local ones to
     * Utilities::MPI::all_gather().
     *
     * @return A map from owner to IndexSet, that contains the local indices
     * of the points that were passed to this function on the calling mpi
     * process, and that falls within the part of triangulation owned by this
     * mpi process.
     */
    std::map<unsigned int, IndexSet>
    insert_global_particles(
      const std::vector<Particle<dim, spacedim>> &particles,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &global_bounding_boxes);

    /**
     * Set the position of the particles by using the values contained in the
     * vector @p input_vector.
     *
     * @tparam VectorType Any of the parallel distributed vectors supported by
     * the library.
     *
     * The vector @p input_vector should have read access to the indices
     * created by extracting the locally relevant ids with
     * locally_owned_particle_ids(), and taking its tensor
     * product with the index set representing the range `[0, spacedim)`, i.e.:
     * @code
     * IndexSet ids = particle_handler.locally_owned_particle_ids().
     *  tensor_product(complete_index_set(spacedim));
     * @endcode
     *
     * The position of the particle with global index `id` is read from
     * spacedim consecutive entries starting from
     * `input_vector[id*spacedim]`.
     *
     * Notice that it is not necessary that the @p input_vector *owns* those
     * indices, however it has to have read access to them (i.e., it can be a
     * distributed vector with ghost entries).
     *
     * If the argument @p displace_particles is set to false, then the new
     * position taken from the values contained in
     * @p input_vector, replacing the previously stored particle position.
     * By default, the particles are displaced by the amount contained in the
     * @p input_vector, i.e., the contents of the vector are considered
     * *offsets* that are added to the previous position.
     *
     * After setting the new position, this function calls internally the method
     * sort_particles_into_subdomains_and_cells(). You should
     * make sure you satisfy the requirements of that function.
     *
     * @param[in] input_vector A parallel distributed vector containing
     * the displacement to apply to each particle, or their new absolute
     * position.
     *
     * @param[in] displace_particles Control if the @p input_vector should
     * be interpreted as a displacement vector, or a vector of absolute
     * positions.
     */
    template <typename VectorType>
    std::enable_if_t<
      std::is_convertible_v<VectorType *, Function<spacedim> *> == false>
    set_particle_positions(const VectorType &input_vector,
                           const bool        displace_particles = true);

    /**
     * Set the position of the particles within the particle handler using a
     * vector of points. The new set of point defined by the
     * vector has to be sufficiently close to the original one to ensure that
     * the sort_particles_into_subdomains_and_cells() function manages to find
     * the new cells in which the particles belong.
     *
     * Points are numbered in the same way they are traversed locally by the
     * ParticleHandler. A typical way to use this method, is to first call the
     * get_particle_positions() function, and then modify the resulting vector.
     *
     * @param [in] new_positions A vector of points of dimension
     * particle_handler.n_locally_owned_particles()
     *
     * @param [in] displace_particles When true, this function adds the value
     * of the vector of points to the
     * current position of the particle, thus displacing them by the
     * amount given by the function. When false, the position of the
     * particle is replaced by the value in the vector.
     */
    void
    set_particle_positions(const std::vector<Point<spacedim>> &new_positions,
                           const bool displace_particles = true);


    /**
     * Set the position of the particles within the particle handler using a
     * function with spacedim components. The new set of point defined by the
     * function has to be sufficiently close to the original one to ensure that
     * the sort_particles_into_subdomains_and_cells algorithm manages to find
     * the new cells in which the particles belong.
     *
     * The function is evaluated at the current location of the particles.
     *
     * @param [in] function A function that has n_components==spacedim that
     * describes either the displacement or the new position of the particles as
     * a function of the current location of the particle.
     *
     * @param [in] displace_particles When true, this function adds the results
     * of the function to the current position of the particle, thus displacing
     * them by the amount given by the function. When false, the position of the
     * particle is replaced by the value of the function.
     */
    void
    set_particle_positions(const Function<spacedim> &function,
                           const bool                displace_particles = true);

    /**
     * Read the position of the particles and store them into the distributed
     * vector @p output_vector. By default the
     * @p output_vector is overwritten by this operation, but you can add to
     * its entries by setting @p add_to_output_vector to `true`.
     *
     * @tparam VectorType Any of the parallel distributed vectors supported by
     * the library.
     *
     * This is the reverse operation of the set_particle_positions() function.
     * The position of the particle with global index `id` is written to
     * spacedim consecutive entries starting from
     * `output_vector[id*spacedim]`.
     *
     * Notice that, if you use a distributed vector type, it is not necessary
     * for the @p output_vector to own the entries corresponding to the indices
     * that will be written. However you should keep in mind that this requires
     * a global communication to distribute the entries above to their
     * respective owners.
     *
     * @param[in, out] output_vector A parallel distributed vector containing
     * the positions of the particles, or updated with the positions of the
     * particles.
     *
     * @param[in] add_to_output_vector Control if the function should set the
     * entries of the @p output_vector or if should add to them.
     */
    template <typename VectorType>
    void
    get_particle_positions(VectorType &output_vector,
                           const bool  add_to_output_vector = false);

    /**
     * Gather the position of the particles within the particle handler in
     * a vector of points. The order of the points is the same on would obtain
     * by iterating over all (local) particles, and querying their locations.
     *
     * @param [in,out] positions A vector preallocated at size
     * `particle_handler.n_locally_owned_articles` and whose points will become
     * the positions of the locally owned particles
     *
     * @param [in] add_to_output_vector When true, the value of the point of
     * the particles is added to the positions vector. When false,
     * the value of the points in the positions vector are replaced by the
     * position of the particles.
     */
    void
    get_particle_positions(std::vector<Point<spacedim>> &positions,
                           const bool add_to_output_vector = false);

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
    register_additional_store_load_functions(
      const std::function<std::size_t()> &size_callback,
      const std::function<void *(const particle_iterator &, void *)>
        &store_callback,
      const std::function<const void *(const particle_iterator &, const void *)>
        &load_callback);

    /**
     * Return the total number of particles that were managed by this class
     * the last time the update_cached_numbers() function was called.
     * The actual number of particles may have changed since then if
     * particles have been added or removed.
     *
     * @return Total number of particles in simulation.
     */
    types::particle_index
    n_global_particles() const;

    /**
     * Return the maximum number of particles per cell the last
     * time the update_cached_numbers() function was called.
     *
     * @return Maximum number of particles in one cell in simulation.
     */
    types::particle_index
    n_global_max_particles_per_cell() const;

    /**
     * Return the number of particles in the local part of the
     * triangulation.
     */
    types::particle_index
    n_locally_owned_particles() const;

    /**
     * Return the next free particle index in the global set
     * of particles the last
     * time the update_cached_numbers() function was called.
     */
    types::particle_index
    get_next_free_particle_index() const;

    /**
     * Extract an IndexSet with global dimensions equal to
     * get_next_free_particle_index(), containing the locally owned
     * particle indices.
     *
     * This function can be used to construct distributed vectors and matrices
     * to manipulate particles using linear algebra operations.
     *
     * Notice that it is the user's responsibility to guarantee that particle
     * indices are unique, and no check is performed to verify that this is the
     * case, nor that the union of all IndexSet objects on each mpi process is
     * complete.
     *
     * @return An IndexSet of size get_next_free_particle_index(), containing
     * n_locally_owned_particle() indices.
     */
    IndexSet
    locally_owned_particle_ids() const;

    /**
     * Return the number of properties each particle has.
     */
    unsigned int
    n_properties_per_particle() const;

    /**
     * Return one past the largest local index (in MPI-local index space)
     * returned by ParticleAccessor::get_local_index(). This number can be
     * larger than locally_owned_particle_ids().n_elements(), because local
     * indices are not necessarily forming a contiguous range and can contain
     * gaps in the range 0 to get_max_local_particle_index(). As a
     * consequence, the number is not updated upon calls to remove_particle()
     * or similar functions, and refreshed only as new particles get added or
     * in sort_particles_into_subdomains_and_cells().
     *
     * This function is appropriate for resizing vectors working with
     * ParticleAccessor::get_local_index().
     */
    types::particle_index
    get_max_local_particle_index() const;

    /**
     * Return a reference to the property pool that owns all particle
     * properties, and organizes them physically.
     */
    PropertyPool<dim, spacedim> &
    get_property_pool() const;

    /**
     * Find and update the cells containing each particle for all locally owned
     * particles. If particles moved out of the local subdomain
     * they will be sent to their new process and inserted there.
     * After this function call every particle is either on its current
     * process and in its current cell, or deleted (if it could not find
     * its new process or cell).
     *
     * The user may attach a function to the signal
     * Particles::ParticleHandler::Signals::particle_lost(). The signal is
     * triggered whenever a particle is deleted, and the connected functions
     * are called passing an iterator to the particle in question, and its last
     * known cell association.
     */
    void
    sort_particles_into_subdomains_and_cells();

    /**
     * Exchange all particles that live in cells that are ghost cells to
     * other processes. Clears and re-populates the ghost_neighbors
     * member variable.
     */
    void
    exchange_ghost_particles(const bool enable_ghost_cache = false);

    /**
     * Update all particles that live in cells that are ghost cells to
     * other processes. In this context, update means to update the
     * location and the properties of the ghost particles assuming that
     * the ghost particles have not changed cells. Consequently, this will
     * not update the reference location of the particles.
     */
    void
    update_ghost_particles();

    /**
     * This function prepares the particle handler for a coarsening and
     * refinement cycle, by storing the necessary information to transfer
     * particles to their new cells. The implementation depends on the
     * triangulation type that is connected to the particle handler and differs
     * between shared and distributed triangulations. This function should be
     * used like the corresponding function with the same name in the
     * SolutionTransfer() class.
     */
    void
    prepare_for_coarsening_and_refinement();

    /**
     * This function unpacks the particle data after a coarsening and
     * refinement cycle, by reading the necessary information to transfer
     * particles to their new cells. This function should be
     * used like the SolutionTransfer::interpolate() function after mesh
     * refinement has finished. Note that this function requires a working
     * mapping, i.e. if you use a mapping class that requires setup after
     * a mesh refinement (e.g. MappingQCache(), or MappingEulerian()), the
     * mapping has to be ready before you can call this function.
     *
     * @note It is important to note, that if you do not call this function
     * after a refinement/coarsening operation (and a repartitioning operation
     * in a distributed triangulation) and the particle handler contained
     * particles before, the particle handler will not be usable any more. Not
     * calling this function after refinement would be similar to trying to
     * access a solution vector after mesh refinement without first using a
     * SolutionTransfer class to transfer the vector to the new mesh.
     */
    void
    unpack_after_coarsening_and_refinement();

    /**
     * This function prepares the particle handler for serialization. This needs
     * to be done before calling Triangulation::save(). This function should be
     * used like the corresponding function with the same name in the
     * SolutionTransfer() class.
     *
     * @note It is important to note, that if you do not call this function
     * before a serialization operation in a distributed triangulation, the
     * particles will not be serialized, and therefore will disappear after
     * deserialization.
     */
    void
    prepare_for_serialization();

    /**
     * Execute the deserialization of the particle data. This needs to be
     * done after calling Triangulation::load(). The data must have been stored
     * before the serialization of the triangulation using the
     * prepare_for_serialization() function. This function should be used like
     * the corresponding function with the same name in the SolutionTransfer
     * class.
     */
    void
    deserialize();

    /**
     * Serialize the contents of this class using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     *
     * This function is used in step-83.
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * A structure that has boost::signal objects for a number of actions that a
     * particle handler can do to itself. How signals can be used in
     * applications is explained in the "Getting notice when a triangulation
     * changes" section in the Triangulation class with more information and
     * examples. In short these signals allow the particle handler to notify
     * applications about certain events inside the particle handler, e.g. when
     * a particle is lost.
     *
     * For documentation on signals, see
     * http://www.boost.org/doc/libs/release/libs/signals2 .
     */
    struct Signals
    {
      /**
       * This signal is triggered whenever the
       * ParticleHandler::sort_particles_into_subdomains_and_cells() function
       * encounters a particle that can not be associated with a cell. This can
       * happen if the particle leaves the domain of the triangulation, or if it
       * leaves the locally known domain in a parallel triangulation (including
       * the ghost cells for a parallel::distributed::triangulation).
       *
       * The connected function receives an iterator to the particle in
       * question, and its last known cell association.
       *
       * This signal is used in step-19.
       */
      boost::signals2::signal<void(
        const typename Particles::ParticleIterator<dim, spacedim> &particle,
        const typename Triangulation<dim, spacedim>::active_cell_iterator
          &cell)>
        particle_lost;
    };

    /**
     * Signals for the events that a particle handler can notify the
     * calling application about.
     */
    mutable Signals signals;

  private:
    /**
     * Insert a particle into the collection of particles from a raw
     * data pointer. This function is used for shipping particles
     * during serialization and refinement and not intended for use outside
     * of this class. Return an iterator to the new position of the particle.
     */
    particle_iterator
    insert_particle(
      const void                                                       *&data,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

    /**
     * Perform the local insertion operation into the particle container. This
     * function is used in the higher-level functions inserting particles.
     */
    particle_iterator
    insert_particle(
      const typename PropertyPool<dim, spacedim>::Handle
        tria_attached_data_index,
      const typename Triangulation<dim, spacedim>::cell_iterator &cell);

    /**
     * Delete all entries in the particles container, and then set the three
     * anchor entries used to distinguish between owned and ghost cells: We
     * add one item to the front, one between and one as the last element in
     * the particle container to be able to iterate across particles without
     * `if` statements, solely relying on ParticleAccessor::operator== to
     * terminate operations, and using the `cell_iterator` inside the particle
     * container to check for valid states.
     */
    void
    reset_particle_container(particle_container &particles);

    /**
     * Address of the triangulation to work on.
     */
    ObserverPointer<const Triangulation<dim, spacedim>,
                    ParticleHandler<dim, spacedim>>
      triangulation;

    /**
     * Address of the mapping to work on.
     */
    ObserverPointer<const Mapping<dim, spacedim>,
                    ParticleHandler<dim, spacedim>>
      mapping;

    /**
     * This object owns and organizes the memory for all particle
     * properties. Since particles reference the property pool, the
     * latter has to be destroyed *after* the particles are destroyed.
     * This is achieved by making sure the `property_pool` member variable
     * precedes the declaration of the `particles`
     * member variable.
     */
    std::unique_ptr<PropertyPool<dim, spacedim>> property_pool;

    /**
     * Set of particles currently living in the locally owned or ghost cells.
     */
    particle_container particles;

    /**
     * Iterator to the end of the list elements of particle_container which
     * belong to locally owned elements. Made const to avoid accidental
     * modification.
     */
    const typename particle_container::iterator owned_particles_end;

    /**
     * List from the active cells on the present MPI process to positions in
     * either `owned_particles` or `ghost_particles` for fast $\mathcal O(1)$
     * access to the particles of a cell.
     */
    std::vector<typename particle_container::iterator> cells_to_particle_cache;

    /**
     * This variable stores how many particles are stored globally. It is
     * calculated by update_cached_numbers().
     */
    types::particle_index global_number_of_particles;

    /**
     * This variable stores how many particles are locally owned. It is
     * calculated by update_cached_numbers().
     */
    types::particle_index number_of_locally_owned_particles;

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
    std::function<std::size_t()> size_callback;

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
    std::function<void *(const particle_iterator &, void *)> store_callback;

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
    std::function<const void *(const particle_iterator &, const void *)>
      load_callback;

    /**
     * This variable is set by the register_data_attach()
     * function and used by the notify_ready_to_unpack() function
     * to check where the particle data was registered in the corresponding
     * triangulation object.
     */
    unsigned int tria_attached_data_index;

    /**
     * Tolerance to be used for GeometryInfo<dim>::is_inside_cell().
     */
    const double tolerance_inside_cell = 1e-12;

    /**
     * The GridTools::Cache is used to store the information about the
     * vertex_to_cells set and the vertex_to_cell_centers vectors to prevent
     * recomputing them every time we sort_into_subdomain_and_cells().
     * This cache is automatically updated when the triangulation has
     * changed. This cache is stored within a unique pointer because the
     * particle handler has a constructor that enables it to be constructed
     * without a triangulation. The cache does not have such a constructor.
     */
    std::unique_ptr<GridTools::Cache<dim, spacedim>> triangulation_cache;

#ifdef DEAL_II_WITH_MPI
    /**
     * Transfer particles that have crossed subdomain boundaries to other
     * processors.
     * All received particles and their new cells will be appended to the
     * class variable `particles` at the right slot.
     *
     * @param [in] particles_to_send All particles that should be sent and
     * their new subdomain_ids are in this map.
     *
     * @param [in] new_cells_for_particles Optional vector of cell
     * iterators with the same structure as @p particles_to_send. If this
     * parameter is given it should contain the cell iterator for every
     * particle to be send in which the particle belongs. This parameter
     * is necessary if the cell information of the particle iterator is
     * outdated (e.g. after particle movement).
     *
     * @param [in] enable_cache Optional bool that enables updating
     * the ghost particles without rebuilding them from scratch by
     * building a cache of type GhostParticlePartitioner, which
     * stores the necessary information to update the ghost particles.
     * Once this cache is built, the ghost particles can be updated
     * by a call to send_recv_particles_properties_and_location().
     */
    void
    send_recv_particles(
      const std::map<types::subdomain_id, std::vector<particle_iterator>>
        &particles_to_send,
      const std::map<
        types::subdomain_id,
        std::vector<
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
        &new_cells_for_particles = std::map<
          types::subdomain_id,
          std::vector<
            typename Triangulation<dim, spacedim>::active_cell_iterator>>(),
      const bool enable_cache = false);

    /**
     * Transfer ghost particles' position and properties assuming that the
     * particles have not changed cells. This routine uses the
     * GhostParticlePartitioner as a caching structure to know which particles
     * are ghost to other processes, and where they need to be sent.  It
     * inherently assumes that particles cannot have changed cell, and writes
     * the result back to the `particles` member variable.
     *
     * @param [in] particles_to_send All particles for which information
     * should be sent and their new subdomain_ids are in this map.
     */
    void
    send_recv_particles_properties_and_location(
      const std::map<types::subdomain_id, std::vector<particle_iterator>>
        &particles_to_send);

#endif

    /**
     * Cache structure used to store the elements which are required to
     * exchange the particle information (location and properties) across
     * processors in order to update the ghost particles. This structure
     * is only used to update the ghost particles.
     */
    internal::GhostParticlePartitioner<dim, spacedim> ghost_particles_cache;

    /**
     * Connect the particle handler to the relevant triangulation signals to
     * appropriately react to changes in the underlying triangulation.
     */
    void
    connect_to_triangulation_signals();

    /**
     * A list of connections with which this object connects to the
     * triangulation to get information about when the triangulation changes.
     */
    std::vector<boost::signals2::connection> tria_listeners;

    /**
     * Function that gets called by the triangulation signals if
     * the structure of the mesh has changed. This function is
     * responsible for resizing the particle container if no
     * particles are stored, i.e. if the usual call to
     * prepare_for_..., and transfer_particles_after_... functions
     * is not necessary.
     */
    void
    post_mesh_change_action();

    /**
     * Function that should be called before every refinement
     * and when writing checkpoints. This function is used to
     * register pack_callback() with the triangulation.
     */
    void
    register_data_attach();

    /**
     * Callback function that should be called after every refinement
     * and after resuming from a checkpoint.  This function is used to
     * call the notify_ready_to_unpack() function of the triangulation
     * and hand over the unpack_callback() function, which will
     * unpack the particle data for each cell.
     */
    void
    notify_ready_to_unpack(const bool serialization);

    /**
     * Called by listener functions from Triangulation for every cell
     * before a refinement step. All particles have to be attached to their
     * cell to be sent around to the new cell and owning process.
     */
    std::vector<char>
    pack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellStatus                                            status) const;

    /**
     * Called by listener functions after a refinement step for each cell
     * to unpack the particle data and transfer it to the local container.
     */
    void
    unpack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellStatus                                            status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range);

    /**
     * Internal function returning an iterator to the begin of the container
     * for owned particles.
     */
    typename particle_container::iterator
    particle_container_owned_begin() const;

    /**
     * Internal function returning an iterator to the end of the container for
     * owned particles.
     */
    typename particle_container::iterator
    particle_container_owned_end() const;

    /**
     * Internal function returning an iterator to the begin of the container
     * for ghost particles.
     */
    typename particle_container::iterator
    particle_container_ghost_begin() const;

    /**
     * Internal function returning an iterator to the end of the container for
     * ghost particles.
     */
    typename particle_container::iterator
    particle_container_ghost_end() const;
  };



  /* ---------------------- inline and template functions ------------------
   */

  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->begin();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin()
  {
    return particle_iterator(particle_container_owned_begin(),
                             *property_pool,
                             0);
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->end();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end()
  {
    return particle_iterator(particle_container_owned_end(), *property_pool, 0);
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->begin_ghost();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost()
  {
    return particle_iterator(particle_container_ghost_begin(),
                             *property_pool,
                             0);
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->end_ghost();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost()
  {
    return particle_iterator(particle_container_ghost_end(), *property_pool, 0);
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_container::iterator
  ParticleHandler<dim, spacedim>::particle_container_owned_begin() const
  {
    // We should always have at least the three anchor entries in the list of
    // particles
    Assert(!particles.empty(), ExcInternalError());
    typename particle_container::iterator begin =
      const_cast<particle_container &>(particles).begin();
    return ++begin;
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_container::iterator
  ParticleHandler<dim, spacedim>::particle_container_owned_end() const
  {
    Assert(!particles.empty(), ExcInternalError());
    return owned_particles_end;
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_container::iterator
  ParticleHandler<dim, spacedim>::particle_container_ghost_begin() const
  {
    Assert(!particles.empty(), ExcInternalError());
    typename particle_container::iterator begin = owned_particles_end;
    return ++begin;
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_container::iterator
  ParticleHandler<dim, spacedim>::particle_container_ghost_end() const
  {
    Assert(!particles.empty(), ExcInternalError());
    typename particle_container::iterator end =
      const_cast<particle_container &>(particles).end();
    return --end;
  }



  template <int dim, int spacedim>
  template <class Archive>
  inline void
  ParticleHandler<dim, spacedim>::serialize(Archive &ar, const unsigned int)
  {
    // Note that we do not serialize the particle data itself (i.e., the
    // 'particles' member variable). Instead we
    // use the serialization functionality of the triangulation class, because
    // this guarantees that data is immediately shipped to new processes if
    // the domain is distributed differently after resuming from a checkpoint.
    //
    // See step-83 for how to serialize ParticleHandler objects.
    ar &global_number_of_particles &global_max_particles_per_cell
      &next_free_particle_index;
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  inline std::enable_if_t<
    std::is_convertible_v<VectorType *, Function<spacedim> *> == false>
  ParticleHandler<dim, spacedim>::set_particle_positions(
    const VectorType &input_vector,
    const bool        displace_particles)
  {
    AssertDimension(input_vector.size(),
                    get_next_free_particle_index() * spacedim);
    for (auto &p : *this)
      {
        Point<spacedim> &position = p.get_location();
        const auto       id       = p.get_id();

        if (displace_particles)
          for (unsigned int i = 0; i < spacedim; ++i)
            position[i] += input_vector[id * spacedim + i];
        else
          for (unsigned int i = 0; i < spacedim; ++i)
            position[i] = input_vector[id * spacedim + i];
      }
    sort_particles_into_subdomains_and_cells();
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  inline void
  ParticleHandler<dim, spacedim>::get_particle_positions(
    VectorType &output_vector,
    const bool  add_to_output_vector)
  {
    AssertDimension(output_vector.size(),
                    get_next_free_particle_index() * spacedim);
    for (const auto &p : *this)
      {
        auto       point = p.get_location();
        const auto id    = p.get_id();
        if (add_to_output_vector)
          for (unsigned int i = 0; i < spacedim; ++i)
            output_vector[id * spacedim + i] += point[i];
        else
          for (unsigned int i = 0; i < spacedim; ++i)
            output_vector[id * spacedim + i] = point[i];
      }
    if (add_to_output_vector)
      output_vector.compress(VectorOperation::add);
    else
      output_vector.compress(VectorOperation::insert);
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
