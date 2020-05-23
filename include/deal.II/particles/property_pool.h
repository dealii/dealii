// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

#ifndef dealii_particles_property_pool_h
#define dealii_particles_property_pool_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * This class manages a memory space in which particles store their
   * properties. Because this is dynamic memory and often every particle
   * needs the same amount, it is more efficient to let this be handled by a
   * central manager that does not need to allocate/deallocate memory every
   * time a particle is constructed/destroyed.
   * The current implementation uses simple new/delete allocation for every
   * block. Additionally, the current implementation
   * assumes the same number of properties per particle, but of course the
   * PropertyType could contain a pointer to dynamically allocated memory
   * with varying sizes per particle (this memory would not be managed by this
   * class).
   * Because PropertyPool only returns handles it could be enhanced internally
   * (e.g. to allow for varying number of properties per handle) without
   * affecting its interface.
   */
  class PropertyPool
  {
  public:
    /**
     * Typedef for the handle that is returned to the particles, and that
     * uniquely identifies the slot of memory that is reserved for this
     * particle.
     */
    using Handle = double *;

    /**
     * Define a default (invalid) value for handles.
     */
    static const Handle invalid_handle;

    /**
     * Constructor. Stores the number of properties per reserved slot.
     */
    PropertyPool(const unsigned int n_properties_per_slot);

    /**
     * Return a new handle that allows accessing the reserved block
     * of memory. If the number of properties is zero this will return an
     * invalid handle.
     */
    Handle
    allocate_properties_array();

    /**
     * Mark the properties corresponding to the handle @p handle as
     * deleted. Calling this function more than once for the same
     * handle causes undefined behavior.
     */
    void
    deallocate_properties_array(const Handle handle);

    /**
     * Return an ArrayView to the properties that correspond to the given
     * handle @p handle.
     */
    ArrayView<double>
    get_properties(const Handle handle);

    /**
     * Reserve the dynamic memory needed for storing the properties of
     * @p size particles.
     */
    void
    reserve(const std::size_t size);

    /**
     * Return how many properties are stored per slot in the pool.
     */
    unsigned int
    n_properties_per_slot() const;

  private:
    /**
     * The number of properties that are reserved per particle.
     */
    const unsigned int n_properties;
  };


} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
