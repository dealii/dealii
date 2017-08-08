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

#ifndef dealii__particles_property_pool_h
#define dealii__particles_property_pool_h

#include <deal.II/base/array_view.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * This class manages the memory space in which particles store their
   * properties. Because this is dynamic memory and every particle needs the
   * same amount it is more efficient to let this be handled by a central
   * manager that does not need to allocate/deallocate memory every time a
   * particle is constructed/destroyed.
   */
  class PropertyPool
  {
  public:
    /**
     * Typedef for the handle that is returned to the particles, and that
     * uniquely identifies the slot of memory that is reserved for this
     * particle.
     */
    typedef double *Handle;

    /**
     * Define a default (invalid) value for handles.
     */
    static const Handle invalid_handle;

    /**
     * Constructor. Stores the number of properties per reserved slot.
     */
    PropertyPool (const unsigned int n_properties_per_slot);

    /**
     * Returns a new handle that allows accessing the reserved block
     * of memory.
     */
    Handle allocate_properties_array ();

    /**
     * Mark the properties corresponding to the handle @p handle as
     * deleted. After calling this function calling get_properties() with
     * the freed handle causes undefined behavior.
     */
    void deallocate_properties_array (const Handle handle);

    /**
     * Return an ArrayView to the properties that correspond to the given
     * handle @p handle.
     */
    ArrayView<double> get_properties (const Handle handle);

    /**
     * Reserves the dynamic memory needed for storing the properties of
     * @p size particles.
     */
    void reserve(const std::size_t size);

  private:
    /**
     * The number of properties that are reserved per particle.
     */
    const unsigned int n_properties;
  };

}

DEAL_II_NAMESPACE_CLOSE

#endif
