// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_objects_orientations_h
#define dealii_tria_objects_orientations_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/reference_cell.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * Class storing orientation information for various objects in a
     * Triangulation.
     *
     * @ingroup reordering
     */
    class TriaObjectsOrientations
    {
    public:
      /**
       * Constructor.
       */
      TriaObjectsOrientations();

      /**
       * Constructor. Sets up objects in the default orientation (orientation
       * = `true`).
       */
      TriaObjectsOrientations(const unsigned int n_objects);

      /**
       * Return number of geometric objects stored by this class.
       */
      unsigned int
      n_objects() const;

      /**
       * Reset the object to a default state.
       */
      void
      reinit(const unsigned int n_objects);

      /**
       * Change the number of stored objects. New objects are constructed in
       * the default orientation (true, false, false).
       */
      void
      resize(const unsigned int n_objects);

      /**
       * Return the size of objects of this kind.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Get the combined orientation of the object, as described in the class
       * documentation.
       */
      types::geometric_orientation
      get_combined_orientation(const unsigned int object) const;

      /**
       * Get the orientation bit of the object.
       */
      bool
      get_orientation(const unsigned int object) const;

      /**
       * Get the rotation bit of the object.
       */
      bool
      get_rotation(const unsigned int object) const;

      /**
       * Get the flip bit of the object.
       */
      bool
      get_flip(const unsigned int object) const;

      /**
       * Set the combined orientation of the object, as described in the class
       * documentation.
       */
      void
      set_combined_orientation(const unsigned int                 object,
                               const types::geometric_orientation value);

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);

    private:
      /**
       * Number of objects.
       */
      unsigned int n_stored_objects;

      /**
       * Orientations.
       */
      std::vector<types::geometric_orientation> object_orientations;
    };

    //----------------------------------------------------------------------//

    inline TriaObjectsOrientations::TriaObjectsOrientations()
    {
      reinit(0);
    }



    inline TriaObjectsOrientations::TriaObjectsOrientations(
      const unsigned int n_objects)
    {
      reinit(n_objects);
    }



    inline void
    TriaObjectsOrientations::reinit(const unsigned int n_objects)
    {
      n_stored_objects = n_objects;
      object_orientations.assign(n_objects,
                                 numbers::default_geometric_orientation);
    }



    inline void
    TriaObjectsOrientations::resize(const unsigned int n_objects)
    {
      object_orientations.resize(n_objects,
                                 numbers::default_geometric_orientation);
      n_stored_objects = n_objects;
    }



    inline std::size_t
    TriaObjectsOrientations::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(n_stored_objects) +
             MemoryConsumption::memory_consumption(object_orientations);
    }



    inline unsigned int
    TriaObjectsOrientations::n_objects() const
    {
      return n_stored_objects;
    }



    inline types::geometric_orientation
    TriaObjectsOrientations::get_combined_orientation(
      const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return object_orientations[object];
    }



    inline bool
    TriaObjectsOrientations::get_orientation(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return !Utilities::get_bit(object_orientations[object], 0);
    }



    inline bool
    TriaObjectsOrientations::get_rotation(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return Utilities::get_bit(object_orientations[object], 1);
    }



    inline bool
    TriaObjectsOrientations::get_flip(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return Utilities::get_bit(object_orientations[object], 2);
    }



    inline void
    TriaObjectsOrientations::set_combined_orientation(
      const unsigned int                 object,
      const types::geometric_orientation value)
    {
      AssertIndexRange(object, n_stored_objects);
      object_orientations[object] = value;
    }



    template <class Archive>
    void
    TriaObjectsOrientations::serialize(Archive &ar, const unsigned int)
    {
      ar &n_stored_objects &object_orientations;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
