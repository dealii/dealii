// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2022 by the deal.II authors
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
     * In deal.II, we express the orientation of an object with three Booleans:
     * orientation, rotate, and flip. The default values for these are true,
     * false, and false. These are represented either as individual booleans or
     * as a 'combined' orientation: the combined format places orientation in
     * the least significant bit, then rotate, then flip.
     *
     * For a quadrilateral, these values correspond to
     * - *orientation* : `true` is the default orientation and `false` means
     *   vertices 1 and 2 are swapped.
     * - *rotation* : all vertices are rotated by 90 degrees clockwise.
     * - *flip* : all vertices are rotated by 180 degrees clockwise.
     *
     * For a triangle, these values correspond to
     * - *orientation* : `true` is the default orientation and `false` means
     *   vertices 1 and 2 are swapped.
     * - *rotate* : all vertices are rotated by 60 degrees clockwise.
     * - *flip* : all vertices are rotated by 120 degrees clockwise.
     *
     * Here, 'clockwise' is relative to the vector defined by the cross
     * product of two lines in their standard orientation (which, e.g., points
     * into the hexahedron for face 0 but out of the hexahedron for face 1).
     *
     * For triangles, to enable indexing from the combined orientation, we do
     * not consider flip-rotate or flip-orient-rotate as those cases are
     * equivalent, respectively, to the identity operation or the orientation =
     * `true` case as flip-rotate is equal to the identity operation. This
     * choice ensures that the integer value of the combined orientation is in
     * $[0, 5]$.
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
      unsigned char
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
      set_combined_orientation(const unsigned int  object,
                               const unsigned char value);

      /**
       * Set the orientation bit of the object.
       */
      void
      set_orientation(const unsigned int object, const bool value);

      /**
       * Set the rotate bit of the object.
       */
      void
      set_rotation(const unsigned int object, const bool value);

      /**
       * Set the flip bit of the object.
       */
      void
      set_flip(const unsigned int object, const bool value);

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
       * Flags.
       */
      std::vector<unsigned char> flags;
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
      // Assign to the default orientation
      flags.assign(n_objects,
                   ReferenceCell::default_combined_face_orientation());
    }



    inline void
    TriaObjectsOrientations::resize(const unsigned int n_objects)
    {
      flags.resize(n_objects,
                   ReferenceCell::default_combined_face_orientation());
      n_stored_objects = n_objects;
    }



    inline std::size_t
    TriaObjectsOrientations::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(n_stored_objects) +
             MemoryConsumption::memory_consumption(flags);
    }



    inline unsigned int
    TriaObjectsOrientations::n_objects() const
    {
      return n_stored_objects;
    }



    inline unsigned char
    TriaObjectsOrientations::get_combined_orientation(
      const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return flags[object];
    }



    inline bool
    TriaObjectsOrientations::get_orientation(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return Utilities::get_bit(flags[object], 0);
    }



    inline bool
    TriaObjectsOrientations::get_rotation(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return Utilities::get_bit(flags[object], 1);
    }



    inline bool
    TriaObjectsOrientations::get_flip(const unsigned int object) const
    {
      AssertIndexRange(object, n_stored_objects);
      return Utilities::get_bit(flags[object], 2);
    }



    inline void
    TriaObjectsOrientations::set_combined_orientation(const unsigned int object,
                                                      const unsigned char value)
    {
      AssertIndexRange(object, n_stored_objects);
      flags[object] = value;
    }



    inline void
    TriaObjectsOrientations::set_orientation(const unsigned int object,
                                             const bool         value)
    {
      AssertIndexRange(object, n_stored_objects);
      Utilities::set_bit(flags[object], 0, value);
    }



    inline void
    TriaObjectsOrientations::set_rotation(const unsigned int object,
                                          const bool         value)
    {
      AssertIndexRange(object, n_stored_objects);
      Utilities::set_bit(flags[object], 1, value);
    }



    inline void
    TriaObjectsOrientations::set_flip(const unsigned int object,
                                      const bool         value)
    {
      AssertIndexRange(object, n_stored_objects);
      Utilities::set_bit(flags[object], 2, value);
    }



    template <class Archive>
    void
    TriaObjectsOrientations::serialize(Archive &ar, const unsigned int)
    {
      ar &n_stored_objects &flags;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
