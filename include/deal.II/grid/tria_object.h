// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_tria_object_h
#define dealii_tria_object_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * View onto the current geometric object and its faces.
     *
     * @note The geometric objects are not saved as separate instances of
     * TriaObject but in a single vector in TriaObjects.
     */
    class TriaObject
    {
    public:
      /**
       * Constructor.
       *
       * @param faces ArrayView inside of the vector in TriaObjects.
       */
      TriaObject(const ArrayView<int> &faces)
        : faces(faces)
      {}

      /**
       * Store the content of @p other in the vector of TriaObjects.
       */
      TriaObject &
      operator=(const std::initializer_list<int> &other)
      {
        AssertDimension(faces.size(), other.size());

        const std::vector<int> other_v = other;

        for (unsigned int i = 0; i < faces.size(); ++i)
          faces[i] = other_v[i];

        return *this;
      }

      /**
       * The same as above but for `unsigned int`.
       */
      TriaObject &
      operator=(const std::initializer_list<unsigned int> &other)
      {
        AssertDimension(faces.size(), other.size());

        const std::vector<unsigned int> other_v = other;

        for (unsigned int i = 0; i < faces.size(); ++i)
          faces[i] = other_v[i];

        return *this;
      }

      /**
       * Return index of the @p i-th face.
       */
      int
      face(const unsigned int i) const
      {
        return faces[i];
      }

      /**
       * Set index of the @p i-th face.
       */
      void
      set_face(const unsigned int i, const int index)
      {
        faces[i] = index;
      }

    private:
      /**
       * List of face this object is made up.
       */
      const ArrayView<int> faces;
    };
  } // namespace TriangulationImplementation
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
