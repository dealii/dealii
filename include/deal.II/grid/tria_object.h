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

#ifndef DOXYGEN
namespace internal
{
  namespace TriangulationImplementation
  {
    class TriaObjectView;
  }
} // namespace internal
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * Class template for the <tt>structdim</tt>-dimensional cells
     * constituting a dealii::Triangulation of dimension <tt>structdim</tt> or
     * lower dimensional objects of higher dimensions.  They are characterized
     * by the (global) indices of their faces, which are cells of dimension
     * <tt>structdim-1</tt> or vertices if <tt>structdim=1</tt>.
     *
     * @note This class is only used during setup of the Triangulation and its
     *   content is saved into a single long vector inside of `TriaObjects`.
     *
     * @author Guido Kanschat, 2007
     */
    class TriaObject
    {
    public:
      /**
       * Constructor.
       */
      TriaObject(const std::initializer_list<int> &faces);

      /**
       * Constructor.
       */
      TriaObject(const std::initializer_list<unsigned int> &faces);

      /**
       * Return the index of the ith face object.
       */
      int
      face(const unsigned int i) const;

      /**
       * Set the index of the ith face object.
       */
      void
      set_face(const unsigned int i, const int index);

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption();

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);

    private:
      /**
       * Global indices of the face iterators bounding this cell if dim@>1,
       * and the two vertex indices in 1d.
       */
      std::vector<int> faces;

      friend TriaObjectView;
    };

    /**
     * View onto the current geometric object and its faces.
     *
     * @note The geometric objects are not saved as separate instances of
     * TriaObject but in a single vector in TriaObjects.
     */
    class TriaObjectView
    {
    public:
      /**
       * Constructor.
       *
       * @param faces ArrayView inside of the vector in TriaObjects.
       */
      TriaObjectView(const ArrayView<int> &faces)
        : faces(faces)
      {}

      /**
       * Store the content of @p other in the vector of TriaObjects.
       */
      TriaObjectView &
      operator=(const TriaObject &other)
      {
        AssertDimension(faces.size(), other.faces.size());

        for (unsigned int i = 0; i < faces.size(); ++i)
          faces[i] = other.faces[i];

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

    //----------------------------------------------------------------------//


    inline TriaObject::TriaObject(const std::initializer_list<int> &faces)
    {
      this->faces = faces;
    }



    inline TriaObject::TriaObject(
      const std::initializer_list<unsigned int> &faces_in)
    {
      this->faces.reserve(faces_in.size());
      for (const auto face : faces_in)
        this->faces.push_back(face);
    }



    inline int
    TriaObject::face(const unsigned int i) const
    {
      AssertIndexRange(i, faces.size());
      return faces[i];
    }



    inline void
    TriaObject::set_face(const unsigned int i, const int index)
    {
      AssertIndexRange(i, faces.size());
      faces[i] = index;
    }



    inline std::size_t
    TriaObject::memory_consumption()
    {
      return MemoryConsumption::memory_consumption(this->faces);
    }


    template <class Archive>
    void
    TriaObject::serialize(Archive &ar, const unsigned int)
    {
      ar &faces;
    }


  } // namespace TriangulationImplementation
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
