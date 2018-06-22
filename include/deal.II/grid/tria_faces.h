// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

#ifndef dealii_tria_faces_h
#define dealii_tria_faces_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria_object.h>
#include <deal.II/grid/tria_objects.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * General template for information belonging to the faces of a
     * triangulation. These classes are similar to the TriaLevel classes. As
     * cells are organised in a hierarchical structure of levels, each
     * triangulation consists of several such TriaLevels. However the faces of
     * a triangulation, lower dimensional objects like lines in 2D or lines
     * and quads in 3D, do not have to be based on such a hierarchical
     * structure. In fact we have to organise them in only one object if we
     * want to enable anisotropic refinement. Therefore the TriaFaces classes
     * store the information belonging to the faces of a triangulation
     * separately from the TriaLevel classes.
     *
     * This general template is only provided to enable a specialization
     * below.
     *
     * @author Tobias Leicht, 2006
     */

    template <int dim>
    class TriaFaces
    {
    public:
      /**
       * Constructor. This constructor is deleted to prevent the use of this template,
       * as only the specializations should be used
       */
      TriaFaces() = delete;
    };



    /**
     * Faces only have a meaning in <tt>dim@>=1</tt>. In <tt>dim=1</tt> they
     * are vertices, which are handled differently, so only for
     * <tt>dim@>=2</tt> the use of TriaFaces is reasonable, for <tt>dim=1</tt>
     * the class is empty.
     */
    template <>
    class TriaFaces<1>
    {
    public:
      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object. Of course this returns 0.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };

    /**
     * In <tt>dim=2</tt> the cells are quads, the faces accordingly are lines.
     */
    template <>
    class TriaFaces<2>
    {
    public:
      /**
       * The TriaObject containing the data of lines.
       */
      TriaObjects<TriaObject<1>> lines;

    public:
      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };

    /**
     * In <tt>dim=3</tt> the cells are hexes, the faces accordingly are quads.
     * In addition to that we also have to enable the storage of lines.
     */
    template <>
    class TriaFaces<3>
    {
    public:
      /**
       * The TriaObject containing the data of quads.
       */

      TriaObjectsQuad3D quads;

      /**
       * The TriaObject containing the data of lines.
       */
      TriaObjects<TriaObject<1>> lines;

    public:
      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    template <class Archive>
    void
    TriaFaces<1>::serialize(Archive &, const unsigned int)
    {}



    template <class Archive>
    void
    TriaFaces<2>::serialize(Archive &ar, const unsigned int)
    {
      ar &lines;
    }



    template <class Archive>
    void
    TriaFaces<3>::serialize(Archive &ar, const unsigned int)
    {
      ar &quads &lines;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
