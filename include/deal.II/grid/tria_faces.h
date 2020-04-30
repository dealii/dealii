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
     * This class contains information belonging to the faces of a
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
     * @author Tobias Leicht, 2006
     * @author Peter Munch, 2020
     */
    class TriaFaces
    {
    public:
      /**
       * Constructor.
       */
      TriaFaces(const unsigned int dim);

      /**
       * Default constructor for Boost::serialization.
       */
      TriaFaces() = default;

      /**
       * Dimension of the underlying triangulation.
       */
      unsigned int dim;

      /**
       * The TriaObject containing the data of quads.
       *
       * @note Used only for dim=3.
       */
      TriaObjects quads;

      /**
       * Orientation of each line of each quad.
       *
       * @note Used only for dim=3.
       */
      std::vector<unsigned char> quads_line_orientations;

      /**
       * The TriaObject containing the data of lines.
       *
       * @note Used only for dim>1.
       */
      TriaObjects lines;

      /**
       * Reserve space for line_orientations.
       *
       * @note Used only for dim=3.
       */
      void
      reserve_space(const unsigned int new_quads_in_pairs,
                    const unsigned int new_quads_single = 0);

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
    TriaFaces::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      if (dim == 2)
        ar &lines;

      if (dim == 3)
        ar &quads &lines &quads_line_orientations;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
