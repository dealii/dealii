// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_faces_h
#define dealii_tria_faces_h

#include <deal.II/base/config.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_objects.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * This class contains information belonging to the faces of a
     * triangulation. These classes are similar to the TriaLevel classes. As
     * cells are organized in a hierarchical structure of levels, each
     * triangulation consists of several such TriaLevels. However the faces of
     * a triangulation, lower dimensional objects like lines in 2d or lines
     * and quads in 3d, do not have to be based on such a hierarchical
     * structure. In fact we have to organise them in only one object if we
     * want to enable anisotropic refinement. Therefore the TriaFaces classes
     * store the information belonging to the faces of a triangulation
     * separately from the TriaLevel classes.
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
       * Orientation of each line of each quad. Like elsewhere, `true` refers to
       * the standard orientation and `false` refers to the reverse orientation.
       *
       * @note Used only for dim=3.
       */
      std::vector<bool> quads_line_orientations;

      /**
       * Whether or not each quad is a Quadrilateral. Since, if dim = 3, faces
       * are either Triangles or Quadrilaterals, it suffices to store a
       * boolean.
       *
       * @note Used only for dim=3.
       */
      std::vector<bool> quad_is_quadrilateral;

      /**
       * Helper accessor function for quad_is_quadrilateral
       */
      ReferenceCell
      get_quad_type(const std::size_t index) const;

      /**
       * Helper accessor function for quad_is_quadrilateral
       */
      void
      set_quad_type(const std::size_t index, const ReferenceCell face_type);

      /**
       * The TriaObject containing the data of lines.
       *
       * @note Used only for dim>1.
       */
      TriaObjects lines;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    inline ReferenceCell
    TriaFaces::get_quad_type(const std::size_t index) const
    {
      AssertIndexRange(index, quad_is_quadrilateral.size());
      return quad_is_quadrilateral[index] ? ReferenceCells::Quadrilateral :
                                            ReferenceCells::Triangle;
    }



    inline void
    TriaFaces::set_quad_type(const std::size_t   index,
                             const ReferenceCell face_type)
    {
      AssertIndexRange(index, quad_is_quadrilateral.size());
      Assert(face_type == ReferenceCells::Quadrilateral ||
               face_type == ReferenceCells::Triangle,
             ExcInternalError());
      if (face_type == ReferenceCells::Quadrilateral)
        quad_is_quadrilateral[index] = true;
      else
        quad_is_quadrilateral[index] = false;
    }



    template <class Archive>
    void
    TriaFaces::serialize(Archive &ar, const unsigned int)
    {
      ar                                        &dim;
      ar &quads &lines &quads_line_orientations &quad_is_quadrilateral;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
