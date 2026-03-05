// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
    template <int dim>
    class TriaFaces
    {
    public:
      /**
       * Constructor.
       *
       * @param[in] max_children_per_quad Maximum number of children (across all
       *            relevant ReferenceCell types) a quad (i.e., a face of a cell
       *            in 3d) in the present Triangulation may have. Only used in
       *            3d.
       *
       * @param[in] max_lines_per_quad Maximum number of faces (i.e., neighbors)
       *            per cell. Like @p max_children_per_quad, this is the maximum
       *            over all relevant ReferenceCell types and is only used in 2d
       *            and 3d.
       */
      TriaFaces(const unsigned int max_children_per_quad,
                const unsigned int max_lines_per_quad);

      /**
       * Default constructor for Boost::serialization.
       */
      TriaFaces() = default;

      /**
       * Resize all internal arrays and populate with default values.
       */
      void
      allocate(const std::size_t n_lines, const std::size_t n_quads);

      /**
       * Maximum number of lines for each quad (i.e., face of a 3d cell).
       * Only used for `dim == 3`.
       */
      unsigned int max_lines_per_quad;

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
      ReferenceCell<dim - 1>
      get_quad_type(const std::size_t index) const;

      /**
       * Helper accessor function for quad_is_quadrilateral
       */
      void
      set_quad_type(const std::size_t            index,
                    const ReferenceCell<dim - 1> face_type);

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



    template <int dim>
    inline ReferenceCell<dim - 1>
    TriaFaces<dim>::get_quad_type(const std::size_t index) const
    {
      AssertIndexRange(index, quad_is_quadrilateral.size());
      if constexpr (dim == 3)
        return quad_is_quadrilateral[index] ? ReferenceCells::Quadrilateral :
                                              ReferenceCells::Triangle;
      else
        {
          DEAL_II_ASSERT_UNREACHABLE();
          return ReferenceCells::Invalid<dim - 1>;
        }
    }



    template <int dim>
    inline void
    TriaFaces<dim>::set_quad_type(const std::size_t            index,
                                  const ReferenceCell<dim - 1> face_type)
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



    template <int dim>
    template <class Archive>
    void
    TriaFaces<dim>::serialize(Archive &ar, const unsigned int)
    {
      ar &quads &lines &quads_line_orientations &quad_is_quadrilateral;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
