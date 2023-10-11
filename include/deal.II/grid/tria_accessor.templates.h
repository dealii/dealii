// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2023 by the deal.II authors
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

#ifndef dealii_tria_accessor_templates_h
#define dealii_tria_accessor_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_faces.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/grid/tria_levels.h>

#include <boost/container/small_vector.hpp>

#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriaAccessorImplementation
  {
    /**
     * Compute the diameter for a given set of vertices. The vertices are
     * interpreted, depending on their count, as the vertices of a particular
     * reference cell.
     */
    template <int dim, int spacedim>
    inline double
    diameter(
      const boost::container::small_vector<Point<spacedim>,
                                           GeometryInfo<dim>::vertices_per_cell>
        vertices)
    {
      const ReferenceCell reference_cell =
        ReferenceCell::n_vertices_to_type(dim, vertices.size());

      if (reference_cell == ReferenceCells::Line)
        // Return the distance between the two vertices
        return (vertices[1] - vertices[0]).norm();
      else if (reference_cell == ReferenceCells::Triangle)
        // Return the longest of the three edges
        return std::max({(vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[2] - vertices[0]).norm()});
      else if (reference_cell == ReferenceCells::Quadrilateral)
        // Return the longer one of the two diagonals of the quadrilateral
        return std::max({(vertices[3] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm()});
      else if (reference_cell == ReferenceCells::Tetrahedron)
        // Return the longest of the six edges of the tetrahedron
        return std::max({(vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[3] - vertices[0]).norm(),
                         (vertices[3] - vertices[1]).norm(),
                         (vertices[3] - vertices[2]).norm()});
      else if (reference_cell == ReferenceCells::Pyramid)
        // Return ...
        return std::max({// the longest diagonal of the quadrilateral base
                         // of the pyramid or ...
                         (vertices[3] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         // the longest edge connected with the apex of the
                         // pyramid
                         (vertices[4] - vertices[0]).norm(),
                         (vertices[4] - vertices[1]).norm(),
                         (vertices[4] - vertices[2]).norm(),
                         (vertices[4] - vertices[3]).norm()});
      else if (reference_cell == ReferenceCells::Wedge)
        // Return ...
        return std::max({// the longest of the 2*3=6 diagonals of the three
                         // quadrilateral sides of the wedge or ...
                         (vertices[4] - vertices[0]).norm(),
                         (vertices[3] - vertices[1]).norm(),
                         (vertices[5] - vertices[1]).norm(),
                         (vertices[4] - vertices[2]).norm(),
                         (vertices[5] - vertices[0]).norm(),
                         (vertices[3] - vertices[2]).norm(),
                         // the longest of the 3*2=6 edges of the two
                         // triangular faces of the wedge
                         (vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[2] - vertices[0]).norm(),
                         (vertices[4] - vertices[3]).norm(),
                         (vertices[5] - vertices[4]).norm(),
                         (vertices[5] - vertices[3]).norm()});
      else if (reference_cell == ReferenceCells::Hexahedron)
        // Return the longest of the four diagonals of the hexahedron
        return std::max({(vertices[7] - vertices[0]).norm(),
                         (vertices[6] - vertices[1]).norm(),
                         (vertices[2] - vertices[5]).norm(),
                         (vertices[3] - vertices[4]).norm()});

      Assert(false, ExcNotImplemented());
      return -1e10;
    }
  } // namespace TriaAccessorImplementation
} // namespace internal


/*--------------------- Functions: TriaAccessorBase -------------------------*/

template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim>::TriaAccessorBase(
  const Triangulation<dim, spacedim> *tria,
  const int                           level,
  const int                           index,
  const AccessorData *)
  : present_level((structdim == dim) ? level : 0)
  , present_index(index)
  , tria(tria)
{
  // non-cells have no level, so a 0
  // should have been passed, or a -1
  // for an end-iterator, or -2 for
  // an invalid (default constructed)
  // iterator
  if (structdim != dim)
    {
      Assert((level == 0) || (level == -1) || (level == -2),
             ExcInternalError());
    }
}


template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim>::TriaAccessorBase(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
  : present_level(a.present_level)
  , present_index(a.present_index)
  , tria(a.tria)
{}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::copy_from(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria          = a.tria;

  if (structdim != dim)
    {
      Assert((present_level == 0) || (present_level == -1) ||
               (present_level == -2),
             ExcInternalError());
    }
}



template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim> &
TriaAccessorBase<structdim, dim, spacedim>::operator=(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria          = a.tria;

  if (structdim != dim)
    {
      Assert((present_level == 0) || (present_level == -1) ||
               (present_level == -2),
             ExcInternalError());
    }
  return *this;
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator==(
  const TriaAccessorBase<structdim, dim, spacedim> &a) const
{
  Assert(tria == a.tria || tria == nullptr || a.tria == nullptr,
         TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria == a.tria) && (present_level == a.present_level) &&
          (present_index == a.present_index));
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator!=(
  const TriaAccessorBase<structdim, dim, spacedim> &a) const
{
  Assert(tria == a.tria || tria == nullptr || a.tria == nullptr,
         TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria != a.tria) || (present_level != a.present_level) ||
          (present_index != a.present_index));
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator<(
  const TriaAccessorBase<structdim, dim, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  if (present_level != other.present_level)
    return (present_level < other.present_level);

  return (present_index < other.present_index);
}



template <int structdim, int dim, int spacedim>
inline int
TriaAccessorBase<structdim, dim, spacedim>::level() const
{
  // This is always zero or invalid
  // if the object is not a cell
  return present_level;
}



template <int structdim, int dim, int spacedim>
inline int
TriaAccessorBase<structdim, dim, spacedim>::index() const
{
  return present_index;
}



template <int structdim, int dim, int spacedim>
inline IteratorState::IteratorStates
TriaAccessorBase<structdim, dim, spacedim>::state() const
{
  if ((present_level >= 0) && (present_index >= 0))
    return IteratorState::valid;
  else if (present_index == -1)
    return IteratorState::past_the_end;
  else
    return IteratorState::invalid;
}



template <int structdim, int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
TriaAccessorBase<structdim, dim, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::operator++()
{
  // this iterator is used for
  // objects without level
  ++this->present_index;

  if (structdim != dim)
    {
      // is index still in the range of
      // the vector? (note that we don't
      // have to set the level, since
      // dim!=1 and the object therefore
      // has no level)
      if (this->present_index >= static_cast<int>(objects().n_objects()))
        this->present_index = -1;
    }
  else
    {
      while (this->present_index >=
             static_cast<int>(
               this->tria->levels[this->present_level]->cells.n_objects()))
        {
          // no -> go one level up until we find
          // one with more than zero cells
          ++this->present_level;
          this->present_index = 0;
          // highest level reached?
          if (this->present_level >=
              static_cast<int>(this->tria->levels.size()))
            {
              // return with past the end pointer
              this->present_level = this->present_index = -1;
              return;
            }
        }
    }
}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::operator--()
{
  // same as operator++
  --this->present_index;

  if (structdim != dim)
    {
      if (this->present_index < 0)
        this->present_index = -1;
    }
  else
    {
      while (this->present_index < 0)
        {
          // no -> go one level down
          --this->present_level;
          // lowest level reached?
          if (this->present_level == -1)
            {
              // return with past the end pointer
              this->present_level = this->present_index = -1;
              return;
            }
          // else
          this->present_index =
            this->tria->levels[this->present_level]->cells.n_objects() - 1;
        }
    }
}



template <int structdim, int dim, int spacedim>
inline dealii::internal::TriangulationImplementation::TriaObjects &
TriaAccessorBase<structdim, dim, spacedim>::objects() const
{
  if (structdim == dim)
    return this->tria->levels[this->present_level]->cells;

  if (structdim == 1 && dim > 1)
    return this->tria->faces->lines;

  if (structdim == 2 && dim > 2)
    return this->tria->faces->quads;

  Assert(false, ExcInternalError());

  return this->tria->levels[this->present_level]->cells;
}



/*---------------------- Functions: InvalidAccessor -------------------------*/

template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(const void *,
                                                           const int,
                                                           const int,
                                                           const AccessorData *)
{
  Assert(false,
         ExcMessage("You are attempting an invalid conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(
  const InvalidAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an invalid conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::copy_from(const InvalidAccessor &)
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::operator==(
  const InvalidAccessor &) const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::operator!=(
  const InvalidAccessor &) const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return true;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::used() const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::has_children() const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::operator++() const
{}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::operator--() const
{}



template <int structdim, int dim, int spacedim>
types::manifold_id
InvalidAccessor<structdim, dim, spacedim>::manifold_id() const
{
  return numbers::flat_manifold_id;
}



template <int structdim, int dim, int spacedim>
unsigned int
InvalidAccessor<structdim, dim, spacedim>::user_index() const
{
  return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::set_user_index(
  const unsigned int) const
{
  Assert(false,
         ExcMessage("You are trying to set the user index of an "
                    "invalid object."));
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::set_manifold_id(
  const types::manifold_id) const
{
  Assert(false,
         ExcMessage("You are trying to set the manifold id of an "
                    "invalid object."));
}



template <int structdim, int dim, int spacedim>
inline Point<spacedim> &
InvalidAccessor<structdim, dim, spacedim>::vertex(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  static Point<spacedim> invalid_vertex;
  return invalid_vertex;
}


template <int structdim, int dim, int spacedim>
inline void *
InvalidAccessor<structdim, dim, spacedim>::line(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return nullptr;
}



template <int structdim, int dim, int spacedim>
inline void *
InvalidAccessor<structdim, dim, spacedim>::quad(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return nullptr;
}


/*------------------------ Functions: TriaAccessor ---------------------------*/


namespace internal
{
  namespace TriaAccessorImplementation
  {
    // make sure that if in the following we
    // write TriaAccessor
    // we mean the *class*
    // dealii::TriaAccessor, not the
    // enclosing namespace
    // dealii::internal::TriaAccessor
    using dealii::TriaAccessor;

    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int dim, int spacedim>
      inline static unsigned int
      line_index(const TriaAccessor<1, dim, spacedim> &, const unsigned int)
      {
        Assert(false,
               ExcMessage("You can't ask for the index of a line bounding "
                          "a one-dimensional cell because it is not "
                          "bounded by lines."));
        return numbers::invalid_unsigned_int;
      }


      template <int dim, int spacedim>
      inline static unsigned int
      line_index(const TriaAccessor<2, dim, spacedim> &accessor,
                 const unsigned int                    i)
      {
        constexpr unsigned int max_faces_per_cell = 4;
        return accessor.objects()
          .cells[accessor.present_index * max_faces_per_cell + i];
      }


      inline static unsigned int
      line_index(const TriaAccessor<3, 3, 3> &accessor, const unsigned int i)
      {
        const auto [face_index, line_index] =
          accessor.reference_cell().standard_line_to_face_and_line_index(i);
        const auto line_within_face_index =
          accessor.reference_cell().standard_to_real_face_line(
            line_index,
            face_index,
            combined_face_orientation(accessor, face_index));

        return accessor.quad(face_index)->line_index(line_within_face_index);
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline static unsigned int
      quad_index(const TriaAccessor<structdim, dim, spacedim> &,
                 const unsigned int)
      {
        Assert(false,
               ExcMessage("You can't ask for the index of a quad bounding "
                          "a one- or two-dimensional cell because it is not "
                          "bounded by quads."));
        return numbers::invalid_unsigned_int;
      }


      inline static unsigned int
      quad_index(const TriaAccessor<3, 3, 3> &accessor, const unsigned int i)
      {
        constexpr unsigned int max_faces_per_cell = 6;
        return accessor.tria->levels[accessor.present_level]
          ->cells.cells[accessor.present_index * max_faces_per_cell + i];
      }



      /**
       * Implementation of the function of some name in the mother class
       */
      template <int structdim, int dim, int spacedim>
      inline static bool
      face_orientation(const TriaAccessor<structdim, dim, spacedim> &,
                       const unsigned int)
      {
        /*
         * Default implementation used in 1d
         *
         * In 1d, face_orientation is always true
         */

        return true;
      }


      template <int spacedim>
      inline static bool
      face_orientation(const TriaAccessor<2, 2, spacedim> &accessor,
                       const unsigned int                  face)
      {
        return line_orientation(accessor, face);
      }


      inline static bool
      face_orientation(const TriaAccessor<3, 3, 3> &accessor,
                       const unsigned int           face)
      {
        return accessor.tria->levels[accessor.present_level]
          ->face_orientations.get_orientation(
            accessor.present_index * GeometryInfo<3>::faces_per_cell + face);
      }



      template <int dim, int spacedim>
      inline static unsigned char
      combined_face_orientation(
        const TriaAccessor<1, dim, spacedim> & /*accessor*/,
        const unsigned int /*face*/)
      {
        // There is only one way to orient a vertex
        return ReferenceCell::default_combined_face_orientation();
      }



      template <int dim, int spacedim>
      inline static unsigned char
      combined_face_orientation(const TriaAccessor<2, dim, spacedim> &accessor,
                                const unsigned int                    face)
      {
        return line_orientation(accessor, face) == true ?
                 ReferenceCell::default_combined_face_orientation() :
                 ReferenceCell::reversed_combined_line_orientation();
      }



      inline static unsigned char
      combined_face_orientation(const TriaAccessor<3, 3, 3> &accessor,
                                const unsigned int           face)
      {
        AssertIndexRange(face, accessor.n_faces());
        return accessor.tria->levels[accessor.present_level]
          ->face_orientations.get_combined_orientation(
            accessor.present_index * GeometryInfo<3>::faces_per_cell + face);
      }


      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline static bool
      face_flip(const TriaAccessor<structdim, dim, spacedim> &,
                const unsigned int)
      {
        /*
         * Default implementation used in 1d and 2d
         *
         * In 1d, face_flip is always false as there is no such concept as
         * "flipped" faces in 1d.
         *
         * In 2d, we currently only support meshes where all faces are in
         * standard orientation, so the result is also false. This also
         * matches the fact that one can *always* orient faces in 2d in such a
         * way that the don't need to be flipped
         */
        return false;
      }



      inline static bool
      face_flip(const TriaAccessor<3, 3, 3> &accessor, const unsigned int face)
      {
        AssertIndexRange(face, accessor.n_faces());
        return accessor.tria->levels[accessor.present_level]
          ->face_orientations.get_flip(
            accessor.present_index * GeometryInfo<3>::faces_per_cell + face);
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline static bool
      face_rotation(const TriaAccessor<structdim, dim, spacedim> &,
                    const unsigned int)
      {
        /*
         * Default implementation used in 1d and 2d
         *
         * In 1d and 2d, face_rotation is always false as there is no such
         * concept as "rotated" faces in 1d and 2d.
         */
        return false;
      }


      inline static bool
      face_rotation(const TriaAccessor<3, 3, 3> &accessor,
                    const unsigned int           face)
      {
        AssertIndexRange(face, accessor.n_faces());

        return accessor.tria->levels[accessor.present_level]
          ->face_orientations.get_rotation(
            accessor.present_index * GeometryInfo<3>::faces_per_cell + face);
      }

      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int dim, int spacedim>
      inline static bool
      line_orientation(const TriaAccessor<1, dim, spacedim> &,
                       const unsigned int)
      {
        return true;
      }


      template <int spacedim>
      inline static bool
      line_orientation(const TriaAccessor<2, 2, spacedim> &accessor,
                       const unsigned int                  line)
      {
        AssertIndexRange(line, accessor.n_lines());
        if (accessor.tria->levels[accessor.present_level]
              ->face_orientations.n_objects() == 0)
          return true; // quads in 2d have no non-standard orientation and
                       // face_orientations is left empty
        else
          return accessor.tria->levels[accessor.present_level]
            ->face_orientations.get_orientation(
              accessor.present_index * GeometryInfo<2>::faces_per_cell + line);
      }


      template <int spacedim>
      inline static bool
      line_orientation(const TriaAccessor<2, 3, spacedim> &accessor,
                       const unsigned int                  line)
      {
        Assert(accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        Assert(accessor.present_index * GeometryInfo<3>::lines_per_face + line <
                 accessor.tria->faces->quads_line_orientations.size(),
               ExcInternalError());

        // quads as part of 3d hexes can have non-standard orientation
        return accessor.tria->faces->quads_line_orientations
          [accessor.present_index * GeometryInfo<3>::lines_per_face + line];
      }


      inline static bool
      line_orientation(const TriaAccessor<3, 3, 3> &accessor,
                       const unsigned int           line)
      {
        Assert(accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        AssertIndexRange(line, accessor.n_lines());

        // First pick a face on which this line is a part of, and the
        // index of the line within.
        const auto [face_index, line_index] =
          accessor.reference_cell().standard_line_to_face_and_line_index(line);
        const auto line_within_face_index =
          accessor.reference_cell().standard_to_real_face_line(
            line_index,
            face_index,
            combined_face_orientation(accessor, face_index));

        // Then query how that line is oriented within that face:
        return accessor.reference_cell().standard_vs_true_line_orientation(
          line_index,
          face_index,
          combined_face_orientation(accessor, face_index),
          accessor.quad(face_index)->line_orientation(line_within_face_index));
      }

      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline static void
      set_combined_face_orientation(
        const TriaAccessor<structdim, dim, spacedim> &,
        const unsigned int,
        const unsigned char)
      {
        Assert(false, ExcInternalError());
      }



      inline static void
      set_combined_face_orientation(const TriaAccessor<3, 3, 3> &accessor,
                                    const unsigned int           face,
                                    const unsigned char combined_orientation)
      {
        AssertIndexRange(face, accessor.n_faces());
        accessor.tria->levels[accessor.present_level]
          ->face_orientations.set_combined_orientation(
            accessor.present_index * GeometryInfo<3>::faces_per_cell + face,
            combined_orientation);
      }

      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int dim, int spacedim>
      inline static void
      set_line_orientation(const TriaAccessor<1, dim, spacedim> &,
                           const unsigned int,
                           const bool)
      {
        Assert(false, ExcInternalError());
      }


      template <int spacedim>
      inline static void
      set_line_orientation(const TriaAccessor<2, 2, spacedim> &,
                           const unsigned int,
                           const bool)
      {
        // quads in 2d have no
        // non-standard orientation
        Assert(false, ExcInternalError());
      }


      template <int spacedim>
      inline static void
      set_line_orientation(const TriaAccessor<2, 3, spacedim> &accessor,
                           const unsigned int                  line,
                           const bool                          value)
      {
        Assert(accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        AssertIndexRange(line, accessor.n_lines());
        Assert(accessor.present_index * GeometryInfo<3>::lines_per_face + line <
                 accessor.tria->faces->quads_line_orientations.size(),
               ExcInternalError());

        // quads as part of 3d hexes can have non-standard orientation
        accessor.tria->faces->quads_line_orientations
          [accessor.present_index * GeometryInfo<3>::lines_per_face + line] =
          value;
      }


      inline static void
      set_line_orientation(const TriaAccessor<3, 3, 3> &,
                           const unsigned int,
                           const bool)
      {
        // it seems like we don't need this
        // one
        Assert(false, ExcNotImplemented());
      }


      /**
       * Implementation of the function of same name in the enclosing class.
       */
      template <int dim, int spacedim>
      inline static unsigned int
      vertex_index(const TriaAccessor<1, dim, spacedim> &accessor,
                   const unsigned int                    corner)
      {
        constexpr unsigned int max_faces_per_cell = 2;
        return accessor.objects()
          .cells[accessor.present_index * max_faces_per_cell + corner];
      }


      template <int dim, int spacedim>
      inline static unsigned int
      vertex_index(const TriaAccessor<2, dim, spacedim> &accessor,
                   const unsigned int                    corner)
      {
        const auto [line_index, vertex_index] =
          accessor.reference_cell().standard_vertex_to_face_and_vertex_index(
            corner);
        const auto vertex_within_line_index =
          accessor.reference_cell().standard_to_real_face_vertex(
            vertex_index, line_index, accessor.line_orientation(line_index));

        return accessor.line(line_index)
          ->vertex_index(vertex_within_line_index);
      }



      inline static unsigned int
      vertex_index(const TriaAccessor<3, 3, 3> &accessor,
                   const unsigned int           corner)
      {
        const auto [face_index, vertex_index] =
          accessor.reference_cell().standard_vertex_to_face_and_vertex_index(
            corner);
        const auto vertex_within_face_index =
          accessor.reference_cell().standard_to_real_face_vertex(
            vertex_index,
            face_index,
            combined_face_orientation(accessor, face_index));

        return accessor.quad(face_index)
          ->vertex_index(vertex_within_face_index);
      }



      template <int dim, int spacedim>
      static std::array<unsigned int, 1>
      get_line_indices_of_cell(const TriaAccessor<1, dim, spacedim> &)
      {
        Assert(false, ExcInternalError());
        return {};
      }



      template <int structdim, int dim, int spacedim>
      static std::array<unsigned int, 4>
      get_line_indices_of_cell(const TriaAccessor<2, dim, spacedim> &cell)
      {
        // For 2d cells the access cell->line_orientation() is already
        // efficient
        std::array<unsigned int, 4> line_indices = {};
        for (const unsigned int line : cell.line_indices())
          line_indices[line] = cell.line_index(line);
        return line_indices;
      }

      /**
       * A helper function to provide faster access to cell->line_index() in
       * 3d
       */
      template <int structdim, int dim, int spacedim>
      static std::array<unsigned int, 12>
      get_line_indices_of_cell(
        const TriaAccessor<structdim, dim, spacedim> &cell)
      {
        std::array<unsigned int, 12> line_indices = {};

        // For hexahedra, the classical access via quads -> lines is too
        // inefficient. Unroll this code here to allow the compiler to inline
        // the necessary functions.
        const auto ref_cell = cell.reference_cell();
        if (ref_cell == ReferenceCells::Hexahedron)
          {
            for (unsigned int f = 4; f < 6; ++f)
              {
                const unsigned char orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * GeometryInfo<3>::faces_per_cell + f);

                // It might seem superfluous to spell out the four indices
                // that get later consumed by a for loop over these four
                // elements; however, for the compiler it is easier to inline
                // the statement of standard_to_real_face_line() when next to
                // each other, as opposed to be interleaved with a
                // line_index() call.
                const std::array<unsigned int, 4> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation),
                   ref_cell.standard_to_real_face_line(2, f, orientation),
                   ref_cell.standard_to_real_face_line(3, f, orientation)}};
                const auto quad = cell.quad(f);
                for (unsigned int l = 0; l < 4; ++l)
                  line_indices[4 * (f - 4) + l] =
                    quad->line_index(my_indices[l]);
              }
            for (unsigned int f = 0; f < 2; ++f)
              {
                const unsigned char orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * GeometryInfo<3>::faces_per_cell + f);
                const std::array<unsigned int, 2> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation)}};
                const auto quad      = cell.quad(f);
                line_indices[8 + f]  = quad->line_index(my_indices[0]);
                line_indices[10 + f] = quad->line_index(my_indices[1]);
              }
          }
        else if (ref_cell == ReferenceCells::Tetrahedron)
          {
            std::array<unsigned int, 3> orientations{
              {combined_face_orientation(cell, 0),
               combined_face_orientation(cell, 1),
               combined_face_orientation(cell, 2)}};
            const std::array<unsigned int, 6> my_indices{
              {ref_cell.standard_to_real_face_line(0, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(2, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(2, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(1, 2, orientations[2])}};
            line_indices[0] = cell.quad(0)->line_index(my_indices[0]);
            line_indices[1] = cell.quad(0)->line_index(my_indices[1]);
            line_indices[2] = cell.quad(0)->line_index(my_indices[2]);
            line_indices[3] = cell.quad(1)->line_index(my_indices[3]);
            line_indices[4] = cell.quad(1)->line_index(my_indices[4]);
            line_indices[5] = cell.quad(2)->line_index(my_indices[5]);
          }
        else
          // For other shapes (wedges, pyramids), we do not currently
          // implement an optimized function.
          for (unsigned int l = 0; l < std::min(12U, cell.n_lines()); ++l)
            line_indices[l] = cell.line_index(l);

        return line_indices;
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 1d specialization
       */
      template <int dim, int spacedim>
      static std::array<unsigned int, 1>
      get_line_orientations_of_cell(const TriaAccessor<1, dim, spacedim> &)
      {
        Assert(false, ExcInternalError());
        return {};
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 2d specialization
       */
      template <int dim, int spacedim>
      static std::array<bool, 4>
      get_line_orientations_of_cell(const TriaAccessor<2, dim, spacedim> &cell)
      {
        // For 2d cells the access cell->line_orientation() is already
        // efficient
        std::array<bool, 4> line_orientations = {};
        for (const unsigned int line : cell.line_indices())
          line_orientations[line] = cell.line_orientation(line);
        return line_orientations;
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 3d specialization
       */
      template <int dim, int spacedim>
      static std::array<bool, 12>
      get_line_orientations_of_cell(const TriaAccessor<3, dim, spacedim> &cell)
      {
        std::array<bool, 12> line_orientations = {};

        // For hexahedra, the classical access via quads -> lines is too
        // inefficient. Unroll this code here to allow the compiler to inline
        // the necessary functions.
        const auto ref_cell = cell.reference_cell();
        if (ref_cell == ReferenceCells::Hexahedron)
          {
            for (unsigned int f = 4; f < 6; ++f)
              {
                const unsigned char orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * GeometryInfo<3>::faces_per_cell + f);

                // It might seem superfluous to spell out the four indices and
                // orientations that get later consumed by a for loop over
                // these four elements; however, for the compiler it is easier
                // to inline the statement of standard_to_real_face_line()
                // when next to each other, as opposed to be interleaved with
                // a line_index() call.
                const std::array<unsigned int, 4> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation),
                   ref_cell.standard_to_real_face_line(2, f, orientation),
                   ref_cell.standard_to_real_face_line(3, f, orientation)}};
                const auto                quad = cell.quad(f);
                const std::array<bool, 4> my_orientations{
                  {ref_cell.standard_vs_true_line_orientation(
                     0, f, orientation, quad->line_orientation(my_indices[0])),
                   ref_cell.standard_vs_true_line_orientation(
                     1, f, orientation, quad->line_orientation(my_indices[1])),
                   ref_cell.standard_vs_true_line_orientation(
                     2, f, orientation, quad->line_orientation(my_indices[2])),
                   ref_cell.standard_vs_true_line_orientation(
                     3,
                     f,
                     orientation,
                     quad->line_orientation(my_indices[3]))}};
                for (unsigned int l = 0; l < 4; ++l)
                  line_orientations[4 * (f - 4) + l] = my_orientations[l];
              }
            for (unsigned int f = 0; f < 2; ++f)
              {
                const unsigned char orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * GeometryInfo<3>::faces_per_cell + f);
                const std::array<unsigned int, 2> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation)}};
                const auto                quad = cell.quad(f);
                const std::array<bool, 2> my_orientations{
                  {ref_cell.standard_vs_true_line_orientation(
                     0, f, orientation, quad->line_orientation(my_indices[0])),
                   ref_cell.standard_vs_true_line_orientation(
                     1,
                     f,
                     orientation,
                     quad->line_orientation(my_indices[1]))}};
                line_orientations[8 + f]  = my_orientations[0];
                line_orientations[10 + f] = my_orientations[1];
              }
          }
        else if (ref_cell == ReferenceCells::Tetrahedron)
          {
            std::array<unsigned int, 3> orientations{
              {combined_face_orientation(cell, 0),
               combined_face_orientation(cell, 1),
               combined_face_orientation(cell, 2)}};
            const std::array<unsigned int, 6> my_indices{
              {ref_cell.standard_to_real_face_line(0, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(2, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(2, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(1, 2, orientations[2])}};
            line_orientations[0] = ref_cell.standard_vs_true_line_orientation(
              0,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[0]));
            line_orientations[1] = ref_cell.standard_vs_true_line_orientation(
              1,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[1]));
            line_orientations[2] = ref_cell.standard_vs_true_line_orientation(
              2,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[2]));
            line_orientations[3] = ref_cell.standard_vs_true_line_orientation(
              1,
              1,
              orientations[1],
              cell.quad(1)->line_orientation(my_indices[3]));
            line_orientations[4] = ref_cell.standard_vs_true_line_orientation(
              2,
              1,
              orientations[1],
              cell.quad(1)->line_orientation(my_indices[4]));
            line_orientations[5] = ref_cell.standard_vs_true_line_orientation(
              1,
              2,
              orientations[2],
              cell.quad(2)->line_orientation(my_indices[5]));
          }
        else
          // For other shapes (wedges, pyramids), we do not currently
          // implement an optimized function
          for (unsigned int l = 0; l < std::min(12U, cell.n_lines()); ++l)
            line_orientations[l] = cell.line_orientation(l);

        return line_orientations;
      }
    };
  } // namespace TriaAccessorImplementation
} // namespace internal



template <int structdim, int dim, int spacedim>
inline TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *parent,
  const int                           level,
  const int                           index,
  const AccessorData                 *local_data)
  : TriaAccessorBase<structdim, dim, spacedim>(parent, level, index, local_data)
{}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::used() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  return this->objects().used[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::vertex_iterator(
  const unsigned int i) const
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>(this->tria,
                                                      0,
                                                      vertex_index(i));
}



template <int structdim, int dim, int spacedim>
inline ReferenceCell
TriaAccessor<structdim, dim, spacedim>::reference_cell() const
{
  if (structdim == 0)
    return ReferenceCells::Vertex;
  else if (structdim == 1)
    return ReferenceCells::Line;
  else if (structdim == dim)
    return this->tria->levels[this->present_level]
      ->reference_cell[this->present_index];
  else
    return this->tria->faces->get_quad_type(this->present_index);
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::vertex_index(
  const unsigned int corner) const
{
  AssertIndexRange(corner, this->n_vertices());

  if (structdim == dim)
    {
      constexpr unsigned int max_vertices_per_cell = 1 << dim;
      const std::size_t      my_index =
        static_cast<std::size_t>(this->present_index) * max_vertices_per_cell;
      AssertIndexRange(my_index + corner,
                       this->tria->levels[this->present_level]
                         ->cell_vertex_indices_cache.size());
      const unsigned int vertex_index =
        this->tria->levels[this->present_level]
          ->cell_vertex_indices_cache[my_index + corner];
      Assert(vertex_index != numbers::invalid_unsigned_int, ExcInternalError());
      return vertex_index;
    }
  else
    return dealii::internal::TriaAccessorImplementation::Implementation::
      vertex_index(*this, corner);
}



template <int structdim, int dim, int spacedim>
inline Point<spacedim> &
TriaAccessor<structdim, dim, spacedim>::vertex(const unsigned int i) const
{
  return const_cast<Point<spacedim> &>(this->tria->vertices[vertex_index(i)]);
}



template <int structdim, int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::line_iterator
  TriaAccessor<structdim, dim, spacedim>::line(const unsigned int i) const
{
  // checks happen in line_index
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator(this->tria, 0, line_index(i));
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::line_index(const unsigned int i) const
{
  AssertIndexRange(i, this->n_lines());

  return dealii::internal::TriaAccessorImplementation::Implementation::
    line_index(*this, i);
}



template <int structdim, int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::quad_iterator
  TriaAccessor<structdim, dim, spacedim>::quad(const unsigned int i) const
{
  // checks happen in quad_index
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator(this->tria, 0, quad_index(i));
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::quad_index(const unsigned int i) const
{
  return dealii::internal::TriaAccessorImplementation::Implementation::
    quad_index(*this, i);
}



template <int structdim, int dim, int spacedim>
inline unsigned char
TriaAccessor<structdim, dim, spacedim>::combined_face_orientation(
  const unsigned int face) const
{
  return dealii::internal::TriaAccessorImplementation::Implementation::
    combined_face_orientation(*this, face);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_orientation(
  const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::
    face_orientation(*this, face);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_flip(const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::
    face_flip(*this, face);
}


template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_rotation(
  const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::
    face_rotation(*this, face);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::line_orientation(
  const unsigned int line) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(line, this->n_lines());

  return dealii::internal::TriaAccessorImplementation::Implementation::
    line_orientation(*this, line);
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_line_orientation(
  const unsigned int line,
  const bool         value) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(line, this->n_lines());

  dealii::internal::TriaAccessorImplementation::Implementation::
    set_line_orientation(*this, line, value);
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_combined_face_orientation(
  const unsigned int  face,
  const unsigned char combined_orientation) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(face, this->n_faces());

  dealii::internal::TriaAccessorImplementation::Implementation::
    set_combined_face_orientation(*this, face, combined_orientation);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_used_flag() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  this->objects().used[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_used_flag() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  this->objects().used[this->present_index] = false;
}


template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::child_index(const unsigned int i) const
{
  Assert(has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  AssertIndexRange(i, n_children());

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;
  return this->objects().children[n_sets_of_two * this->present_index + i / 2] +
         i % 2;
}



template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::isotropic_child_index(
  const unsigned int i) const
{
  AssertIndexRange(i, GeometryInfo<structdim>::max_children_per_cell);

  switch (structdim)
    {
      case 1:
        return child_index(i);
      case 2:
        {
          const RefinementCase<2> this_refinement_case(
            static_cast<std::uint8_t>(refinement_case()));

          Assert(this_refinement_case != RefinementCase<2>::no_refinement,
                 TriaAccessorExceptions::ExcCellHasNoChildren());

          if (this_refinement_case == RefinementCase<2>::cut_xy)
            return child_index(i);
          else if ((this_refinement_case == RefinementCase<2>::cut_x) &&
                   (child(i % 2)->refinement_case() ==
                    RefinementCase<2>::cut_y))
            return child(i % 2)->child_index(i / 2);
          else if ((this_refinement_case == RefinementCase<2>::cut_y) &&
                   (child(i / 2)->refinement_case() ==
                    RefinementCase<2>::cut_x))
            return child(i / 2)->child_index(i % 2);
          else
            Assert(
              false,
              ExcMessage(
                "This cell has no grandchildren equivalent to isotropic refinement"));
          break;
        }

      case 3:
        Assert(false, ExcNotImplemented());
    }
  return -1;
}



template <int structdim, int dim, int spacedim>
RefinementCase<structdim>
TriaAccessor<structdim, dim, spacedim>::refinement_case() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));

  switch (structdim)
    {
      case 1:
        return (RefinementCase<structdim>(
          this->objects().children[this->present_index] != -1 ?
            // cast the branches
            // here first to uchar
            // and then (above) to
            // RefinementCase<structdim>
            // so that the
            // conversion is valid
            // even for the case
            // structdim>1 (for
            // which this part of
            // the code is dead
            // anyway)
            static_cast<std::uint8_t>(RefinementCase<1>::cut_x) :
            static_cast<std::uint8_t>(RefinementCase<1>::no_refinement)));

      default:
        Assert(static_cast<unsigned int>(this->present_index) <
                 this->objects().refinement_cases.size(),
               ExcIndexRange(this->present_index,
                             0,
                             this->objects().refinement_cases.size()));

        return (static_cast<RefinementCase<structdim>>(
          this->objects().refinement_cases[this->present_index]));
    }
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<structdim, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::child(const unsigned int i) const

{
  // checking of 'i' happens in child_index
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> q(
    this->tria, (dim == structdim ? this->level() + 1 : 0), child_index(i));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &child) const
{
  const auto n_children = this->n_children();
  for (unsigned int child_n = 0; child_n < n_children; ++child_n)
    if (this->child(child_n) == child)
      return child_n;

  Assert(false,
         ExcMessage("The given child is not a child of the current object."));
  return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<structdim, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::isotropic_child(
  const unsigned int i) const
{
  // checking of 'i' happens in child() or
  // child_index() called below
  switch (structdim)
    {
      case 1:
        // no anisotropic refinement in 1d
        return child(i);

      case 2:
        {
          const RefinementCase<2> this_refinement_case(
            static_cast<std::uint8_t>(refinement_case()));

          Assert(this_refinement_case != RefinementCase<2>::no_refinement,
                 TriaAccessorExceptions::ExcCellHasNoChildren());

          if (this_refinement_case == RefinementCase<2>::cut_xy)
            return child(i);
          else if ((this_refinement_case == RefinementCase<2>::cut_x) &&
                   (child(i % 2)->refinement_case() ==
                    RefinementCase<2>::cut_y))
            return child(i % 2)->child(i / 2);
          else if ((this_refinement_case == RefinementCase<2>::cut_y) &&
                   (child(i / 2)->refinement_case() ==
                    RefinementCase<2>::cut_x))
            return child(i / 2)->child(i % 2);
          else
            Assert(
              false,
              ExcMessage(
                "This cell has no grandchildren equivalent to isotropic refinement"));
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
  // we don't get here but have to return
  // something...
  return child(0);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::has_children() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;
  return (this->objects().children[n_sets_of_two * this->present_index] != -1);
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::n_children() const
{
  return GeometryInfo<structdim>::n_children(refinement_case());
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_refinement_case(
  const RefinementCase<structdim> &refinement_case) const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  Assert(static_cast<unsigned int>(this->present_index) <
           this->objects().refinement_cases.size(),
         ExcIndexRange(this->present_index,
                       0,
                       this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index] = refinement_case;
}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::clear_refinement_case() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  Assert(static_cast<unsigned int>(this->present_index) <
           this->objects().refinement_cases.size(),
         ExcIndexRange(this->present_index,
                       0,
                       this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index] =
    RefinementCase<structdim>::no_refinement;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_children(const unsigned int i,
                                                     const int index) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(i % 2 == 0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;

  Assert(
    // clearing the child index for a cell
    (index == -1) ||
      // if setting the child index for the i'th child (with i==0),
      // then the index must be a non-negative number
      (i == 0 && !this->has_children() && (index >= 0)) ||
      // if setting the child index for the i'th child (with i>0),
      // then the previously stored index must be the invalid
      // index
      (i > 0 && this->has_children() && (index >= 0) &&
       this->objects().children[n_sets_of_two * this->present_index + i / 2] ==
         -1),
    TriaAccessorExceptions::ExcCantSetChildren(index));

  this->objects().children[n_sets_of_two * this->present_index + i / 2] = index;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_children() const
{
  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;

  for (unsigned int i = 0; i < n_sets_of_two; ++i)
    set_children(2 * i, -1);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::user_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_flags[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_user_flag() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::clear_user_flag() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = false;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_flag() const
{
  set_user_flag();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_flag();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_flag() const
{
  clear_user_flag();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_flag();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_data() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().clear_user_data(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_user_pointer(void *p) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_pointer() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = nullptr;
}



template <int structdim, int dim, int spacedim>
void *
TriaAccessor<structdim, dim, spacedim>::user_pointer() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_pointer(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_pointer(
  void *p) const
{
  set_user_pointer(p);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_pointer(p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_pointer() const
{
  clear_user_pointer();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_pointer();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_user_index(
  const unsigned int p) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = 0;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::user_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_index(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_index(
  const unsigned int p) const
{
  set_user_index(p);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_index(p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_index() const
{
  clear_user_index();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_index();
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::max_refinement_depth() const
{
  if (!this->has_children())
    return 0;

  unsigned int max_depth = 1;
  for (unsigned int c = 0; c < n_children(); ++c)
    max_depth = std::max(max_depth, child(c)->max_refinement_depth() + 1);
  return max_depth;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_active_descendants() const
{
  if (!this->has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c = 0; c < n_children(); ++c)
        sum += this->child(c)->n_active_descendants();
      return sum;
    }
}



template <int structdim, int dim, int spacedim>
types::boundary_id
TriaAccessor<structdim, dim, spacedim>::boundary_id() const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects()
    .boundary_or_material_id[this->present_index]
    .boundary_id;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_boundary_id(
  const types::boundary_id boundary_ind) const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(boundary_ind != numbers::internal_face_boundary_id,
         ExcMessage("You are trying to set the boundary_id to an invalid "
                    "value (numbers::internal_face_boundary_id is reserved)."));
  Assert(this->at_boundary(),
         ExcMessage("You are trying to set the boundary_id of an "
                    "internal object, which is not allowed!"));

  this->objects().boundary_or_material_id[this->present_index].boundary_id =
    boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_boundary_id_internal(
  const types::boundary_id boundary_ind) const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().boundary_or_material_id[this->present_index].boundary_id =
    boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_all_boundary_ids(
  const types::boundary_id boundary_ind) const
{
  set_boundary_id(boundary_ind);

  switch (structdim)
    {
      case 1:
        // 1d objects have no sub-objects
        // where we have to do anything
        break;

      case 2:
        // for boundary quads also set
        // boundary_id of bounding lines
        for (unsigned int i = 0; i < this->n_lines(); ++i)
          this->line(i)->set_boundary_id(boundary_ind);
        break;

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::at_boundary() const
{
  // error checking is done
  // in boundary_id()
  return (boundary_id() != numbers::internal_face_boundary_id);
}



template <int structdim, int dim, int spacedim>
const Manifold<dim, spacedim> &
TriaAccessor<structdim, dim, spacedim>::get_manifold() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->get_manifold(this->manifold_id());
}


template <int structdim, int dim, int spacedim>
types::manifold_id
TriaAccessor<structdim, dim, spacedim>::manifold_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects().manifold_id[this->present_index];
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_manifold_id(
  const types::manifold_id manifold_ind) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().manifold_id[this->present_index] = manifold_ind;
}


template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_all_manifold_ids(
  const types::manifold_id manifold_ind) const
{
  set_manifold_id(manifold_ind);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->set_all_manifold_ids(manifold_ind);

  switch (structdim)
    {
      case 1:
        if (dim == 1)
          {
            (*this->tria->vertex_to_manifold_id_map_1d)[vertex_index(0)] =
              manifold_ind;
            (*this->tria->vertex_to_manifold_id_map_1d)[vertex_index(1)] =
              manifold_ind;
          }
        break;

      case 2:
        // for quads/simplices also set manifold_id of bounding lines
        for (unsigned int i = 0; i < this->n_lines(); ++i)
          this->line(i)->set_manifold_id(manifold_ind);
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::diameter() const
{
  boost::container::small_vector<Point<spacedim>,
                                 GeometryInfo<structdim>::vertices_per_cell>
    vertices(this->n_vertices());

  for (unsigned int v = 0; v < vertices.size(); ++v)
    vertices[v] = this->vertex(v);

  return internal::TriaAccessorImplementation::diameter<structdim, spacedim>(
    vertices);
}



template <int dim, int spacedim>
double
CellAccessor<dim, spacedim>::diameter(
  const Mapping<dim, spacedim> &mapping) const
{
  return internal::TriaAccessorImplementation::diameter<dim, spacedim>(
    mapping.get_vertices(typename Triangulation<dim, spacedim>::cell_iterator(
      this->tria, this->level(), this->index())));
}



template <int structdim, int dim, int spacedim>
std::pair<Point<spacedim>, double>
TriaAccessor<structdim, dim, spacedim>::enclosing_ball() const
{
  // If the object is one dimensional,
  // the enclosing ball is the initial iterate
  // i.e., the ball's center and diameter are
  // the center and the diameter of the object.
  if (structdim == 1)
    return std::make_pair((this->vertex(1) + this->vertex(0)) * 0.5,
                          (this->vertex(1) - this->vertex(0)).norm() * 0.5);

  // The list is_initial_guess_vertex contains bool values and has
  // the same size as the number of vertices per object.
  // The entries of is_initial_guess_vertex are set true only for those
  // two vertices corresponding to the largest diagonal which is being used
  // to construct the initial ball.
  // We employ this mask to skip these two vertices while enlarging the ball.
  std::vector<bool> is_initial_guess_vertex(this->n_vertices());

  // First let all the vertices be outside
  std::fill(is_initial_guess_vertex.begin(),
            is_initial_guess_vertex.end(),
            false);

  // Get an initial guess by looking at the largest diagonal
  Point<spacedim> center;
  double          radius = 0;

  switch (structdim)
    {
      case 2:
        {
          const Point<spacedim> p30(this->vertex(3) - this->vertex(0));
          const Point<spacedim> p21(this->vertex(2) - this->vertex(1));
          if (p30.norm() > p21.norm())
            {
              center                     = this->vertex(0) + 0.5 * p30;
              radius                     = p30.norm() / 2.;
              is_initial_guess_vertex[3] = true;
              is_initial_guess_vertex[0] = true;
            }
          else
            {
              center                     = this->vertex(1) + 0.5 * p21;
              radius                     = p21.norm() / 2.;
              is_initial_guess_vertex[2] = true;
              is_initial_guess_vertex[1] = true;
            }
          break;
        }
      case 3:
        {
          const Point<spacedim>     p70(this->vertex(7) - this->vertex(0));
          const Point<spacedim>     p61(this->vertex(6) - this->vertex(1));
          const Point<spacedim>     p25(this->vertex(2) - this->vertex(5));
          const Point<spacedim>     p34(this->vertex(3) - this->vertex(4));
          const std::vector<double> diagonals = {p70.norm(),
                                                 p61.norm(),
                                                 p25.norm(),
                                                 p34.norm()};
          const std::vector<double>::const_iterator it =
            std::max_element(diagonals.begin(), diagonals.end());
          if (it == diagonals.begin())
            {
              center                     = this->vertex(0) + 0.5 * p70;
              is_initial_guess_vertex[7] = true;
              is_initial_guess_vertex[0] = true;
            }
          else if (it == diagonals.begin() + 1)
            {
              center                     = this->vertex(1) + 0.5 * p61;
              is_initial_guess_vertex[6] = true;
              is_initial_guess_vertex[1] = true;
            }
          else if (it == diagonals.begin() + 2)
            {
              center                     = this->vertex(5) + 0.5 * p25;
              is_initial_guess_vertex[2] = true;
              is_initial_guess_vertex[5] = true;
            }
          else
            {
              center                     = this->vertex(4) + 0.5 * p34;
              is_initial_guess_vertex[3] = true;
              is_initial_guess_vertex[4] = true;
            }
          radius = *it * 0.5;
          break;
        }
      default:
        Assert(false, ExcNotImplemented());
        return std::pair<Point<spacedim>, double>();
    }

  // For each vertex that is found to be geometrically outside the ball
  // enlarge the ball  so that the new ball contains both the previous ball
  // and the given vertex.
  for (const unsigned int v : this->vertex_indices())
    if (!is_initial_guess_vertex[v])
      {
        const double distance = center.distance(this->vertex(v));
        if (distance > radius)
          {
            // we found a vertex which is outside of the ball
            // extend it (move center and change radius)
            const Point<spacedim> pCV(center - this->vertex(v));
            radius = (distance + radius) * 0.5;
            center = this->vertex(v) + pCV * (radius / distance);

            // Now the new ball constructed in this block
            // encloses the vertex (v) that was found to be geometrically
            // outside the old ball.
          }
      }
#ifdef DEBUG
  bool all_vertices_within_ball = true;

  // Set all_vertices_within_ball false if any of the vertices of the object
  // are geometrically outside the ball
  for (const unsigned int v : this->vertex_indices())
    if (center.distance(this->vertex(v)) >
        radius + 100. * std::numeric_limits<double>::epsilon())
      {
        all_vertices_within_ball = false;
        break;
      }
  // If all the vertices are not within the ball throw error
  Assert(all_vertices_within_ball, ExcInternalError());
#endif
  return std::make_pair(center, radius);
}


template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::minimum_vertex_distance() const
{
  switch (structdim)
    {
      case 1:
        return (this->vertex(1) - this->vertex(0)).norm();
      case 2:
      case 3:
        {
          double min = std::numeric_limits<double>::max();
          for (const unsigned int i : this->vertex_indices())
            for (unsigned int j = i + 1; j < this->n_vertices(); ++j)
              min = std::min(min,
                             (this->vertex(i) - this->vertex(j)) *
                               (this->vertex(i) - this->vertex(j)));
          return std::sqrt(min);
        }
      default:
        Assert(false, ExcNotImplemented());
        return -1e10;
    }
}


template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::is_translation_of(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &o) const
{
  // go through the vertices and check... The
  // cell is a translation of the previous
  // one in case the distance between the
  // individual vertices in the two cell is
  // the same for all the vertices. So do the
  // check by first getting the distance on
  // the first vertex, and then checking
  // whether all others have the same down to
  // rounding errors (we have to be careful
  // here because the calculation of the
  // displacement between one cell and the
  // next can already result in the loss of
  // one or two digits), so we choose 1e-12
  // times the distance between the zeroth
  // vertices here.
  bool                      is_translation = true;
  const Tensor<1, spacedim> dist           = o->vertex(0) - this->vertex(0);
  const double              tol_square     = 1e-24 * dist.norm_square();
  for (unsigned int i = 1; i < this->n_vertices(); ++i)
    {
      const Tensor<1, spacedim> dist_new =
        (o->vertex(i) - this->vertex(i)) - dist;
      if (dist_new.norm_square() > tol_square)
        {
          is_translation = false;
          break;
        }
    }
  return is_translation;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_vertices() const
{
  return this->reference_cell().n_vertices();
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_lines() const
{
  return this->reference_cell().n_lines();
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_faces() const
{
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));

  return this->reference_cell().n_faces();
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::vertex_indices() const
{
  return {0U, n_vertices()};
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::line_indices() const
{
  return {0U, n_lines()};
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::face_indices() const
{
  return {0U, n_faces()};
}



/*----------------- Functions: TriaAccessor<0,dim,spacedim> -----------------*/

template <int dim, int spacedim>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *tria,
  const unsigned int                  vertex_index)
  : tria(tria)
  , global_vertex_index(vertex_index)
{}



template <int dim, int spacedim>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *tria,
  const int /*level*/,
  const int index,
  const AccessorData *)
  : tria(tria)
  , global_vertex_index(index)
{}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::copy_from(const TriaAccessor &t)
{
  tria                = t.tria;
  global_vertex_index = t.global_vertex_index;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator<(
  const TriaAccessor<0, dim, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  return (global_vertex_index < other.global_vertex_index);
}



template <int dim, int spacedim>
inline IteratorState::IteratorStates
TriaAccessor<0, dim, spacedim>::state() const
{
  if (global_vertex_index != numbers::invalid_unsigned_int)
    return IteratorState::valid;
  else
    return IteratorState::past_the_end;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::level()
{
  return 0;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::index() const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
TriaAccessor<0, dim, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::operator++()
{
  ++global_vertex_index;
  if (global_vertex_index >= tria->n_vertices())
    global_vertex_index = numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::operator--()
{
  if (global_vertex_index != numbers::invalid_unsigned_int)
    {
      if (global_vertex_index != 0)
        --global_vertex_index;
      else
        global_vertex_index = numbers::invalid_unsigned_int;
    }
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator==(const TriaAccessor &t) const
{
  const bool result =
    ((tria == t.tria) && (global_vertex_index == t.global_vertex_index));

  return result;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator!=(const TriaAccessor &t) const
{
  return !(*this == t);
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::vertex_index(const unsigned int) const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline Point<spacedim> &
TriaAccessor<0, dim, spacedim>::vertex(const unsigned int) const
{
  return const_cast<Point<spacedim> &>(
    this->tria->vertices[global_vertex_index]);
}



template <int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::line_iterator
  TriaAccessor<0, dim, spacedim>::line(const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator();
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::line_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::quad_iterator
  TriaAccessor<0, dim, spacedim>::quad(const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator();
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::quad_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::diameter() const
{
  return 0.;
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::extent_in_direction(const unsigned int) const
{
  return 0.;
}



template <int dim, int spacedim>
inline Point<spacedim>
TriaAccessor<0, dim, spacedim>::center(const bool, const bool) const
{
  return this->tria->vertices[global_vertex_index];
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::measure() const
{
  return 0.;
}



template <int dim, int spacedim>
inline unsigned char
TriaAccessor<0, dim, spacedim>::combined_face_orientation(
  const unsigned int /*face*/)
{
  return 0;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_orientation(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_flip(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_rotation(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::line_orientation(const unsigned int /*line*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::has_children()
{
  return false;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::n_children()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::n_active_descendants()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::max_refinement_depth()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<0, dim, spacedim>> &)
{
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<0, dim, spacedim>::child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>();
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<0, dim, spacedim>::isotropic_child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>();
}



template <int dim, int spacedim>
inline RefinementCase<0>
TriaAccessor<0, dim, spacedim>::refinement_case()
{
  return {RefinementPossibilities<0>::no_refinement};
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::child_index(const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::isotropic_child_index(const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::used() const
{
  return tria->vertex_used(global_vertex_index);
}



/*------------------- Functions: TriaAccessor<0,1,spacedim> -----------------*/

template <int spacedim>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const Triangulation<1, spacedim> *tria,
  const VertexKind                  vertex_kind,
  const unsigned int                vertex_index)
  : tria(tria)
  , vertex_kind(vertex_kind)
  , global_vertex_index(vertex_index)
{}



template <int spacedim>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const Triangulation<1, spacedim> *tria,
  const int                         level,
  const int                         index,
  const AccessorData *)
  : tria(tria)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  // in general, calling this constructor should yield an error -- users should
  // instead call the one immediately above. however, if you create something
  // like Triangulation<1>::face_iterator() then this calls the default
  // constructor of the iterator which calls the accessor with argument list
  // (0,-2,-2,0), so in this particular case accept this call and create an
  // object that corresponds to the default constructed (invalid) vertex
  // accessor
  (void)level;
  (void)index;
  Assert((level == -2) && (index == -2), ExcInternalError());
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::copy_from(const TriaAccessor &t)
{
  tria                = t.tria;
  vertex_kind         = t.vertex_kind;
  global_vertex_index = t.global_vertex_index;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::copy_from(
  const TriaAccessorBase<0, 1, spacedim> &)
{
  // We cannot convert from TriaAccessorBase to
  // TriaAccessor<0,1,spacedim> because the latter is not derived from
  // the former. We should never get here.
  Assert(false, ExcInternalError());
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator<(
  const TriaAccessor<0, 1, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  return (global_vertex_index < other.global_vertex_index);
}



template <int spacedim>
inline IteratorState::IteratorStates
TriaAccessor<0, 1, spacedim>::state()
{
  return IteratorState::valid;
}


template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::level()
{
  return 0;
}



template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::index() const
{
  return global_vertex_index;
}



template <int spacedim>
inline const Triangulation<1, spacedim> &
TriaAccessor<0, 1, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::operator++() const
{
  Assert(false, ExcNotImplemented());
}


template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::operator--() const
{
  Assert(false, ExcNotImplemented());
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator==(const TriaAccessor &t) const
{
  const bool result =
    ((tria == t.tria) && (global_vertex_index == t.global_vertex_index));
  // if we point to the same vertex,
  // make sure we know the same about
  // it
  if (result == true)
    Assert(vertex_kind == t.vertex_kind, ExcInternalError());

  return result;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator!=(const TriaAccessor &t) const
{
  return !(*this == t);
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::vertex_index(const unsigned int i) const
{
  AssertIndexRange(i, 1);
  (void)i;
  return global_vertex_index;
}



template <int spacedim>
inline Point<spacedim> &
TriaAccessor<0, 1, spacedim>::vertex(const unsigned int i) const
{
  AssertIndexRange(i, 1);
  (void)i;
  return const_cast<Point<spacedim> &>(
    this->tria->vertices[global_vertex_index]);
}



template <int spacedim>
inline Point<spacedim>
TriaAccessor<0, 1, spacedim>::center() const
{
  return this->tria->vertices[global_vertex_index];
}



template <int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<1, spacedim>::line_iterator
  TriaAccessor<0, 1, spacedim>::line(const unsigned int)
{
  return {};
}


template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::line_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}


template <int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<1, spacedim>::quad_iterator
  TriaAccessor<0, 1, spacedim>::quad(const unsigned int)
{
  return {};
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::quad_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}


template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::at_boundary() const
{
  return vertex_kind != interior_vertex;
}


template <int spacedim>
inline types::boundary_id
TriaAccessor<0, 1, spacedim>::boundary_id() const
{
  switch (vertex_kind)
    {
      case left_vertex:
      case right_vertex:
        {
          Assert(tria->vertex_to_boundary_id_map_1d->find(
                   this->vertex_index()) !=
                   tria->vertex_to_boundary_id_map_1d->end(),
                 ExcInternalError());

          return (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()];
        }

      default:
        return numbers::internal_face_boundary_id;
    }
}



template <int spacedim>
inline const Manifold<1, spacedim> &
TriaAccessor<0, 1, spacedim>::get_manifold() const
{
  return this->tria->get_manifold(this->manifold_id());
}



template <int spacedim>
inline types::manifold_id
TriaAccessor<0, 1, spacedim>::manifold_id() const
{
  if (tria->vertex_to_manifold_id_map_1d->find(this->vertex_index()) !=
      tria->vertex_to_manifold_id_map_1d->end())
    return (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()];
  else
    return numbers::flat_manifold_id;
}


template <int spacedim>
inline unsigned char
TriaAccessor<0, 1, spacedim>::combined_face_orientation(
  const unsigned int /*face*/)
{
  return 0;
}


template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_orientation(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_flip(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_rotation(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::line_orientation(const unsigned int /*line*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::has_children()
{
  return false;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::n_children()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::n_active_descendants()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::max_refinement_depth()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<0, 1, spacedim>> &)
{
  return numbers::invalid_unsigned_int;
}



template <int spacedim>
inline TriaIterator<TriaAccessor<0, 1, spacedim>>
TriaAccessor<0, 1, spacedim>::child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, 1, spacedim>>();
}


template <int spacedim>
inline TriaIterator<TriaAccessor<0, 1, spacedim>>
TriaAccessor<0, 1, spacedim>::isotropic_child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, 1, spacedim>>();
}


template <int spacedim>
inline RefinementCase<0>
TriaAccessor<0, 1, spacedim>::refinement_case()
{
  return {RefinementPossibilities<0>::no_refinement};
}

template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::child_index(const unsigned int)
{
  return -1;
}


template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::isotropic_child_index(const unsigned int)
{
  return -1;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_boundary_id(const types::boundary_id b) const
{
  Assert(tria->vertex_to_boundary_id_map_1d->find(this->vertex_index()) !=
           tria->vertex_to_boundary_id_map_1d->end(),
         ExcMessage("You can't set the boundary_id of a face of a cell that is "
                    "not actually at the boundary."));

  (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_manifold_id(const types::manifold_id b)
{
  (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_all_boundary_ids(
  const types::boundary_id b) const
{
  set_boundary_id(b);
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_all_manifold_ids(const types::manifold_id b)
{
  set_manifold_id(b);
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::used() const
{
  return tria->vertex_used(global_vertex_index);
}



template <int spacedim>
inline ReferenceCell
TriaAccessor<0, 1, spacedim>::reference_cell() const
{
  return ReferenceCells::Vertex;
}



template <int spacedim>
unsigned int
TriaAccessor<0, 1, spacedim>::n_vertices() const
{
  return 1;
}



template <int spacedim>
unsigned int
TriaAccessor<0, 1, spacedim>::n_lines() const
{
  return 0;
}



template <int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<0, 1, spacedim>::vertex_indices() const
{
  return {0U, n_vertices()};
}



template <int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<0, 1, spacedim>::line_indices() const
{
  return {0U, n_lines()};
}

/*------------------ Functions: CellAccessor<dim,spacedim> ------------------*/


template <int dim, int spacedim>
inline CellAccessor<dim, spacedim>::CellAccessor(
  const Triangulation<dim, spacedim> *parent,
  const int                           level,
  const int                           index,
  const AccessorData                 *local_data)
  : TriaAccessor<dim, dim, spacedim>(parent, level, index, local_data)
{}



template <int dim, int spacedim>
inline CellAccessor<dim, spacedim>::CellAccessor(
  const TriaAccessor<dim, dim, spacedim> &cell_accessor)
  : TriaAccessor<dim, dim, spacedim>(
      static_cast<const TriaAccessor<dim, dim, spacedim> &>(cell_accessor))
{}



namespace internal
{
  namespace CellAccessorImplementation
  {
    template <int spacedim>
    inline dealii::TriaIterator<dealii::TriaAccessor<0, 1, spacedim>>
    get_face(const dealii::CellAccessor<1, spacedim> &cell,
             const unsigned int                       i)
    {
      dealii::TriaAccessor<0, 1, spacedim> a(
        &cell.get_triangulation(),
        ((i == 0) && cell.at_boundary(0) ?
           dealii::TriaAccessor<0, 1, spacedim>::left_vertex :
           ((i == 1) && cell.at_boundary(1) ?
              dealii::TriaAccessor<0, 1, spacedim>::right_vertex :
              dealii::TriaAccessor<0, 1, spacedim>::interior_vertex)),
        dealii::internal::TriaAccessorImplementation::Implementation::
          vertex_index(cell, i));
      return dealii::TriaIterator<dealii::TriaAccessor<0, 1, spacedim>>(a);
    }


    template <int spacedim>
    inline dealii::TriaIterator<dealii::TriaAccessor<1, 2, spacedim>>
    get_face(const dealii::CellAccessor<2, spacedim> &cell,
             const unsigned int                       i)
    {
      return cell.line(i);
    }


    template <int spacedim>
    inline dealii::TriaIterator<dealii::TriaAccessor<2, 3, spacedim>>
    get_face(const dealii::CellAccessor<3, spacedim> &cell,
             const unsigned int                       i)
    {
      return cell.quad(i);
    }
  } // namespace CellAccessorImplementation
} // namespace internal



template <int dim, int spacedim>
inline TriaIterator<CellAccessor<dim, spacedim>>
CellAccessor<dim, spacedim>::child(const unsigned int i) const
{
  TriaIterator<CellAccessor<dim, spacedim>> q(this->tria,
                                              this->present_level + 1,
                                              this->child_index(i));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                      GeometryInfo<dim>::max_children_per_cell>
CellAccessor<dim, spacedim>::child_iterators() const
{
  boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                 GeometryInfo<dim>::max_children_per_cell>
    child_iterators(this->n_children());

  for (unsigned int i = 0; i < this->n_children(); ++i)
    child_iterators[i] = this->child(i);

  return child_iterators;
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>
CellAccessor<dim, spacedim>::face(const unsigned int i) const
{
  AssertIndexRange(i, this->n_faces());
  return dealii::internal::CellAccessorImplementation::get_face(*this, i);
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::face_iterator_to_index(
  const TriaIterator<TriaAccessor<dim - 1, dim, spacedim>> &face) const
{
  for (const unsigned int face_n : this->face_indices())
    if (this->face(face_n) == face)
      return face_n;

  Assert(false,
         ExcMessage("The given face is not a face of the current cell."));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline boost::container::small_vector<
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
  GeometryInfo<dim>::faces_per_cell>
CellAccessor<dim, spacedim>::face_iterators() const
{
  boost::container::small_vector<
    TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
    GeometryInfo<dim>::faces_per_cell>
    face_iterators(this->n_faces());

  for (const unsigned int i : this->face_indices())
    face_iterators[i] =
      dealii::internal::CellAccessorImplementation::get_face(*this, i);

  return face_iterators;
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::face_index(const unsigned int i) const
{
  switch (dim)
    {
      case 1:
        {
          return this->vertex_index(i);
        }

      case 2:
        return this->line_index(i);

      case 3:
        return this->quad_index(i);

      default:
        return numbers::invalid_unsigned_int;
    }
}



template <int dim, int spacedim>
inline int
CellAccessor<dim, spacedim>::neighbor_index(const unsigned int face_no) const
{
  AssertIndexRange(face_no, this->n_faces());
  return this->tria->levels[this->present_level]
    ->neighbors[this->present_index * GeometryInfo<dim>::faces_per_cell +
                face_no]
    .second;
}



template <int dim, int spacedim>
inline int
CellAccessor<dim, spacedim>::neighbor_level(const unsigned int face_no) const
{
  AssertIndexRange(face_no, this->n_faces());
  return this->tria->levels[this->present_level]
    ->neighbors[this->present_index * GeometryInfo<dim>::faces_per_cell +
                face_no]
    .first;
}



template <int dim, int spacedim>
inline RefinementCase<dim>
CellAccessor<dim, spacedim>::refine_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for refinement must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert(this->is_active() || !this->tria->levels[this->present_level]
                                 ->refine_flags[this->present_index],
         ExcRefineCellNotActive());
  return RefinementCase<dim>(
    this->tria->levels[this->present_level]->refine_flags[this->present_index]);
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::set_refine_flag(
  const RefinementCase<dim> refinement_case) const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  Assert(!coarsen_flag_set(), ExcCellFlaggedForCoarsening());

  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    refinement_case;
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::clear_refine_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    RefinementCase<dim>::no_refinement;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::flag_for_face_refinement(
  const unsigned int             face_no,
  const RefinementCase<dim - 1> &face_refinement_case) const
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(face_no, this->n_faces());
  AssertIndexRange(face_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);

  // the new refinement case is a combination
  // of the minimum required one for the given
  // face refinement and the already existing
  // flagged refinement case
  RefinementCase<dim> old_ref_case = refine_flag_set();
  RefinementCase<dim> new_ref_case =
    (old_ref_case |
     GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(
       face_refinement_case,
       face_no,
       this->face_orientation(face_no),
       this->face_flip(face_no),
       this->face_rotation(face_no)));
  set_refine_flag(new_ref_case);
  // return, whether we had to change the
  // refinement flag
  return new_ref_case != old_ref_case;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::flag_for_line_refinement(
  const unsigned int line_no) const
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(line_no, this->n_lines());

  // the new refinement case is a combination
  // of the minimum required one for the given
  // line refinement and the already existing
  // flagged refinement case
  RefinementCase<dim>
    old_ref_case = refine_flag_set(),
    new_ref_case =
      old_ref_case |
      GeometryInfo<dim>::min_cell_refinement_case_for_line_refinement(line_no);
  set_refine_flag(new_ref_case);
  // return, whether we had to change the
  // refinement flag
  return new_ref_case != old_ref_case;
}



template <>
inline dealii::internal::SubfaceCase<1>
CellAccessor<1>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}

template <>
inline dealii::internal::SubfaceCase<1>
CellAccessor<1, 2>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}


template <>
inline dealii::internal::SubfaceCase<1>
CellAccessor<1, 3>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}


template <>
inline dealii::internal::SubfaceCase<2>
CellAccessor<2>::subface_case(const unsigned int face_no) const
{
  Assert(is_active(), TriaAccessorExceptions::ExcCellNotActive());
  AssertIndexRange(face_no, this->n_faces());
  return ((face(face_no)->has_children()) ?
            dealii::internal::SubfaceCase<2>::case_x :
            dealii::internal::SubfaceCase<2>::case_none);
}

template <>
inline dealii::internal::SubfaceCase<2>
CellAccessor<2, 3>::subface_case(const unsigned int face_no) const
{
  Assert(is_active(), TriaAccessorExceptions::ExcCellNotActive());
  AssertIndexRange(face_no, this->n_faces());
  return ((face(face_no)->has_children()) ?
            dealii::internal::SubfaceCase<2>::case_x :
            dealii::internal::SubfaceCase<2>::case_none);
}


template <>
inline dealii::internal::SubfaceCase<3>
CellAccessor<3>::subface_case(const unsigned int face_no) const
{
  Assert(is_active(), TriaAccessorExceptions::ExcCellNotActive());
  AssertIndexRange(face_no, this->n_faces());
  switch (static_cast<std::uint8_t>(face(face_no)->refinement_case()))
    {
      case RefinementCase<3>::no_refinement:
        return dealii::internal::SubfaceCase<3>::case_none;
      case RefinementCase<3>::cut_x:
        if (face(face_no)->child(0)->has_children())
          {
            Assert(face(face_no)->child(0)->refinement_case() ==
                     RefinementCase<2>::cut_y,
                   ExcInternalError());
            if (face(face_no)->child(1)->has_children())
              {
                Assert(face(face_no)->child(1)->refinement_case() ==
                         RefinementCase<2>::cut_y,
                       ExcInternalError());
                return dealii::internal::SubfaceCase<3>::case_x1y2y;
              }
            else
              return dealii::internal::SubfaceCase<3>::case_x1y;
          }
        else
          {
            if (face(face_no)->child(1)->has_children())
              {
                Assert(face(face_no)->child(1)->refinement_case() ==
                         RefinementCase<2>::cut_y,
                       ExcInternalError());
                return dealii::internal::SubfaceCase<3>::case_x2y;
              }
            else
              return dealii::internal::SubfaceCase<3>::case_x;
          }
      case RefinementCase<3>::cut_y:
        if (face(face_no)->child(0)->has_children())
          {
            Assert(face(face_no)->child(0)->refinement_case() ==
                     RefinementCase<2>::cut_x,
                   ExcInternalError());
            if (face(face_no)->child(1)->has_children())
              {
                Assert(face(face_no)->child(1)->refinement_case() ==
                         RefinementCase<2>::cut_x,
                       ExcInternalError());
                return dealii::internal::SubfaceCase<3>::case_y1x2x;
              }
            else
              return dealii::internal::SubfaceCase<3>::case_y1x;
          }
        else
          {
            if (face(face_no)->child(1)->has_children())
              {
                Assert(face(face_no)->child(1)->refinement_case() ==
                         RefinementCase<2>::cut_x,
                       ExcInternalError());
                return dealii::internal::SubfaceCase<3>::case_y2x;
              }
            else
              return dealii::internal::SubfaceCase<3>::case_y;
          }
      case RefinementCase<3>::cut_xy:
        return dealii::internal::SubfaceCase<3>::case_xy;
      default:
        Assert(false, ExcInternalError());
    }
  // we should never get here
  return dealii::internal::SubfaceCase<3>::case_none;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::coarsen_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for coarsening must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert(this->is_active() || !this->tria->levels[this->present_level]
                                 ->coarsen_flags[this->present_index],
         ExcRefineCellNotActive());
  return this->tria->levels[this->present_level]
    ->coarsen_flags[this->present_index];
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::set_coarsen_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  Assert(!refine_flag_set(), ExcCellFlaggedForRefinement());

  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] =
    true;
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::clear_coarsen_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] =
    false;
}



template <int dim, int spacedim>
inline TriaIterator<CellAccessor<dim, spacedim>>
CellAccessor<dim, spacedim>::neighbor(const unsigned int face_no) const
{
  TriaIterator<CellAccessor<dim, spacedim>> q(this->tria,
                                              neighbor_level(face_no),
                                              neighbor_index(face_no));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_active() const
{
  return !this->has_children();
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_locally_owned() const
{
  Assert(this->is_active(),
         ExcMessage("is_locally_owned() can only be called on active cells!"));
#ifndef DEAL_II_WITH_MPI
  return true;
#else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  return (locally_owned_subdomain == numbers::invalid_subdomain_id ||
          this->subdomain_id() == locally_owned_subdomain);

#endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_locally_owned_on_level() const
{
#ifndef DEAL_II_WITH_MPI
  return true;
#else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  return (locally_owned_subdomain == numbers::invalid_subdomain_id ||
          this->level_subdomain_id() == locally_owned_subdomain);

#endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_ghost() const
{
  Assert(this->is_active(),
         ExcMessage("is_ghost() can only be called on active cells!"));
  if (this->has_children())
    return false;

#ifndef DEAL_II_WITH_MPI
  return false;
#else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition rules out that case as all cells to a
  // serial triangulation are locally owned and none is ghosted. The second
  // and third conditions check whether the cell's subdomain is not the
  // locally owned one and not artificial.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  const types::subdomain_id subdomain_id = this->subdomain_id();
  return (locally_owned_subdomain != numbers::invalid_subdomain_id &&
          subdomain_id != locally_owned_subdomain &&
          subdomain_id != numbers::artificial_subdomain_id);

#endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_ghost_on_level() const
{
#ifndef DEAL_II_WITH_MPI
  return false;
#else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  const types::subdomain_id subdomain_id = this->level_subdomain_id();
  return (locally_owned_subdomain != numbers::invalid_subdomain_id &&
          subdomain_id != locally_owned_subdomain &&
          subdomain_id != numbers::artificial_subdomain_id);

#endif
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_artificial() const
{
  Assert(this->is_active(),
         ExcMessage("is_artificial() can only be called on active cells!"));
#ifndef DEAL_II_WITH_MPI
  return false;
#else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition rules out that case as all cells to a
  // serial triangulation are locally owned and none is artificial.
  return (this->tria->locally_owned_subdomain() !=
            numbers::invalid_subdomain_id &&
          this->subdomain_id() == numbers::artificial_subdomain_id);

#endif
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_artificial_on_level() const
{
#ifndef DEAL_II_WITH_MPI
  return false;
#else
  return (this->tria->locally_owned_subdomain() !=
            numbers::invalid_subdomain_id &&
          this->level_subdomain_id() == numbers::artificial_subdomain_id);
#endif
}



template <int dim, int spacedim>
inline types::subdomain_id
CellAccessor<dim, spacedim>::subdomain_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(this->is_active(),
         ExcMessage("subdomain_id() can only be called on active cells!"));
  return this->tria->levels[this->present_level]
    ->subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
inline types::subdomain_id
CellAccessor<dim, spacedim>::level_subdomain_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]
    ->level_subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::neighbor_face_no(const unsigned int neighbor) const
{
  const unsigned int n2 = neighbor_of_neighbor_internal(neighbor);
  if (n2 != numbers::invalid_unsigned_int)
    // return this value as the
    // neighbor is not coarser
    return n2;
  else
    // the neighbor is coarser
    return neighbor_of_coarser_neighbor(neighbor).first;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_level_cell()
{
  return false;
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::active_cell_index() const
{
  Assert(this->is_active(), TriaAccessorExceptions::ExcCellNotActive());
  return this->tria->levels[this->present_level]
    ->active_cell_indices[this->present_index];
}



template <int dim, int spacedim>
inline types::global_cell_index
CellAccessor<dim, spacedim>::global_active_cell_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(this->is_active(),
         ExcMessage(
           "global_active_cell_index() can only be called on active cells!"));

  return this->tria->levels[this->present_level]
    ->global_active_cell_indices[this->present_index];
}



template <int dim, int spacedim>
inline types::global_cell_index
CellAccessor<dim, spacedim>::global_level_cell_index() const
{
  return this->tria->levels[this->present_level]
    ->global_level_cell_indices[this->present_index];
}


DEAL_II_NAMESPACE_CLOSE

#endif
