// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2017 by the deal.II authors
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

#ifndef dealii_tria_accessor_templates_h
#define dealii_tria_accessor_templates_h


#include <deal.II/base/config.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_levels.h>
#include <deal.II/grid/tria_faces.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/distributed/tria_base.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  template <int, int> class Triangulation;
}


/*------------------------ Functions: TriaAccessorBase ---------------------------*/

template <int structdim, int dim, int spacedim>
inline
TriaAccessorBase<structdim,dim,spacedim>::TriaAccessorBase (
  const Triangulation<dim,spacedim> *tria,
  const int                          level,
  const int                          index,
  const AccessorData *)
  :
  present_level((structdim==dim) ? level : 0),
  present_index (index),
  tria (tria)
{

  // non-cells have no level, so a 0
  // should have been passed, or a -1
  // for an end-iterator, or -2 for
  // an invalid (default constructed)
  // iterator
  if (structdim != dim)
    {
      Assert ((level == 0) || (level == -1) || (level == -2),
              ExcInternalError());
    }
}


template <int structdim, int dim, int spacedim>
inline
TriaAccessorBase<structdim,dim,spacedim>::TriaAccessorBase (const TriaAccessorBase<structdim,dim,spacedim> &a)
  :
  present_level(a.present_level),
  present_index(a.present_index),
  tria(a.tria)
{}


template <int structdim, int dim, int spacedim>
inline
void
TriaAccessorBase<structdim,dim,spacedim>::copy_from (const TriaAccessorBase<structdim,dim,spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria = a.tria;

  if (structdim != dim)
    {
      Assert ((present_level == 0) || (present_level == -1) || (present_level == -2),
              ExcInternalError());
    }
}



template <int structdim, int dim, int spacedim>
inline
TriaAccessorBase<structdim,dim,spacedim> &
TriaAccessorBase<structdim,dim,spacedim>::operator= (const TriaAccessorBase<structdim,dim,spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria = a.tria;

  if (structdim != dim)
    {
      Assert ((present_level == 0) || (present_level == -1) || (present_level == -2),
              ExcInternalError());
    }
  return *this;
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessorBase<structdim,dim,spacedim>::operator == (const TriaAccessorBase<structdim,dim,spacedim> &a) const
{
  Assert (tria == a.tria || tria == nullptr || a.tria == nullptr,
          TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria == a.tria) &&
          (present_level == a.present_level) &&
          (present_index == a.present_index));
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessorBase<structdim,dim,spacedim>::operator != (const TriaAccessorBase<structdim,dim,spacedim> &a) const
{
  Assert (tria == a.tria || tria == nullptr || a.tria == nullptr,
          TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria != a.tria) ||
          (present_level != a.present_level) ||
          (present_index != a.present_index));
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessorBase<structdim,dim,spacedim>::operator < (const TriaAccessorBase<structdim,dim,spacedim> &other) const
{
  Assert (tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  if (present_level != other.present_level)
    return (present_level < other.present_level);

  return (present_index < other.present_index);

}



template <int structdim, int dim, int spacedim>
inline
int
TriaAccessorBase<structdim,dim,spacedim>::level () const
{
  // This is always zero or invalid
  // if the object is not a cell
  return present_level;
}



template <int structdim, int dim, int spacedim>
inline
int
TriaAccessorBase<structdim,dim,spacedim>::index () const
{
  return present_index;
}



template <int structdim, int dim, int spacedim>
inline
IteratorState::IteratorStates
TriaAccessorBase<structdim,dim,spacedim>::state () const
{
  if ((present_level>=0) && (present_index>=0))
    return IteratorState::valid;
  else if (present_index==-1)
    return IteratorState::past_the_end;
  else
    return IteratorState::invalid;
}



template <int structdim, int dim, int spacedim>
inline
const Triangulation<dim,spacedim> &
TriaAccessorBase<structdim,dim,spacedim>::get_triangulation () const
{
  return *tria;
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessorBase<structdim,dim,spacedim>::operator ++ ()
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
      if (this->present_index
          >=
          static_cast<int>(objects().cells.size()))
        this->present_index = -1;
    }
  else
    {
      while (this->present_index
             >=
             static_cast<int>(this->tria->levels[this->present_level]->cells.cells.size()))
        {
          // no -> go one level up until we find
          // one with more than zero cells
          ++this->present_level;
          this->present_index = 0;
          // highest level reached?
          if (this->present_level >= static_cast<int>(this->tria->levels.size()))
            {
              // return with past the end pointer
              this->present_level = this->present_index = -1;
              return;
            }
        }
    }
}


template <int structdim, int dim, int spacedim>
inline
void
TriaAccessorBase<structdim,dim,spacedim>::operator -- ()
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
          this->present_index = this->tria->levels[this->present_level]->cells.cells.size()-1;
        }
    }
}


namespace internal
{
  namespace TriaAccessorBaseImplementation
  {
    /**
     * Out of a face object, get the sub-objects of dimensionality given by
     * the last argument.
     */
    template <int dim>
    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<1> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<dim> *faces,
                 const std::integral_constant<int, 1>)
    {
      return &faces->lines;
    }


    template <int dim>
    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<2> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<dim> *faces,
                 const std::integral_constant<int, 2>)
    {
      return &faces->quads;
    }

    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<1> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<1> *,
                 const std::integral_constant<int, 1>)
    {
      Assert (false, ExcInternalError());
      return nullptr;
    }

    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<2> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<2> *,
                 const std::integral_constant<int, 2>)
    {
      Assert (false, ExcInternalError());
      return nullptr;
    }

    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<3> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<3> *,
                 const std::integral_constant<int, 3>)
    {
      Assert (false, ExcInternalError());
      return nullptr;
    }

    /**
     * This function should never be used, but we need it for the template
     * instantiation of TriaAccessorBase<dim,dim,spacedim>::objects() const
     */
    template <int dim>
    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<3> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaFaces<dim> *,
                 const std::integral_constant<int, 3>)
    {
      Assert (false, ExcInternalError());
      return nullptr;
    }

    /**
     * Copy the above functions for cell objects.
     */
    template <int structdim, int dim>
    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<structdim> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<dim> > *,
                 const std::integral_constant<int, structdim>)
    {
      Assert (false, ExcInternalError());
      return nullptr;
    }

    template <int dim>
    inline
    dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<dim> > *
    get_objects (dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<dim> > *cells,
                 const std::integral_constant<int, dim>)
    {
      return cells;
    }
  }
}



template <int structdim, int dim, int spacedim>
inline
dealii::internal::TriangulationImplementation::TriaObjects<dealii::internal::TriangulationImplementation::TriaObject<structdim> > &
TriaAccessorBase<structdim,dim,spacedim>::objects() const
{
  if (structdim != dim)
    // get sub-objects. note that the
    // current class is only used for
    // objects that are *not* cells
    return *dealii::internal::TriaAccessorBaseImplementation::get_objects (this->tria->faces.get(),
           std::integral_constant<int, structdim> ());
  else
    return *dealii::internal::TriaAccessorBaseImplementation::get_objects (&this->tria->levels[this->present_level]->cells,
           std::integral_constant<int, structdim> ());
}



/*------------------------ Functions: InvalidAccessor ---------------------------*/

template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::
InvalidAccessor (const Triangulation<dim,spacedim> *,
                 const int,
                 const int,
                 const AccessorData *)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::
InvalidAccessor (const InvalidAccessor &i)
  :
  TriaAccessorBase<structdim,dim,spacedim> (static_cast<const TriaAccessorBase<structdim,dim,spacedim>&>(i))
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::
copy_from (const InvalidAccessor &)
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::
operator == (const InvalidAccessor &) const
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
InvalidAccessor<structdim, dim, spacedim>::
operator != (const InvalidAccessor &) const
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
InvalidAccessor<structdim, dim, spacedim>::used () const
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
InvalidAccessor<structdim, dim, spacedim>::has_children () const
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
InvalidAccessor<structdim, dim, spacedim>::operator ++ () const
{}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::operator -- () const
{}



template <int structdim, int dim, int spacedim>
types::manifold_id
InvalidAccessor<structdim, dim, spacedim>::manifold_id () const
{
  return numbers::invalid_manifold_id;
}


template <int structdim, int dim, int spacedim>
inline
Point<spacedim> &
InvalidAccessor<structdim, dim, spacedim>::vertex (const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  static Point<spacedim> invalid_vertex;
  return invalid_vertex;
}


template <int structdim, int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator
InvalidAccessor<structdim, dim, spacedim>::line (const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator();
}



template <int structdim, int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator
InvalidAccessor<structdim, dim, spacedim>::quad (const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator();
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
      inline
      static
      unsigned int
      line_index (const TriaAccessor<1, dim, spacedim> &,
                  const unsigned int)
      {
        Assert (false,
                ExcMessage ("You can't ask for the index of a line bounding "
                            "a one-dimensional cell because it is not "
                            "bounded by lines."));
        return numbers::invalid_unsigned_int;
      }


      template <int dim, int spacedim>
      inline
      static
      unsigned int
      line_index (const TriaAccessor<2, dim, spacedim> &accessor,
                  const unsigned int i)
      {
        return accessor.objects().cells[accessor.present_index].face(i);
      }


      template <int dim, int spacedim>
      inline
      static
      unsigned int
      line_index (const TriaAccessor<3, dim, spacedim> &accessor,
                  const unsigned int i)
      {
        // get the line index by asking the
        // quads. first assume standard orientation
        //
        // set up a table that for each
        // line describes a) from which
        // quad to take it, b) which line
        // therein it is if the face is
        // oriented correctly
        static const unsigned int lookup_table[12][2] =
        {
          { 4, 0 }, // take first four lines from bottom face
          { 4, 1 },
          { 4, 2 },
          { 4, 3 },

          { 5, 0 }, // second four lines from top face
          { 5, 1 },
          { 5, 2 },
          { 5, 3 },

          { 0, 0 }, // the rest randomly
          { 1, 0 },
          { 0, 1 },
          { 1, 1 }
        };

        // respect non-standard faces by calling the
        // reordering function from GeometryInfo

        const unsigned int quad_index=lookup_table[i][0];
        const unsigned int std_line_index=lookup_table[i][1];

        const unsigned int line_index=GeometryInfo<dim>::standard_to_real_face_line(
                                        std_line_index,
                                        accessor.face_orientation(quad_index),
                                        accessor.face_flip(quad_index),
                                        accessor.face_rotation(quad_index));

        return (accessor.quad(quad_index)->line_index(line_index));
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      unsigned int
      quad_index (const TriaAccessor<structdim, dim, spacedim> &,
                  const unsigned int)
      {
        Assert (false,
                ExcMessage ("You can't ask for the index of a quad bounding "
                            "a one- or two-dimensional cell because it is not "
                            "bounded by quads."));
        return numbers::invalid_unsigned_int;
      }


      template <int dim, int spacedim>
      inline
      static
      unsigned int
      quad_index (const TriaAccessor<3, dim, spacedim> &accessor,
                  const unsigned int i)
      {
        Assert (i<GeometryInfo<3>::quads_per_cell,
                ExcIndexRange(i,0,GeometryInfo<3>::quads_per_cell));
        return accessor.tria->levels[accessor.present_level]
               ->cells.cells[accessor.present_index].face(i);
      }



      /**
       * Implementation of the function of some name in the mother class
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      bool
      face_orientation (const TriaAccessor<structdim, dim, spacedim> &,
                        const unsigned int)
      {
        /*
         * Default implementation used in 1d and 2d
         *
         * In 1d and 2d, face_orientation is always true
         */

        return true;
      }


      template <int dim, int spacedim>
      inline
      static
      bool
      face_orientation (const TriaAccessor<3, dim, spacedim> &accessor,
                        const unsigned int face)
      {
        return (accessor.tria->levels[accessor.present_level]
                ->cells.face_orientation(accessor.present_index, face));
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      bool
      face_flip (const TriaAccessor<structdim, dim, spacedim> &,
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


      template <int dim, int spacedim>
      inline
      static
      bool
      face_flip (const TriaAccessor<3, dim, spacedim> &accessor,
                 const unsigned int face)
      {
        Assert (face<GeometryInfo<3>::faces_per_cell,
                ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
        Assert (accessor.present_index * GeometryInfo<3>::faces_per_cell + face
                < accessor.tria->levels[accessor.present_level]
                ->cells.face_flips.size(),
                ExcInternalError());

        return (accessor.tria->levels[accessor.present_level]
                ->cells.face_flips[accessor.present_index *
                                   GeometryInfo<3>::faces_per_cell
                                   + face]);
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      bool
      face_rotation (const TriaAccessor<structdim, dim, spacedim> &,
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


      template <int dim, int spacedim>
      inline
      static
      bool
      face_rotation (const TriaAccessor<3, dim, spacedim> &accessor,
                     const unsigned int face)
      {
        Assert (face<GeometryInfo<3>::faces_per_cell,
                ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
        Assert (accessor.present_index * GeometryInfo<3>::faces_per_cell + face
                < accessor.tria->levels[accessor.present_level]
                ->cells.face_rotations.size(),
                ExcInternalError());

        return (accessor.tria->levels[accessor.present_level]
                ->cells.face_rotations[accessor.present_index *
                                       GeometryInfo<3>::faces_per_cell
                                       + face]);
      }

      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int dim, int spacedim>
      inline
      static
      bool
      line_orientation (const TriaAccessor<1, dim, spacedim> &,
                        const unsigned int)
      {
        return true;
      }


      template <int spacedim>
      inline
      static
      bool
      line_orientation (const TriaAccessor<2, 2, spacedim> &,
                        const unsigned int)
      {
        // quads in 2d have no non-standard orientation
        return true;
      }


      template <int spacedim>
      inline
      static
      bool
      line_orientation (const TriaAccessor<2, 3, spacedim> &accessor,
                        const unsigned int line)
      {
        // quads as part of 3d hexes can have non-standard orientation

        //TODO: why is this face_orientation, not line_orientation as in the setter function?
        return accessor.tria->faces->quads.face_orientation(accessor.present_index, line);
      }


      template <int dim, int spacedim>
      inline
      static
      bool
      line_orientation (const TriaAccessor<3, dim, spacedim> &accessor,
                        const unsigned int line)
      {
        Assert (accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        Assert (line<GeometryInfo<3>::lines_per_cell,
                ExcIndexRange (line, 0, GeometryInfo<3>::lines_per_cell));

        // get the line index by asking the
        // quads. first assume standard orientation
        //
        // set up a table that for each
        // line describes a) from which
        // quad to take it, b) which line
        // therein it is if the face is
        // oriented correctly
        static const unsigned int lookup_table[12][2] =
        {
          { 4, 0 }, // take first four lines from bottom face
          { 4, 1 },
          { 4, 2 },
          { 4, 3 },

          { 5, 0 }, // second four lines from top face
          { 5, 1 },
          { 5, 2 },
          { 5, 3 },

          { 0, 0 }, // the rest randomly
          { 1, 0 },
          { 0, 1 },
          { 1, 1 }
        };

        const unsigned int quad_index=lookup_table[line][0];
        const unsigned int std_line_index=lookup_table[line][1];

        const unsigned int line_index=GeometryInfo<dim>::standard_to_real_face_line(
                                        std_line_index,
                                        accessor.face_orientation(quad_index),
                                        accessor.face_flip(quad_index),
                                        accessor.face_rotation(quad_index));

        // now we got to the correct line and ask
        // the quad for its line_orientation. however, if
        // the face is rotated, it might be possible,
        // that a standard orientation of the line
        // with respect to the face corresponds to a
        // non-standard orientation for the line with
        // respect to the cell.
        //
        // set up a table indicating if the two
        // standard orientations coincide
        //
        // first index: two pairs of lines 0(lines
        // 0/1) and 1(lines 2/3)
        //
        // second index: face_orientation; 0:
        // opposite normal, 1: standard
        //
        // third index: face_flip; 0: standard, 1:
        // face rotated by 180 degrees
        //
        // forth index: face_rotation: 0: standard,
        // 1: face rotated by 90 degrees

        static const bool bool_table[2][2][2][2] =
        {
          { { { true, false },   // lines 0/1, face_orientation=false, face_flip=false, face_rotation=false and true
              { false, true }
            },  // lines 0/1, face_orientation=false, face_flip=true, face_rotation=false and true
            { { true, true },    // lines 0/1, face_orientation=true, face_flip=false, face_rotation=false and true
              { false, false }
            }
          },// lines 0/1, face_orientation=true, face_flip=true, face_rotation=false and true

          { { { true, true },    // lines 2/3 ...
              { false, false }
            },
            { { true, false },
              { false, true }
            }
          }
        };


        return (accessor.quad(quad_index)
                ->line_orientation(line_index)
                == bool_table[std_line_index/2]
                [accessor.face_orientation(quad_index)]
                [accessor.face_flip(quad_index)]
                [accessor.face_rotation(quad_index)]);
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      void
      set_face_orientation (const TriaAccessor<structdim, dim, spacedim> &,
                            const unsigned int,
                            const bool)
      {
        Assert (false, ExcInternalError());
      }


      template <int dim, int spacedim>
      inline
      static
      void
      set_face_orientation (const TriaAccessor<3, dim, spacedim> &accessor,
                            const unsigned int face,
                            const bool value)
      {
        Assert (accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        Assert (face<GeometryInfo<3>::faces_per_cell,
                ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
        Assert (accessor.present_index * GeometryInfo<3>::faces_per_cell + face
                < accessor.tria->levels[accessor.present_level]
                ->cells.face_orientations.size(),
                ExcInternalError());
        accessor.tria->levels[accessor.present_level]
        ->cells.face_orientations[accessor.present_index *
                                  GeometryInfo<3>::faces_per_cell
                                  +
                                  face] = value;
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      void
      set_face_flip (const TriaAccessor<structdim, dim, spacedim> &,
                     const unsigned int,
                     const bool)
      {
        Assert (false, ExcInternalError());
      }


      template <int dim, int spacedim>
      inline
      static
      void
      set_face_flip (const TriaAccessor<3, dim, spacedim> &accessor,
                     const unsigned int face,
                     const bool value)
      {
        Assert (face<GeometryInfo<3>::faces_per_cell,
                ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
        Assert (accessor.present_index * GeometryInfo<3>::faces_per_cell + face
                < accessor.tria->levels[accessor.present_level]
                ->cells.face_flips.size(),
                ExcInternalError());

        accessor.tria->levels[accessor.present_level]
        ->cells.face_flips[accessor.present_index *
                           GeometryInfo<3>::faces_per_cell
                           + face] = value;
      }



      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int structdim, int dim, int spacedim>
      inline
      static
      void
      set_face_rotation (const TriaAccessor<structdim, dim, spacedim> &,
                         const unsigned int,
                         const bool)
      {
        Assert (false, ExcInternalError());
      }


      template <int dim, int spacedim>
      inline
      static
      void
      set_face_rotation (const TriaAccessor<3, dim, spacedim> &accessor,
                         const unsigned int face,
                         const bool value)
      {
        Assert (face<GeometryInfo<3>::faces_per_cell,
                ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
        Assert (accessor.present_index * GeometryInfo<3>::faces_per_cell + face
                < accessor.tria->levels[accessor.present_level]
                ->cells.face_rotations.size(),
                ExcInternalError());

        accessor.tria->levels[accessor.present_level]
        ->cells.face_rotations[accessor.present_index *
                               GeometryInfo<3>::faces_per_cell
                               + face] = value;
      }

      /**
       * Implementation of the function of some name in the mother class.
       */
      template <int dim, int spacedim>
      inline
      static
      void
      set_line_orientation (const TriaAccessor<1, dim, spacedim> &,
                            const unsigned int,
                            const bool)
      {
        Assert (false, ExcInternalError());
      }


      template <int spacedim>
      inline
      static
      void
      set_line_orientation (const TriaAccessor<2, 2, spacedim> &,
                            const unsigned int,
                            const bool)
      {
        // quads in 2d have no
        // non-standard orientation
        Assert (false, ExcInternalError());
      }


      template <int spacedim>
      inline
      static
      void
      set_line_orientation (const TriaAccessor<2, 3, spacedim> &accessor,
                            const unsigned int line,
                            const bool value)
      {
        Assert (accessor.used(), TriaAccessorExceptions::ExcCellNotUsed());
        Assert (line<GeometryInfo<3>::lines_per_face,
                ExcIndexRange (line, 0, GeometryInfo<3>::lines_per_face));
        Assert (accessor.present_index * GeometryInfo<3>::lines_per_face + line
                < accessor.tria->faces->quads.line_orientations.size(),
                ExcInternalError());
        // quads as part of 3d hexes
        // can have non-standard
        // orientation
        accessor.tria->faces->quads.line_orientations[accessor.present_index *
                                                      GeometryInfo<3>::lines_per_face
                                                      + line]
          = value;
      }


      template <int dim, int spacedim>
      inline
      static
      void
      set_line_orientation (const TriaAccessor<3, dim, spacedim> &,
                            const unsigned int,
                            const bool)
      {
        // it seems like we don't need this
        // one
        Assert (false, ExcNotImplemented());
      }


      /**
       * Implementation of the function of same name in the enclosing class.
       */
      template <int dim, int spacedim>
      inline
      static
      unsigned int
      vertex_index (const TriaAccessor<1,dim,spacedim> &accessor,
                    const unsigned int corner)
      {
        return accessor.objects().cells[accessor.present_index].face (corner);
      }


      template <int dim, int spacedim>
      inline
      static
      unsigned int
      vertex_index (const TriaAccessor<2,dim,spacedim> &accessor,
                    const unsigned int corner)
      {
        // table used to switch the vertices, if the
        // line orientation is wrong,
        //
        // first index: line orientation 0: false or
        // 1: true=standard
        //
        // second index: vertex index to be switched
        // (or not)

        static const unsigned int switch_table[2][2]= {{1,0},{0,1}};

        return accessor.line(corner%2)
               ->vertex_index(switch_table[accessor.line_orientation(corner%2)][corner/2]);
      }



      template <int dim, int spacedim>
      inline
      static
      unsigned int
      vertex_index (const TriaAccessor<3,dim,spacedim> &accessor,
                    const unsigned int corner)
      {
        // get the corner indices by asking either
        // the bottom or the top face for its
        // vertices. handle non-standard faces by
        // calling the vertex reordering function
        // from GeometryInfo

        // bottom face (4) for first four vertices,
        // top face (5) for the rest
        const unsigned int face_index=4+corner/4;

        return accessor.quad(face_index)
               ->vertex_index(GeometryInfo<dim>
                              ::standard_to_real_face_vertex(corner%4,
                                                             accessor.face_orientation(face_index),
                                                             accessor.face_flip(face_index),
                                                             accessor.face_rotation(face_index)));
      }
    };
  }
}




template <int structdim, int dim, int spacedim>
inline
TriaAccessor<structdim, dim, spacedim>::
TriaAccessor (const Triangulation<dim,spacedim> *parent,
              const int                 level,
              const int                 index,
              const AccessorData       *local_data)
  :
  TriaAccessorBase<structdim,dim,spacedim> (parent, level, index, local_data)
{}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::used () const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));
  return this->objects().used[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline
TriaIterator<TriaAccessor<0,dim,spacedim> >
TriaAccessor<structdim,dim,spacedim>::vertex_iterator (const unsigned int i) const
{
  return TriaIterator<TriaAccessor<0,dim,spacedim> >(this->tria, 0, vertex_index (i));
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim, dim, spacedim>::
vertex_index (const unsigned int corner) const
{
  Assert (corner<GeometryInfo<structdim>::vertices_per_cell,
          ExcIndexRange(corner,0,GeometryInfo<structdim>::vertices_per_cell));

  return dealii::internal::TriaAccessorImplementation::Implementation::vertex_index (*this, corner);
}



template <int structdim, int dim, int spacedim>
inline
Point<spacedim> &
TriaAccessor<structdim, dim, spacedim>::vertex (const unsigned int i) const
{
  return const_cast<Point<spacedim> &> (this->tria->vertices[vertex_index(i)]);
}



template <int structdim, int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator
TriaAccessor<structdim,dim,spacedim>::line (const unsigned int i) const
{
  // checks happen in line_index
  return typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator
         (this->tria, 0, line_index (i));
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::line_index (const unsigned int i) const
{
  Assert (i < GeometryInfo<structdim>::lines_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<structdim>::lines_per_cell));

  return dealii::internal::TriaAccessorImplementation::Implementation::line_index (*this, i);
}




template <int structdim, int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator
TriaAccessor<structdim,dim,spacedim>::quad (const unsigned int i) const
{
  // checks happen in quad_index
  return typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator
         (this->tria, 0, quad_index (i));
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::quad_index (const unsigned int i) const
{
  return dealii::internal::TriaAccessorImplementation::Implementation::quad_index (*this, i);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_orientation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::face_orientation (*this, face);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_flip (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::face_flip (*this, face);
}


template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_rotation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return dealii::internal::TriaAccessorImplementation::Implementation::face_rotation (*this, face);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::line_orientation (const unsigned int line) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (line<GeometryInfo<structdim>::lines_per_cell,
          ExcIndexRange (line, 0, GeometryInfo<structdim>::lines_per_cell));

  return dealii::internal::TriaAccessorImplementation::Implementation::line_orientation (*this, line);
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_orientation (const unsigned int face,
                                                            const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  dealii::internal::TriaAccessorImplementation::Implementation::set_face_orientation (*this, face, value);
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_flip (const unsigned int face,
                                                     const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  dealii::internal::TriaAccessorImplementation::Implementation::set_face_flip (*this, face, value);
}


template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_rotation (const unsigned int face,
                                                         const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  dealii::internal::TriaAccessorImplementation::Implementation::set_face_rotation (*this, face, value);
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_line_orientation (const unsigned int line,
                                                            const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (line<GeometryInfo<structdim>::lines_per_cell,
          ExcIndexRange (line, 0, GeometryInfo<structdim>::lines_per_cell));

  dealii::internal::TriaAccessorImplementation::Implementation::set_line_orientation (*this, line, value);
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim, dim, spacedim>::set_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));
  this->objects().used[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));
  this->objects().used[this->present_index] = false;
}


template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::
child_index (const unsigned int i) const
{
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  Assert (i<n_children(),
          ExcIndexRange (i, 0, n_children()));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two = GeometryInfo<structdim>::max_children_per_cell/2;
  return this->objects().children[n_sets_of_two*this->present_index+i/2]+i%2;
}



template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::
isotropic_child_index (const unsigned int i) const
{
  Assert (i<GeometryInfo<structdim>::max_children_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<structdim>::max_children_per_cell));

  switch (structdim)
    {
    case 1:
      return child_index (i);
    case 2:
    {
      const RefinementCase<2>
      this_refinement_case (static_cast<std::uint8_t>(refinement_case()));

      Assert (this_refinement_case != RefinementCase<2>::no_refinement,
              TriaAccessorExceptions::ExcCellHasNoChildren());

      if (this_refinement_case == RefinementCase<2>::cut_xy)
        return child_index(i);
      else if ((this_refinement_case == RefinementCase<2>::cut_x)
               &&
               (child(i%2)->refinement_case()==RefinementCase<2>::cut_y))
        return child(i%2)->child_index(i/2);
      else if ((this_refinement_case == RefinementCase<2>::cut_y)
               &&
               (child(i/2)->refinement_case()==RefinementCase<2>::cut_x))
        return child(i/2)->child_index(i%2);
      else
        Assert(false,
               ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
      break;
    }

    case 3:
      Assert (false, ExcNotImplemented());
    }
  return -1;
}



template <int structdim, int dim, int spacedim>
RefinementCase<structdim>
TriaAccessor<structdim, dim, spacedim>::refinement_case() const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));

  switch (structdim)
    {
    case 1:
      return (RefinementCase<structdim>
              (this->objects().children[this->present_index] != -1
               ?
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
      Assert (static_cast<unsigned int> (this->present_index) <
              this->objects().refinement_cases.size(),
              ExcIndexRange(this->present_index, 0,
                            this->objects().refinement_cases.size()));

      return (static_cast<RefinementCase<structdim> >
              (this->objects().refinement_cases[this->present_index]));
    }
}




template <int structdim, int dim, int spacedim>
inline
TriaIterator<TriaAccessor<structdim,dim,spacedim> >
TriaAccessor<structdim,dim,spacedim>::child (const unsigned int i) const

{
  // checking of 'i' happens in child_index
  const TriaIterator<TriaAccessor<structdim,dim,spacedim> >
  q (this->tria,
     (dim == structdim ? this->level() + 1 : 0),
     child_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
          ExcInternalError());

  return q;
}



template <int structdim, int dim, int spacedim>
inline
TriaIterator<TriaAccessor<structdim,dim,spacedim> >
TriaAccessor<structdim,dim,spacedim>::
isotropic_child (const unsigned int i) const
{
  // checking of 'i' happens in child() or
  // child_index() called below
  switch (structdim)
    {
    case 1:
      // no anisotropic refinement in 1D
      return child(i);

    case 2:
    {
      const RefinementCase<2>
      this_refinement_case (static_cast<std::uint8_t>(refinement_case()));

      Assert (this_refinement_case != RefinementCase<2>::no_refinement,
              TriaAccessorExceptions::ExcCellHasNoChildren());

      if (this_refinement_case == RefinementCase<2>::cut_xy)
        return child(i);
      else if ((this_refinement_case == RefinementCase<2>::cut_x)
               &&
               (child(i%2)->refinement_case()==RefinementCase<2>::cut_y))
        return child(i%2)->child(i/2);
      else if ((this_refinement_case == RefinementCase<2>::cut_y)
               &&
               (child(i/2)->refinement_case()==RefinementCase<2>::cut_x))
        return child(i/2)->child(i%2);
      else
        Assert(false,
               ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }
  // we don't get here but have to return
  // something...
  return child(0);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two = GeometryInfo<structdim>::max_children_per_cell/2;
  return (this->objects().children[n_sets_of_two * this->present_index] != -1);
}




template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::n_children () const
{
  return GeometryInfo<structdim>::n_children(refinement_case());
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim, dim, spacedim>::
set_refinement_case (const RefinementCase<structdim> &refinement_case) const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));
  Assert (static_cast<unsigned int> (this->present_index) <
          this->objects().refinement_cases.size(),
          ExcIndexRange(this->present_index, 0,
                        this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index] = refinement_case;
}


template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim, dim, spacedim>::clear_refinement_case () const
{
  Assert (this->state() == IteratorState::valid,
          TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(*this));
  Assert (static_cast<unsigned int> (this->present_index) <
          this->objects().refinement_cases.size(),
          ExcIndexRange(this->present_index, 0,
                        this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index]
    = RefinementCase<structdim>::no_refinement;
}





template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_children (const unsigned int i,
                                                      const int index) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i%2==0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two = GeometryInfo<structdim>::max_children_per_cell/2;

  Assert (
    // clearing the child index for a cell
    (index==-1)
    ||
    // if setting the child index for the i'th child (with i==0),
    // then the index must be a non-negative number
    (i==0 && !this->has_children() && (index>=0))
    ||
    // if setting the child index for the i'th child (with i>0),
    // then the previously stored index must be the invalid
    // index
    (i>0  &&  this->has_children() && (index>=0) &&
     this->objects().children[n_sets_of_two*this->present_index+i/2] == -1),
    TriaAccessorExceptions::ExcCantSetChildren(index));

  this->objects().children[n_sets_of_two*this->present_index+i/2] = index;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_children () const
{
  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two = GeometryInfo<structdim>::max_children_per_cell/2;

  for (unsigned int i=0; i<n_sets_of_two; ++i)
    set_children (2*i,-1);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::user_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_flags[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::clear_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = false;
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::recursively_set_user_flag () const
{
  set_user_flag ();

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_set_user_flag ();
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::recursively_clear_user_flag () const
{
  clear_user_flag ();

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_clear_user_flag ();
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::clear_user_data () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().clear_user_data(this->present_index);
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::set_user_pointer (void *p) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::clear_user_pointer () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = nullptr;
}



template <int structdim, int dim, int spacedim>
void *TriaAccessor<structdim,dim,spacedim>::user_pointer () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_pointer(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim,dim,spacedim>::recursively_set_user_pointer (void *p) const
{
  set_user_pointer (p);

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_set_user_pointer (p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim,dim,spacedim>::recursively_clear_user_pointer () const
{
  clear_user_pointer ();

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_clear_user_pointer ();
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::set_user_index (unsigned int p) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim,dim,spacedim>::clear_user_index () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = 0;
}



template <int structdim, int dim, int spacedim>
unsigned int TriaAccessor<structdim,dim,spacedim>::user_index () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_index(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim,dim,spacedim>::recursively_set_user_index (unsigned int p) const
{
  set_user_index (p);

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_set_user_index (p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim,dim,spacedim>::recursively_clear_user_index () const
{
  clear_user_index ();

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_clear_user_index ();
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::max_refinement_depth () const
{
  if (!this->has_children())
    return 0;

  unsigned int max_depth = 1;
  for (unsigned int c=0; c<n_children(); ++c)
    max_depth = std::max (max_depth,
                          child(c)->max_refinement_depth() + 1);
  return max_depth;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::number_of_children () const
{
  if (!this->has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<n_children(); ++c)
        sum += this->child(c)->number_of_children();
      return sum;
    }
}



template <int structdim, int dim, int spacedim>
types::boundary_id
TriaAccessor<structdim, dim, spacedim>::boundary_id () const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects().boundary_or_material_id[this->present_index].boundary_id;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_boundary_id (const types::boundary_id boundary_ind) const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (boundary_ind != numbers::internal_face_boundary_id,
          ExcMessage("You are trying to set the boundary_id to an illegal value (numbers::internal_face_boundary_id is reserved)."));
  Assert (this->at_boundary(),
          ExcMessage("You are trying to set the boundary_id of an internal object, which is illegal!"));

  this->objects().boundary_or_material_id[this->present_index].boundary_id = boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_boundary_id_internal (const types::boundary_id boundary_ind) const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().boundary_or_material_id[this->present_index].boundary_id = boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_all_boundary_ids (const types::boundary_id boundary_ind) const
{
  set_boundary_id (boundary_ind);

  switch (structdim)
    {
    case 1:
      // 1d objects have no sub-objects
      // where we have to do anything
      break;

    case 2:
      // for boundary quads also set
      // boundary_id of bounding lines
      for (unsigned int i=0; i<4; ++i)
        this->line(i)->set_boundary_id (boundary_ind);
      break;

    default:
      Assert (false, ExcNotImplemented());
    }
}



template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::at_boundary () const
{
  // error checking is done
  // in boundary_id()
  return (boundary_id() != numbers::internal_face_boundary_id);
}



template <int structdim, int dim, int spacedim>
const Boundary<dim,spacedim> &
TriaAccessor<structdim, dim, spacedim>::get_boundary () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  // Get the default (manifold_id)
  const types::manifold_id mi = this->objects().manifold_id[this->present_index];

  // In case this is not valid, check
  // the boundary id, after having
  // casted it to a manifold id
  if (mi == numbers::invalid_manifold_id)
    return this->tria->get_boundary(structdim < dim ?
                                    this->objects().boundary_or_material_id[this->present_index].boundary_id:
                                    dim < spacedim ?
                                    this->objects().boundary_or_material_id[this->present_index].material_id:
                                    numbers::invalid_manifold_id);
  else
    return this->tria->get_boundary(mi);
}


template <int structdim, int dim, int spacedim>
const Manifold<dim,spacedim> &
TriaAccessor<structdim, dim, spacedim>::get_manifold () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->get_manifold(this->manifold_id());
}


template <int structdim, int dim, int spacedim>
types::manifold_id
TriaAccessor<structdim, dim, spacedim>::manifold_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects().manifold_id[this->present_index];
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_manifold_id (const types::manifold_id manifold_ind) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().manifold_id[this->present_index] = manifold_ind;
}


template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_all_manifold_ids (const types::manifold_id manifold_ind) const
{
  set_manifold_id (manifold_ind);

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->set_all_manifold_ids (manifold_ind);

  switch (structdim)
    {
    case 1:
      if (dim == 1)
        {
          (*this->tria->vertex_to_manifold_id_map_1d)
          [vertex_index(0)] = manifold_ind;
          (*this->tria->vertex_to_manifold_id_map_1d)
          [vertex_index(1)] = manifold_ind;
        }
      break;

    case 2:
      // for quads also set manifold_id of bounding lines
      for (unsigned int i=0; i<4; ++i)
        this->line(i)->set_manifold_id (manifold_ind);
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
}



template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::diameter () const
{
  switch (structdim)
    {
    case 1:
      return (this->vertex(1)-this->vertex(0)).norm();
    case 2:
      return std::max((this->vertex(3)-this->vertex(0)).norm(),
                      (this->vertex(2)-this->vertex(1)).norm());
    case 3:
      return std::max( std::max((this->vertex(7)-this->vertex(0)).norm(),
                                (this->vertex(6)-this->vertex(1)).norm()),
                       std::max((this->vertex(2)-this->vertex(5)).norm(),
                                (this->vertex(3)-this->vertex(4)).norm()) );
    default:
      Assert (false, ExcNotImplemented());
      return -1e10;
    }
}



template <int structdim, int dim, int spacedim>
std::pair<Point<spacedim>,double>
TriaAccessor<structdim, dim, spacedim>::enclosing_ball () const
{
  // If the object is one dimensional,
  // the enclosing ball is the initial iterate
  // i.e., the ball's center and diameter are
  // the center and the diameter of the object.
  if (structdim==1)
    return std::make_pair( (this->vertex(1)+this->vertex(0))*0.5,
                           (this->vertex(1)-this->vertex(0)).norm()*0.5);

  // The list is_initial_guess_vertex contains bool values and has
  // the same size as the number of vertices per object.
  // The entries of is_initial_guess_vertex are set true only for those
  // two vertices corresponding to the largest diagonal which is being used
  // to construct the initial ball.
  // We employ this mask to skip these two vertices while enlarging the ball.
  std::array<bool, GeometryInfo<structdim>::vertices_per_cell> is_initial_guess_vertex;

  //First let all the vertices be outside
  std::fill(is_initial_guess_vertex.begin(),
            is_initial_guess_vertex.end(),
            false);

  // Get an initial guess by looking at the largest diagonal
  Point<spacedim> center;
  double radius = 0;

  switch (structdim)
    {
    case 2:
    {
      const Point<spacedim> p30( this->vertex(3)-this->vertex(0));
      const Point<spacedim> p21( this->vertex(2)-this->vertex(1));
      if (p30.norm() > p21.norm())
        {
          center = this->vertex(0) + 0.5* p30;
          radius = p30.norm()/2.;
          is_initial_guess_vertex[3] = true;
          is_initial_guess_vertex[0] = true;
        }
      else
        {
          center = this->vertex(1) + 0.5* p21;
          radius = p21.norm()/2.;
          is_initial_guess_vertex[2] = true;
          is_initial_guess_vertex[1] = true;
        }
      break;
    }
    case 3:
    {
      const Point<spacedim> p70( this->vertex(7)-this->vertex(0));
      const Point<spacedim> p61( this->vertex(6)-this->vertex(1));
      const Point<spacedim> p25( this->vertex(2)-this->vertex(5));
      const Point<spacedim> p34( this->vertex(3)-this->vertex(4));
      const std::vector<double> diagonals= { p70.norm(),
                                             p61.norm(),
                                             p25.norm(),
                                             p34.norm()
                                           };
      const std::vector<double>::const_iterator
      it = std::max_element( diagonals.begin(), diagonals.end());
      if (it == diagonals.begin())
        {
          center = this->vertex(0) + 0.5* p70;
          is_initial_guess_vertex[7] = true;
          is_initial_guess_vertex[0] = true;
        }
      else if (it == diagonals.begin()+1)
        {
          center = this->vertex(1) + 0.5* p61;
          is_initial_guess_vertex[6] = true;
          is_initial_guess_vertex[1] = true;
        }
      else if (it == diagonals.begin()+2)
        {
          center = this->vertex(5) + 0.5* p25;
          is_initial_guess_vertex[2] = true;
          is_initial_guess_vertex[5] = true;
        }
      else
        {
          center = this->vertex(4) + 0.5* p34;
          is_initial_guess_vertex[3] = true;
          is_initial_guess_vertex[4] = true;
        }
      radius = *it * 0.5;
      break;
    }
    default:
      Assert (false, ExcNotImplemented());
      return std::pair<Point<spacedim>,double>();
    }

  // For each vertex that is found to be geometrically outside the ball
  // enlarge the ball  so that the new ball contains both the previous ball
  // and the given vertex.
  for (unsigned int v = 0; v < GeometryInfo<structdim>::vertices_per_cell; ++v)
    if (!is_initial_guess_vertex[v])
      {
        const double distance = center.distance(this->vertex(v));
        if (distance > radius)
          {
            // we found a vertex which is outside of the ball
            // extend it (move center and change radius)
            const Point<spacedim> pCV (center - this->vertex(v));
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
  for (unsigned int v = 0; v < GeometryInfo<structdim>::vertices_per_cell; ++v)
    if (center.distance(this->vertex(v)) > radius + 100. *std::numeric_limits<double>::epsilon())
      {
        all_vertices_within_ball = false;
        break;
      }
  // If all the vertices are not within the ball throw error
  Assert (all_vertices_within_ball, ExcInternalError());
#endif
  return std::make_pair(center, radius);
}


template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::minimum_vertex_distance () const
{
  switch (structdim)
    {
    case 1:
      return (this->vertex(1)-this->vertex(0)).norm();
    case 2:
    case 3:
    {
      double min = std::numeric_limits<double>::max();
      for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
        for (unsigned int j=i+1; j<GeometryInfo<structdim>::vertices_per_cell; ++j)
          min = std::min(min, (this->vertex(i)-this->vertex(j)) * (this->vertex(i)-this->vertex(j)));
      return std::sqrt(min);
    }
    default:
      Assert (false, ExcNotImplemented());
      return -1e10;
    }
}


template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::
is_translation_of (const TriaIterator<TriaAccessor<structdim,dim,spacedim> > &o) const
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
  bool is_translation = true;
  const Tensor<1,spacedim> dist = o->vertex(0) - this->vertex(0);
  const double tol_square = 1e-24 * dist.norm_square();
  for (unsigned int i=1; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
    {
      const Tensor<1,spacedim> dist_new = (o->vertex(i) - this->vertex(i)) - dist;
      if (dist_new.norm_square() > tol_square)
        {
          is_translation = false;
          break;
        }
    }
  return is_translation;
}


/*------------------------ Functions: TriaAccessor<0,dim,spacedim> -----------------------*/

template <int dim, int spacedim>
inline
TriaAccessor<0, dim, spacedim>::
TriaAccessor (const Triangulation<dim, spacedim> *tria,
              const unsigned int                  vertex_index)
  :
  tria (tria),
  global_vertex_index (vertex_index)
{}



template <int dim, int spacedim>
inline
TriaAccessor<0, dim, spacedim>::
TriaAccessor (const Triangulation<dim,spacedim> *tria,
              const int /*level*/,
              const int index,
              const AccessorData *)
  :
  tria (tria),
  global_vertex_index (index)
{}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline
TriaAccessor<0, dim, spacedim>::
TriaAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &)
  :
  tria (nullptr),
  global_vertex_index (numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline
TriaAccessor<0, dim, spacedim>::
TriaAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
  :
  tria (nullptr),
  global_vertex_index (numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
inline
void
TriaAccessor<0, dim, spacedim>::copy_from (const TriaAccessor &t)
{
  tria = t.tria;
  global_vertex_index = t.global_vertex_index;
}



template <int dim, int spacedim>
inline
IteratorState::IteratorStates
TriaAccessor<0, dim, spacedim>::state () const
{
  if (global_vertex_index != numbers::invalid_unsigned_int)
    return IteratorState::valid;
  else
    return IteratorState::past_the_end;
}



template <int dim, int spacedim>
inline
int
TriaAccessor<0, dim, spacedim>::level ()
{
  return 0;
}



template <int dim, int spacedim>
inline
int
TriaAccessor<0, dim, spacedim>::index () const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline
void
TriaAccessor<0, dim, spacedim>::operator ++ ()
{
  ++global_vertex_index;
  if (global_vertex_index >= tria->n_vertices())
    global_vertex_index = numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline
void
TriaAccessor<0, dim, spacedim>::operator -- ()
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
inline
bool
TriaAccessor<0, dim, spacedim>::operator == (const TriaAccessor &t) const
{
  const bool result = ((tria == t.tria)
                       &&
                       (global_vertex_index == t.global_vertex_index));

  return result;
}



template <int dim, int spacedim>
inline
bool
TriaAccessor<0, dim, spacedim>::operator != (const TriaAccessor &t) const
{
  return !(*this==t);
}



template <int dim, int spacedim>
inline
unsigned int
TriaAccessor<0, dim, spacedim>::vertex_index (const unsigned int) const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline
Point<spacedim> &
TriaAccessor<0, dim, spacedim>::vertex (const unsigned int) const
{
  return const_cast<Point<spacedim> &> (this->tria->vertices[global_vertex_index]);
}



template <int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator
TriaAccessor<0, dim, spacedim>::line (const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::line_iterator();
}



template <int dim, int spacedim>
inline
unsigned int
TriaAccessor<0, dim, spacedim>::line_index (const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator
TriaAccessor<0, dim, spacedim>::quad (const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::Iterators<dim,spacedim>::quad_iterator();
}



template <int dim, int spacedim>
inline
unsigned int
TriaAccessor<0, dim, spacedim>::quad_index (const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline
double
TriaAccessor<0, dim, spacedim>::diameter () const
{
  return 0.;
}



template <int dim, int spacedim>
inline
double
TriaAccessor<0, dim, spacedim>::extent_in_direction (const unsigned int) const
{
  return 0.;
}



template <int dim, int spacedim>
inline
Point<spacedim>
TriaAccessor<0, dim, spacedim>::center (const bool,
                                        const bool) const
{
  return this->tria->vertices[global_vertex_index];
}



template <int dim, int spacedim>
inline
double
TriaAccessor<0, dim, spacedim>::measure () const
{
  return 0.;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::face_orientation (const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::face_flip (const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::face_rotation (const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::line_orientation (const unsigned int /*line*/)
{
  return false;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::has_children ()
{
  return false;
}



template <int dim, int spacedim>
inline
unsigned int TriaAccessor<0, dim, spacedim>::n_children()
{
  return 0;
}



template <int dim, int spacedim>
inline
unsigned int TriaAccessor<0, dim, spacedim>::number_of_children ()
{
  return 0;
}



template <int dim, int spacedim>
inline
unsigned int TriaAccessor<0, dim, spacedim>::max_refinement_depth ()
{
  return 0;
}



template <int dim, int spacedim>
inline
TriaIterator<TriaAccessor<0,dim,spacedim> >
TriaAccessor<0, dim, spacedim>::child (const unsigned int)
{
  return TriaIterator<TriaAccessor<0,dim,spacedim> >();
}



template <int dim, int spacedim>
inline
TriaIterator<TriaAccessor<0,dim,spacedim> >
TriaAccessor<0, dim, spacedim>::isotropic_child (const unsigned int)
{
  return TriaIterator<TriaAccessor<0,dim,spacedim> >();
}



template <int dim, int spacedim>
inline
RefinementCase<0> TriaAccessor<0, dim, spacedim>::refinement_case ()
{
  return RefinementCase<0>(RefinementPossibilities<0>::no_refinement);
}



template <int dim, int spacedim>
inline
int TriaAccessor<0, dim, spacedim>::child_index (const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline
int TriaAccessor<0, dim, spacedim>::isotropic_child_index (const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline
bool TriaAccessor<0, dim, spacedim>::used () const
{
  return tria->vertex_used(global_vertex_index);
}



/*------------------------ Functions: TriaAccessor<0,1,spacedim> -----------------------*/

template <int spacedim>
inline
TriaAccessor<0, 1, spacedim>::
TriaAccessor (const Triangulation<1,spacedim> *tria,
              const VertexKind      vertex_kind,
              const unsigned int    vertex_index)
  :
  tria (tria),
  vertex_kind (vertex_kind),
  global_vertex_index (vertex_index)
{}



template <int spacedim>
inline
TriaAccessor<0, 1, spacedim>::
TriaAccessor (const Triangulation<1,spacedim> *tria,
              const int level,
              const int index,
              const AccessorData *)
  :
  tria (tria),
  vertex_kind (interior_vertex),
  global_vertex_index (numbers::invalid_unsigned_int)
{
  // in general, calling this constructor should yield an error -- users should
  // instead call the one immediately above. however, if you create something
  // like Triangulation<1>::face_iterator() then this calls the default constructor
  // of the iterator which calls the accessor with argument list (0,-2,-2,0), so
  // in this particular case accept this call and create an object that corresponds
  // to the default constructed (invalid) vertex accessor
  (void)level;
  (void)index;
  Assert ((level == -2) && (index == -2), ExcInternalError());
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline
TriaAccessor<0, 1, spacedim>::
TriaAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &)
  :
  tria (nullptr),
  vertex_kind (interior_vertex),
  global_vertex_index (numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline
TriaAccessor<0, 1, spacedim>::
TriaAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
  :
  tria (nullptr),
  vertex_kind (interior_vertex),
  global_vertex_index (numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
inline
void
TriaAccessor<0, 1, spacedim>::copy_from (const TriaAccessor &t)
{
  tria = t.tria;
  vertex_kind = t.vertex_kind;
  global_vertex_index = t.global_vertex_index;
}



template <int spacedim>
inline
IteratorState::IteratorStates
TriaAccessor<0, 1, spacedim>::state ()
{
  return IteratorState::valid;
}


template <int spacedim>
inline
int
TriaAccessor<0, 1, spacedim>::level ()
{
  return 0;
}



template <int spacedim>
inline
int
TriaAccessor<0, 1, spacedim>::index () const
{
  return global_vertex_index;
}



template <int spacedim>
inline
void
TriaAccessor<0, 1, spacedim>::operator ++ () const
{
  Assert (false, ExcNotImplemented());
}


template <int spacedim>
inline
void
TriaAccessor<0, 1, spacedim>::operator -- () const
{
  Assert (false, ExcNotImplemented());
}



template <int spacedim>
inline
bool
TriaAccessor<0, 1, spacedim>::operator == (const TriaAccessor &t) const
{
  const bool result = ((tria == t.tria)
                       &&
                       (global_vertex_index == t.global_vertex_index));
  // if we point to the same vertex,
  // make sure we know the same about
  // it
  if (result == true)
    Assert (vertex_kind == t.vertex_kind, ExcInternalError());

  return result;
}



template <int spacedim>
inline
bool
TriaAccessor<0, 1, spacedim>::operator != (const TriaAccessor &t) const
{
  return !(*this==t);
}



template <int spacedim>
inline
unsigned int
TriaAccessor<0, 1, spacedim>::vertex_index (const unsigned int i) const
{
  Assert(i==0, ExcIndexRange(i, 0, 1));
  (void)i;
  return global_vertex_index;
}



template <int spacedim>
inline
Point<spacedim> &
TriaAccessor<0, 1, spacedim>::vertex (const unsigned int i) const
{
  Assert(i==0, ExcIndexRange(i, 0, 1));
  (void)i;
  return const_cast<Point<spacedim> &> (this->tria->vertices[global_vertex_index]);
}



template <int spacedim>
inline
Point<spacedim>
TriaAccessor<0, 1, spacedim>::center () const
{
  return this->tria->vertices[global_vertex_index];
}



template <int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<1,spacedim>::line_iterator
TriaAccessor<0, 1, spacedim>::line (const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::Iterators<1,spacedim>::line_iterator();
}


template <int spacedim>
inline
unsigned int
TriaAccessor<0, 1, spacedim>::line_index (const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}


template <int spacedim>
inline
typename dealii::internal::TriangulationImplementation::Iterators<1,spacedim>::quad_iterator
TriaAccessor<0, 1, spacedim>::quad (const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::Iterators<1,spacedim>::quad_iterator();
}



template <int spacedim>
inline
unsigned int
TriaAccessor<0, 1, spacedim>::quad_index (const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}


template <int spacedim>
inline
bool
TriaAccessor<0, 1, spacedim>::at_boundary () const
{
  return vertex_kind != interior_vertex;
}


template <int spacedim>
inline
types::boundary_id
TriaAccessor<0, 1, spacedim>::boundary_id () const
{
  switch (vertex_kind)
    {
    case left_vertex:
    case right_vertex:
    {
      Assert (tria->vertex_to_boundary_id_map_1d->find (this->vertex_index())
              != tria->vertex_to_boundary_id_map_1d->end(),
              ExcInternalError());

      return (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()];
    }

    default:
      return numbers::internal_face_boundary_id;
    }
}



template <int spacedim>
inline
const Manifold<1, spacedim> &
TriaAccessor<0, 1, spacedim>::get_manifold () const
{
  return this->tria->get_manifold(this->manifold_id());
}



template <int spacedim>
inline
types::manifold_id
TriaAccessor<0, 1, spacedim>::manifold_id () const
{
  if ( tria->vertex_to_manifold_id_map_1d->find (this->vertex_index())
       != tria->vertex_to_manifold_id_map_1d->end())
    return (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()];
  else
    return numbers::invalid_manifold_id;
}


template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::face_orientation (const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::face_flip (const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::face_rotation (const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::line_orientation (const unsigned int /*line*/)
{
  return false;
}



template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::has_children ()
{
  return false;
}



template <int spacedim>
inline
unsigned int TriaAccessor<0, 1, spacedim>::n_children()
{
  return 0;
}



template <int spacedim>
inline
unsigned int TriaAccessor<0, 1, spacedim>::number_of_children ()
{
  return 0;
}



template <int spacedim>
inline
unsigned int TriaAccessor<0, 1, spacedim>::max_refinement_depth ()
{
  return 0;
}


template <int spacedim>
inline
TriaIterator<TriaAccessor<0,1,spacedim> >
TriaAccessor<0, 1, spacedim>::child (const unsigned int)
{
  return TriaIterator<TriaAccessor<0,1,spacedim> >();
}


template <int spacedim>
inline
TriaIterator<TriaAccessor<0,1,spacedim> >
TriaAccessor<0, 1, spacedim>::isotropic_child (const unsigned int)
{
  return TriaIterator<TriaAccessor<0,1,spacedim> >();
}


template <int spacedim>
inline
RefinementCase<0> TriaAccessor<0, 1, spacedim>::refinement_case ()
{
  return RefinementCase<0>(RefinementPossibilities<0>::no_refinement);
}

template <int spacedim>
inline
int TriaAccessor<0, 1, spacedim>::child_index (const unsigned int)
{
  return -1;
}


template <int spacedim>
inline
int TriaAccessor<0, 1, spacedim>::isotropic_child_index (const unsigned int)
{
  return -1;
}



template <int spacedim>
inline
void
TriaAccessor<0, 1, spacedim>::set_boundary_id (const types::boundary_id b)
{
  Assert (tria->vertex_to_boundary_id_map_1d->find (this->vertex_index())
          != tria->vertex_to_boundary_id_map_1d->end(),
          ExcInternalError());

  (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline
void
TriaAccessor<0, 1, spacedim>::set_manifold_id (const types::manifold_id b)
{
  (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline
void TriaAccessor<0, 1, spacedim>::set_all_boundary_ids (const types::boundary_id b)
{
  set_boundary_id (b);
}



template <int spacedim>
inline
void TriaAccessor<0, 1, spacedim>::set_all_manifold_ids (const types::manifold_id b)
{
  set_manifold_id (b);
}



template <int spacedim>
inline
bool TriaAccessor<0, 1, spacedim>::used () const
{
  return tria->vertex_used(global_vertex_index);
}

/*------------------------ Functions: CellAccessor<dim,spacedim> -----------------------*/


template <int dim, int spacedim>
inline
CellAccessor<dim,spacedim>::
CellAccessor (const Triangulation<dim,spacedim> *parent,
              const int                 level,
              const int                 index,
              const AccessorData       *local_data)
  :
  TriaAccessor<dim, dim, spacedim> (parent, level, index, local_data)
{}



template <int dim, int spacedim>
inline
CellAccessor<dim,spacedim>::CellAccessor (const TriaAccessor<dim,dim,spacedim> &cell_accessor)
  :
  TriaAccessor<dim, dim, spacedim> (static_cast<const TriaAccessor<dim, dim, spacedim>&>(cell_accessor))
{}



namespace internal
{
  namespace CellAccessorImplementation
  {
    template <int spacedim>
    inline
    dealii::TriaIterator<dealii::TriaAccessor<0, 1, spacedim> >
    get_face (const dealii::CellAccessor<1,spacedim> &cell,
              const unsigned int i)
    {
      dealii::TriaAccessor<0, 1, spacedim>
      a (&cell.get_triangulation(),
         ((i == 0) && cell.at_boundary(0)
          ?
          dealii::TriaAccessor<0, 1, spacedim>::left_vertex
          :
          ((i == 1) && cell.at_boundary(1)
           ?
           dealii::TriaAccessor<0, 1, spacedim>::right_vertex
           :
           dealii::TriaAccessor<0, 1, spacedim>::interior_vertex)),
         cell.vertex_index(i));
      return dealii::TriaIterator<dealii::TriaAccessor<0, 1, spacedim> >(a);
    }


    template <int spacedim>
    inline
    dealii::TriaIterator<dealii::TriaAccessor<1, 2, spacedim> >
    get_face (const dealii::CellAccessor<2,spacedim> &cell,
              const unsigned int i)
    {
      return cell.line(i);
    }


    template <int spacedim>
    inline
    dealii::TriaIterator<dealii::TriaAccessor<2, 3, spacedim> >
    get_face (const dealii::CellAccessor<3,spacedim> &cell,
              const unsigned int i)
    {
      return cell.quad(i);
    }
  }
}



template <int dim, int spacedim>
inline
TriaIterator<TriaAccessor<dim-1, dim, spacedim> >
CellAccessor<dim,spacedim>::face (const unsigned int i) const
{
  return dealii::internal::CellAccessorImplementation::get_face (*this, i);
}



template <int dim, int spacedim>
inline
unsigned int
CellAccessor<dim,spacedim>::face_index (const unsigned int i) const
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
inline
int
CellAccessor<dim,spacedim>::neighbor_index (const unsigned int i) const
{
  AssertIndexRange (i,GeometryInfo<dim>::faces_per_cell);
  return this->tria->levels[this->present_level]->
         neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].second;
}



template <int dim, int spacedim>
inline
int
CellAccessor<dim,spacedim>::neighbor_level (const unsigned int i) const
{
  AssertIndexRange (i, GeometryInfo<dim>::faces_per_cell);
  return this->tria->levels[this->present_level]->
         neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].first;
}



template <int dim, int spacedim>
inline
RefinementCase<dim>
CellAccessor<dim,spacedim>::refine_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for refinement must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert (this->active() ||  !this->tria->levels[this->present_level]->refine_flags[this->present_index],
          ExcRefineCellNotActive());
  return RefinementCase<dim>(this->tria->levels[this->present_level]->refine_flags[this->present_index]);
}



template <int dim, int spacedim>
inline
void
CellAccessor<dim,spacedim>::set_refine_flag (const RefinementCase<dim> refinement_case) const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  Assert (!coarsen_flag_set(),
          ExcCellFlaggedForCoarsening());

  this->tria->levels[this->present_level]->refine_flags[this->present_index] = refinement_case;
}



template <int dim, int spacedim>
inline
void
CellAccessor<dim,spacedim>::clear_refine_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    RefinementCase<dim>::no_refinement;
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::flag_for_face_refinement (const unsigned int face_no,
                                                      const RefinementCase<dim-1> &face_refinement_case) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));
  Assert (face_no<GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange(face_no,0,GeometryInfo<dim>::faces_per_cell));
  Assert (face_refinement_case < RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(face_refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));

  // the new refinement case is a combination
  // of the minimum required one for the given
  // face refinement and the already existing
  // flagged refinement case
  RefinementCase<dim> old_ref_case = refine_flag_set();
  RefinementCase<dim>
  new_ref_case = (old_ref_case
                  | GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(face_refinement_case,
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
inline
bool
CellAccessor<dim,spacedim>::flag_for_line_refinement (const unsigned int line_no) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));
  Assert (line_no<GeometryInfo<dim>::lines_per_cell,
          ExcIndexRange(line_no,0,GeometryInfo<dim>::lines_per_cell));

  // the new refinement case is a combination
  // of the minimum required one for the given
  // line refinement and the already existing
  // flagged refinement case
  RefinementCase<dim> old_ref_case=refine_flag_set(),
                      new_ref_case=old_ref_case
                                   | GeometryInfo<dim>::min_cell_refinement_case_for_line_refinement(line_no);
  set_refine_flag(new_ref_case);
  // return, whether we had to change the
  // refinement flag
  return new_ref_case != old_ref_case;
}



template <>
inline
dealii::internal::SubfaceCase<1>
CellAccessor<1>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}

template <>
inline
dealii::internal::SubfaceCase<1>
CellAccessor<1,2>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}


template <>
inline
dealii::internal::SubfaceCase<1>
CellAccessor<1,3>::subface_case(const unsigned int) const
{
  return dealii::internal::SubfaceCase<1>::case_none;
}


template <>
inline
dealii::internal::SubfaceCase<2>
CellAccessor<2>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<2>::faces_per_cell,
         ExcIndexRange(face_no,0,GeometryInfo<2>::faces_per_cell));
  return ((face(face_no)->has_children()) ?
          dealii::internal::SubfaceCase<2>::case_x :
          dealii::internal::SubfaceCase<2>::case_none);
}

template <>
inline
dealii::internal::SubfaceCase<2>
CellAccessor<2,3>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<2>::faces_per_cell,
         ExcIndexRange(face_no,0,GeometryInfo<2>::faces_per_cell));
  return ((face(face_no)->has_children()) ?
          dealii::internal::SubfaceCase<2>::case_x :
          dealii::internal::SubfaceCase<2>::case_none);
}


template <>
inline
dealii::internal::SubfaceCase<3>
CellAccessor<3>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<3>::faces_per_cell,
         ExcIndexRange(face_no,0,GeometryInfo<3>::faces_per_cell));
  switch (static_cast<std::uint8_t> (face(face_no)->refinement_case()))
    {
    case RefinementCase<3>::no_refinement:
      return dealii::internal::SubfaceCase<3>::case_none;
    case RefinementCase<3>::cut_x:
      if (face(face_no)->child(0)->has_children())
        {
          Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<2>::cut_y,
                 ExcInternalError());
          if (face(face_no)->child(1)->has_children())
            {
              Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_y,
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
              Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_y,
                     ExcInternalError());
              return dealii::internal::SubfaceCase<3>::case_x2y;
            }
          else
            return dealii::internal::SubfaceCase<3>::case_x;
        }
    case RefinementCase<3>::cut_y:
      if (face(face_no)->child(0)->has_children())
        {
          Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<2>::cut_x,
                 ExcInternalError());
          if (face(face_no)->child(1)->has_children())
            {
              Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_x,
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
              Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_x,
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
inline
bool
CellAccessor<dim,spacedim>::coarsen_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for coarsening must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert (this->active() ||  !this->tria->levels[this->present_level]->coarsen_flags[this->present_index],
          ExcRefineCellNotActive());
  return this->tria->levels[this->present_level]->coarsen_flags[this->present_index];
}



template <int dim, int spacedim>
inline
void
CellAccessor<dim,spacedim>::set_coarsen_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  Assert (!refine_flag_set(), ExcCellFlaggedForRefinement());

  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] = true;
}



template <int dim, int spacedim>
inline
void
CellAccessor<dim,spacedim>::clear_coarsen_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] = false;
}



template <int dim, int spacedim>
inline
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim,spacedim>::neighbor (const unsigned int i) const
{
  TriaIterator<CellAccessor<dim,spacedim> >
  q (this->tria, neighbor_level (i), neighbor_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
          ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim,spacedim>::child (const unsigned int i) const
{
  TriaIterator<CellAccessor<dim,spacedim> >
  q (this->tria, this->present_level+1, this->child_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
          ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::active () const
{
  return !this->has_children();
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_locally_owned () const
{
  Assert (this->active(),
          ExcMessage("is_locally_owned() can only be called on active cells!"));
#ifndef DEAL_II_WITH_MPI
  return true;
#else
  if (is_artificial())
    return false;

  const parallel::Triangulation<dim,spacedim> *pt
    = dynamic_cast<const parallel::Triangulation<dim,spacedim> *>(this->tria);

  if (pt == nullptr)
    return true;
  else
    return (this->subdomain_id() == pt->locally_owned_subdomain());

#endif
}


template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_locally_owned_on_level () const
{

#ifndef DEAL_II_WITH_MPI
  return true;
#else

  const parallel::Triangulation<dim,spacedim> *pt
    = dynamic_cast<const parallel::Triangulation<dim,spacedim> *>(this->tria);

  if (pt == nullptr)
    return true;
  else
    return (this->level_subdomain_id() == pt->locally_owned_subdomain());

#endif
}


template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_ghost () const
{
  Assert (this->active(),
          ExcMessage("is_ghost() can only be called on active cells!"));
  if (is_artificial() || this->has_children())
    return false;

#ifndef DEAL_II_WITH_MPI
  return false;
#else

  const parallel::Triangulation<dim,spacedim> *pt
    = dynamic_cast<const parallel::Triangulation<dim,spacedim> *>(this->tria);

  if (pt == nullptr)
    return false;
  else
    return (this->subdomain_id() != pt->locally_owned_subdomain());

#endif
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_artificial () const
{
  Assert (this->active(),
          ExcMessage("is_artificial() can only be called on active cells!"));
#ifndef DEAL_II_WITH_MPI
  return false;
#else

  const parallel::Triangulation<dim,spacedim> *pt
    = dynamic_cast<const parallel::Triangulation<dim,spacedim> *>(this->tria);

  if (pt == nullptr)
    return false;
  else
    return this->subdomain_id() == numbers::artificial_subdomain_id;

#endif
}



template <int dim, int spacedim>
inline
types::subdomain_id
CellAccessor<dim, spacedim>::subdomain_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (this->active(),
          ExcMessage("subdomain_id() can only be called on active cells!"));
  return this->tria->levels[this->present_level]->subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
inline
unsigned int
CellAccessor<dim,spacedim>::neighbor_face_no (const unsigned int neighbor) const
{
  const unsigned int n2=neighbor_of_neighbor_internal(neighbor);
  if (n2!=numbers::invalid_unsigned_int)
    // return this value as the
    // neighbor is not coarser
    return n2;
  else
    // the neighbor is coarser
    return neighbor_of_coarser_neighbor(neighbor).first;
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_level_cell()
{
  return false;
}


DEAL_II_NAMESPACE_CLOSE

#endif
