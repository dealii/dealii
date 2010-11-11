//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_accessor_templates_h
#define __deal2__tria_accessor_templates_h


#include <base/config.h>
#include <base/geometry_info.h>
#include <base/template_constraints.h>
#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_faces.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.templates.h>
#include <distributed/tria.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


/*------------------------ Functions: TriaAccessorBase ---------------------------*/

template <int structdim, int dim, int spacedim>
inline
TriaAccessorBase<structdim,dim,spacedim>::TriaAccessorBase (
  const Triangulation<dim,spacedim>* tria,
  const int                          level,
  const int                          index,
  const AccessorData*)
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
TriaAccessorBase<structdim,dim,spacedim>&
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
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return ((present_level == a.present_level) &&
	  (present_index == a.present_index));
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessorBase<structdim,dim,spacedim>::operator != (const TriaAccessorBase<structdim,dim,spacedim> &a) const
{
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return ((present_level != a.present_level) ||
	  (present_index != a.present_index));
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
  else
    if (present_index==-1)
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
  namespace TriaAccessorBase
  {
				     /**
				      * Out of a face object, get the
				      * sub-objects of dimensionality
				      * given by the last argument.
				      */
    template <int dim>
    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> >*
    get_objects (internal::Triangulation::TriaFaces<dim> *faces,
		 const internal::int2type<1>)
    {
      return &faces->lines;
    }


    template <int dim>
    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> >*
    get_objects (internal::Triangulation::TriaFaces<dim> *faces,
		 const internal::int2type<2>)
    {
      return &faces->quads;
    }

    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> >*
    get_objects (internal::Triangulation::TriaFaces<1>*,
		 const internal::int2type<1>)
    {
      Assert (false, ExcInternalError());
      return 0;
    }

    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> >*
    get_objects (internal::Triangulation::TriaFaces<2>*,
		 const internal::int2type<2>)
    {
      Assert (false, ExcInternalError());
      return 0;
    }

    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<3> >*
    get_objects (internal::Triangulation::TriaFaces<3>*,
		 const internal::int2type<3>)
    {
      Assert (false, ExcInternalError());
      return 0;
    }

				     /**
				      * This function should never be
				      * used, but we need it for the
				      * template instantiation of TriaAccessorBase<dim,dim,spacedim>::objects() const
				      */
    template <int dim>
    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<3> >*
    get_objects (internal::Triangulation::TriaFaces<dim> *faces,
		 const internal::int2type<3>)
    {
      Assert (false, ExcInternalError());
      return 0;
    }

				     /**
				      * Copy the above functions for
				      * cell objects.
				      */
    template <int structdim, int dim>
    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<structdim> >*
    get_objects (internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<dim> >*,
		 const internal::int2type<structdim>)
    {
      Assert (false, ExcInternalError());
      return 0;
    }

    template <int dim>
    inline
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<dim> >*
    get_objects (internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<dim> >* cells,
		 const internal::int2type<dim>)
    {
      return cells;
    }
  }
}



template <int structdim, int dim, int spacedim>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<structdim> > &
TriaAccessorBase<structdim,dim,spacedim>::objects() const
{
  if (structdim != dim)
				     // get sub-objects. note that the
				     // current class is only used for
				     // objects that are *not* cells
    return *internal::TriaAccessorBase::get_objects (this->tria->faces,
						     internal::int2type<structdim> ());
  else
    return *internal::TriaAccessorBase::get_objects (&this->tria->levels[this->present_level]->cells,
						     internal::int2type<structdim> ());
}



/*------------------------ Functions: InvalidAccessor ---------------------------*/

template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::
InvalidAccessor (const Triangulation<dim,spacedim> *,
		 const int                 ,
		 const int                 ,
		 const AccessorData       *)
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



/*------------------------ Functions: TriaAccessor ---------------------------*/


namespace internal
{
  namespace TriaAccessor
  {
				     // make sure that if in the following we
				     // write TriaAccessor
				     // we mean the *class*
				     // dealii::TriaAccessor, not the
				     // enclosing namespace
				     // internal::TriaAccessor
    using dealii::TriaAccessor;

/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
    struct Implementation
    {
					 /**
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int dim, int spacedim>
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
	static
	unsigned int
	line_index (const TriaAccessor<2, dim, spacedim> &accessor,
		    const unsigned int i)
	  {
	    return accessor.objects().cells[accessor.present_index].face(i);
	  }


	template <int dim, int spacedim>
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
	      { { 4, 0 }, // take first four lines from bottom face
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
		{ 1, 1 }};

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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	bool
	face_orientation (const TriaAccessor<structdim, dim, spacedim> &,
			  const unsigned int)
	  {
	    return true;
	  }


	template <int dim, int spacedim>
	static
	bool
	face_orientation (const TriaAccessor<3, dim, spacedim> &accessor,
			  const unsigned int face)
	  {
	    return (accessor.tria->levels[accessor.present_level]
		    ->cells.face_orientation(accessor.present_index, face));
	  }



					 /**
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	bool
	face_flip (const TriaAccessor<structdim, dim, spacedim> &,
		   const unsigned int)
	  {
	    return true;
	  }


	template <int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	bool
	face_rotation (const TriaAccessor<structdim, dim, spacedim> &,
		       const unsigned int)
	  {
	    return true;
	  }


	template <int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int dim, int spacedim>
	static
	bool
	line_orientation (const TriaAccessor<1, dim, spacedim> &,
			  const unsigned int)
	  {
	    return true;
	  }


	template <int spacedim>
	static
	bool
	line_orientation (const TriaAccessor<2, 2, spacedim> &,
			  const unsigned int)
	  {
					     // quads in 2d have no
					     // non-standard orientation
	    return true;
	  }


	template <int spacedim>
	static
	bool
	line_orientation (const TriaAccessor<2, 3, spacedim> &accessor,
			  const unsigned int line)
	  {
					     // quads as part of 3d hexes
					     // can have non-standard
					     // orientation
//TODO: why is this face_orientation, not line_orientation as in the setter function?
	    return accessor.tria->faces->quads.face_orientation(accessor.present_index, line);
	  }


	template <int dim, int spacedim>
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
	      { { 4, 0 }, // take first four lines from bottom face
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
		{ 1, 1 }};

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
					     // with respect to the face corrsponds to a
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
	      { { { { true, false },   // lines 0/1, face_orientation=false, face_flip=false, face_rotation=false and true
		    { false, true }},  // lines 0/1, face_orientation=false, face_flip=true, face_rotation=false and true
		  { { true, true },    // lines 0/1, face_orientation=true, face_flip=false, face_rotation=false and true
		    { false, false }}},// linea 0/1, face_orientation=true, face_flip=true, face_rotation=false and true

		{ { { true, true },    // lines 2/3 ...
		    { false, false }},
		  { { true, false },
		    { false, true }}}};


	    return (accessor.quad(quad_index)
		    ->line_orientation(line_index)
		    == bool_table[std_line_index/2]
		    [accessor.face_orientation(quad_index)]
		    [accessor.face_flip(quad_index)]
		    [accessor.face_rotation(quad_index)]);
	  }



					 /**
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	void
	set_face_orientation (const TriaAccessor<structdim, dim, spacedim> &,
			      const unsigned int,
			      const bool)
	  {
	    Assert (false, ExcInternalError());
	  }


	template <int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	void
	set_face_flip (const TriaAccessor<structdim, dim, spacedim> &,
		       const unsigned int,
		       const bool)
	  {
	    Assert (false, ExcInternalError());
	  }


	template <int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int structdim, int dim, int spacedim>
	static
	void
	set_face_rotation (const TriaAccessor<structdim, dim, spacedim> &,
			   const unsigned int,
			   const bool)
	  {
	    Assert (false, ExcInternalError());
	  }


	template <int dim, int spacedim>
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
					  * Implementation of the function
					  * of some name in the mother
					  * class.
					  */
	template <int dim, int spacedim>
	static
	void
	set_line_orientation (const TriaAccessor<1, dim, spacedim> &,
			      const unsigned int,
			      const bool)
	  {
	    Assert (false, ExcInternalError());
	  }


	template <int spacedim>
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
					  * Implementation of the function of same
					  * name in the enclosing class.
					  */
	template <int dim, int spacedim>
	static
	unsigned int
	vertex_index (const TriaAccessor<1,dim,spacedim> &accessor,
		      const unsigned int corner)
	  {
	    return accessor.objects().cells[accessor.present_index].face (corner);
	  }


	template <int dim, int spacedim>
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

	    static const unsigned int switch_table[2][2]={{1,0},{0,1}};

	    return accessor.line(corner%2)
	      ->vertex_index(switch_table[accessor.line_orientation(corner%2)][corner/2]);
	  }



	template <int dim, int spacedim>
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
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return this->objects().used[this->present_index];
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::
vertex_index (const unsigned int corner) const
{
  Assert (corner<GeometryInfo<structdim>::vertices_per_cell,
	  ExcIndexRange(corner,0,GeometryInfo<structdim>::vertices_per_cell));

  return internal::TriaAccessor::Implementation::vertex_index (*this, corner);
}



template <int structdim, int dim, int spacedim>
Point<spacedim> &
TriaAccessor<structdim, dim, spacedim>::vertex (const unsigned int i) const
{
  return const_cast<Point<spacedim> &> (this->tria->vertices[vertex_index(i)]);
}



template <int structdim, int dim, int spacedim>
inline
typename internal::Triangulation::Iterators<dim,spacedim>::line_iterator
TriaAccessor<structdim,dim,spacedim>::line (const unsigned int i) const
{
				   // checks happen in line_index
  return typename internal::Triangulation::Iterators<dim,spacedim>::line_iterator
    (this->tria, 0, line_index (i));
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::line_index (const unsigned int i) const
{
  Assert (i < GeometryInfo<structdim>::lines_per_cell,
	  ExcIndexRange (i, 0, GeometryInfo<structdim>::lines_per_cell));

  return internal::TriaAccessor::Implementation::line_index (*this, i);
}




template <int structdim, int dim, int spacedim>
inline
typename internal::Triangulation::Iterators<dim,spacedim>::quad_iterator
TriaAccessor<structdim,dim,spacedim>::quad (const unsigned int i) const
{
				   // checks happen in quad_index
  return typename internal::Triangulation::Iterators<dim,spacedim>::quad_iterator
    (this->tria, 0, quad_index (i));
}



template <int structdim, int dim, int spacedim>
inline
unsigned int
TriaAccessor<structdim,dim,spacedim>::quad_index (const unsigned int i) const
{
  return internal::TriaAccessor::Implementation::quad_index (*this, i);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_orientation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return internal::TriaAccessor::Implementation::face_orientation (*this, face);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_flip (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return internal::TriaAccessor::Implementation::face_flip (*this, face);
}


template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::face_rotation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return internal::TriaAccessor::Implementation::face_rotation (*this, face);
}



template <int structdim, int dim, int spacedim>
inline
bool
TriaAccessor<structdim,dim,spacedim>::line_orientation (const unsigned int line) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (line<GeometryInfo<structdim>::lines_per_cell,
	  ExcIndexRange (line, 0, GeometryInfo<structdim>::lines_per_cell));

  return internal::TriaAccessor::Implementation::line_orientation (*this, line);
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_orientation (const unsigned int face,
							     const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  internal::TriaAccessor::Implementation::set_face_orientation (*this, face, value);
}



template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_flip (const unsigned int face,
						      const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  internal::TriaAccessor::Implementation::set_face_flip (*this, face, value);
}


template <int structdim, int dim, int spacedim>
inline
void
TriaAccessor<structdim,dim,spacedim>::set_face_rotation (const unsigned int face,
							  const bool value) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  internal::TriaAccessor::Implementation::set_face_rotation (*this, face, value);
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

  internal::TriaAccessor::Implementation::set_line_orientation (*this, line, value);
}



template <int structdim, int dim, int spacedim>
void TriaAccessor<structdim, dim, spacedim>::set_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  this->objects().used[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  this->objects().used[this->present_index] = false;
}


template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::
parent_index () const
{
  Assert (this->present_level > 0, TriaAccessorExceptions::ExcCellHasNoParent ());

  					// the parent of two consecutive cells
  					// is stored only once, since it is
  					// the same
  return this->tria->levels[this->present_level]->parents[this->present_index / 2];
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
	  this_refinement_case (static_cast<unsigned char>(refinement_case()));

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
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());

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
		     static_cast<unsigned char>(RefinementCase<1>::cut_x) :
		     static_cast<unsigned char>(RefinementCase<1>::no_refinement)));

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
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

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
	  this_refinement_case (static_cast<unsigned char>(refinement_case()));

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
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());

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
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
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
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
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

  Assert ((index==-1) ||
	  (i==0 && !this->has_children() && (index>=0)) ||
	  (i>0  &&  this->has_children() && (index>=0) &&
	   this->objects().children[n_sets_of_two*this->present_index+i/2] == -1),
	  TriaAccessorExceptions::ExcCantSetChildren(index));

  this->objects().children[n_sets_of_two*this->present_index+i/2] = index;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_parent (const unsigned int parent_index)
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (this->present_level > 0, TriaAccessorExceptions::ExcCellHasNoParent ());
  this->tria->levels[this->present_level]->parents[this->present_index / 2]
    = parent_index;
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
  this->objects().user_pointer(this->present_index) = 0;
}



template <int structdim, int dim, int spacedim>
void * TriaAccessor<structdim,dim,spacedim>::user_pointer () const
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
unsigned char
TriaAccessor<structdim, dim, spacedim>::boundary_indicator () const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects().material_id[this->present_index];
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_boundary_indicator (const unsigned char boundary_ind) const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().material_id[this->present_index] = boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set_all_boundary_indicators (const unsigned char boundary_ind) const
{
  set_boundary_indicator (boundary_ind);
  set_boundary_indicator (boundary_ind);

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
	      this->line(i)->set_boundary_indicator (boundary_ind);
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
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
}



template <int structdim, int dim, int spacedim>
const Boundary<dim,spacedim> &
TriaAccessor<structdim, dim, spacedim>::get_boundary () const
{
  Assert (structdim<dim, ExcImpossibleInDim(dim));
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->tria->get_boundary(this->objects()
				  .material_id[this->present_index]);
}



template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::diameter () const
{
  switch (structdim)
    {
      case 1:
	    return std::sqrt((this->vertex(1)-this->vertex(0)).square());
      case 2:
	    return std::sqrt(std::max((this->vertex(3)-this->vertex(0)).square(),
				      (this->vertex(2)-this->vertex(1)).square()));
      case 3:
	    return std::sqrt(std::max( std::max((this->vertex(7)-this->vertex(0)).square(),
						(this->vertex(6)-this->vertex(1)).square()),
				       std::max((this->vertex(2)-this->vertex(5)).square(),
						(this->vertex(3)-this->vertex(4)).square()) ));
      default:
	    Assert (false, ExcNotImplemented());
	    return -1e10;
    }
}


template <int structdim, int dim, int spacedim>
Point<spacedim>
TriaAccessor<structdim, dim, spacedim>::center () const
{
  Point<spacedim> p;
  for (unsigned int v=0; v<GeometryInfo<structdim>::vertices_per_cell; ++v)
    p += vertex(v);
  return p/GeometryInfo<structdim>::vertices_per_cell;
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
  const Point<spacedim> dist = o->vertex(0) - this->vertex(0);
  const double tol_square = 1e-24 * dist.norm_square();
  for (unsigned int i=1; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
    {
      const Point<spacedim> dist_new = (o->vertex(i) - this->vertex(i)) - dist;
      if (dist_new.norm_square() > tol_square)
	{
	  is_translation = false;
	  break;
	}
    }
  return is_translation;
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
CellAccessor<dim,spacedim>::
CellAccessor (const CellAccessor &cell_accessor)
                :
		TriaAccessor<dim, dim, spacedim> (static_cast<const TriaAccessor<dim, dim, spacedim>&>(cell_accessor))
{}



namespace internal
{
  namespace CellAccessor
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
    dealii::TriaIterator<dealii::TriaAccessor<1, 2, spacedim> >
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
  return internal::CellAccessor::get_face (*this, i);
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
	Assert (false, ExcImpossibleInDim(1));
	return numbers::invalid_unsigned_int;
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
internal::SubfaceCase<1>
CellAccessor<1>::subface_case(const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
  return internal::SubfaceCase<1>::case_none;
}

template <>
inline
internal::SubfaceCase<1>
CellAccessor<1,2>::subface_case(const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
  return internal::SubfaceCase<1>::case_none;
}


template <>
inline
internal::SubfaceCase<2>
CellAccessor<2>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<2>::faces_per_cell,
	 ExcIndexRange(face_no,0,GeometryInfo<2>::faces_per_cell));
  return (face(face_no)->has_children()) ? internal::SubfaceCase<2>::case_x : internal::SubfaceCase<2>::case_none;
}

template <>
inline
internal::SubfaceCase<2>
CellAccessor<2,3>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<2>::faces_per_cell,
	 ExcIndexRange(face_no,0,GeometryInfo<2>::faces_per_cell));
  return (face(face_no)->has_children()) ? internal::SubfaceCase<2>::case_x : internal::SubfaceCase<2>::case_none;
}


template <>
inline
internal::SubfaceCase<3>
CellAccessor<3>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<3>::faces_per_cell,
	 ExcIndexRange(face_no,0,GeometryInfo<3>::faces_per_cell));
  switch (static_cast<unsigned char> (face(face_no)->refinement_case()))
    {
      case RefinementCase<3>::no_refinement:
	    return internal::SubfaceCase<3>::case_none;
	    break;
      case RefinementCase<3>::cut_x:
	    if (face(face_no)->child(0)->has_children())
	      {
 		Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<2>::cut_y,
 		       ExcInternalError());
		if (face(face_no)->child(1)->has_children())
		  {
 		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_y,
 			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_x1y2y;
		  }
		else
		  return internal::SubfaceCase<3>::case_x1y;
	      }
	    else
	      {
		if (face(face_no)->child(1)->has_children())
		  {
  		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_y,
  			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_x2y;
		  }
		else
		  return internal::SubfaceCase<3>::case_x;
	      }
	    break;
      case RefinementCase<3>::cut_y:
	    if (face(face_no)->child(0)->has_children())
	      {
 		Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<2>::cut_x,
 		       ExcInternalError());
		if (face(face_no)->child(1)->has_children())
		  {
 		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_x,
 			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_y1x2x;
		  }
		else
		  return internal::SubfaceCase<3>::case_y1x;
	      }
	    else
	      {
		if (face(face_no)->child(1)->has_children())
		  {
  		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<2>::cut_x,
  			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_y2x;
		  }
		else
		  return internal::SubfaceCase<3>::case_y;
	      }
	    break;
      case RefinementCase<3>::cut_xy:
	    return internal::SubfaceCase<3>::case_xy;
	    break;
      default:
	    Assert(false, ExcInternalError());
    }
				   // we should never get here
  return internal::SubfaceCase<3>::case_none;
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
	  TriaAccessorExceptions::ExcUnusedCellAsNeighbor());

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
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

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
CellAccessor<dim,spacedim>::is_ghost () const
{
#ifndef DEAL_II_USE_P4EST
  return false;
#else
  const types::subdomain_id_t subdomain = this->subdomain_id();
  if (subdomain == types::artificial_subdomain_id)
    return false;

  const parallel::distributed::Triangulation<dim,spacedim> *pdt
    = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim> *>(this->tria);

  if (pdt == 0)
    return false;
  else
    return (subdomain != pdt->locally_owned_subdomain());
#endif
}



template <int dim, int spacedim>
inline
bool
CellAccessor<dim,spacedim>::is_artificial () const
{
#ifndef DEAL_II_USE_P4EST
  return false;
#else
  return (this->subdomain_id() == types::artificial_subdomain_id);
#endif
}



template <int dim, int spacedim>
inline
unsigned int
CellAccessor<dim,spacedim>::neighbor_face_no (const unsigned int neighbor) const
{
  if (dim==1)
    return neighbor_of_neighbor(neighbor);
  else
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
}

DEAL_II_NAMESPACE_CLOSE

#endif
