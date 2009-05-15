//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_accessor_templates_h
#define __deal2__dof_accessor_templates_h


#include <base/config.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_levels.h>
#include <dofs/dof_faces.h>
#include <hp/dof_levels.h>
#include <hp/dof_faces.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, class DH>
DoFAccessor<structdim,DH>::
DoFAccessor ()
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
inline
DoFAccessor<structdim,DH>::
DoFAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
	     const int                 level,
	     const int                 index,
	     const DH                 *dof_handler)
		:
		internal::DoFAccessor::Inheritance<structdim,DH::dimension,DH::space_dimension>::BaseClass (tria,
													    level,
													    index),
                dof_handler(const_cast<DH*>(dof_handler))
{}



template <int structdim, class DH>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim,DH>::
DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
template <int dim2, class DH2>
DoFAccessor<structdim,DH>::
DoFAccessor (const DoFAccessor<dim2, DH2> &)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::set_dof_handler (DH *dh)
{
  Assert (dh != 0, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, class DH>
inline
const DH &
DoFAccessor<structdim,DH>::get_dof_handler () const
{
  return *this->dof_handler;
}


template <int structdim, class DH>
inline
DoFAccessor<structdim,DH> &
DoFAccessor<structdim,DH>::operator = (const DoFAccessor<structdim,DH> &da)
{
  this->set_dof_handler (da.dof_handler);
  return *this;
}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::copy_from (const DoFAccessor<structdim,DH> &a)
{
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
}



template <int structdim, class DH>
inline
bool
DoFAccessor<structdim,DH>::operator == (const DoFAccessor<structdim,DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator == (a));
}



template <int structdim, class DH>
inline
bool
DoFAccessor<structdim,DH>::operator != (const DoFAccessor<structdim,DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator != (a));
}



template <int structdim, class DH>
inline
TriaIterator<DoFAccessor<structdim,DH> >
DoFAccessor<structdim,DH>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  int next_level;
  if (DH::dimension==structdim)
    next_level = this->level()+1;
  else
    next_level = 0;
  
  TriaIterator<DoFAccessor<structdim,DH> > q (this->tria,
					      next_level,
					      this->child_index (i),
					      this->dof_handler);
  
				   // make sure that we either created
				   // a past-the-end iterator or one
				   // pointing to a used cell
  Assert ((q.state() == IteratorState::past_the_end)
	  ||
	  q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

  return q;
}



namespace internal
{
  namespace DoFAccessor
  {
				     /**
				      * A class like the one with same
				      * name in tria.cc. See there for
				      * more information.
				      */
    struct Implementation
    {
					 /**
					  * Implementations of the
					  * get_dof_index/set_dof_index functions.
					  */
	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<1,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>)
	  {
	    return dof_handler.levels[obj_level]->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<1,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>,
		       const unsigned int              global_index)
	  {
	    dof_handler.levels[obj_level]->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }
    

	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    return dof_handler.faces->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>,
		       const unsigned int              global_index)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    dof_handler.faces->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }
    
    
	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<2>)
	  {
	    return dof_handler.levels[obj_level]->quads.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<2>,
		       const unsigned int              global_index)
	  {
	    dof_handler.levels[obj_level]->quads.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    return dof_handler.faces->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<1>,
		       const unsigned int              global_index)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    dof_handler.faces->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }



	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<2>)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    return dof_handler.faces->quads.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<2>,
		       const unsigned int              global_index)
	  {
					     // faces have no levels
	    Assert (obj_level == 0, ExcInternalError());
	    dof_handler.faces->quads.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }



	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<3>)
	  {
	    return dof_handler.levels[obj_level]->hexes.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index);	
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
		       const unsigned int              obj_level,
		       const unsigned int              obj_index,
		       const unsigned int              fe_index,
		       const unsigned int              local_index,
		       internal::int2type<3>,
		       const unsigned int              global_index)
	  {
	    dof_handler.levels[obj_level]->hexes.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &)
	  {
	    return dof_handler.levels[obj_level]->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.levels[obj_level]->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.faces->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);  
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<2> &)
	  {
	    return dof_handler.levels[obj_level]->quads.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<2> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.levels[obj_level]->quads.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<1> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.faces->lines.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);  
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<2> &)
	  {
	    return dof_handler.faces->quads.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<2> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.faces->quads.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<3> &)
	  {
	    return dof_handler.levels[obj_level]->hexes.
	      get_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     obj_level);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
		       const unsigned int       obj_level,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const internal::int2type<3> &,
		       const unsigned int       global_index)
	  {
	    dof_handler.levels[obj_level]->hexes.
	      set_dof_index (dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index,
			     obj_level);  
	  }


	template <int structdim, int dim, int spacedim>
	static
	bool
	fe_index_is_active (const dealii::DoFHandler<dim,spacedim> &,
			    const unsigned int,
			    const unsigned int,
			    const unsigned int fe_index,
			    const internal::int2type<structdim> &)
	  {
	    return (fe_index == 0);
	  }



	template <int structdim, int dim, int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::DoFHandler<dim,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const internal::int2type<structdim> &)
	  {
					     // check that the object we look
					     // at is in fact active. the
					     // problem is that we have
					     // templatized on the
					     // dimensionality of the object,
					     // so it may be a cell, a face,
					     // or a line. we have a bit of
					     // trouble doing this all in the
					     // generic case, so only check if
					     // it is either a cell or a
					     // line. the only case this
					     // leaves out is faces in 3d --
					     // let's hope that this never is
					     // a problem
	    Assert ((dim==structdim
		     ?
		     typename
		     dealii::DoFHandler<dim,spacedim>::
		     raw_cell_iterator (&dof_handler.get_tria(),
					obj_level,
					obj_index,
					&dof_handler)->used()
		     :
		     (structdim==1
		      ?
		      typename
		      dealii::DoFHandler<dim,spacedim>::
		      raw_line_iterator (&dof_handler.get_tria(),
					 obj_level,
					 obj_index,
					 &dof_handler)->used()
		      :
		      true))
		    == true,
		    ExcMessage ("This cell is not active and therefore can't be "
				"queried for its active FE indices"));
	    return 1;
	  }



	template <int structdim, int dim, int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::DoFHandler<dim,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<structdim> &)
	  {
					     // check that the object we look
					     // at is in fact active. the
					     // problem is that we have
					     // templatized on the
					     // dimensionality of the object,
					     // so it may be a cell, a face,
					     // or a line. we have a bit of
					     // trouble doing this all in the
					     // generic case, so only check if
					     // it is either a cell or a
					     // line. the only case this
					     // leaves out is faces in 3d --
					     // let's hope that this never is
					     // a problem
	    Assert ((dim==structdim
		     ?
		     typename
		     dealii::DoFHandler<dim,spacedim>::
		     raw_cell_iterator (&dof_handler.get_tria(),
					obj_level,
					obj_index,
					&dof_handler)->used()
		     :
		     (structdim==1
		      ?
		      typename
		      dealii::DoFHandler<dim,spacedim>::
		      raw_line_iterator (&dof_handler.get_tria(),
					 obj_level,
					 obj_index,
					 &dof_handler)->used()
		      :
		      true))
		    == true,
		    ExcMessage ("This cell is not active and therefore can't be "
				"queried for its active FE indices"));
	    Assert (n == 0, ExcIndexRange (n, 0, 1));
  
	    return dealii::DoFHandler<dim,spacedim>::default_fe_index;
	  }
	

	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<1> &)
	  {
	    return dof_handler.levels[obj_level]->lines.fe_index_is_active(dof_handler,
									   obj_index,
									   fe_index,
									   obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.levels[obj_level]->lines.n_active_fe_indices (dof_handler,
									     obj_index);
	  }



	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.levels[obj_level]->lines.nth_active_fe_index (dof_handler,
									     obj_level,
									     obj_index,
									     n);
	  }


	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.fe_index_is_active(dof_handler,
							       obj_index,
							       fe_index,
							       obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			     const unsigned int ,
			     const unsigned int obj_index,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.n_active_fe_indices (dof_handler,
								 obj_index);
	  }


	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.nth_active_fe_index (dof_handler,
								 obj_level,
								 obj_index,
								 n);
	  }

  

	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<2> &)
	  {
	    return dof_handler.levels[obj_level]->quads.fe_index_is_active(dof_handler,
									   obj_index,
									   fe_index,
									   obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const internal::int2type<2> &)
	  {
	    return dof_handler.levels[obj_level]->quads.n_active_fe_indices (dof_handler,
									     obj_index);
	  }



	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<2> &)
	  {
	    return dof_handler.levels[obj_level]->quads.nth_active_fe_index (dof_handler,
									     obj_level,
									     obj_index,
									     n);
	  }



	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.fe_index_is_active(dof_handler,
							       obj_index,
							       fe_index,
							       obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int,
			     const unsigned int obj_index,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.n_active_fe_indices (dof_handler,
								 obj_index);
	  }



	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<1> &)
	  {
	    return dof_handler.faces->lines.nth_active_fe_index (dof_handler,
								 obj_level,
								 obj_index,
								 n);
	  }



	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<2> &)
	  {
	    return dof_handler.faces->quads.fe_index_is_active(dof_handler,
							       obj_index,
							       fe_index,
							       obj_level);
	  }

	template <int spacedim>
	static
	bool
	fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			    const unsigned int obj_level,
			    const unsigned int obj_index,
			    const unsigned int fe_index,
			    const internal::int2type<3> &)
	  {
	    return dof_handler.levels[obj_level]->hexes.fe_index_is_active(dof_handler,
									   obj_index,
									   fe_index,
									   obj_level);
	  }


	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int ,
			     const unsigned int obj_index,
			     const internal::int2type<2> &)
	  {
	    return dof_handler.faces->quads.n_active_fe_indices (dof_handler,
								 obj_index);
	  }



	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<2> &)
	  {
	    return dof_handler.faces->quads.nth_active_fe_index (dof_handler,
								 obj_level,
								 obj_index,
								 n);
	  }



	template <int spacedim>
	static
	unsigned int
	n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const internal::int2type<3> &)
	  {
	    return dof_handler.levels[obj_level]->hexes.n_active_fe_indices (dof_handler,
									     obj_index);
	  }



	template <int spacedim>
	static
	unsigned int
	nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
			     const unsigned int obj_level,
			     const unsigned int obj_index,
			     const unsigned int n,
			     const internal::int2type<3> &)
	  {
	    return dof_handler.levels[obj_level]->hexes.nth_active_fe_index (dof_handler,
									     obj_level,
									     obj_index,
									     n);
	  }

					 /**
					  * Set the @p local_index-th
					  * degree of freedom
					  * corresponding to the finite
					  * element specified by @p
					  * fe_index on the vertex with
					  * global number @p
					  * vertex_index to @p
					  * global_index.
					  */
	template <int dim, int spacedim>
	static
	void
	set_vertex_dof_index (dealii::DoFHandler<dim,spacedim> &dof_handler,
			      const unsigned int vertex_index,
			      const unsigned int fe_index,
			      const unsigned int local_index,
			      const unsigned int global_index)
	  {
	    Assert ((fe_index == dealii::DoFHandler<dim,spacedim>::default_fe_index),
		    ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
	    Assert (dof_handler.selected_fe != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (local_index < dof_handler.selected_fe->dofs_per_vertex,
		    ExcIndexRange(local_index, 0,
				  dof_handler.selected_fe->dofs_per_vertex));

	    dof_handler.vertex_dofs[vertex_index *
				    dof_handler.selected_fe->dofs_per_vertex
				    + local_index]
	      = global_index;
	  }


	template <int dim, int spacedim>
	static
	void
	set_vertex_dof_index (dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
			      const unsigned int vertex_index,
			      const unsigned int fe_index,
			      const unsigned int local_index,
			      const unsigned int global_index)
	  {
	    Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
		     ExcMessage ("You need to specify a FE index when working "
				 "with hp DoFHandlers"));
	    Assert (dof_handler.finite_elements != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (local_index < (*dof_handler.finite_elements)[fe_index].dofs_per_vertex,
		    ExcIndexRange(local_index, 0,
				  (*dof_handler.finite_elements)[fe_index].dofs_per_vertex));
	    Assert (fe_index < dof_handler.finite_elements->size(),
		    ExcInternalError());
	    Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
		    numbers::invalid_unsigned_int,
		    ExcMessage ("This vertex is unused and has no DoFs associated with it"));

					     // hop along the list of index
					     // sets until we find the one
					     // with the correct fe_index, and
					     // then poke into that
					     // part. trigger an exception if
					     // we can't find a set for this
					     // particular fe_index
	    const unsigned int starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
	    unsigned int *pointer              = &dof_handler.vertex_dofs[starting_offset];
	    while (true)
	      {
		Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

		const unsigned int this_fe_index = *pointer;
	
		Assert (this_fe_index != numbers::invalid_unsigned_int,
			ExcInternalError());
		Assert (this_fe_index < dof_handler.finite_elements->size(),
			ExcInternalError());

		if (this_fe_index == fe_index)
		  {
		    *(pointer + 1 + local_index) = global_index;
		    return;
		  }
		else
		  pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
	      }
	  }
	
	
					 /**
					  * Get the @p local_index-th
					  * degree of freedom
					  * corresponding to the finite
					  * element specified by @p
					  * fe_index on the vertex with
					  * global number @p
					  * vertex_index to @p
					  * global_index.
					  */

	template <int dim, int spacedim>
	static
	unsigned int
	get_vertex_dof_index (const dealii::DoFHandler<dim,spacedim> &dof_handler,
			      const unsigned int vertex_index,
			      const unsigned int fe_index,
			      const unsigned int local_index)
	  {
	    Assert ((fe_index == dealii::DoFHandler<dim,spacedim>::default_fe_index),
		    ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
	    Assert (dof_handler.selected_fe != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (local_index < dof_handler.selected_fe->dofs_per_vertex,
		    ExcIndexRange(local_index, 0,
				  dof_handler.selected_fe->dofs_per_vertex));

	    return
	      dof_handler.vertex_dofs[vertex_index *
				      dof_handler.selected_fe->dofs_per_vertex
				      + local_index];
	  }  


	template<int dim, int spacedim>
	static
	unsigned int
	get_vertex_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
			      const unsigned int vertex_index,
			      const unsigned int fe_index,
			      const unsigned int local_index)
	  {
	    Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
		     ExcMessage ("You need to specify a FE index when working "
				 "with hp DoFHandlers"));
	    Assert (dof_handler.finite_elements != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (local_index < (*dof_handler.finite_elements)[fe_index].dofs_per_vertex,
		    ExcIndexRange(local_index, 0,
				  (*dof_handler.finite_elements)[fe_index].dofs_per_vertex));
	    Assert (vertex_index < dof_handler.vertex_dofs_offsets.size(),
		    ExcIndexRange (vertex_index, 0,
				   dof_handler.vertex_dofs_offsets.size()));
	    Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
		    numbers::invalid_unsigned_int,
		    ExcMessage ("This vertex is unused and has no DoFs associated with it"));

					     // hop along the list of index
					     // sets until we find the one
					     // with the correct fe_index, and
					     // then poke into that
					     // part. trigger an exception if
					     // we can't find a set for this
					     // particular fe_index
	    const unsigned int starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
	    const unsigned int *pointer        = &dof_handler.vertex_dofs[starting_offset];
	    while (true)
	      {
		Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

		const unsigned int this_fe_index = *pointer;
	
		Assert (this_fe_index != numbers::invalid_unsigned_int,
			ExcInternalError());
		Assert (this_fe_index < dof_handler.finite_elements->size(),
			ExcInternalError());

		if (this_fe_index == fe_index)
		  return *(pointer + 1 + local_index);
		else
		  pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
	      }
	  }


					 /**
					  * Return the number of
					  * different finite elements
					  * that are active on a given
					  * vertex.
					  */
	template<int dim, int spacedim>
	static
	unsigned int
	n_active_vertex_fe_indices (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
				    const unsigned int vertex_index)
	  {
	    Assert (dof_handler.finite_elements != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));

					     // if this vertex is unused, return 0
	    if (dof_handler.vertex_dofs_offsets[vertex_index] == numbers::invalid_unsigned_int)
	      return 0;

					     // hop along the list of index
					     // sets and count the number of
					     // hops
	    const unsigned int starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
	    const unsigned int *pointer        = &dof_handler.vertex_dofs[starting_offset];

	    Assert (*pointer != numbers::invalid_unsigned_int,
		    ExcInternalError());
    
	    unsigned int counter = 0;
	    while (true)
	      {
		Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

		const unsigned int this_fe_index = *pointer;
	
		if (this_fe_index == numbers::invalid_unsigned_int)
		  return counter;
		else
		  {
		    pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
		    ++counter;
		  }
	      }
	  }



					 /**
					  * Return the fe index of the
					  * n-th finite element active
					  * on a given vertex.
					  */
	template<int dim, int spacedim>
	static
	unsigned int
	nth_active_vertex_fe_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
				    const unsigned int vertex_index,
				    const unsigned int n)
	  {
	    Assert (dof_handler.finite_elements != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (n < n_active_vertex_fe_indices(dof_handler, vertex_index),
		    ExcIndexRange (n, 0, n_active_vertex_fe_indices(dof_handler,
								    vertex_index)));
					     // make sure we don't ask on
					     // unused vertices
	    Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
		    numbers::invalid_unsigned_int,
		    ExcInternalError());

					     // hop along the list of index
					     // sets and count the number of
					     // hops
	    const unsigned int starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
	    const unsigned int *pointer        = &dof_handler.vertex_dofs[starting_offset];

	    Assert (*pointer != numbers::invalid_unsigned_int,
		    ExcInternalError());
    
	    unsigned int counter = 0;
	    while (true)
	      {
		Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

		const unsigned int this_fe_index = *pointer;
	
		Assert (this_fe_index < dof_handler.finite_elements->size(),
			ExcInternalError());

		if (counter == n)
		  return this_fe_index;

		Assert (this_fe_index != numbers::invalid_unsigned_int,
			ExcInternalError());
	
		pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
		++counter;	
	      }
	  }
  


					 /**
					  * Return whether a particular
					  * finite element index is
					  * active on the specified
					  * vertex.
					  */
	template<int dim, int spacedim>
	static
	bool
	fe_is_active_on_vertex (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
				const unsigned int vertex_index,
				const unsigned int fe_index)
	  {
	    Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
		     ExcMessage ("You need to specify a FE index when working "
				 "with hp DoFHandlers"));
	    Assert (dof_handler.finite_elements != 0,
		    ExcMessage ("No finite element collection is associated with "
				"this DoFHandler"));
	    Assert (fe_index < dof_handler.finite_elements->size(),
		    ExcInternalError());

					     // make sure we don't ask on
					     // unused vertices
	    Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
		    numbers::invalid_unsigned_int,
		    ExcInternalError());

					     // hop along the list of index
					     // sets and see whether we find
					     // the given index
	    const unsigned int starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
	    const unsigned int *pointer        = &dof_handler.vertex_dofs[starting_offset];

	    Assert (*pointer != numbers::invalid_unsigned_int,
		    ExcInternalError());
    
	    while (true)
	      {
		Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

		const unsigned int this_fe_index = *pointer;
	
		Assert (this_fe_index < dof_handler.finite_elements->size(),
			ExcInternalError());

		if (this_fe_index == numbers::invalid_unsigned_int)
		  return false;
		else
		  if (this_fe_index == fe_index)
		    return true;
		  else
		    pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
	      }
	  }
  
    };
  }
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::dof_index (const unsigned int i,
				const unsigned int fe_index) const
{
				   // access the respective DoF 
  return internal::DoFAccessor::Implementation::get_dof_index (*this->dof_handler,
							       this->level(),
							       this->present_index,
							       fe_index,
							       i,
							       internal::int2type<dim>());
}



template <int dim, class DH>
inline
void
DoFAccessor<dim,DH>::set_dof_index (const unsigned int i,
				    const unsigned int index,
				    const unsigned int fe_index) const
{
				   // access the respective DoF
  internal::DoFAccessor::Implementation::set_dof_index (*this->dof_handler,
							this->level(),
							this->present_index,
							fe_index,
							i,
							internal::int2type<dim>(),
							index);
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::n_active_fe_indices () const
{
				   // access the respective DoF
  return
    internal::DoFAccessor::Implementation::
    n_active_fe_indices (*this->dof_handler,
			      this->level(),
			      this->present_index,
			      internal::int2type<dim>());
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::nth_active_fe_index (const unsigned int n) const
{
				   // access the respective DoF
  return
    internal::DoFAccessor::Implementation::
    nth_active_fe_index (*this->dof_handler,
			      this->level(),
			      this->present_index,
			      n,
			      internal::int2type<dim>());
}



template <int dim, class DH>
inline
bool
DoFAccessor<dim,DH>::fe_index_is_active (const unsigned int fe_index) const
{
				   // access the respective DoF
  return
    internal::DoFAccessor::Implementation::
    fe_index_is_active (*this->dof_handler,
			this->level(),
			this->present_index,
			fe_index,
			internal::int2type<dim>());
}



template <int structdim, class DH>
inline
unsigned int
DoFAccessor<structdim, DH>::vertex_dof_index (const unsigned int vertex,
					      const unsigned int i,
					      const unsigned int fe_index) const
{
  return
    internal::DoFAccessor::Implementation::get_vertex_dof_index
    (*this->dof_handler,
     this->vertex_index(vertex),
     fe_index,
     i);
}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim, DH>::set_vertex_dof_index (const unsigned int vertex,
						  const unsigned int i,
						  const unsigned int index,
						  const unsigned int fe_index) const
{
  internal::DoFAccessor::Implementation::set_vertex_dof_index
    (*this->dof_handler,
     this->vertex_index(vertex),
     fe_index,
     i,
     index);
}






namespace internal
{
  namespace DoFAccessor
  {
    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const FiniteElement<dim,spacedim> &fe,
	    const unsigned int)
    {
      return fe;
    }


  
    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const dealii::hp::FECollection<dim,spacedim> &fe,
	    const unsigned int                            index)
    {
      return fe[index];
    }
  }
}


template <int dim, class DH>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFAccessor<dim,DH>::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));
  
  return internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), fe_index);
}



namespace internal
{
  namespace DoFAccessor
  {
    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<1,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<2; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
      for (unsigned int d=0; d<dofs_per_line; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }


    
    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<2,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
			 dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<4; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
				       // now copy dof numbers from the line. for
				       // lines with the wrong orientation (which
				       // might occur in 3d), we have already made
				       // sure that we're ok by picking the correct
				       // vertices (this happens automatically in
				       // the vertex() function). however, if the
				       // line is in wrong orientation, we look at
				       // it in flipped orientation and we will have
				       // to adjust the shape function indices that
				       // we see to correspond to the correct
				       // (face-local) ordering.
      for (unsigned int line=0; line<4; ++line)
	for (unsigned int d=0; d<dofs_per_line; ++d)
	  *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
						   adjust_line_dof_index_for_line_orientation(d,
											      accessor.line_orientation(line)),
						   fe_index);
      for (unsigned int d=0; d<dofs_per_quad; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }

    

    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<3,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
			 dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad,
			 dofs_per_hex    = accessor.get_fe(fe_index).dofs_per_hex;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<8; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
				       // now copy dof numbers from the line. for
				       // lines with the wrong orientation, we have
				       // already made sure that we're ok by picking
				       // the correct vertices (this happens
				       // automatically in the vertex()
				       // function). however, if the line is in
				       // wrong orientation, we look at it in
				       // flipped orientation and we will have to
				       // adjust the shape function indices that we
				       // see to correspond to the correct
				       // (cell-local) ordering.
      for (unsigned int line=0; line<12; ++line)
	for (unsigned int d=0; d<dofs_per_line; ++d)
	  *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
						   adjust_line_dof_index_for_line_orientation(d,
											      accessor.line_orientation(line)),fe_index);
				       // now copy dof numbers from the face. for
				       // faces with the wrong orientation, we
				       // have already made sure that we're ok by
				       // picking the correct lines and vertices
				       // (this happens automatically in the
				       // line() and vertex() functions). however,
				       // if the face is in wrong orientation, we
				       // look at it in flipped orientation and we
				       // will have to adjust the shape function
				       // indices that we see to correspond to the
				       // correct (cell-local) ordering. The same
				       // applies, if the face_rotation or
				       // face_orientation is non-standard
      for (unsigned int quad=0; quad<6; ++quad)
	for (unsigned int d=0; d<dofs_per_quad; ++d)
	  *next++ = accessor.quad(quad)->dof_index(accessor.get_fe(fe_index).
						   adjust_quad_dof_index_for_face_orientation(d,
											      accessor.face_orientation(quad),
											      accessor.face_flip(quad),
											      accessor.face_rotation(quad)),
						   fe_index);
      for (unsigned int d=0; d<dofs_per_hex; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }
  }
}


template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::get_dof_indices (std::vector<unsigned int> &dof_indices,
					    const unsigned int         fe_index) const
{
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  switch (structdim)
    {
      case 1:
	    Assert (dof_indices.size() ==
		    (2*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     this->dof_handler->get_fe()[fe_index].dofs_per_line),
		    ExcVectorDoesNotMatch());
	    break;
      case 2:
	    Assert (dof_indices.size() ==
		    (4*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     4*this->dof_handler->get_fe()[fe_index].dofs_per_line +
		     this->dof_handler->get_fe()[fe_index].dofs_per_quad),
		    ExcVectorDoesNotMatch());
	    break;
      case 3:
	    Assert (dof_indices.size() ==
		    (8*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     12*this->dof_handler->get_fe()[fe_index].dofs_per_line +
		     6*this->dof_handler->get_fe()[fe_index].dofs_per_quad +
		     this->dof_handler->get_fe()[fe_index].dofs_per_hex),
		    ExcVectorDoesNotMatch());
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  

				   // this function really only makes
				   // sense if either a) there are
				   // degrees of freedom defined on
				   // the present object, or b) the
				   // object is non-active objects but
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (this->fe_index_is_active (fe_index)
	  ||
	  (this->dof_handler->get_fe()[fe_index].dofs_per_cell ==
	   GeometryInfo<structdim>::vertices_per_cell *
	   this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
	  ExcInternalError());

				   // now do the actual work
  internal::DoFAccessor::get_dof_indices (*this, dof_indices, fe_index);
}





template <int structdim, class DH>
inline
typename internal::DoFHandler::Iterators<DH>::line_iterator
DoFAccessor<structdim,DH>::line (const unsigned int i) const
{
  Assert (structdim > 1, ExcImpossibleInDim(structdim));
				   // checking of 'i' happens in
				   // line_index(i)
  
  return typename internal::DoFHandler::Iterators<DH>::line_iterator
    (
      this->tria,
      0,  // only sub-objects are allowed, which have no level
      this->line_index(i),
      this->dof_handler
    );
}


template <int structdim, class DH>
inline
typename internal::DoFHandler::Iterators<DH>::quad_iterator
DoFAccessor<structdim,DH>::quad (const unsigned int i) const
{
  Assert (structdim > 2, ExcImpossibleInDim(structdim));
				   // checking of 'i' happens in
				   // quad_index(i)

  return typename internal::DoFHandler::Iterators<DH>::quad_iterator
    (
      this->tria,
      0,  // only sub-objects are allowed, which have no level
      this->quad_index (i),
      this->dof_handler
    );
}





/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessor
  {
				     // make sure we refer to class
				     // dealii::DoFCellAccessor, not
				     // namespace
				     // internal::DoFCellAccessor
    using dealii::DoFCellAccessor;
    using dealii::DoFHandler;
    
/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
    struct Implementation
    {
					 /**
					  * Implement the updating of the
					  * cache. Currently not
					  * implemented for hp::DoFHandler
					  * objects.
					  */
	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<1,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<1>::vertices_per_cell))
	      return;
	
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;
	
					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());
	
	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);
	
	    for (unsigned int vertex=0; vertex<2; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
	    for (unsigned int d=0; d<dofs_per_line; ++d)
	      *next++ = accessor.dof_index(d);
	  }



	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<2,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<2>::vertices_per_cell))
	      return;
  
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_quad   = accessor.get_fe().dofs_per_quad,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;

					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());

	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

	    for (unsigned int vertex=0; vertex<4; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
	    for (unsigned int line=0; line<4; ++line)
	      for (unsigned int d=0; d<dofs_per_line; ++d)
		*next++ = accessor.line(line)->dof_index(d);
	    for (unsigned int d=0; d<dofs_per_quad; ++d)
	      *next++ = accessor.dof_index(d);
	  }



	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<3,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<3>::vertices_per_cell))
	      return;
  
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_quad   = accessor.get_fe().dofs_per_quad,
			       dofs_per_hex    = accessor.get_fe().dofs_per_hex,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;

					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());

	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

	    for (unsigned int vertex=0; vertex<8; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
					     // now copy dof numbers from the line. for
					     // lines with the wrong orientation, we have
					     // already made sure that we're ok by picking
					     // the correct vertices (this happens
					     // automatically in the vertex()
					     // function). however, if the line is in
					     // wrong orientation, we look at it in
					     // flipped orientation and we will have to
					     // adjust the shape function indices that we
					     // see to correspond to the correct
					     // (cell-local) ordering.
	    for (unsigned int line=0; line<12; ++line)
	      for (unsigned int d=0; d<dofs_per_line; ++d)
		*next++ = accessor.line(line)->dof_index(accessor.dof_handler->get_fe().
							 adjust_line_dof_index_for_line_orientation(d,
												    accessor.line_orientation(line)));
					     // now copy dof numbers from the face. for
					     // faces with the wrong orientation, we
					     // have already made sure that we're ok by
					     // picking the correct lines and vertices
					     // (this happens automatically in the
					     // line() and vertex() functions). however,
					     // if the face is in wrong orientation, we
					     // look at it in flipped orientation and we
					     // will have to adjust the shape function
					     // indices that we see to correspond to the
					     // correct (cell-local) ordering. The same
					     // applies, if the face_rotation or
					     // face_orientation is non-standard
	    for (unsigned int quad=0; quad<6; ++quad)
	      for (unsigned int d=0; d<dofs_per_quad; ++d)
		*next++ = accessor.quad(quad)->dof_index(accessor.dof_handler->get_fe().
							 adjust_quad_dof_index_for_face_orientation(d,
												    accessor.face_orientation(quad),
												    accessor.face_flip(quad),
												    accessor.face_rotation(quad)));
	    for (unsigned int d=0; d<dofs_per_hex; ++d)
	      *next++ = accessor.dof_index(d);
	  }


					 // implementation for the case of
					 // hp::DoFHandler objects. it's
					 // not implemented there, for no
					 // space dimension
	template <int dim, int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &)
	  {
//TODO[WB]: should implement a dof indices cache for hp as well
  
					     // not implemented, but should also
					     // not be called
	    Assert (false, ExcNotImplemented());
	  }


					 /**
					  * A function that collects the
					  * global indices of degrees of
					  * freedom. This function works
					  * for ::DoFHandler and all
					  * template arguments and copies
					  * the data out of the cache that
					  * we hold for each cell.
					  */
	template <int dim, int spacedim>
	static
	void
	get_dof_indices (const DoFCellAccessor<DoFHandler<dim,spacedim> > &accessor,
			 std::vector<unsigned int>                        &dof_indices)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (dof_indices.size() == accessor.get_fe().dofs_per_cell,
		    typename BaseClass::ExcVectorDoesNotMatch());

					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices
	    Assert (!accessor.has_children()
		    ||
		    (accessor.get_fe().dofs_per_cell ==
		     accessor.get_fe().dofs_per_vertex * GeometryInfo<dim>::vertices_per_cell),
		    ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
	    unsigned int *cache = &accessor.dof_handler->levels[accessor.level()]
				  ->cell_dof_indices_cache[accessor.present_index * accessor.get_fe().dofs_per_cell]; 
	    for (unsigned int i=0; i<accessor.get_fe().dofs_per_cell; ++i, ++cache)
	      dof_indices[i] = *cache;
	  }      

					 /**
					  * Same function as above except
					  * that it works for
					  * hp::DoFHandler objects that do
					  * not have a cache for the local
					  * DoF indices.
					  */
	template <int dim, int spacedim>
	static
	void
	get_dof_indices (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			 std::vector<unsigned int>                            &dof_indices)
	  {
					     // no caching for hp::DoFHandler implemented 
	    accessor.dealii::DoFAccessor<dim,dealii::hp::DoFHandler<dim,spacedim> >::get_dof_indices (dof_indices,
										      accessor.active_fe_index());
	  }


					 /**
					  * A function that collects the
					  * values of degrees of freedom. This
					  * function works for ::DoFHandler
					  * and all template arguments and
					  * uses the data from the cache of
					  * indices that we hold for each
					  * cell.
					  */
	template <int dim, int spacedim, class InputVector, typename number>
	static
	void
	get_dof_values (const DoFCellAccessor<DoFHandler<dim,spacedim> > &accessor,
			const InputVector &values,
			dealii::Vector<number> &local_values)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (local_values.size() == accessor.get_fe().dofs_per_cell,
		    typename BaseClass::ExcVectorDoesNotMatch());
	    Assert (values.size() == accessor.get_dof_handler().n_dofs(),
		    typename BaseClass::ExcVectorDoesNotMatch());

					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices
	    Assert (!accessor.has_children()
		    ||
		    (accessor.get_fe().dofs_per_cell ==
		     accessor.get_fe().dofs_per_vertex * GeometryInfo<dim>::vertices_per_cell),
		    ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
	    unsigned int *cache = &accessor.dof_handler->levels[accessor.level()]
				  ->cell_dof_indices_cache[accessor.present_index * accessor.get_fe().dofs_per_cell]; 
	    for (unsigned int i=0; i<accessor.get_fe().dofs_per_cell; ++i, ++cache)
	      local_values(i) = values(*cache);
	  }      

					 /**
					  * Same function as above except
					  * that it works for
					  * hp::DoFHandler objects that do
					  * not have a cache for the local
					  * DoF indices.
					  */
	template <int dim, int spacedim, class InputVector, typename number>
	static
	void
	get_dof_values (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			const InputVector &values,
			dealii::Vector<number>    &local_values)
	  {
					     // no caching for hp::DoFHandler
					     // implemented
	    const unsigned int dofs_per_cell = accessor.get_fe().dofs_per_cell;
  
	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    get_dof_indices (accessor, local_dof_indices);

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      local_values(i) = values(local_dof_indices[i]);
	  }


					 /**
					  * Same set of functions as above
					  * except that it sets rather than
					  * gets values
					  */
	template <int dim, int spacedim, class OutputVector, typename number>
	static
	void
	set_dof_values (const DoFCellAccessor<DoFHandler<dim,spacedim> > &accessor,
			const dealii::Vector<number> &local_values,
			OutputVector &values)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (local_values.size() == accessor.get_fe().dofs_per_cell,
		    typename BaseClass::ExcVectorDoesNotMatch());
	    Assert (values.size() == accessor.get_dof_handler().n_dofs(),
		    typename BaseClass::ExcVectorDoesNotMatch());

					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices
	    Assert (!accessor.has_children()
		    ||
		    (accessor.get_fe().dofs_per_cell ==
		     accessor.get_fe().dofs_per_vertex * GeometryInfo<dim>::vertices_per_cell),
		    ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
	    unsigned int *cache = &accessor.dof_handler->levels[accessor.level()]
				  ->cell_dof_indices_cache[accessor.present_index * accessor.get_fe().dofs_per_cell]; 
	    for (unsigned int i=0; i<accessor.get_fe().dofs_per_cell; ++i, ++cache)
	      values(*cache) = local_values(i);
	  }      

	template <int dim, int spacedim, class OutputVector, typename number>
	static
	void
	set_dof_values (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			const dealii::Vector<number> &local_values,
			OutputVector &values)
	  {
					     // no caching for hp::DoFHandler
					     // implemented
	    const unsigned int dofs_per_cell = accessor.get_fe().dofs_per_cell;
  
	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    get_dof_indices (accessor, local_dof_indices);

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      values(local_dof_indices[i]) = local_values(i);
	  }
	

					 /**
					  * Do what the active_fe_index
					  * function in the parent class
					  * is supposed to do.
					  */
	template <int dim, int spacedim>
	static
	unsigned int
	active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim> > &)
	  {
					     // ::DoFHandler only supports a
					     // single active fe with index
					     // zero
	    return 0;
	  }



	template <int dim, int spacedim>
	static
	unsigned int
	active_fe_index (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor)
	  {
	    Assert (static_cast<unsigned int>(accessor.level()) < accessor.dof_handler->levels.size(),
		    ExcMessage ("DoFHandler not initialized"));
	    Assert (static_cast<std::vector<unsigned int>::size_type>(accessor.present_index) <
		    accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size (),
		    ExcIndexRange (accessor.present_index, 0,
				   accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size ()));
	    return accessor.dof_handler->levels[accessor.level()]
	      ->active_fe_indices[accessor.present_index];
	  }



					 /**
					  * Do what the
					  * set_active_fe_index function
					  * in the parent class is
					  * supposed to do.
					  */
	template <int dim, int spacedim>
	static
	void
	set_active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim> > &,
			     const unsigned int                                i)
	  {
					     // ::DoFHandler only supports a
					     // single active fe with index
					     // zero
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (i == 0, typename BaseClass::ExcInvalidObject());
	  }



	template <int dim, int spacedim>
	static
	void
	set_active_fe_index (DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			     const unsigned int                              i)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (accessor.dof_handler != 0,
		    typename BaseClass::ExcInvalidObject());
	    Assert (static_cast<unsigned int>(accessor.level()) <
		    accessor.dof_handler->levels.size(),
		    ExcMessage ("DoFHandler not initialized"));
	    Assert (static_cast<std::vector<unsigned int>::size_type>(accessor.present_index) <
		    accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size (),
		    ExcIndexRange (accessor.present_index, 0,
				   accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size ()));
	    accessor.dof_handler->levels[accessor.level()]
	      ->active_fe_indices[accessor.present_index] = i;
	  }
    };
  }
}




template <class DH>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
                 const int                 level,
                 const int                 index,
                 const AccessorData       *local_data)
                :
                DoFAccessor<DH::dimension,DH> (tria,level,index,local_data)
{}



template <class DH>
template <int structdim2, int dim2, int spacedim2>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, typename BaseClass::ExcInvalidObject());
}



template <class DH>
template <int dim2, class DH2>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const DoFAccessor<dim2,DH2> &)
{
  Assert (false, typename BaseClass::ExcInvalidObject());
}


template <class DH>
inline
typename internal::DoFHandler::Iterators<DH>::cell_iterator
DoFCellAccessor<DH>::neighbor (const unsigned int i) const
{
  typename internal::DoFHandler::Iterators<DH>::cell_iterator
    q (this->tria,
       this->neighbor_level (i),
       this->neighbor_index (i),
       this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsNeighbor());
#endif
  return q;
}


template <class DH>
inline
typename internal::DoFHandler::Iterators<DH>::cell_iterator
DoFCellAccessor<DH>::child (const unsigned int i) const
{
  typename internal::DoFHandler::Iterators<DH>::cell_iterator
    q (this->tria,
       this->level()+1,
       this->child_index (i),
       this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsChild());
#endif
  return q;
}



template <class DH>
typename internal::DoFHandler::Iterators<DH>::face_iterator
DoFCellAccessor<DH>::face (const unsigned int i) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell, ExcIndexRange (i, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dim = DH::dimension;  
  Assert (dim > 1, ExcImpossibleInDim(1));

  switch (dim)
    {
      case 2:
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator
	      (this->tria,
	       0,
	       this->line_index (i),
	       this->dof_handler);
	    
      case 3:
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator
	      (this->tria,
	       0,
	       this->quad_index (i),
	       this->dof_handler);
	    
      default:
	    Assert (false, ExcNotImplemented());
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator();
    }
}



template <class DH>
inline
void
DoFCellAccessor<DH>::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  internal::DoFCellAccessor::Implementation::get_dof_indices (*this, dof_indices);
}



template <class DH>
template <class InputVector, typename number>
void
DoFCellAccessor<DH>::get_dof_values (const InputVector &values,
				     Vector<number>    &local_values) const
{
  internal::DoFCellAccessor::Implementation
    ::get_dof_values (*this, values, local_values);
}



template <class DH>
template <class OutputVector, typename number>
void
DoFCellAccessor<DH>::set_dof_values (const Vector<number> &local_values,
				     OutputVector         &values) const
{
  internal::DoFCellAccessor::Implementation
    ::set_dof_values (*this, local_values, values);
}



template <class DH>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFCellAccessor<DH>::get_fe () const
{
  return internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), active_fe_index());
}



template <class DH>
inline
unsigned int
DoFCellAccessor<DH>::active_fe_index () const
{
  return internal::DoFCellAccessor::Implementation::active_fe_index (*this);
}



template <class DH>
inline
void
DoFCellAccessor<DH>::set_active_fe_index (const unsigned int i)
{
  internal::DoFCellAccessor::Implementation::set_active_fe_index (*this, i);
}




template <class DH>
template <typename number, typename OutputVector>
inline
void
DoFCellAccessor<DH>::
distribute_local_to_global (const Vector<number> &local_source,
			    OutputVector         &global_destination) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.size() == this->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.size(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

//TODO[WB]: This function could me made more efficient. First, it allocates memory, which could be avoided by passing in another argument as a scratch array (or by using the dof_indices cache for non-hp DoFHandlers, see the functions in the Implementation class above). second, the elementwise access is really slow if we use PETSc vectors/matrices. This should be fixed eventually
  
				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  this->get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
}



template <class DH>
template <typename number, typename OutputMatrix>
inline
void
DoFCellAccessor<DH>::
distribute_local_to_global (const FullMatrix<number> &local_source,
			    OutputMatrix             &global_destination) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.m() == this->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.m(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());

  const unsigned int n_dofs = local_source.m();

//TODO[WB]: This function could me made more efficient. First, it allocates memory, which could be avoided by passing in another argument as a scratch array (or by using the dof_indices cache for non-hp DoFHandlers, see the functions in the Implementation class above). second, the elementwise access is really slow if we use PETSc vectors/matrices. This should be fixed eventually

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  this->get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
}


DEAL_II_NAMESPACE_CLOSE

#endif
