//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <dofs/dof_levels.h>
#include <dofs/dof_faces.h>
#include <dofs/dof_constraints.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_levels.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <base/geometry_info.h>
#include <fe/fe.h>
#include <lac/sparse_matrix.h>
#include <base/exceptions.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN



/* ------------------------ MGVertexDoFs ----------------------------------- */

template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::MGVertexDoFs ()
		:
		coarsest_level (deal_II_numbers::invalid_unsigned_int),
		finest_level (0),
		indices (0)
{}


template <int dim>
void MGDoFHandler<dim>::MGVertexDoFs::init (const unsigned int cl,
					    const unsigned int fl,
					    const unsigned int dofs_per_vertex)
{
  if (indices != 0)
    {
      delete[] indices;
      indices = 0;
    }

  coarsest_level = cl;
  finest_level   = fl;

				   // in the case where an invalid
				   // entry is requested, just leave
				   // everything else alone
  if (cl > fl)
    return;
  
  const unsigned int n_levels = finest_level-coarsest_level + 1;
  
  indices = new unsigned int[n_levels * dofs_per_vertex];
  Assert (indices != 0, ExcNoMemory ());

  for (unsigned int i=0; i<n_levels*dofs_per_vertex; ++i)
    indices[i] = DoFHandler<dim>::invalid_dof_index;
}


template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::~MGVertexDoFs ()
{
  delete[] indices;
}


template <int dim>
typename MGDoFHandler<dim>::MGVertexDoFs &
MGDoFHandler<dim>::MGVertexDoFs::operator = (const MGVertexDoFs &)
{
  Assert (false,
	  ExcMessage ("We don't know how many dofs per vertex there are, so there's "
		      "no way we can copy the 'indices' array"));
  return *this;
}



template <int dim>
unsigned int MGDoFHandler<dim>::MGVertexDoFs::get_coarsest_level () const {
  return coarsest_level;
}


template <int dim>
unsigned int MGDoFHandler<dim>::MGVertexDoFs::get_finest_level () const {
  return finest_level;
}


/* ------------------------ MGDoFHandler ------------------------------------- */

template <int dim>
const unsigned int MGDoFHandler<dim>::dimension;


template <int dim>
MGDoFHandler<dim>::MGDoFHandler (const Triangulation<dim> &tria) :
		DoFHandler<dim> (tria),
		mg_faces (NULL)
{}


template <int dim>
MGDoFHandler<dim>::~MGDoFHandler ()
{
  clear ();
}


#if deal_II_dimension == 1

template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::begin_raw (const unsigned int level) const
{
  return begin_raw_line (level);
}


template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::begin (const unsigned int level) const
{
  return begin_line (level);
}


template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::begin_active (const unsigned int level) const
{
  return begin_active_line (level);
}


template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::end () const
{
  return end_line ();
}


template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw () const
{
  return last_raw_line ();
}


template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw (const unsigned int level) const
{
  return last_raw_line (level);
}


template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last () const
{
  return last_line ();
}


template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last (const unsigned int level) const
{
  return last_line (level);
}


template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active () const
{
  return last_active_line ();
}


template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active (const unsigned int level) const
{
  return last_active_line (level);
}


template <>
internal::MGDoFHandler::Iterators<1>::raw_face_iterator
MGDoFHandler<1>::begin_raw_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::face_iterator
MGDoFHandler<1>::begin_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::active_face_iterator
MGDoFHandler<1>::begin_active_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::raw_face_iterator
MGDoFHandler<1>::end_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::raw_face_iterator
MGDoFHandler<1>::end_raw_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::active_face_iterator
MGDoFHandler<1>::end_active_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::raw_face_iterator
MGDoFHandler<1>::last_raw_face () const 
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::face_iterator
MGDoFHandler<1>::last_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
internal::MGDoFHandler::Iterators<1>::active_face_iterator
MGDoFHandler<1>::last_active_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::begin_raw_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::begin_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::begin_active_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::end_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::end_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::end_raw_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::end_active_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::begin_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::begin_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::begin_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::end_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::end_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::end_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::end_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::last_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::last_raw_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::last_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::last_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::last_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}


template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::last_active_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}

#endif


#if deal_II_dimension == 2

template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::begin_raw (const unsigned int level) const
{
  return begin_raw_quad (level);
}


template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::begin (const unsigned int level) const
{
  return begin_quad (level);
}


template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::begin_active (const unsigned int level) const
{
  return begin_active_quad (level);
}


template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::end () const
{
  return end_quad ();
}


template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw () const
{
  return last_raw_quad ();
}


template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw (const unsigned int level) const
{
  return last_raw_quad (level);
}


template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last () const
{
  return last_quad ();
}


template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last (const unsigned int level) const
{
  return last_quad (level);
}


template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active () const
{
  return last_active_quad ();
}


template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active (const unsigned int level) const
{
  return last_active_quad (level);
}


template <>
internal::MGDoFHandler::Iterators<2>::raw_face_iterator
MGDoFHandler<2>::begin_raw_face () const
{
  return begin_raw_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::face_iterator
MGDoFHandler<2>::begin_face () const
{
  return begin_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::active_face_iterator
MGDoFHandler<2>::begin_active_face () const
{
  return begin_active_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::raw_face_iterator
MGDoFHandler<2>::end_face () const
{
  return end_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::raw_face_iterator
MGDoFHandler<2>::end_raw_face () const
{
  return end_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::active_face_iterator
MGDoFHandler<2>::end_active_face () const
{
  return end_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::raw_face_iterator
MGDoFHandler<2>::last_raw_face () const
{
  return last_raw_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::face_iterator
MGDoFHandler<2>::last_face () const
{
  return last_line ();
}


template <>
internal::MGDoFHandler::Iterators<2>::active_face_iterator
MGDoFHandler<2>::last_active_face () const
{
  return last_active_line ();
}


template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::begin_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::begin_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::begin_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::end_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::end_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::end_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::end_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::last_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::last_raw_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::last_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::last_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::last_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::last_active_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}


#endif


#if deal_II_dimension == 3

template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::begin_raw (const unsigned int level) const
{
  return begin_raw_hex (level);
}


template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::begin (const unsigned int level) const
{
  return begin_hex (level);
}


template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::begin_active (const unsigned int level) const
{
  return begin_active_hex (level);
}


template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::end () const
{
  return end_hex ();
}


template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::last_raw () const
{
  return last_raw_hex ();
}


template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::last_raw (const unsigned int level) const
{
  return last_raw_hex (level);
}


template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::last () const
{
  return last_hex ();
}


template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::last (const unsigned int level) const
{
  return last_hex (level);
}


template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::last_active () const
{
  return last_active_hex ();
}


template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::last_active (const unsigned int level) const
{
  return last_active_hex (level);
}


template <>
internal::MGDoFHandler::Iterators<3>::raw_face_iterator
MGDoFHandler<3>::begin_raw_face () const
{
  return begin_raw_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::face_iterator
MGDoFHandler<3>::begin_face () const
{
  return begin_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::active_face_iterator
MGDoFHandler<3>::begin_active_face () const
{
  return begin_active_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::raw_face_iterator
MGDoFHandler<3>::end_face () const
{
  return end_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::raw_face_iterator
MGDoFHandler<3>::end_raw_face () const
{
  return end_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::active_face_iterator
MGDoFHandler<3>::end_active_face () const
{
  return end_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::raw_face_iterator
MGDoFHandler<3>::last_raw_face () const
{
  return last_raw_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::face_iterator
MGDoFHandler<3>::last_face () const
{
  return last_quad ();
}


template <>
internal::MGDoFHandler::Iterators<3>::active_face_iterator
MGDoFHandler<3>::last_active_face () const
{
  return last_active_quad ();
}


#endif


template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::begin_raw_line (const unsigned int level) const
{
  Assert (dim==1 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());
  return raw_line_iterator (this->tria,
			    dim==1 ? this->tria->begin_raw_line(level)->level() : 0,
			    this->tria->begin_raw_line(level)->index(),
			    this);
}


template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::begin_line (const unsigned int level) const
{
  Assert (dim==1 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());
  return line_iterator (this->tria,
			dim==1 ? this->tria->begin_line(level)->level() : 0,
			this->tria->begin_line(level)->index(),
			this);
}


template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::begin_active_line (const unsigned int level) const
{
  Assert (dim==1 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());
  return active_line_iterator (this->tria,
			       dim==1 ? this->tria->begin_active_line(level)->level() : 0,
			       this->tria->begin_active_line(level)->index(),
			       this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::begin_raw_quad (const unsigned int level) const
{
  Assert (dim==2 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());  
  return raw_quad_iterator (this->tria,
			    dim==2 ? this->tria->begin_raw_quad(level)->level() : 0,
			    this->tria->begin_raw_quad(level)->index(),
			    this);
}


template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::begin_quad (const unsigned int level) const
{
  Assert (dim==2 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());  
  return quad_iterator (this->tria,
			dim==2 ? this->tria->begin_quad(level)->level() : 0,
			this->tria->begin_quad(level)->index(),
			this);
}


template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::begin_active_quad (const unsigned int level) const
{
  Assert (dim==2 || level == 0,
	  typename Triangulation<dim>::ExcFacesHaveNoLevel());  
  return active_quad_iterator (this->tria,
			       dim==2 ? this->tria->begin_active_quad(level)->level() : 0,
			       this->tria->begin_active_quad(level)->index(),
			       this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_hex_iterator
MGDoFHandler<dim>::begin_raw_hex (const unsigned int level) const
{
  return raw_hex_iterator (this->tria,
			   this->tria->begin_raw_hex(level)->level(),
			   this->tria->begin_raw_hex(level)->index(),
			   this);
}


template <int dim>
typename MGDoFHandler<dim>::hex_iterator
MGDoFHandler<dim>::begin_hex (const unsigned int level) const
{
  return hex_iterator (this->tria,
		       this->tria->begin_hex(level)->level(),
		       this->tria->begin_hex(level)->index(),
		       this);
}


template <int dim>
typename MGDoFHandler<dim>::active_hex_iterator
MGDoFHandler<dim>::begin_active_hex (const unsigned int level) const
{
  return active_hex_iterator (this->tria,
			      this->tria->begin_active_hex(level)->level(),
			      this->tria->begin_active_hex(level)->index(),
			      this);
}



				 // We can't ask the respective tria->-function here, as
				 // that would include dereferencing a past-the-end iterator
				 // which is not allowed. Therefore we have to repeat the
				 // code from tria.cc
				 // This is valid for all functions whisch return past the
				 // end iterators, namely all functions end_*()
template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::end_line () const
{
  return raw_line_iterator (this->tria, -1, -1, this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::end_quad () const
{
  return raw_quad_iterator (this->tria, -1, -1, this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_hex_iterator
MGDoFHandler<dim>::end_hex () const
{
  return raw_hex_iterator (this->tria, -1, -1, this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::raw_cell_iterator
MGDoFHandler<dim>::end_raw (const unsigned int level) const
{
  return (level == mg_levels.size()-1 ?
	  end() :
	  begin_raw (level+1));
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::cell_iterator
MGDoFHandler<dim>::end (const unsigned int level) const
{
  return (level == mg_levels.size()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::active_cell_iterator
MGDoFHandler<dim>::end_active (const unsigned int level) const
{
  return (level == mg_levels.size()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::raw_line_iterator
MGDoFHandler<dim>::end_raw_line (const unsigned int level) const
{
  return raw_line_iterator(this->tria,
			   this->tria->end_raw_line(level)->level(),
			   this->tria->end_raw_line(level)->index(),
			   this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::line_iterator
MGDoFHandler<dim>::end_line (const unsigned int level) const
{
  return line_iterator(this->tria,
		       this->tria->end_line(level)->level(),
		       this->tria->end_line(level)->index(),
		       this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::active_line_iterator
MGDoFHandler<dim>::end_active_line (const unsigned int level) const
{
  return active_line_iterator(this->tria,
			      this->tria->end_active_line(level)->level(),
			      this->tria->end_active_line(level)->index(),
			      this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::raw_quad_iterator
MGDoFHandler<dim>::end_raw_quad (const unsigned int level) const
{
  return raw_quad_iterator(this->tria,
			   this->tria->end_raw_quad(level)->level(),
			   this->tria->end_raw_quad(level)->index(),
			   this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::quad_iterator
MGDoFHandler<dim>::end_quad (const unsigned int level) const
{
  return quad_iterator(this->tria,
		       this->tria->end_quad(level)->level(),
		       this->tria->end_quad(level)->index(),
		       this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::active_quad_iterator
MGDoFHandler<dim>::end_active_quad (const unsigned int level) const
{
  return active_quad_iterator(this->tria,
			      this->tria->end_active_quad(level)->level(),
			      this->tria->end_active_quad(level)->index(),
			      this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::raw_hex_iterator
MGDoFHandler<dim>::end_raw_hex (const unsigned int level) const
{
  return raw_hex_iterator(this->tria,
			  this->tria->end_raw_hex(level)->level(),
			  this->tria->end_raw_hex(level)->index(),
			  this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::hex_iterator
MGDoFHandler<dim>::end_hex (const unsigned int level) const
{
  return hex_iterator(this->tria,
		      this->tria->end_hex(level)->level(),
		      this->tria->end_hex(level)->index(),
		      this);
}


template <int dim>
typename internal::MGDoFHandler::Iterators<dim>::active_hex_iterator
MGDoFHandler<dim>::end_active_hex (const unsigned int level) const
{
  return active_hex_iterator(this->tria,
			     this->tria->end_active_hex(level)->level(),
			     this->tria->end_active_hex(level)->index(),
			     this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::last_raw_line (const unsigned int level) const
{
  return raw_line_iterator (this->tria,
			    this->tria->last_raw_line(level)->level(),
			    this->tria->last_raw_line(level)->index(),
			    this);
}


template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::last_line (const unsigned int level) const
{
  return line_iterator (this->tria,
			this->tria->last_line(level)->level(),
			this->tria->last_line(level)->index(),
			this);
}


template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::last_active_line (const unsigned int level) const
{
  return active_line_iterator (this->tria,
			       this->tria->last_active_line(level)->level(),
			       this->tria->last_active_line(level)->index(),
			       this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::last_raw_quad (const unsigned int level) const
{
  return raw_quad_iterator (this->tria,
			    this->tria->last_raw_quad(level)->level(),
			    this->tria->last_raw_quad(level)->index(),
			    this);
}


template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::last_quad (const unsigned int level) const
{
  return quad_iterator (this->tria,
			this->tria->last_quad(level)->level(),
			this->tria->last_quad(level)->index(),
			this);
}


template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::last_active_quad (const unsigned int level) const
{
  return active_quad_iterator (this->tria,
			       this->tria->last_active_quad(level)->level(),
			       this->tria->last_active_quad(level)->index(),
			       this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_hex_iterator
MGDoFHandler<dim>::last_raw_hex (const unsigned int level) const
{
  return raw_hex_iterator (this->tria,
			   this->tria->last_raw_hex(level)->level(),
			   this->tria->last_raw_hex(level)->index(),
			   this);
}


template <int dim>
typename MGDoFHandler<dim>::hex_iterator
MGDoFHandler<dim>::last_hex (const unsigned int level) const
{
  return hex_iterator (this->tria,
		       this->tria->last_hex(level)->level(),
		       this->tria->last_hex(level)->index(),
		       this);
}


template <int dim>
typename MGDoFHandler<dim>::active_hex_iterator
MGDoFHandler<dim>::last_active_hex (const unsigned int level) const
{
  return active_hex_iterator (this->tria,
			      this->tria->last_active_hex(level)->level(),
			      this->tria->last_active_hex(level)->index(),
			      this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::last_raw_line () const
{
  if (dim==1)
        return last_raw_line (mg_levels.size()-1);
  else
    return raw_line_iterator (this->tria,
			      this->tria->last_raw_line()->level(),
			      this->tria->last_raw_line()->index(),
			      this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::last_raw_quad () const
{
  if (dim == 2)
    return last_raw_quad (mg_levels.size()-1);
  else
    return raw_quad_iterator (this->tria,
			      this->tria->last_raw_quad()->level(),
			      this->tria->last_raw_quad()->index(),
			      this);
}


template <int dim>
typename MGDoFHandler<dim>::raw_hex_iterator
MGDoFHandler<dim>::last_raw_hex () const
{
  if (dim==3)
        return last_raw_hex (mg_levels.size()-1);
  else
    return raw_hex_iterator (this->tria,
			     this->tria->last_raw_hex()->level(),
			     this->tria->last_raw_hex()->index(),
			     this);
}


template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::last_line () const
{
  if (dim==1)
            return last_line (mg_levels.size()-1);
  else
    return line_iterator (this->tria,
			  this->tria->last_line()->level(),
			  this->tria->last_line()->index(),
			  this);
}


template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::last_quad () const
{
  if (dim == 2)
    return last_quad (mg_levels.size()-1);
  else
    return quad_iterator (this->tria,
			  this->tria->last_quad()->level(),
			  this->tria->last_quad()->index(),
			  this);
}


template <int dim>
typename MGDoFHandler<dim>::hex_iterator
MGDoFHandler<dim>::last_hex () const
{
  if (dim==3)
        return last_hex (mg_levels.size()-1);
  else
    return hex_iterator (this->tria,
			 this->tria->last_hex()->level(),
			 this->tria->last_hex()->index(),
			 this);
}


template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::last_active_line () const
{
  if (dim==1)
            return last_active_line (mg_levels.size()-1);
  else
    return active_line_iterator (this->tria,
				 this->tria->last_active_line()->level(),
				 this->tria->last_active_line()->index(),
				 this);
}


template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::last_active_quad () const
{
  if (dim==2)
    return last_active_quad (mg_levels.size()-1);
  else
    return active_quad_iterator (this->tria,
				 this->tria->last_active_quad()->level(),
				 this->tria->last_active_quad()->index(),
				 this);
}


template <int dim>
typename MGDoFHandler<dim>::active_hex_iterator
MGDoFHandler<dim>::last_active_hex () const
{
  if (dim == 3)
    return last_active_hex (mg_levels.size()-1);
  else
    return active_hex_iterator (this->tria,
				this->tria->last_active_hex()->level(),
				this->tria->last_active_hex()->index(),
				this);
}


//---------------------------------------------------------------------------


#if deal_II_dimension == 1 
template <>
template <>
unsigned int
MGDoFHandler<1>::get_dof_index<1> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index) const
{
  return this->mg_levels[obj_level]->lines.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
void
MGDoFHandler<1>::set_dof_index<1> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index,
				 const unsigned int       global_index) const
{
  this->mg_levels[obj_level]->lines.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}
#endif

#if deal_II_dimension == 2
template <>
template <>
unsigned int
MGDoFHandler<2>::get_dof_index<1> (const unsigned int       ,
				   const unsigned int       obj_index,
				   const unsigned int       fe_index,
				   const unsigned int       local_index) const
{
  return this->mg_faces->lines.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
void
MGDoFHandler<2>::set_dof_index<1> (const unsigned int       ,
				   const unsigned int       obj_index,
				   const unsigned int       fe_index,
				   const unsigned int       local_index,
				   const unsigned int       global_index) const
{
  this->mg_faces->lines.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}


template <>
template <>
unsigned int
MGDoFHandler<2>::get_dof_index<2> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index) const
{
  return this->mg_levels[obj_level]->quads.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
void
MGDoFHandler<2>::set_dof_index<2> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index,
				 const unsigned int       global_index) const
{
  this->mg_levels[obj_level]->quads.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}
#endif

#if deal_II_dimension == 3
template <>
template <>
unsigned int
MGDoFHandler<3>::get_dof_index<1> (const unsigned int       ,
				   const unsigned int       obj_index,
				   const unsigned int       fe_index,
				   const unsigned int       local_index) const
{
  return this->mg_faces->lines.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
void
MGDoFHandler<3>::set_dof_index<1> (const unsigned int       ,
				   const unsigned int       obj_index,
				   const unsigned int       fe_index,
				   const unsigned int       local_index,
				   const unsigned int       global_index) const
{
  this->mg_faces->lines.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}


template <>
template <>
unsigned int
MGDoFHandler<3>::get_dof_index<2> (const unsigned int       ,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index) const
{
  return this->mg_faces->quads.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
unsigned int
MGDoFHandler<3>::get_dof_index<3> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index) const
{
  return this->mg_levels[obj_level]->hexes.
    get_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index);
}


template <>
template <>
void
MGDoFHandler<3>::set_dof_index<2> (const unsigned int       ,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index,
				 const unsigned int       global_index) const
{
  this->mg_faces->quads.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}


template <>
template <>
void
MGDoFHandler<3>::set_dof_index<3> (const unsigned int       obj_level,
				 const unsigned int       obj_index,
				 const unsigned int       fe_index,
				 const unsigned int       local_index,
				 const unsigned int       global_index) const
{
  this->mg_levels[obj_level]->hexes.
    set_dof_index (*this,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index);
  
}
#endif


template <int dim>
void MGDoFHandler<dim>::distribute_dofs (const FiniteElement<dim> &fe,
					 const unsigned int        offset)
{
				   // first distribute global dofs
  DoFHandler<dim>::distribute_dofs (fe);


				   // reserve space for the MG dof numbers
  reserve_space ();
  mg_used_dofs.resize (this->tria->n_levels(), 0);

				   // Clear user flags because we will
				   // need them. But first we save
				   // them and make sure that we
				   // restore them later such that at
				   // the end of this function the
				   // Triangulation will be in the
				   // same state as it was at the
				   // beginning of this function.
  std::vector<bool> user_flags;
  this->tria->save_user_flags(user_flags);
  const_cast<Triangulation<dim> &>(*(this->tria)).clear_user_flags ();

				   // now distribute indices on each level
				   // separately
  for (unsigned int level=0; level<this->tria->n_levels(); ++level)
    {
      unsigned int next_free_dof = offset;
      cell_iterator cell = begin(level),
		    endc = end(level);

      for (; cell != endc; ++cell) 
	next_free_dof = distribute_dofs_on_cell (cell, next_free_dof);
  
      mg_used_dofs[level] = next_free_dof;
    };
  
				   // finally restore the user flags
  const_cast<Triangulation<dim> &>(*(this->tria)).load_user_flags(user_flags);
}


#if deal_II_dimension == 1

template <>
unsigned int
MGDoFHandler<1>::distribute_dofs_on_cell (cell_iterator &cell,
					  unsigned int   next_free_dof)
{

				   // distribute dofs of vertices
  if (this->selected_fe->dofs_per_vertex > 0)
    for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
      {
	cell_iterator neighbor = cell->neighbor(v);
	
	if (neighbor.state() == IteratorState::valid)
	  {
					     // has neighbor already been processed?
	    if (neighbor->user_flag_set() &&
		(neighbor->level() == cell->level()))
					       // copy dofs if the neighbor is on
					       // the same level (only then are
					       // mg dofs the same)
	      {
		if (v==0) 
		  for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
		    cell->set_mg_vertex_dof_index (cell->level(), 0, d,
						   neighbor->mg_vertex_dof_index (cell->level(), 1, d));
		else
		  for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
		    cell->set_mg_vertex_dof_index (cell->level(), 1, d,
						   neighbor->mg_vertex_dof_index (cell->level(), 0, d));
		
						 // next neighbor
		continue;
	      };
	  };
	
					 // otherwise: create dofs newly
	for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	  cell->set_mg_vertex_dof_index (cell->level(), v, d, next_free_dof++);
      };
  
				   // dofs of line
  if (this->selected_fe->dofs_per_line > 0)
    for (unsigned int d=0; d<this->selected_fe->dofs_per_line; ++d)
      cell->set_mg_dof_index (cell->level(), d, next_free_dof++);

				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
}

#endif


#if deal_II_dimension == 2

template <>
unsigned int
MGDoFHandler<2>::distribute_dofs_on_cell (cell_iterator &cell,
					  unsigned int   next_free_dof) {
  if (this->selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<2>::invalid_dof_index)
	for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);
    
  				   // for the four sides
  if (this->selected_fe->dofs_per_line > 0)
    for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
      {
	line_iterator line = cell->line(side);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->mg_dof_index(cell->level(), 0) == DoFHandler<2>::invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<this->selected_fe->dofs_per_line; ++d)
	    line->set_mg_dof_index (cell->level(), d, next_free_dof++);	    
      };


				   // dofs of quad
  if (this->selected_fe->dofs_per_quad > 0)
    for (unsigned int d=0; d<this->selected_fe->dofs_per_quad; ++d)
      cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
}

#endif


#if deal_II_dimension == 3

template <>
unsigned int
MGDoFHandler<3>::distribute_dofs_on_cell (cell_iterator &cell,
					  unsigned int   next_free_dof) {
  if (this->selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<3>::invalid_dof_index)
	for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);
    
  				   // for the lines
  if (this->selected_fe->dofs_per_line > 0)
    for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
      {
	line_iterator line = cell->line(l);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<this->selected_fe->dofs_per_line; ++d)
	    line->set_mg_dof_index (cell->level(), d, next_free_dof++);	    
      };

  				   // for the quads
  if (this->selected_fe->dofs_per_quad > 0)
    for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
      {
	quad_iterator quad = cell->quad(q);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (quad->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<this->selected_fe->dofs_per_quad; ++d)
	    quad->set_mg_dof_index (cell->level(), d, next_free_dof++);	    
      };


				   // dofs of cell
  if (this->selected_fe->dofs_per_hex > 0)
    for (unsigned int d=0; d<this->selected_fe->dofs_per_hex; ++d)
      cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
}

#endif


template <int dim>
void
MGDoFHandler<dim>::clear ()
{
				   // release own memory
  clear_space ();

				   // let base class release its mem
				   // as well
  DoFHandler<dim>::clear ();  
}


template <int dim>
unsigned int MGDoFHandler<dim>::n_dofs (const unsigned int level) const {
  Assert (level < mg_used_dofs.size(), ExcInvalidLevel(level));
  
  return mg_used_dofs[level];
}


#if deal_II_dimension == 1

template <>
void MGDoFHandler<1>::renumber_dofs (const unsigned int level,
				     const std::vector<unsigned int> &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level), DoFHandler<1>::ExcRenumberingIncomplete());
  
				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // ways it's not only more correct but also
				   // faster
  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	i->set_index (level, d, this->selected_fe->dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->selected_fe->dofs_per_vertex)]);

  for (std::vector<unsigned int>::iterator i=mg_levels[level]->lines.dofs.begin();
       i!=mg_levels[level]->lines.dofs.end(); ++i) 
    {
      if (*i != DoFHandler<1>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}

#endif


#if deal_II_dimension == 2

template <>
void MGDoFHandler<2>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level), 
	  DoFHandler<2>::ExcRenumberingIncomplete());
  
  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	i->set_index (level, d, this->selected_fe->dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->selected_fe->dofs_per_vertex)]);

  if (this->selected_fe->dofs_per_line > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->tria->save_user_flags(user_flags);
      const_cast<Triangulation<2> &>(*(this->tria)).clear_user_flags ();
      
				       // flag all lines adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement
      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell; ++cell)
	for (unsigned int line=0; line < GeometryInfo<2>::faces_per_cell; ++line)
	  cell->face(line)->set_user_flag();
      
      line_iterator line = begin_line(),
		 endline = end_line();
      
      for( ; line != endline; ++line)
	if (line->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->selected_fe->dofs_per_line; ++d)
	      line->set_mg_dof_index (level, d, new_numbers[line->mg_dof_index(level, d)]);
	    line->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<2> &>(*(this->tria)).load_user_flags (user_flags);
    }

  for (std::vector<unsigned int>::iterator i=mg_levels[level]->quads.dofs.begin();
       i!=mg_levels[level]->quads.dofs.end(); ++i)
    {
      if (*i != DoFHandler<2>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}

#endif


#if deal_II_dimension == 3

template <>
void MGDoFHandler<3>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level),
	  DoFHandler<3>::ExcRenumberingIncomplete());
  
  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->selected_fe->dofs_per_vertex; ++d)
	i->set_index (level, d, this->selected_fe->dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->selected_fe->dofs_per_vertex)]);

				   // LINE DoFs
  if (this->selected_fe->dofs_per_line > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->tria->save_user_flags(user_flags);
      const_cast<Triangulation<3> &>(*(this->tria)).clear_user_flags ();
      
				       // flag all lines adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement
      
      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell ; ++cell)
	for (unsigned int line=0; line < GeometryInfo<3>::lines_per_cell; ++line)
	  cell->line(line)->set_user_flag();
      

      line_iterator line = begin_line(),
		 endline = end_line();
      
      for( ; line != endline; ++line)
	if (line->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->selected_fe->dofs_per_line; ++d)
	      line->set_mg_dof_index (level, d, new_numbers[line->mg_dof_index(level, d)]);
	    line->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<3> &>(*(this->tria)).load_user_flags (user_flags);
    }

				   // QUAD DoFs
  if (this->selected_fe->dofs_per_quad > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->tria->save_user_flags(user_flags);
      const_cast<Triangulation<3> &>(*(this->tria)).clear_user_flags ();
      
				       // flag all quads adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement
      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell ; ++cell)
	for (unsigned int quad=0; quad < GeometryInfo<3>::faces_per_cell; ++quad)
	  cell->face(quad)->set_user_flag();
      
      quad_iterator quad = begin_quad(),
		 endline = end_quad();
      
      for( ; quad != endline; ++quad)
	if (quad->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->selected_fe->dofs_per_quad; ++d)
	      quad->set_mg_dof_index (level, d, new_numbers[quad->mg_dof_index(level, d)]);
	    quad->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<3> &>(*(this->tria)).load_user_flags (user_flags);
    }

				   //HEX DoFs
  for (std::vector<unsigned int>::iterator i=mg_levels[level]->hexes.dofs.begin();
       i!=mg_levels[level]->hexes.dofs.end(); ++i)
    {
      if (*i != DoFHandler<3>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}

#endif


#if deal_II_dimension == 1

template <>
void MGDoFHandler<1>::reserve_space () {
  const unsigned int dim = 1;
  
  Assert (this->selected_fe != 0, DoFHandler<dim>::ExcNoFESelected());
  Assert (this->tria->n_levels() > 0, DoFHandler<dim>::ExcInvalidTriangulation());

				   //////////////////////////
				   // DESTRUCTION
  clear_space ();

				   ////////////////////////////
				   // CONSTRUCTION
  
				   // first allocate space for the
				   // lines on each level
  for (unsigned int i=0; i<this->tria->n_levels(); ++i) 
    {
      mg_levels.push_back (new internal::DoFHandler::DoFLevel<1>);

      mg_levels.back()->lines.dofs = std::vector<unsigned int>(this->tria->n_raw_lines(i) *
							 this->selected_fe->dofs_per_line,
							 DoFHandler<1>::invalid_dof_index);
    };

				   // now allocate space for the
				   // vertices. To this end, we need
				   // to construct as many objects as
				   // there are vertices and let them
				   // allocate enough space for their
				   // vertex indices on the levels they
				   // live on. We need therefore to
				   // count to how many levels a cell
				   // belongs to, which we do by looping
				   // over all cells and storing the
				   // maximum and minimum level each
				   // vertex we pass by  belongs to
  mg_vertex_dofs.resize (this->tria->n_vertices());

  std::vector<unsigned int> min_level (this->tria->n_vertices(), this->tria->n_levels());
  std::vector<unsigned int> max_level (this->tria->n_vertices(), 0);

  Triangulation<dim>::cell_iterator cell = this->tria->begin(),
				    endc = this->tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	 ++vertex)
      {
	const unsigned int vertex_index = cell->vertex_index(vertex);
	if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
	  min_level[vertex_index] = cell->level();
	if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
	  max_level[vertex_index] = cell->level();
      };

				   // now allocate the needed space
  for (unsigned int vertex=0; vertex<this->tria->n_vertices(); ++vertex)
    if (this->tria->vertex_used(vertex))
      {
	Assert (min_level[vertex] < this->tria->n_levels(),   ExcInternalError());
	Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());
	
	mg_vertex_dofs[vertex].init (min_level[vertex],
				     max_level[vertex],
				     this->get_fe().dofs_per_vertex);
      }
    else
      {
	Assert (min_level[vertex] == this->tria->n_levels(),   ExcInternalError());
	Assert (max_level[vertex] == 0, ExcInternalError());
	
					 // reset to original state
	mg_vertex_dofs[vertex].init (1, 0, 0);
      }
}

#endif


#if deal_II_dimension == 2

template <>
void MGDoFHandler<2>::reserve_space () {
  const unsigned int dim = 2;
  
  Assert (this->selected_fe != 0, DoFHandler<dim>::ExcNoFESelected());
  Assert (this->tria->n_levels() > 0, DoFHandler<2>::ExcInvalidTriangulation());
  
				   ////////////////////////////
				   // DESTRUCTION
  clear_space ();

				   ////////////////////////////
				   // CONSTRUCTION
  
				   // first allocate space for the
				   // lines and quads on each level
  for (unsigned int i=0; i<this->tria->n_levels(); ++i) 
    {
      mg_levels.push_back (new internal::DoFHandler::DoFLevel<2>);

      mg_levels.back()->quads.dofs = std::vector<unsigned int> (this->tria->n_raw_quads(i) *
							  this->selected_fe->dofs_per_quad,
							  DoFHandler<2>::invalid_dof_index);
    };
  
  mg_faces = new internal::DoFHandler::DoFFaces<2>;
  mg_faces->lines.dofs = std::vector<unsigned int> (this->tria->n_raw_lines() *
						    this->selected_fe->dofs_per_line,
						    DoFHandler<2>::invalid_dof_index);


				   // now allocate space for the
				   // vertices. To this end, we need
				   // to construct as many objects as
				   // there are vertices and let them
				   // allocate enough space for their
				   // vertex indices on the levels they
				   // live on. We need therefore to
				   // count to how many levels a cell
				   // belongs to, which we do by looping
				   // over all cells and storing the
				   // maximum and minimum level each
				   // vertex we pass by  belongs to
  mg_vertex_dofs.resize (this->tria->n_vertices());

				   // initialize these arrays with
				   // invalid values (min>max)
  std::vector<unsigned int> min_level (this->tria->n_vertices(),
				       this->tria->n_levels());
  std::vector<unsigned int> max_level (this->tria->n_vertices(), 0);

  Triangulation<dim>::cell_iterator cell = this->tria->begin(),
				    endc = this->tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	 ++vertex)
      {
	const unsigned int vertex_index = cell->vertex_index(vertex);
	if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
	  min_level[vertex_index] = cell->level();
	if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
	  max_level[vertex_index] = cell->level();
      }


                                   // now allocate the needed space
  for (unsigned int vertex=0; vertex<this->tria->n_vertices(); ++vertex)
    if (this->tria->vertex_used(vertex))
      {
	Assert (min_level[vertex] < this->tria->n_levels(), ExcInternalError());
	Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

	mg_vertex_dofs[vertex].init (min_level[vertex],
				     max_level[vertex],
				     this->get_fe().dofs_per_vertex);
      }
    else
      {
	Assert (min_level[vertex] == this->tria->n_levels(),   ExcInternalError());
	Assert (max_level[vertex] == 0, ExcInternalError());

					 // reset to original state
	mg_vertex_dofs[vertex].init (1, 0, 0);
      }
}

#endif


#if deal_II_dimension == 3

template <>
void MGDoFHandler<3>::reserve_space () {
  const unsigned int dim = 3;
  
  Assert (this->selected_fe != 0, DoFHandler<3>::ExcNoFESelected());
  Assert (this->tria->n_levels() > 0, DoFHandler<3>::ExcInvalidTriangulation());
  
				   ////////////////////////////
				   // DESTRUCTION
  clear_space ();

				   ////////////////////////////
				   // CONSTRUCTION
  
				   // first allocate space for the
				   // lines and quads on each level
  for (unsigned int i=0; i<this->tria->n_levels(); ++i) 
    {
      mg_levels.push_back (new internal::DoFHandler::DoFLevel<3>);

      mg_levels.back()->hexes.dofs
	= std::vector<unsigned int> (this->tria->n_raw_hexs(i) *
				     this->selected_fe->dofs_per_hex,
				     DoFHandler<3>::invalid_dof_index);
    };
  mg_faces = new internal::DoFHandler::DoFFaces<3>;
  mg_faces->lines.dofs = std::vector<unsigned int> (this->tria->n_raw_lines() *
						    this->selected_fe->dofs_per_line,
						    DoFHandler<3>::invalid_dof_index);
  mg_faces->quads.dofs = std::vector<unsigned int> (this->tria->n_raw_quads() *
						    this->selected_fe->dofs_per_quad,
						    DoFHandler<3>::invalid_dof_index);


				   // now allocate space for the
				   // vertices. To this end, we need
				   // to construct as many objects as
				   // there are vertices and let them
				   // allocate enough space for their
				   // vertex indices on the levels they
				   // live on. We need therefore to
				   // count to how many levels a cell
				   // belongs to, which we do by looping
				   // over all cells and storing the
				   // maximum and minimum level each
				   // vertex we pass by  belongs to
  mg_vertex_dofs.resize (this->tria->n_vertices());

  std::vector<unsigned int> min_level (this->tria->n_vertices(), this->tria->n_levels());
  std::vector<unsigned int> max_level (this->tria->n_vertices(), 0);

  Triangulation<dim>::cell_iterator cell = this->tria->begin(),
				    endc = this->tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	 ++vertex)
      {
	const unsigned int vertex_index = cell->vertex_index(vertex);
	if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
	  min_level[vertex_index] = cell->level();
	if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
	  max_level[vertex_index] = cell->level();
      };


				   // now allocate the needed space
  for (unsigned int vertex=0; vertex<this->tria->n_vertices(); ++vertex)
    if (this->tria->vertex_used(vertex))
      {
	Assert (min_level[vertex] < this->tria->n_levels(), ExcInternalError());
	Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());
	
	mg_vertex_dofs[vertex].init (min_level[vertex],
				     max_level[vertex],
				     this->get_fe().dofs_per_vertex);
      }
    else
      {
	Assert (min_level[vertex] == this->tria->n_levels(), ExcInternalError());
	Assert (max_level[vertex] == 0, ExcInternalError());
	
					 // reset to original state
	mg_vertex_dofs[vertex].init (1, 0, 0);
      }
}

#endif


template <int dim>
void MGDoFHandler<dim>::clear_space ()
{
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<mg_levels.size(); ++i)
    delete mg_levels[i];
  mg_levels.clear ();
  delete mg_faces;
  mg_faces=NULL;

				   // also delete vector of vertex
				   // indices
  std::vector<MGVertexDoFs> tmp;
  std::swap (mg_vertex_dofs, tmp);
}



// explicit instantiations
template class MGDoFHandler<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE
