//----------------------------  dof_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_handler.cc  ---------------------------


#include <base/memory_consumption.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_levels.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_levels.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>

#include <set>
#include <algorithm>


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif




template <int dim>
const unsigned int DoFHandler<dim>::invalid_dof_index;



template <int dim>
DoFHandler<dim>::DoFHandler (Triangulation<dim> &tria) :
		tria(&tria),
		used_dofs (0)
{};


template <int dim>
DoFHandler<dim>::~DoFHandler ()
{
				   // release allocated memory
  clear ();
};


#if deal_II_dimension == 1

template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};


template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::begin (const unsigned int level) const {
  return begin_line (level);
};


template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};


template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::end () const {
  return end_line ();
};


template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw () const {
  return last_raw_line ();
};


template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};


template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::last () const {
  return last_line ();
};


template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::last (const unsigned int level) const {
  return last_line (level);
};


template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active () const {
  return last_active_line ();
};


template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};


template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::begin_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::begin_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::begin_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::end_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::last_raw_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::last_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::last_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::last_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::last_active_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::last_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::begin_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::begin_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::begin_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::end_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::last_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::last_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::last_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::last_raw_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::last_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::last_active_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_hex_iterator
DoFHandler<1>::begin_raw_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::hex_iterator
DoFHandler<1>::begin_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_hex_iterator
DoFHandler<1>::begin_active_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_hex_iterator
DoFHandler<1>::end_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_hex_iterator
DoFHandler<1>::last_raw_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::raw_hex_iterator
DoFHandler<1>::last_raw_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::hex_iterator
DoFHandler<1>::last_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::hex_iterator
DoFHandler<1>::last_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_hex_iterator
DoFHandler<1>::last_active_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::active_hex_iterator
DoFHandler<1>::last_active_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};

#endif


#if deal_II_dimension == 2

template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};


template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};


template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};


template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::end () const {
  return end_quad ();
};


template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw () const {
  return last_raw_quad ();
};


template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};


template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::last () const {
  return last_quad ();
};


template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::last (const unsigned int level) const {
  return last_quad (level);
};


template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active () const {
  return last_active_quad ();
};


template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};


template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::begin_raw_face (const unsigned int level) const {
  return begin_raw_line (level);
};


template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::begin_face (const unsigned int level) const {
  return begin_line (level);
};


template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::begin_active_face (const unsigned int level) const {
  return begin_active_line (level);
};


template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::end_face () const {
  return end_line ();
};


template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::last_raw_face () const {
  return last_raw_line ();
};


template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::last_raw_face (const unsigned int level) const {
  return last_raw_line (level);
};


template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::last_face () const {
  return last_line ();
};


template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::last_face (const unsigned int level) const {
  return last_line (level);
};


template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::last_active_face () const {
  return last_active_line ();
};


template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::last_active_face (const unsigned int level) const {
  return last_active_line (level);
};


template <>
DoFHandler<2>::raw_hex_iterator
DoFHandler<2>::begin_raw_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::hex_iterator
DoFHandler<2>::begin_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::active_hex_iterator
DoFHandler<2>::begin_active_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::raw_hex_iterator
DoFHandler<2>::end_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::raw_hex_iterator
DoFHandler<2>::last_raw_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::raw_hex_iterator
DoFHandler<2>::last_raw_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::hex_iterator
DoFHandler<2>::last_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::hex_iterator
DoFHandler<2>::last_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::active_hex_iterator
DoFHandler<2>::last_active_hex (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<2>::active_hex_iterator
DoFHandler<2>::last_active_hex () const {
  Assert (false, ExcNotImplemented());
  return 0;
};


#endif


#if deal_II_dimension == 3

template <>
DoFHandler<3>::raw_cell_iterator
DoFHandler<3>::begin_raw (const unsigned int level) const {
  return begin_raw_hex (level);
};


template <>
DoFHandler<3>::cell_iterator
DoFHandler<3>::begin (const unsigned int level) const {
  return begin_hex (level);
};


template <>
DoFHandler<3>::active_cell_iterator
DoFHandler<3>::begin_active (const unsigned int level) const {
  return begin_active_hex (level);
};


template <>
DoFHandler<3>::raw_cell_iterator
DoFHandler<3>::end () const {
  return end_hex ();
};


template <>
DoFHandler<3>::raw_cell_iterator
DoFHandler<3>::last_raw () const {
  return last_raw_hex ();
};


template <>
DoFHandler<3>::raw_cell_iterator
DoFHandler<3>::last_raw (const unsigned int level) const {
  return last_raw_hex (level);
};


template <>
DoFHandler<3>::cell_iterator
DoFHandler<3>::last () const {
  return last_hex ();
};


template <>
DoFHandler<3>::cell_iterator
DoFHandler<3>::last (const unsigned int level) const {
  return last_hex (level);
};


template <>
DoFHandler<3>::active_cell_iterator
DoFHandler<3>::last_active () const {
  return last_active_hex ();
};


template <>
DoFHandler<3>::active_cell_iterator
DoFHandler<3>::last_active (const unsigned int level) const {
  return last_active_hex (level);
};


template <>
DoFHandler<3>::raw_face_iterator
DoFHandler<3>::begin_raw_face (const unsigned int level) const {
  return begin_raw_quad (level);
};


template <>
DoFHandler<3>::face_iterator
DoFHandler<3>::begin_face (const unsigned int level) const {
  return begin_quad (level);
};


template <>
DoFHandler<3>::active_face_iterator
DoFHandler<3>::begin_active_face (const unsigned int level) const {
  return begin_active_quad (level);
};


template <>
DoFHandler<3>::raw_face_iterator
DoFHandler<3>::end_face () const {
  return end_quad ();
};


template <>
DoFHandler<3>::raw_face_iterator
DoFHandler<3>::last_raw_face () const {
  return last_raw_quad ();
};


template <>
DoFHandler<3>::raw_face_iterator
DoFHandler<3>::last_raw_face (const unsigned int level) const {
  return last_raw_quad (level);
};


template <>
DoFHandler<3>::face_iterator
DoFHandler<3>::last_face () const {
  return last_quad ();
};


template <>
DoFHandler<3>::face_iterator
DoFHandler<3>::last_face (const unsigned int level) const {
  return last_quad (level);
};


template <>
DoFHandler<3>::active_face_iterator
DoFHandler<3>::last_active_face () const {
  return last_active_quad ();
};


template <>
DoFHandler<3>::active_face_iterator
DoFHandler<3>::last_active_face (const unsigned int level) const {
  return last_active_quad (level);
};


#endif


template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::begin_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->begin_raw_line(level)->level(),
			    tria->begin_raw_line(level)->index(),
			    this);
};


template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::begin_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->begin_line(level)->level(),
			tria->begin_line(level)->index(),
			this);
};


template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::begin_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->begin_active_line(level)->level(),
			       tria->begin_active_line(level)->index(),
			       this);
};


template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::begin_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->begin_raw_quad(level)->level(),
			    tria->begin_raw_quad(level)->index(),
			    this);
};


template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::begin_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->begin_quad(level)->level(),
			tria->begin_quad(level)->index(),
			this);
};


template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::begin_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->begin_active_quad(level)->level(),
			       tria->begin_active_quad(level)->index(),
			       this);
};


template <int dim>
typename DoFHandler<dim>::raw_hex_iterator
DoFHandler<dim>::begin_raw_hex (const unsigned int level) const {
  return raw_hex_iterator (tria,
			   tria->begin_raw_hex(level)->level(),
			   tria->begin_raw_hex(level)->index(),
			   this);
};


template <int dim>
typename DoFHandler<dim>::hex_iterator
DoFHandler<dim>::begin_hex (const unsigned int level) const {
  return hex_iterator (tria,
		       tria->begin_hex(level)->level(),
		       tria->begin_hex(level)->index(),
		       this);
};


template <int dim>
typename DoFHandler<dim>::active_hex_iterator
DoFHandler<dim>::begin_active_hex (const unsigned int level) const {
  return active_hex_iterator (tria,
			      tria->begin_active_hex(level)->level(),
			      tria->begin_active_hex(level)->index(),
			      this);
};


template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::end_line () const {
  return raw_line_iterator (tria, -1, -1, this);
};


template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::end_quad () const {
  return raw_quad_iterator (tria, -1, -1, this);
};


template <int dim>
typename DoFHandler<dim>::raw_hex_iterator
DoFHandler<dim>::end_hex () const {
  return raw_hex_iterator (tria, -1, -1, this);
};


template <int dim>
typename DoFDimensionInfo<dim>::raw_cell_iterator
DoFHandler<dim>::end_raw (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end() :
	  begin_raw (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::cell_iterator
DoFHandler<dim>::end (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::active_cell_iterator
DoFHandler<dim>::end_active (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::raw_face_iterator
DoFHandler<dim>::end_raw_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_face() :
	  begin_raw_face (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::face_iterator
DoFHandler<dim>::end_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  face_iterator(end_face()) :
	  begin_face (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::active_face_iterator
DoFHandler<dim>::end_active_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_face_iterator(end_face()) :
	  begin_active_face (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::raw_line_iterator
DoFHandler<dim>::end_raw_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_line() :
	  begin_raw_line (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::line_iterator
DoFHandler<dim>::end_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  line_iterator(end_line()) :
	  begin_line (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::active_line_iterator
DoFHandler<dim>::end_active_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_line_iterator(end_line()) :
	  begin_active_line (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::raw_quad_iterator
DoFHandler<dim>::end_raw_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_quad() :
	  begin_raw_quad (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::quad_iterator
DoFHandler<dim>::end_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  quad_iterator(end_quad()) :
	  begin_quad (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::active_quad_iterator
DoFHandler<dim>::end_active_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_quad_iterator(end_quad()) :
	  begin_active_quad (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::raw_hex_iterator
DoFHandler<dim>::end_raw_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_hex() :
	  begin_raw_hex (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::hex_iterator
DoFHandler<dim>::end_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  hex_iterator(end_hex()) :
	  begin_hex (level+1));
};


template <int dim>
typename DoFDimensionInfo<dim>::active_hex_iterator
DoFHandler<dim>::end_active_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_hex_iterator(end_hex()) :
	  begin_active_hex (level+1));
};


template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::last_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->last_raw_line(level)->level(),
			    tria->last_raw_line(level)->index(),
			    this);
};


template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::last_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->last_line(level)->level(),
			tria->last_line(level)->index(),
			this);
};


template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::last_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->last_active_line(level)->level(),
			       tria->last_active_line(level)->index(),
			       this);
};


template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::last_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->last_raw_quad(level)->level(),
			    tria->last_raw_quad(level)->index(),
			    this);
};


template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::last_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->last_quad(level)->level(),
			tria->last_quad(level)->index(),
			this);
};


template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::last_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->last_active_quad(level)->level(),
			       tria->last_active_quad(level)->index(),
			       this);
};


template <int dim>
typename DoFHandler<dim>::raw_hex_iterator
DoFHandler<dim>::last_raw_hex (const unsigned int level) const {
  return raw_hex_iterator (tria,
			   tria->last_raw_hex(level)->level(),
			   tria->last_raw_hex(level)->index(),
			   this);
};


template <int dim>
typename DoFHandler<dim>::hex_iterator
DoFHandler<dim>::last_hex (const unsigned int level) const {
  return hex_iterator (tria,
		       tria->last_hex(level)->level(),
		       tria->last_hex(level)->index(),
		       this);
};


template <int dim>
typename DoFHandler<dim>::active_hex_iterator
DoFHandler<dim>::last_active_hex (const unsigned int level) const {
  return active_hex_iterator (tria,
			      tria->last_active_hex(level)->level(),
			      tria->last_active_hex(level)->index(),
			      this);
};


template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::last_raw_quad () const {
  return last_raw_quad (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::raw_hex_iterator
DoFHandler<dim>::last_raw_hex () const {
  return last_raw_hex (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::last_line () const {
  return last_line (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::last_quad () const {
  return last_quad (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::hex_iterator
DoFHandler<dim>::last_hex () const {
  return last_hex (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::last_active_line () const {
  return last_active_line (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::last_active_quad () const {
  return last_active_quad (levels.size()-1);
};


template <int dim>
typename DoFHandler<dim>::active_hex_iterator
DoFHandler<dim>::last_active_hex () const {
  return last_active_hex (levels.size()-1);
};


//------------------------------------------------------------------


#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::n_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return 2*selected_fe->dofs_per_vertex;
};



template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (FunctionMap::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((i->first == 0) || (i->first == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
};



template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (std::set<unsigned char>::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((*i == 0) || (*i == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
};

#endif


template <int dim>
unsigned int DoFHandler<dim>::n_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  
  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);

				   // loop over all faces to check
				   // whether they are at a
				   // boundary. note that we need not
				   // take special care of single
				   // lines (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every boundary line is
				   // also part of a boundary face.
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      };
  return boundary_dofs.size();
};    



template <int dim>
unsigned int
DoFHandler<dim>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());
  
  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      };
  return boundary_dofs.size();
};    



template <int dim>
unsigned int
DoFHandler<dim>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());
  
  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (find (boundary_indicators.begin(),
	      boundary_indicators.end(),
	      face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      };
  return boundary_dofs.size();
};    



template <int dim>
const Triangulation<dim> & DoFHandler<dim>::get_tria () const
{
  return *tria;
};



template <int dim>
unsigned int
DoFHandler<dim>::memory_consumption () const
{
  unsigned int mem = (MemoryConsumption::memory_consumption (tria) +
		      MemoryConsumption::memory_consumption (selected_fe) +
		      MemoryConsumption::memory_consumption (tria) +
		      MemoryConsumption::memory_consumption (levels) +
		      MemoryConsumption::memory_consumption (used_dofs) +
		      MemoryConsumption::memory_consumption (vertex_dofs));
  for (unsigned int i=0; i<levels.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*levels[i]);
  
  return mem;
};



template <int dim>
void DoFHandler<dim>::distribute_dofs (const FiniteElement<dim> &ff,
				       const unsigned int        offset)
{
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
  selected_fe = &ff;
  
  reserve_space ();

				   // clear user flags because we will
				   // need them
  tria->clear_user_flags ();
  
  unsigned int next_free_dof = offset;   
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell) 
    next_free_dof = distribute_dofs_on_cell (cell, next_free_dof);
  
  used_dofs = next_free_dof;
};


template <int dim>
void DoFHandler<dim>::clear () {
				   // release lock to old fe
  selected_fe = 0;

				   // release memory
  clear_space ();
};


#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator &cell,
						     unsigned int          next_free_dof) {

				   // distribute dofs of vertices
  for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
    {
      cell_iterator neighbor = cell->neighbor(v);

      if (neighbor.state() == IteratorState::valid)
	{
					   // find true neighbor; may be its
					   // a child of @p{neighbor}
	  while (neighbor->has_children())
	    neighbor = neighbor->child(v==0 ? 1 : 0);

					   // has neighbor already been processed?
	  if (neighbor->user_flag_set())
					     // copy dofs
	    {
	      if (v==0) 
		for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (0, d,
					      neighbor->vertex_dof_index (1, d));
	      else
		for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (1, d,
					      neighbor->vertex_dof_index (0, d));

					       // next neighbor
	      continue;
	    };
	};
            
				       // otherwise: create dofs newly
      for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	cell->set_vertex_dof_index (v, d, next_free_dof++);
    };
  
				   // dofs of line
  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
    cell->set_dof_index (d, next_free_dof++);

				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif


#if deal_II_dimension == 2

template <>
unsigned int DoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator &cell,
						     unsigned int          next_free_dof) {
  if (selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->vertex_dof_index(vertex, 0) == invalid_dof_index)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  cell->set_vertex_dof_index (vertex, d, next_free_dof++);
    
  				   // for the four sides
  if (selected_fe->dofs_per_line > 0)
    for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
      {
	line_iterator line = cell->line(side);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->dof_index(0) == invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	    line->set_dof_index (d, next_free_dof++);	    
      };


				   // dofs of quad
  if (selected_fe->dofs_per_quad > 0)
    for (unsigned int d=0; d<selected_fe->dofs_per_quad; ++d)
      cell->set_dof_index (d, next_free_dof++);


				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif


#if deal_II_dimension == 3

template <>
unsigned int DoFHandler<3>::distribute_dofs_on_cell (active_cell_iterator &cell,
						     unsigned int          next_free_dof) {
  if (selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->vertex_dof_index(vertex, 0) == invalid_dof_index)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  cell->set_vertex_dof_index (vertex, d, next_free_dof++);
    
  				   // for the lines
  if (selected_fe->dofs_per_line > 0)
    for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
      {
	line_iterator line = cell->line(l);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->dof_index(0) == invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	    line->set_dof_index (d, next_free_dof++);	    
      };

    				   // for the quads
  if (selected_fe->dofs_per_quad > 0)
    for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
      {
	quad_iterator quad = cell->quad(q);
	
					 // distribute dofs if necessary:
					 // check whether quad dof is already
					 // numbered (check only first dof)
	if (quad->dof_index(0) == invalid_dof_index)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<selected_fe->dofs_per_quad; ++d)
	    quad->set_dof_index (d, next_free_dof++);	    
      };


				   // dofs of hex
  if (selected_fe->dofs_per_hex > 0)
    for (unsigned int d=0; d<selected_fe->dofs_per_hex; ++d)
      cell->set_dof_index (d, next_free_dof++);


				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif


#if deal_II_dimension == 1

template <>
void DoFHandler<1>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
{
  Assert (new_numbers.size() == n_dofs(), ExcRenumberingIncomplete());
#ifdef DEBUG
				   // assert that the new indices are
				   // consecutively numbered
  if (true)
    {
      std::vector<unsigned int> tmp(new_numbers);
      std::sort (tmp.begin(), tmp.end());
      std::vector<unsigned int>::const_iterator p = tmp.begin();
      unsigned int                         i = 0;
      for (; p!=tmp.end(); ++p, ++i)
	Assert (*p == i, ExcNewNumbersNotConsecutive(i));
    };
#endif

				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // way it's not only more correct but also
				   // faster; note, however, that dof numbers
				   // may be invalid_dof_index, namely when the appropriate
				   // vertex/line/etc is unused
  for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
    if (*i != invalid_dof_index)
      *i = new_numbers[*i];
    else
				       // if index is invalid_dof_index: check if this one
				       // really is unused
      Assert (tria->vertices_used[(i-vertex_dofs.begin()) /
				 selected_fe->dofs_per_vertex] == false,
	      ExcInternalError ());
  
  for (unsigned int level=0; level<levels.size(); ++level) 
    for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
	 i!=levels[level]->line_dofs.end(); ++i)
      if (*i != invalid_dof_index)
	*i = new_numbers[*i];
};

#endif


#if deal_II_dimension == 2

template <>
void DoFHandler<2>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
{
  Assert (new_numbers.size() == n_dofs(), ExcRenumberingIncomplete());
#ifdef DEBUG
				   // assert that the new indices are
				   // consecutively numbered
  if (true)
    {
      std::vector<unsigned int> tmp(new_numbers);
      std::sort (tmp.begin(), tmp.end());
      std::vector<unsigned int>::const_iterator p = tmp.begin();
      unsigned int                         i = 0;
      for (; p!=tmp.end(); ++p, ++i)
	Assert (*p == i, ExcNewNumbersNotConsecutive(i));
    };
#endif

				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // way it's not only more correct but also
				   // faster; note, however, that dof numbers
				   // may be invalid_dof_index, namely when the appropriate
				   // vertex/line/etc is unused
  for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
    if (*i != invalid_dof_index)
      *i = new_numbers[*i];
    else
				       // if index is invalid_dof_index: check if this one
				       // really is unused
      Assert (tria->vertices_used[(i-vertex_dofs.begin()) /
				 selected_fe->dofs_per_vertex] == false,
	      ExcInternalError ());
  
  for (unsigned int level=0; level<levels.size(); ++level) 
    {
      for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
	   i!=levels[level]->line_dofs.end(); ++i)
	if (*i != invalid_dof_index)
	  *i = new_numbers[*i];
      for (std::vector<unsigned int>::iterator i=levels[level]->quad_dofs.begin();
	   i!=levels[level]->quad_dofs.end(); ++i)
	if (*i != invalid_dof_index)
	  *i = new_numbers[*i];
    };
};

#endif


#if deal_II_dimension == 3

template <>
void DoFHandler<3>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
{
  Assert (new_numbers.size() == n_dofs(), ExcRenumberingIncomplete());
#ifdef DEBUG
				   // assert that the new indices are
				   // consecutively numbered
  if (true)
    {
      std::vector<unsigned int> tmp(new_numbers);
      std::sort (tmp.begin(), tmp.end());
      std::vector<unsigned int>::const_iterator p = tmp.begin();
      unsigned int                              i = 0;
      for (; p!=tmp.end(); ++p, ++i)
	Assert (*p == i, ExcNewNumbersNotConsecutive(i));
    };
#endif

				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // way it's not only more correct but also
				   // faster; note, however, that dof numbers
				   // may be invalid_dof_index, namely when the appropriate
				   // vertex/line/etc is unused
  for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
    if (*i != invalid_dof_index)
      *i = new_numbers[*i];
    else
				       // if index is invalid_dof_index: check if this one
				       // really is unused
      Assert (tria->vertices_used[(i-vertex_dofs.begin()) /
				 selected_fe->dofs_per_vertex] == false,
	      ExcInternalError ());
  
  for (unsigned int level=0; level<levels.size(); ++level) 
    {
      for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
	   i!=levels[level]->line_dofs.end(); ++i)
	if (*i != invalid_dof_index)
	  *i = new_numbers[*i];
      for (std::vector<unsigned int>::iterator i=levels[level]->quad_dofs.begin();
	   i!=levels[level]->quad_dofs.end(); ++i)
	if (*i != invalid_dof_index)
	  *i = new_numbers[*i];
      for (std::vector<unsigned int>::iterator i=levels[level]->hex_dofs.begin();
	   i!=levels[level]->hex_dofs.end(); ++i)
	if (*i != invalid_dof_index)
	  *i = new_numbers[*i];
    };
};

#endif


#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::max_couplings_between_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return std::min(3*selected_fe->dofs_per_vertex +
		  2*selected_fe->dofs_per_line, n_dofs());
};



template <>
unsigned int DoFHandler<1>::max_couplings_between_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return selected_fe->dofs_per_vertex;
};

#endif


#if deal_II_dimension == 2

template <>
unsigned int DoFHandler<2>::max_couplings_between_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // get these numbers by drawing pictures
				   // and counting...
				   // example:
				   //   |     |     |
				   // --x-----x--x--X--
				   //   |     |  |  |
				   //   |     x--x--x
				   //   |     |  |  |
				   // --x--x--*--x--x--
				   //   |  |  |     |
				   //   x--x--x     |
				   //   |  |  |     |
				   // --X--x--x-----x--
				   //   |     |     |
				   // x = vertices connected with center vertex *;
				   //   = total of 19
				   // (the X vertices are connected with * if
				   // the vertices adjacent to X are hanging
				   // nodes)
				   // count lines -> 28 (don't forget to count
				   // mother and children separately!)
  unsigned int max_couplings;
  switch (tria->max_adjacent_cells())
    {
      case 4:
	    max_couplings=19*selected_fe->dofs_per_vertex +
			  28*selected_fe->dofs_per_line +
			  8*selected_fe->dofs_per_quad;
	    break;
      case 5:
	    max_couplings=21*selected_fe->dofs_per_vertex +
			  31*selected_fe->dofs_per_line +
			  9*selected_fe->dofs_per_quad;
	    break;
      case 6:
	    max_couplings=28*selected_fe->dofs_per_vertex +
			  42*selected_fe->dofs_per_line +
			  12*selected_fe->dofs_per_quad;
	    break;
      case 7:
	    max_couplings=30*selected_fe->dofs_per_vertex +
			  45*selected_fe->dofs_per_line +
			  13*selected_fe->dofs_per_quad;
	    break;
      case 8:
	    max_couplings=37*selected_fe->dofs_per_vertex +
			  56*selected_fe->dofs_per_line +
			  16*selected_fe->dofs_per_quad;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
	    max_couplings=0;
    };
  return std::min(max_couplings,n_dofs());
};



template <>
unsigned int DoFHandler<2>::max_couplings_between_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return 3*selected_fe->dofs_per_vertex + 2*selected_fe->dofs_per_line;
};

#endif


#if deal_II_dimension == 3

template <>
unsigned int DoFHandler<3>::max_couplings_between_dofs () const
{
//TODO:[?] Invent significantly better estimates than the ones in this function  
  Assert (selected_fe != 0, ExcNoFESelected());

				   // doing the same thing here is a rather
				   // complicated thing, compared to the 2d
				   // case, since it is hard to draw pictures
				   // with several refined hexahedra :-) so I
				   // presently only give a coarse estimate
				   // for the case that at most 8 hexes meet
				   // at each vertex
				   //
				   // can anyone give better estimate here?
  const unsigned int max_adjacent_cells = tria->max_adjacent_cells();

  unsigned int max_couplings;
  if (max_adjacent_cells <= 8)
    max_couplings=7*7*7*selected_fe->dofs_per_vertex +
		  7*6*7*3*selected_fe->dofs_per_line +
		  9*4*7*3*selected_fe->dofs_per_quad +
		  27*selected_fe->dofs_per_hex;
  else
    {
      Assert (false, ExcNotImplemented());
      max_couplings=0;
    }
  
  return std::min(max_couplings,n_dofs());
};


template <>
unsigned int DoFHandler<3>::max_couplings_between_boundary_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
				   // we need to take refinement of
				   // one boundary face into consideration
				   // here; in fact, this function returns
				   // what #max_coupling_between_dofs<2>
				   // returns
				   //
				   // we assume here, that only four faces
				   // meet at the boundary; this assumption
				   // is not justified and needs to be
				   // fixed some time. fortunately, ommitting
				   // it for now does no harm since the
				   // matrix will cry foul if its requirements
				   // are not satisfied
  return (19*selected_fe->dofs_per_vertex +
	  28*selected_fe->dofs_per_line +
	  8*selected_fe->dofs_per_quad);
};


#endif



#if deal_II_dimension == 1

template <>
void DoFHandler<1>::reserve_space () {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());

                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  clear_space ();
  
  vertex_dofs.resize(tria->vertices.size()*selected_fe->dofs_per_vertex,
		     invalid_dof_index);
    
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<1>);

      levels.back()->line_dofs = std::vector<unsigned int>(tria->levels[i]->lines.lines.size() *
							   selected_fe->dofs_per_line,
							   invalid_dof_index);
    };
};


#endif


#if deal_II_dimension == 2

template <>
void DoFHandler<2>::reserve_space () {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  clear_space ();
  
  vertex_dofs.resize(tria->vertices.size()*selected_fe->dofs_per_vertex,
		     invalid_dof_index);

  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<2>);

      levels.back()->line_dofs = std::vector<unsigned int> (tria->levels[i]->lines.lines.size() *
							    selected_fe->dofs_per_line,
							    invalid_dof_index);
      levels.back()->quad_dofs = std::vector<unsigned int> (tria->levels[i]->quads.quads.size() *
							    selected_fe->dofs_per_quad,
							    invalid_dof_index);
    };
};

#endif


#if deal_II_dimension == 3

template <>
void DoFHandler<3>::reserve_space () {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  clear_space ();
  
  vertex_dofs.resize(tria->vertices.size()*selected_fe->dofs_per_vertex,
		     invalid_dof_index);

  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<3>);

      levels.back()->line_dofs = std::vector<unsigned int> (tria->levels[i]->lines.lines.size() *
							    selected_fe->dofs_per_line,
							    invalid_dof_index);
      levels.back()->quad_dofs = std::vector<unsigned int> (tria->levels[i]->quads.quads.size() *
							    selected_fe->dofs_per_quad,
							    invalid_dof_index);
      levels.back()->hex_dofs = std::vector<unsigned int> (tria->levels[i]->hexes.hexes.size() *
							   selected_fe->dofs_per_hex,
							   invalid_dof_index);
    };
};

#endif


template <int dim>
void DoFHandler<dim>::clear_space () {  
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.resize (0);

  std::vector<unsigned int> tmp;
  std::swap (vertex_dofs, tmp);
};


/*-------------- Explicit Instantiations -------------------------------*/
template class DoFHandler<deal_II_dimension>;
