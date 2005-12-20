//----------------------------  hp_dof_accessor.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_accessor.cc  ------------------------


#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_accessor.templates.h>
#include <dofs/dof_levels.h>
#include <dofs/hp_dof_handler.h>
#include <dofs/hp_dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <vector>


/*------------------------- Functions: DoFCellAccessor -----------------------*/

#if deal_II_dimension == 1

template <>
TriaIterator<1, DoFObjectAccessor<0,1,hpDoFHandler> >
DoFCellAccessor<1,hpDoFHandler>::face (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return TriaIterator<1, DoFObjectAccessor<0,1, hpDoFHandler> >();
}

#endif


#if deal_II_dimension == 2

template <>
TriaIterator<2, DoFObjectAccessor<1,2,hpDoFHandler> >
DoFCellAccessor<2,hpDoFHandler>::face (const unsigned int i) const
{
  return this->line(i);
}

#endif


#if deal_II_dimension == 3

template <>
TriaIterator<3, DoFObjectAccessor<2, 3, hpDoFHandler> >
DoFCellAccessor<3,hpDoFHandler>::face (const unsigned int i) const
{
  return this->quad(i);
}

#endif


/*----------------- Functions: DoFObjectAccessor<1,dim,hpDoFHander> ----------*/

template <>
void DoFObjectAccessor<1, 1, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<1, hpDoFHandler> BaseClass;

  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_line));

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_line_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->line_dofs[offset+i] = index;
}


template <>
void DoFObjectAccessor<1, 2, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<2, hpDoFHandler> BaseClass;

  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_line));

//TODO:[?] In two dimension it could happen that we have different active_fe_indices
// on a line between to cells. Hence we have to differentiate between these two cases.
// Unfortunately, this requires more information then available now.

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_line_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->line_dofs[offset+i] = index;
}


template <>
void DoFObjectAccessor<1, 3, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<3, hpDoFHandler> BaseClass;

  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_line));

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_line_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->line_dofs[offset+i] = index;
}



/*----------------- Functions: DoFObjectAccessor<2,dim,hpDoFHander> ----------*/

template <>
void DoFObjectAccessor<2, 2, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<2, hpDoFHandler> BaseClass;

  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_quad));

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_quad_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->quad_dofs[offset+i] = index;
}


template <>
void DoFObjectAccessor<2, 3, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<3, hpDoFHandler> BaseClass;

  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_quad));

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_quad_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->quad_dofs[offset+i] = index;
}


/*----------------- Functions: DoFObjectAccessor<3,dim,hpDoFHander> ----------*/

template <>
void DoFObjectAccessor<3, 3, hpDoFHandler>::set_dof_index (const unsigned int i,
							   const unsigned int index) const
{
  typedef DoFAccessor<3, hpDoFHandler> BaseClass;
    
  Assert (this->dof_handler != 0,
	  BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->get_fe() != 0,
	  BaseClass::ExcInvalidObject());
  Assert (i<this->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, this->get_fe().dofs_per_hex));

  const unsigned int offset = this->dof_handler->levels[this->present_level]
                              ->dof_hex_index_offset[this->present_index];
  this->dof_handler->levels[this->present_level]
      ->hex_dofs[offset+i] = index;
}
