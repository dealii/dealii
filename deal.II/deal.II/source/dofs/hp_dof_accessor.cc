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

