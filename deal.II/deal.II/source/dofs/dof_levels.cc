//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>


unsigned int
DoFLevel<1>::memory_consumption () const
{
  return MemoryConsumption::memory_consumption (line_dofs);
}



unsigned int
DoFLevel<2>::memory_consumption () const
{
  return (DoFLevel<1>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quad_dofs));
}



unsigned int
DoFLevel<3>::memory_consumption () const
{
  return (DoFLevel<2>::memory_consumption () +
	  MemoryConsumption::memory_consumption (hex_dofs));
}
