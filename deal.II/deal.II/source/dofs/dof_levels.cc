//----------------------------  dof_levels.cc  ---------------------------
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
//----------------------------  dof_levels.cc  ---------------------------


#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>


unsigned int
DoFLevel<1>::memory_consumption () const
{
  return MemoryConsumption::memory_consumption (line_dofs);
};



unsigned int
DoFLevel<2>::memory_consumption () const
{
  return (DoFLevel<1>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quad_dofs));
};



unsigned int
DoFLevel<3>::memory_consumption () const
{
  return (DoFLevel<2>::memory_consumption () +
	  MemoryConsumption::memory_consumption (hex_dofs));
};
