//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <dofs/number_cache.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    NumberCache::NumberCache ()
		    :
		    n_global_dofs (0),
		    n_locally_owned_dofs (0)
    {}



    unsigned int
    NumberCache::memory_consumption () const
    {
      return
	MemoryConsumption::memory_consumption (n_global_dofs) +
	MemoryConsumption::memory_consumption (n_locally_owned_dofs) +
	MemoryConsumption::memory_consumption (locally_owned_dofs) +
	MemoryConsumption::memory_consumption (n_locally_owned_dofs_per_processor) +
	MemoryConsumption::memory_consumption (locally_owned_dofs_per_processor);
    }
  }
}

DEAL_II_NAMESPACE_CLOSE
