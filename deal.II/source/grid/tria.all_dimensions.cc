//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


/*
 * Single out some functions which are needed by all dimensions, but
 * which are not template. They thus have the same name and when we
 * try to link with the libraries for different dimensions at the same
 * time, we get linker errors for functions defined more than once. By
 * putting these functions in a single file, the linker is allowed to
 * use it only once and throw away all other versions of this file in
 * the other libraries.
 */


#include <base/memory_consumption.h>
#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_boundary.h>

#include <numeric>

DEAL_II_NAMESPACE_OPEN



bool
SubCellData::check_consistency (const unsigned int dim) const
{
  switch (dim)
    {
      case 1:
	    return ((boundary_lines.size() == 0) &&
		    (boundary_quads.size() == 0));
      case 2:
	    return (boundary_quads.size() == 0);
    };
  return true;
}


namespace internal
{
  namespace Triangulation
  {

    NumberCache<1>::NumberCache ()
		    :
		    n_levels (0),
                    n_lines (0),
                    n_active_lines (0)
                                       // all other fields are
                                       // default constructed
    {}



    unsigned int
    NumberCache<1>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (n_levels) +
              MemoryConsumption::memory_consumption (n_lines) +
              MemoryConsumption::memory_consumption (n_lines_level) +
              MemoryConsumption::memory_consumption (n_active_lines) +
              MemoryConsumption::memory_consumption (n_active_lines_level));
    }


    NumberCache<2>::NumberCache () :
                    n_quads (0),
                    n_active_quads (0)
                                       // all other fields are
                                       // default constructed
    {}



    unsigned int
    NumberCache<2>::memory_consumption () const
    {
      return (NumberCache<1>::memory_consumption () +
              MemoryConsumption::memory_consumption (n_quads) +
              MemoryConsumption::memory_consumption (n_quads_level) +
              MemoryConsumption::memory_consumption (n_active_quads) +
              MemoryConsumption::memory_consumption (n_active_quads_level));
    }



    NumberCache<3>::NumberCache () :
                    n_hexes (0),
                    n_active_hexes (0)
                                       // all other fields are
                                       // default constructed
    {}



    unsigned int
    NumberCache<3>::memory_consumption () const
    {
      return (NumberCache<2>::memory_consumption () +
              MemoryConsumption::memory_consumption (n_hexes) +
              MemoryConsumption::memory_consumption (n_hexes_level) +
              MemoryConsumption::memory_consumption (n_active_hexes) +
              MemoryConsumption::memory_consumption (n_active_hexes_level));
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

