//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <dofs/dof_handler.h>
#include <hp/dof_handler.h>

#include<numerics/vectors.templates.h>

DEAL_II_NAMESPACE_OPEN


void
VectorTools::subtract_mean_value(Vector<double>     &v,
				 const std::vector<bool> &p_select)
{
  const unsigned int n = v.size();
  Assert(n == p_select.size(),
         ExcDimensionMismatch(n, p_select.size()));

  double       s       = 0;
  unsigned int counter = 0;
  
  for (unsigned int i=0; i<n; ++i)
    if (p_select[i])
      {
	s += v(i);
	++counter;
      }
  Assert (counter > 0, ExcNoComponentSelected());
  
  s /= counter;
  
  for (unsigned int i=0; i<n; ++i)
    if (p_select[i])
      v(i) -= s;  
}



// ---------------------------- explicit instantiations --------------------
#include "vectors.inst"


DEAL_II_NAMESPACE_CLOSE
