//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/config.h>
#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace deal_II_numbers 
{
  bool is_finite (const double x) 
  {
#ifdef DEAL_II_HAVE_ISFINITE
    return std::isfinite (x);
#else
				     // check against infinities. not
				     // that if x is a NaN, then both
				     // comparisons will be false
    return ((x >= -std::numeric_limits<double>::max())
	    &&
DEAL_II_NAMESPACE_CLOSE
	    (x <= std::numeric_limits<double>::max()));
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE
