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
#ifndef __deal2__numbers_h
#define __deal2__numbers_h

#include <base/config.h>
#include <cmath>
#include <limits>


/**
 * Namespace for the declaration of universal constants. Since the
 * availability in <tt>math.h</tt> is not always guaranteed, we put
 * them here.
 *
 * The constants defined here are a subset of the <tt>M_XXX</tt> constants
 * sometimes declared in the system include file <tt>math.h</tt>, but without
 * the prefix <tt>M_</tt>.
 *
 * In addition to that, we declare  <tt>invalid_unsigned_int</tt> to be the
 * largest unsigned integer representable; this value is widely used in
 * the library as a marker for an invalid index, an invalid size of an
 * array, and similar purposes.
 *
 * Most of the members of this namespace are defined in
 * <tt>base/config.h</tt>. Nevertheless, the inline functions are in
 * <tt>base/numbers.h</tt>.
 */
namespace deal_II_numbers 
{
  				             /**
				              * Return @p true if the given
				              * value is a finite floating point
				              * number, i.e. is neither plus or
				              * minus infinity nor NaN (not a
				              * number).
				              */

  inline bool is_finite (const double x) 
  {
#ifdef DEAL_II_HAVE_ISFINITE
    return std::isfinite (x);
#else
				     // check against infinities. not
				     // that if x is a NaN, then both
				     // comparisons will be false
    return ((x >= std::numeric_limits<double>::min())
	    &&
	    (x <= std::numeric_limits<double>::max()));
#endif
  }
}

#endif
