//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2007 by the deal.II authors
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

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the declaration of universal constants. Since the
 * availability in <tt>math.h</tt> is not always guaranteed, we put
 * them here. Since this file is included by <tt>base/config.h</tt>,
 * they are available to the whole library.
 *
 * The constants defined here are a subset of the <tt>M_XXX</tt> constants
 * sometimes declared in the system include file <tt>math.h</tt>, but without
 * the prefix <tt>M_</tt>.
 *
 * In addition to that, we declare  <tt>invalid_unsigned_int</tt> to be the
 * largest unsigned integer representable; this value is widely used in
 * the library as a marker for an invalid index, an invalid size of an
 * array, and similar purposes.
 */
namespace deal_II_numbers {
                                             /**
                                              * Representation of the
                                              * largest number that
                                              * can be put into an
                                              * unsigned integer. This
                                              * value is widely used
                                              * throughout the library
                                              * as a marker for an
                                              * invalid unsigned
                                              * integer value, such as
                                              * an invalid array
                                              * index, an invalid
                                              * array size, and the
                                              * like.
                                              */
  static const unsigned int 
    invalid_unsigned_int = static_cast<unsigned int> (-1);

                                             /**
                                              * e
                                              */
  static const double  E       = 2.7182818284590452354;

                                             /**
                                              * log_2 e
                                              */
  static const double  LOG2E   = 1.4426950408889634074;

                                             /**
                                              * log_10 e
                                              */
  static const double  LOG10E  = 0.43429448190325182765;

                                             /**
                                              * log_e 2
                                              */
  static const double  LN2     = 0.69314718055994530942;

                                             /**
                                              * log_e 10
                                              */
  static const double  LN10    = 2.30258509299404568402;

                                             /**
                                              * pi
                                              */
  static const double  PI      = 3.14159265358979323846;

                                             /**
                                              * pi/2
                                              */
  static const double  PI_2    = 1.57079632679489661923;

                                             /**
                                              * pi/4
                                              */
  static const double  PI_4    = 0.78539816339744830962;

                                             /**
                                              * sqrt(2)
                                              */
  static const double  SQRT2   = 1.41421356237309504880;

                                             /**
                                              * 1/sqrt(2)
                                              */
  static const double  SQRT1_2 = 0.70710678118654752440;

				             /**
				              * Return @p true if the given
				              * value is a finite floating
				              * point number, i.e. is neither
				              * plus or minus infinity nor NaN
				              * (not a number).
					      *
					      * Note that the argument type of
					      * this function is
					      * <code>double</code>. In other
					      * words, if you give a very
					      * large number of type
					      * <code>long double</code>, this
					      * function may return
					      * <code>false</code> even if the
					      * number is finite with respect
					      * to type <code>long
					      * double</code>.
				              */
  bool is_finite (const double x);
}

DEAL_II_NAMESPACE_CLOSE

#endif
