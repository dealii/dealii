//----------------------------  tests.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tests.h  ---------------------------

#ifndef __tests_tests_h
#define __tests_tests_h

// common definitions used in all the tests


#include <base/logstream.h>
#include <cmath>

// overload floating point output operators for LogStream so that small
// numbers below a certain threshold are simply printed as zeros. this removes
// a number of possible places where output may differ depending on platform,
// compiler options, etc, simply because round-off is different.
inline
LogStream & operator << (LogStream &logstream,
                         const double d)
{
  if (std::fabs (d) < 1e-10)
    logstream.operator << <double> (0.);
  else
    logstream.operator << <double> (d);
  return logstream;
}


inline
LogStream & operator << (LogStream &logstream,
                         const float d)
{
  if (std::fabs (d) < 1e-8)
    logstream.operator << <float> (0.);
  else
    logstream.operator << <float> (d);
  return logstream;
}


#endif // __tests_tests_h
