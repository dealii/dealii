//----------------------------  tests.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
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

#include <base/config.h>
#include <base/logstream.h>
#include <base/exceptions.h>
#include <cmath>

// overload floating point output operators for LogStream so that small
// numbers below a certain threshold are simply printed as zeros. this removes
// a number of possible places where output may differ depending on platform,
// compiler options, etc, simply because round-off is different.
inline
LogStream & operator << (LogStream &logstream,
                         const float d)
{
  if (std::fabs (d) < 1e-8)
    logstream.
#ifdef DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG
      template
#endif
      operator << <float> (0.);
  else
    logstream.
#ifdef DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG
      template
#endif
      operator << <float> (d);
  return logstream;
}


struct SwitchOffStacktrace
{
    SwitchOffStacktrace ()
      {
	deal_II_exceptions::suppress_stacktrace_in_exceptions ();
      }
} deal_II_stacktrace_dummy;


#endif // __tests_tests_h
