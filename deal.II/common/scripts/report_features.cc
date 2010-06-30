//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include <base/config.h>

#include <iostream>

#ifdef HAVE_LIBUMFPACK
extern "C" {
#include <umfpack.h>
}
#endif

#ifdef DEAL_II_USE_MUMPS
#  include <base/utilities.h>
#  include <dmumps_c.h>
#endif 


int main()
{
#ifdef HAVE_LIBBLAS
  std::cout << "dealii-feature: BLAS=yes" << std::endl;
#endif
  
#ifdef HAVE_LIBLAPACK
  std::cout << "dealii-feature: LAPACK=yes" << std::endl;
#endif
  
#ifdef HAVE_LIBUMFPACK
  std::cout << "dealii-feature: UMFPACK="
	    << UMFPACK_MAIN_VERSION << '.'
	    << UMFPACK_SUB_VERSION << '.'
	    << UMFPACK_SUBSUB_VERSION << std::endl;
#endif

#if defined(HAVE_HSL_MA27) || defined(HAVE_HSL_MA47)
  std::cout << "dealii-feature: HSL=";
#ifdef HAVE_HSL_MA27
  std::cout << "MA27";
#endif
#ifdef HAVE_HSL_MA47
  std::cout << "MA47";
#endif
  std::cout << std::endl;
#endif
}
