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

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#include <mpi.h>
#endif

#ifdef DEAL_II_USE_MUMPS
#  include <base/utilities.h>
#  include <dmumps_c.h>
#endif 

#ifdef DEAL_II_USE_SLEPC
#  include <slepcversion.h>
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

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#  ifdef OMPI_MAJOR_VERSION
  std::cout << "dealii-feature: MPI=OpenMPI<br>"
	    << OMPI_MAJOR_VERSION << '.'
	    << OMPI_MINOR_VERSION << '.'
	    << OMPI_RELEASE_VERSION << std::endl;
#  else
  std::cout << "dealii-feature: MPI="
	    << MPI__VERSION << '.'
	    << MPI_SUBVERSION << std::endl;  
#  endif
#endif

#ifdef DEAL_II_USE_MUMPS
  std::cout << "dealii-feature: MUMPS=yes" << std::endl;
#endif

#ifdef DEAL_II_USE_SLEPC
  std::cout << "dealii-feature: SLEPc="
	    << SLEPC_VERSION_MAJOR << '.'
	    << SLEPC_VERSION_MINOR << '.'
	    << SLEPC_VERSION_SUBMINOR << 'p'
	    << SLEPC_VERSION_PATCH << std::endl;
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
