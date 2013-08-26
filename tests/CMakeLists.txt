## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Setup and run the testsuite:
#
# TODO: Describe and document the following:
# DEAL_II_TEST_DIFF
# DEAL_II_TEST_TIME_LIMIT
# NUMDIFF_DIR
#

#
# We need perl for the testsuite:
#
FIND_PACKAGE(Perl)
IF(NOT PERL_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find a perl installation which is required for running the test suite"
    )
ENDIF()

#
# And a diff tool, preferably numdiff:
#
FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )
FIND_PROGRAM(DIFF_EXECUTABLE
  NAMES diff
  )
IF( NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND"
    AND DIFF_EXECUTABLE MATCHES "-NOTFOUND" )
  MESSAGE(FATAL_ERROR
    "Could not find diff or numdiff. One of those are required for running the testsuite."
    )
ENDIF()

IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
  SET_IF_EMPTY(DEAL_II_TEST_DIFF numdiff -a 1e-6 -q -s ' \t\n:')
ELSE()
  SET_IF_EMPTY(DEAL_II_TEST_DIFF diff)
ENDIF()

#
# Set a default time limit of 60 seconds:
#
SET_IF_EMPTY(DEAL_II_TEST_TIME_LIMIT 120)


MESSAGE(STATUS "Proceed to test definitions")

#$(if $(DEAL_II_USE_MPI:no=),      $(if $(USE_CONTRIB_P4EST:no=),mpi,),) \
#$(if $(DEAL_II_USE_MPI:no=),      $(if $(USE_CONTRIB_P4EST:no=),gla,),) \
#$(if $(USE_CONTRIB_P4EST:no=),    distributed_grids,) \
#$(if $(USE_CONTRIB_PETSC:no=),    petsc,) \
#$(if $(USE_CONTRIB_SLEPC:no=),    slepc,) \
#$(if $(USE_CONTRIB_TRILINOS:no=), trilinos,) \

#a-framework
#all-headers
ADD_SUBDIRECTORY(aniso)
ADD_SUBDIRECTORY(base)
#benchmarks
#bits
ADD_SUBDIRECTORY(codim_one)
ADD_SUBDIRECTORY(deal.II)
#fail
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
#hp
#integrators
ADD_SUBDIRECTORY(lac)
#lib
ADD_SUBDIRECTORY(matrix_free)
#mesh_converter
ADD_SUBDIRECTORY(multigrid)
#serialization

#distributed_grids
#gla
#mpi
#petsc
#slepc
#trilinos

IF(DEAL_II_WITH_LAPACK)
  ADD_SUBDIRECTORY(lapack)
ENDIF()

IF(DEAL_II_WITH_METIS)
  ADD_SUBDIRECTORY(metis)
ENDIF()

IF(DEAL_II_WITH_UMFPACK)
  ADD_SUBDIRECTORY(umfpack)
ENDIF()
