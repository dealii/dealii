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
# Setup the testsuite:
#
INCLUDE(setup_testsuite)

#
# And include all subdirectories:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Proceed to test definitions now:")

#ADD_SUBDIRECTORY(a-framework) -- no idea what to do with that ...

ADD_SUBDIRECTORY(all-headers)
ADD_SUBDIRECTORY(aniso)
ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(bits)
ADD_SUBDIRECTORY(build_tests)
ADD_SUBDIRECTORY(codim_one)
ADD_SUBDIRECTORY(deal.II)

IF(DEAL_II_WITH_P4EST)
  ADD_SUBDIRECTORY(distributed_grids)
ENDIF()

#ADD_SUBDIRECTORY(fail) -- the good old failing tests

ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(integrators)
ADD_SUBDIRECTORY(lac)

IF(DEAL_II_WITH_LAPACK)
  ADD_SUBDIRECTORY(lapack)
ENDIF()

ADD_SUBDIRECTORY(matrix_free)

IF(DEAL_II_WITH_METIS)
  ADD_SUBDIRECTORY(metis)
ENDIF()

# TODO: Really with P4EST?
IF(DEAL_II_WITH_MPI AND DEAL_II_WITH_P4EST)
  ADD_SUBDIRECTORY(gla)
  ADD_SUBDIRECTORY(mpi)
ENDIF()

ADD_SUBDIRECTORY(multigrid)

IF(DEAL_II_WITH_PETSC)
  ADD_SUBDIRECTORY(petsc)
ENDIF()

ADD_SUBDIRECTORY(serialization)

IF(DEAL_II_WITH_SLEPC)
  ADD_SUBDIRECTORY(slepc)
ENDIF()

IF(DEAL_II_WITH_TRILINOS)
  ADD_SUBDIRECTORY(trilinos)
ENDIF()

IF(DEAL_II_WITH_UMFPACK)
  ADD_SUBDIRECTORY(umfpack)
ENDIF()

