#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for mpi support:
#

MACRO(FEATURE_MPI_FIND_EXTERNAL var)
  FIND_PACKAGE(MPI QUIET)

  IF(MPI_CXX_FOUND)
    # Hide some variables:
    MARK_AS_ADVANCED(
      MPI_EXTRA_LIBRARY
      MPI_LIBRARY
      )
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL)
  ADD_FLAGS(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}")

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

  # The user has to know the location of the mpi headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})

  # TODO:
  SET(DEAL_II_USE_MPICXX ON CACHE BOOL
    "Set the compiler to the detected mpi wrapper ${MPI_CXX_COMPILER}"
    )

  IF(DEAL_II_MPICXX)
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable mpi library!\n"
    "Please ensure that an mpi library is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection, or set the relevant variables by hand in ccmake.\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MPI)
