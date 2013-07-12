#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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
  #
  # Enable Fortran support so that MPI_Fortran_LIBRARIES is set up.
  #
  IF(NOT CMAKE_Fortran_COMPILER_WORKS)
    ENABLE_LANGUAGE(Fortran OPTIONAL)
  ENDIF()

  #
  # Obey a manual user override: If MPI_CXX_FOUND is set to true in the
  # cache, we skip the FIND_PACKAGE calls:
  #
  IF(MPI_CXX_FOUND)
    SET(MPI_FOUND TRUE)
  ENDIF()

  #
  # If CMAKE_CXX_COMPILER is already an MPI wrapper, use it to determine
  # the mpi implementation. If MPI_CXX_COMPILER is defined use the value
  # directly.
  #
  SET_IF_EMPTY(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  SET_IF_EMPTY(MPI_C_COMPILER ${CMAKE_C_COMPILER}) # for good measure
  SET_IF_EMPTY(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER}) # for good measure
  FIND_PACKAGE(MPI)

  IF(NOT MPI_CXX_FOUND)
    #
    # CMAKE_CXX_COMPILER is apparently not an mpi wrapper.
    # So, let's be a bit more aggressive in finding MPI if DEAL_II_WITH_MPI
    # is set.
    #
    IF(DEAL_II_WITH_MPI)
      MESSAGE(STATUS
        "MPI not found but DEAL_II_WITH_MPI is set to TRUE."
        " Try again with more aggressive search paths:"
        )
      SET(MPI_FOUND) # clear this value so that FIND_PACKAGE runs again.
      UNSET(MPI_CXX_COMPILER CACHE)
      UNSET(MPI_C_COMPILER CACHE)
      UNSET(MPI_Fortran_COMPILER CACHE)
      FIND_PACKAGE(MPI)
    ENDIF()
  ENDIF()

  #
  # Manually clean up variables:
  #
  FOREACH(_lang C CXX Fortran)
    IF(MPI_${_lang}_LIBRARIES MATCHES "-NOTFOUND")
      SET(MPI_${_lang}_LIBRARIES)
    ENDIF()
  ENDFOREACH()

  IF(MPI_CXX_FOUND)
    # Hide some variables:
    MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY MPI_LIBRARY)

    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL)
  ADD_FLAGS(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  ADD_FLAGS(DEAL_II_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}")

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

  # The user has to know the location of the mpi headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
ENDMACRO()


MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable mpi library!\n"
    "Please ensure that an mpi library is installed on your computer\n"
    "and set CMAKE_CXX_COMPILER and CMAKE_C_COMPILER to the appropriate mpi\n"
    "wrappers:\n"
    "    $ CC=\".../mpicc\" CXX=\".../mpicxx\" cmake <...>\n"
    "    $ cmake -DCMAKE_C_COMPILER=\".../mpicc\" -DCMAKE_CXX_COMPIER=\".../mpicxx\" <...>\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MPI)
