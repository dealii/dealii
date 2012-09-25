#
# Configuration for mpi support:
#

MACRO(FEATURE_MPI_FIND_EXTERNAL var)
  FIND_PACKAGE(MPI)

  IF(MPI_CXX_FOUND)
    # Hide some variables:
    MARK_AS_ADVANCED(
      MPI_EXTRA_LIBRARY
      MPI_LIBRARY
      )

    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL var)

  ADD_FLAGS(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}")

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

  # The user has to know the location of the mpi headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})


  SET(DEAL_II_SET_MPI_COMPILER ON CACHE BOOL
    "Set compiler to the detected mpi wrapper"
    )
  MARK_AS_ADVANCED(DEAL_II_SET_MPI_COMPILER)

  IF(DEAL_II_SET_MPI_COMPILER)
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    SET(CMAKE_C_COMPILER   ${MPI_C_COMPILER})
  ENDIF()

  SET(DEAL_II_COMPILER_SUPPORTS_MPI TRUE)

  SET(${var} TRUE)
ENDMACRO()

SET(FEATURE_MPI_CUSTOM_ERROR_MESSAGE TRUE)

MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "Could not find any suitable mpi library!\n\n"
    "Please ensure that an mpi library is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection, or set the relevant variables by hand in ccmake.\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MPI)

