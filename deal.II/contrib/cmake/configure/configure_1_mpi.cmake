#
# Configuration for mpi support:
#

MACRO(FEATURE_MPI_FIND_EXTERNAL var)
  FIND_PACKAGE(MPI)

  IF(MPI_CXX_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL var)

  ADD_FLAGS(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}")

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

  # The user has to know the location of the mpi headers as well:
  LIST(APPEND DEAL_II_EXTERNAL_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})

  # TODO: Set up the rest:

  #MPI_CXX_COMPILER        MPI Compiler wrapper for CXX
  #MPIEXEC                 Executable for running MPI programs
  #MPIEXEC_NUMPROC_FLAG    Flag to pass to MPIEXEC before giving it the number of processors to run on
  #MPIEXEC_PREFLAGS        Flags to pass to MPIEXEC directly before the executable to run.
  #MPIEXEC_POSTFLAGS       Flags to pass to MPIEXEC after other flags

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

