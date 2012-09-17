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

  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  LIST(APPEND deal_ii_external_libraries ${MPI_CXX_LIBRARIES})

  SET(DEAL_II_COMPILER_SUPPORTS_MPI TRUE)


  # TODO: Set up the rest:

  #MPI_CXX_COMPILER        MPI Compiler wrapper for CXX
  #MPI_CXX_COMPILE_FLAGS   Compilation flags for MPI programs


  #MPIEXEC                    Executable for running MPI programs
  #MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving
  #                           it the number of processors to run on
  #MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly
  #                           before the executable to run.
  #MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after other flags

  SET(${var} TRUE)
ENDMACRO()

SET(FEATURE_MPI_CUSTOM_ERROR_MESSAGE TRUE)

MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find any suitable mpi library!

Please ensure that an mpi library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")
ENDMACRO()


CONFIGURE_FEATURE(MPI)
