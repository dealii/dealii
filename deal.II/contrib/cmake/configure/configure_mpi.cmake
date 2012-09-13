FIND_PACKAGE(MPI REQUIRED CXX)

# TODO: A deal.II specific error message if mpi is not found

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

LIST(APPEND deal_ii_external_libraries ${MPI_CXX_LIBRARIES})

SET(DEAL_II_COMPILER_SUPPORTS_MPI TRUE)


#   MPI_CXX_COMPILER        MPI Compiler wrapper for CXX
#   MPI_CXX_COMPILE_FLAGS   Compilation flags for MPI programs


#   MPIEXEC                    Executable for running MPI programs
#   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving
#                              it the number of processors to run on
#   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly
#                              before the executable to run.
#   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after other flags
