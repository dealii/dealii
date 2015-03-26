FIND_PACKAGE(deal.II 8.0 REQUIRED HINTS ${DEAL_II_DIR})
SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER} CACHE STRING "CXX Compiler.")

#
# Silence warnings:
#
FOREACH(_var MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_POSTFLAGS MPIEXEC_PREFLAGS)
  SET(${_var} ${${_var}})
ENDFOREACH()

#
# A custom target that does absolutely nothing. It is used in the main
# project to trigger a "make rebuild_cache" if necessary.
#
ADD_CUSTOM_TARGET(regenerate)
