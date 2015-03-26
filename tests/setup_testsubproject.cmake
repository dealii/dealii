FIND_PACKAGE(deal.II 8.0 REQUIRED HINTS ${DEAL_II_DIR})

SET(CMAKE_BUILD_TYPE DebugRelease)
DEAL_II_INITIALIZE_CACHED_VARIABLES()

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
