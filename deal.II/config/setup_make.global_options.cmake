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


###########################################################################
#                                                                         #
#          Configure and install the old Makefile configuration:          #
#                                                                         #
###########################################################################

#
# Transform some cmake lists into a string that the old Makefile
# mechanism actually understands:
#
TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_C}"
  ${DEAL_II_USER_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

FOREACH(build ${DEAL_II_BUILD_TYPES})
  TO_STRING_AND_ADD_PREFIX(MAKEFILE_DEFINITIONS_${build}
    "-D"
    ${DEAL_II_USER_DEFINITIONS}
    ${DEAL_II_USER_DEFINITIONS_${build}}
    )

  #
  # Add an rpath directive in front of each library, so that libraries
  # outside of the default search directories will be found by the runtime
  # dynamic linker:
  #
  FOREACH(lib
      ${MAKEFILE_LIBRARY_${build}}
      ${DEAL_II_EXTERNAL_LIBRARIES}
      ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
      )
    GET_FILENAME_COMPONENT(path ${lib} PATH)
    LIST(APPEND MAKEFILE_LIBS_${build} "-Wl,-rpath -Wl,${path}")
    LIST(APPEND MAKEFILE_LIBS_${build} ${lib})
  ENDFOREACH()
  TO_STRING(MAKEFILE_LIBS_${build} ${MAKEFILE_LIBS_${build}})
  #
  # Put our linker flags directly in front of this string:
  #
  SET(MAKEFILE_LIBS_${build}
    "${CMAKE_SHARED_LINKER_FLAGS} ${DEAL_II_SHARED_LINKER_FLAGS_${build}} ${MAKEFILE_LIBS_${build}}"
    )
ENDFOREACH()

#
# Boilerplate: The Make.global_options expects variables to be set to
# yes, as is common for Makefiles.
#
COND_SET_TO_YES(DEAL_II_WITH_TBB MAKEFILE_enablethreads)
COND_SET_TO_YES(DEAL_II_WITH_FUNCTIONPARSER MAKEFILE_enableparser)
COND_SET_TO_YES(BUILD_SHARED_LIBS MAKEFILE_enableshared)

COND_SET_TO_YES(DEAL_II_WITH_PETSC MAKEFILE_PETSC)
COND_SET_TO_YES(DEAL_II_USE_PETSC_DEV MAKEFILE_PETSC_DEV)
COND_SET_TO_YES(DEAL_II_WITH_TRILINOS MAKEFILE_TRILINOS)
COND_SET_TO_YES(DEAL_II_WITH_BLAS MAKEFILE_BLAS)
COND_SET_TO_YES(DEAL_II_WITH_LAPACK MAKEFILE_LAPACK)
COND_SET_TO_YES(DEAL_II_WITH_ARPACK MAKEFILE_ARPACK)
COND_SET_TO_YES(DEAL_II_WITH_METIS MAKEFILE_METIS)
COND_SET_TO_YES(DEAL_II_WITH_UMFPACK MAKEFILE_UMFPACK)
COND_SET_TO_YES(DEAL_II_WITH_P4EST MAKEFILE_P4EST)
COND_SET_TO_YES(DEAL_II_WITH_MPI MAKEFILE_MPI)

CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/config/Make.global_options.in
  ${CMAKE_BINARY_DIR}/config/Make.global_options
  )

CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/config/Version.in
  ${CMAKE_BINARY_DIR}/config/Version
  )

INSTALL(FILES
  ${CMAKE_BINARY_DIR}/config/template-arguments
  ${CMAKE_BINARY_DIR}/config/Make.global_options
  DESTINATION common
  COMPONENT compat_files
  )

INSTALL(FILES
  ${CMAKE_BINARY_DIR}/config/Version
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  COMPONENT compat_files
  )

