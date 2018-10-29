## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# This file sets up the pkgconfig configuration files consisting of
#
#   deal.II_base       # common part to all the rest
#   deal.II_debug      # generated when buildtype is debug
#   deal.II_release    # generated when buildtype is release
#   deal.II_general    # generated when buildtype is different from above
#
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct setups.
#
# pkgconfig looks for *.pc files in an environmental variable called
# PKG_CONFIG_PATH. So, to use the library in the build directory issue
# export PKG_CONFIG_PATH=/path/to/BUILD_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#
# To use the library in the installed location
# export PKG_CONFIG_PATH=/path/to/INSTALL_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#


#
# For binary dir (AKA build dir):
#

set(DEAL_II_PATH_PC ${CMAKE_BINARY_DIR})

set (CONFIG_INCLUDE_DIRS_PC "")
foreach(_s ${CONFIG_INCLUDE_DIRS_BUILD_PC})
  string(CONCAT CONFIG_INCLUDE_DIRS_PC "${CONFIG_INCLUDE_DIRS_PC} -I${_s} "  )
endforeach()

# populate the variables
# CONFIG_LIBRARIES_DEBUG_PC 
# CONFIG_LIBRARIES_RELEASE_PC
# CONFIG_LIBRARIES_GENERAL_PC
foreach (_build ${DEAL_II_BUILD_TYPES})
  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_build "general")
  else()
    if(_build MATCHES DEBUG)
      set(_build "debug")
    else()
      set(_build "release")
    endif()
  endif()

  string(TOUPPER ${_build} _B)
  set(PC_LIBRARIES ${CONFIG_LIBRARIES_${_B}})
  set(CONFIG_LIBRARIES_PC "")
  foreach(_s ${PC_LIBRARIES})
    # check if _s contains lib
    # if _s is somethink like /path/to/trilinos/lib/libteuchos.so
    # we have  to add /path/to/trilinos/lib to the runpath
    if(${_s} MATCHES "lib")
      string(CONCAT CONFIG_LIBRARIES_PC " ${CONFIG_LIBRARIES_PC} ${_s} ")

      # work for the runpath
      string(FIND ${_s} "/" _i REVERSE)
      string(SUBSTRING ${_s} 0 ${_i} _p)
      list(APPEND PC_RPATH ${_p})
    else()
      # here _s is something like "m" for the libm.so so we add the -l flag
      string(CONCAT CONFIG_LIBRARIES_PC " ${CONFIG_LIBRARIES_PC} -l${_s} "  )
    endif()
  endforeach()

  set(CONFIG_LIBRARIES_${_B}_PC ${CONFIG_LIBRARIES_PC})
  
  list(REMOVE_DUPLICATES PC_RPATH)
  foreach(_r ${PC_RPATH})
    string(CONCAT CONFIG_LIBRARIES_${_B}_PC " ${CONFIG_LIBRARIES_${_B}_PC} -Wl,-rpath,${_r} "  )
  endforeach()

endforeach()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_base.pc.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}/pkgconfig/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  @ONLY
  )

# generate the requested .pc files
foreach (_build ${DEAL_II_BUILD_TYPES})
  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_build "general")
  else()
    if(_build MATCHES DEBUG)
      set(_build "debug")
    else()
      set(_build "release")
    endif()
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_${_build}.pc.in
    ${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}/pkgconfig/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    @ONLY
    )
endforeach()

#############################################################################
#
# For installation:
#

set(DEAL_II_PATH_PC ${CMAKE_INSTALL_PREFIX})

set (CONFIG_INCLUDE_DIRS_PC "")
foreach(_s ${CONFIG_INCLUDE_DIRS_INSTALL_PC})
  string(CONCAT CONFIG_INCLUDE_DIRS_PC "${CONFIG_INCLUDE_DIRS_PC} -I${_s} "  )
endforeach()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_base.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  @ONLY
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  DESTINATION ${DEAL_II_LIBRARY_RELDIR}/pkgconfig
  COMPONENT library
  )


# generate the requested .pc files
foreach (_build ${DEAL_II_BUILD_TYPES})
  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_build "general")
  else()
    if(_build MATCHES DEBUG)
      set(_build "debug")
    else()
      set(_build "release")
    endif()
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_${_build}.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    @ONLY
    )
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    DESTINATION ${DEAL_II_LIBRARY_RELDIR}/pkgconfig
    COMPONENT library
    )
endforeach()

