## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

IF(DEAL_II_COMPONENT_PACKAGE)
  MESSAGE(STATUS "Setting up CPack")
  SET(CPACK_GENERATOR "Bundle")

  CONFIGURE_FILE(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/mac_startup_script.sh.in
    ${CMAKE_BINARY_DIR}/cpack/mac_startup_script.sh
    @ONLY
    )

  CONFIGURE_FILE(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-terminal.in
    ${CMAKE_BINARY_DIR}/cpack/dealii-terminal
    @ONLY
    )

  CONFIGURE_FILE(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii.conf.in
    ${CMAKE_BINARY_DIR}/cpack/dealii.conf
    @ONLY
    )

  CONFIGURE_FILE(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/Info.plist.in
    ${CMAKE_BINARY_DIR}/cpack/Info.plist
    @ONLY
    )

  SET(CPACK_PACKAGE_ICON
    "${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-icon.icns"
    )

  set(CPACK_PACKAGE_FILE_NAME
    "dealii-${DEAL_II_PACKAGE_VERSION}"
    )
  MESSAGE(STATUS "  Disk filename: ${CPACK_PACKAGE_FILE_NAME}.dmg")

  set(CPACK_BUNDLE_NAME
    "${DEAL_II_CPACK_BUNDLE_NAME}"
    )
  MESSAGE(STATUS "  Application: ${DEAL_II_CPACK_BUNDLE_NAME}.app")

  SET(CPACK_BUNDLE_ICON
    "${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-icon.icns"
    )

  SET(CPACK_BUNDLE_PLIST
    "${CMAKE_BINARY_DIR}/cpack/Info.plist"
    )

  SET(CPACK_BUNDLE_STARTUP_COMMAND
    "${CMAKE_BINARY_DIR}/cpack/mac_startup_script.sh"
    )

  INSTALL(FILES
    ${CMAKE_BINARY_DIR}/cpack/dealii.conf
    DESTINATION ${DEAL_II_SHARE_RELDIR}
    )

  INSTALL(PROGRAMS
    ${CMAKE_BINARY_DIR}/cpack/dealii-terminal
    DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    )

  IF(NOT "${DEAL_II_CPACK_EXTERNAL_LIBS}" STREQUAL "")
    SET(_SRC "/Applications/${DEAL_II_CPACK_BUNDLE_NAME}.app/Contents/Resources/${DEAL_II_CPACK_EXTERNAL_LIBS}/")
    IF(IS_DIRECTORY ${_SRC})
       MESSAGE(STATUS "  Will copy ${_SRC} *as is* in the generated package")
       INSTALL(DIRECTORY ${_SRC}
         DESTINATION ${DEAL_II_CPACK_EXTERNAL_LIBS}
         USE_SOURCE_PERMISSIONS
         )
    ENDIF()
  ENDIF()

  INCLUDE(CPack)
  MESSAGE(STATUS "Setting up CPack - Done")
ENDIF()
