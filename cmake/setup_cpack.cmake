## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

if(DEAL_II_COMPONENT_PACKAGE)
  message(STATUS "Setting up CPack")
  set(CPACK_GENERATOR "Bundle")

  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/mac_startup_script.sh.in
    ${CMAKE_BINARY_DIR}/cpack/mac_startup_script.sh
    @ONLY
    )

  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-terminal.in
    ${CMAKE_BINARY_DIR}/cpack/dealii-terminal
    @ONLY
    )

  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii.conf.in
    ${CMAKE_BINARY_DIR}/cpack/dealii.conf
    @ONLY
    )

  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/Info.plist.in
    ${CMAKE_BINARY_DIR}/cpack/Info.plist
    @ONLY
    )

  set(CPACK_PACKAGE_ICON
    "${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-icon.icns"
    )

  set(CPACK_PACKAGE_FILE_NAME
    "dealii-${DEAL_II_PACKAGE_VERSION}"
    )
  message(STATUS "  Disk filename: ${CPACK_PACKAGE_FILE_NAME}.dmg")

  set(CPACK_BUNDLE_NAME
    "${DEAL_II_CPACK_BUNDLE_NAME}"
    )
  message(STATUS "  Application: ${DEAL_II_CPACK_BUNDLE_NAME}.app")

  set(CPACK_BUNDLE_ICON
    "${CMAKE_SOURCE_DIR}/cmake/cpack-mac-bundle/dealii-icon.icns"
    )

  set(CPACK_BUNDLE_PLIST
    "${CMAKE_BINARY_DIR}/cpack/Info.plist"
    )

  set(CPACK_BUNDLE_STARTUP_COMMAND
    "${CMAKE_BINARY_DIR}/cpack/mac_startup_script.sh"
    )

  install(FILES
    ${CMAKE_BINARY_DIR}/cpack/dealii.conf
    DESTINATION ${DEAL_II_SHARE_RELDIR}
    )

  install(PROGRAMS
    ${CMAKE_BINARY_DIR}/cpack/dealii-terminal
    DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    )

  if(NOT "${DEAL_II_CPACK_EXTERNAL_LIBS}" STREQUAL "")
    set(_SRC "/Applications/${DEAL_II_CPACK_BUNDLE_NAME}.app/Contents/Resources/${DEAL_II_CPACK_EXTERNAL_LIBS}/")
    if(IS_DIRECTORY ${_SRC})
       message(STATUS "  Will copy ${_SRC} *as is* in the generated package")
       install(DIRECTORY ${_SRC}
         DESTINATION ${DEAL_II_CPACK_EXTERNAL_LIBS}
         USE_SOURCE_PERMISSIONS
         )
    endif()
  endif()

  include(CPack)
  message(STATUS "Setting up CPack - Done")
endif()
