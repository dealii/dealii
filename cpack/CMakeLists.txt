## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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

  SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/cpack/mac-bundle/dealii-icon.icns")
  set(CPACK_PACKAGE_FILE_NAME "${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}")
  set(CPACK_BUNDLE_NAME "${DEAL_II_PACKAGE_NAME}")
  SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/cpack/mac-bundle/dealii-icon.icns")
  SET(CPACK_BUNDLE_PLIST "${CMAKE_SOURCE_DIR}/cpack/mac-bundle/Info.plist")
  SET(CPACK_BUNDLE_STARTUP_COMMAND "${CMAKE_SOURCE_DIR}/cpack/mac-bundle/mac_startup_script.sh")

  INSTALL(FILES ${CMAKE_SOURCE_DIR}/cpack/mac-bundle/dealii.conf DESTINATION etc/)
  INSTALL(FILES ${CMAKE_SOURCE_DIR}/cpack/mac-bundle/dealii.location DESTINATION etc/)
  INSTALL(PROGRAMS ${CMAKE_SOURCE_DIR}/cpack/mac-bundle/dealii-terminal DESTINATION bin/)

  OPTION(DEAL_II_INSTALLS_EXTERNAL_LIBS
         "Install external libraries together with deal.II."
         OFF
  )

  IF(DEAL_II_INSTALLS_EXTERNAL_LIBS)
    SET(DEAL_II_EXTERNAL_LIBS_TREE "NOT_SET" CACHE PATH
	"Path to tree of external libraries that will be installed in bundle package."
    )
    INSTALL(DIRECTORY ${DEAL_II_EXTERNAL_LIBS_TREE}/ 
            DESTINATION opt
            USE_SOURCE_PERMISSIONS)
  ENDIF()

  INCLUDE(CPack)
  MESSAGE(STATUS "Setting up CPack - Done")
ENDIF()
