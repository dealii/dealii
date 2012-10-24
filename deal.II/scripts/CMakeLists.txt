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

ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Build and install the old dependency_resolution and report_features
  # executables:
  #

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_BINARY_DIR}/include/
    )

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  ADD_EXECUTABLE(report_features report_features.cc)
  TARGET_LINK_LIBRARIES(report_features ${DEAL_II_EXTERNAL_LIBRARIES})

  INSTALL(TARGETS make_dependencies report_features expand_instantiations
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )
ENDIF()

