## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


#
# Recurse into the sub-directories:
#
add_subdirectory(deal.II)

#
# For CMake versions >= 3.28 onwards, we attach the header files as a
# FILE_SET to the target directly, which is then installed in
# cmake/config/CMakeLists.txt alongside the target. For older CMake
# versions we need to install headers by hand:
#
if(CMAKE_VERSION VERSION_LESS 3.28)
  #
  # Add a rule for how to install the header files that are part of the source tree:
  #
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deal.II
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )

  #
  # and don't forget to install all generated header files, too:
  #
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deal.II
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )
endif()
