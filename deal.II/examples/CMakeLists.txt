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

IF(DEAL_II_COMPONENT_EXAMPLES)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
    COMPONENT examples
    PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" EXCLUDE
    PATTERN "doc" EXCLUDE
    PATTERN "doxygen" EXCLUDE
    PATTERN ".svn" EXCLUDE
    PATTERN "TODO" EXCLUDE
    )
ENDIF()
