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
  #
  # Add a rule on how to install the examples:
  #

  file(GLOB _steps ${CMAKE_CURRENT_SOURCE_DIR}/step-*)

  FOREACH(_step ${_steps})
    INSTALL(DIRECTORY ${_step}
      DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
      COMPONENT examples
      PATTERN ".svn" EXCLUDE
      PATTERN "doc" EXCLUDE
      PATTERN "TODO" EXCLUDE
      )
  ENDFOREACH()

  INSTALL(FILES GNUmakefile
    DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
    COMPONENT examples
    )

ENDIF()
