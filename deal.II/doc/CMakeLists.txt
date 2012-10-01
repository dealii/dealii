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

#
# Set up all necessary bits for the documentation
#

IF(DEAL_II_COMPONENT_DOCUMENTATION)

  #
  # Generate the documentation via doxygen:
  #
  ADD_SUBDIRECTORY(doxygen)


  #
  # Install the homepage:
  #

  INSTALL(DIRECTORY .
    DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}
    COMPONENT documentation
    FILES_MATCHING REGEX "^.*(html|css|eps|fig|png)$"  # TODO
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/title.html.in
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    )
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}
    COMPONENT documentation
    )

ENDIF(DEAL_II_COMPONENT_DOCUMENTATION)

