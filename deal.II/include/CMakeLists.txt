#
# Configure config.h
#
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/base/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II/base/config.h
  )


#
#TODO: This perl script has to die...:
#
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/
  COMMAND perl
  ARGS ${CMAKE_SOURCE_DIR}/scripts/lapack_templates.pl
       ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/lac/lapack_templates.h.in
       > ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  )

ADD_CUSTOM_TARGET(lapack_templates ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  )

#
# Add a rule on how to install the header files:
#

INSTALL(DIRECTORY deal.II
  DESTINATION ${DEAL_II_INCLUDE_RELDIR}
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
  )

#
# and don't forget to install all generated header files, too:
#
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_INCLUDE_RELDIR}
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
  )
