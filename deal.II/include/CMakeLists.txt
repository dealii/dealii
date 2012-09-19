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

ADD_DEPENDENCIES(deal_ii_target_dependencies lapack_templates)


#
# Add a rule on how to install the header files:
#

INSTALL(DIRECTORY deal.II
  DESTINATION include
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  )

#
# and don't forget to install all generated header files, too:
#
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION include
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  )
