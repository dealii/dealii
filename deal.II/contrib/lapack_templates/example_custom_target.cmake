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
