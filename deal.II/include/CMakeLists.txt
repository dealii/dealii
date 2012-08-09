CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/lac/lapack_templates.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/base/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II/base/config.h
  )
