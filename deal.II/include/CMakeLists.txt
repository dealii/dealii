CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/lac/lapack_templates.h.in
  ${CMAKE_BINARY_DIR}/configured/deal.II/lac/lapack_templates.h #TODO
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/base/config.h.in
  ${CMAKE_BINARY_DIR}/configured/deal.II/base/config.h #TODO
  )
