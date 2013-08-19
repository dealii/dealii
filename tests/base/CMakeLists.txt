GET_FILENAME_COMPONENT(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)

FILE(GLOB _tests "anisotropic_1.cc")

FOREACH(_test ${_tests})
  GET_FILENAME_COMPONENT(_test ${_test} NAME_WE)
  DEAL_II_ADD_TEST(${_category} ${_test})
ENDFOREACH()
