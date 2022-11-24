cmake_minimum_required(VERSION 3.1.0)
project(TestBoost)

include_directories(${BOOST_INCLUDE_DIRS})

add_library(derived SHARED polymorphic_derived2.cpp)
add_executable(test_boost test_dll_exported.cpp)

target_link_libraries(derived ${BOOST_LIBRARIES})
target_link_libraries(test_boost derived ${BOOST_LIBRARIES})

add_custom_target(run
                  COMMAND test_boost
                  DEPENDS test_boost
                  WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})
