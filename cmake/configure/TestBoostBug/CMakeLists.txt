CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(TestBoost)

SET(BOOST_INCLUDE_DIRS CACHE STRING "The include paths for BOOST")
SET(BOOST_LIBRARIES CACHE STRING "The BOOST libraries to use")

# We assume the BOOST_LIBRARIES to be separated by '|'
# This matches the change in the separator we do in
# ../configure_boost.cmake
STRING(REPLACE "|" ";" BOOST_LIBRARIES_SEPARATED ${BOOST_LIBRARIES})
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIRS})

ADD_LIBRARY(derived SHARED polymorphic_derived2.cpp)
ADD_EXECUTABLE(test_boost test_dll_exported.cpp)

TARGET_LINK_LIBRARIES(derived ${BOOST_LIBRARIES_SEPARATED})
TARGET_LINK_LIBRARIES(test_boost derived ${BOOST_LIBRARIES_SEPARATED})

ADD_CUSTOM_TARGET(run
                  COMMAND test_boost
                  DEPENDS test_boost
                  WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})
