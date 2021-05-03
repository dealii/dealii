CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)
PROJECT(TestBoost)

INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIRS})

ADD_LIBRARY(derived SHARED polymorphic_derived2.cpp)
ADD_EXECUTABLE(test_boost test_dll_exported.cpp)

TARGET_LINK_LIBRARIES(derived ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(test_boost derived ${BOOST_LIBRARIES})

ADD_CUSTOM_TARGET(run
                  COMMAND test_boost
                  DEPENDS test_boost
                  WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})
