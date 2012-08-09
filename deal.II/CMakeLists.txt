PROJECT (deal.II)
set (VERSION "8.0.pre")
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.8)


#TODO:
SET (CMAKE_CXX_FLAGS "-std=c++0x -D_REENTRANT -fPIC -O2 -march=native")
SET (CMAKE_INCLUDE_CURRENT_DIR "true")


#TODO:
INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/include #TODO
  ${CMAKE_SOURCE_DIR}/contrib #TODO

  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )


ADD_SUBDIRECTORY(common/expand_instantiations)
#ADD_SUBDIRECTORY(contrib)

INCLUDE(${CMAKE_SOURCE_DIR}/contrib/cmake/macros/macro_expand_instantiations.cmake)

ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(source)
