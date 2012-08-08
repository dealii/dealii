PROJECT (deal.II)
set (VERSION "8.0.pre")
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.8)

SET (CMAKE_CXX_FLAGS "-std=c++0x -D_REENTRANT -fPIC -O2 -march=native")

SET (CMAKE_INCLUDE_CURRENT_DIR "true")

INCLUDE_DIRECTORIES(
  "${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/include" #TODO
  "${CMAKE_SOURCE_DIR}/contrib" #TODO
  "${CMAKE_SOURCE_DIR}/include/"
  )

ADD_SUBDIRECTORY(common/expand_instantiations)

ADD_SUBDIRECTORY(source)
