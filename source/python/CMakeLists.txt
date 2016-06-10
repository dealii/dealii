## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

FIND_LIBRARY(Boost_PYTHON_LIBRARY boost_python PATH ${Boost_LIBRARY_DIRS} NO_DEFAULT_PATH)
INCLUDE(FindPythonLibs)
INCLUDE(FindPythonInterp)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

SET(_src
  wrappers.cc
  export_point.cc
  export_triangulation.cc
  point_wrapper.cc
  triangulation_wrapper.cc
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/python/*.h
  )


PYTHON_ADD_MODULE(PyDealII ${_src})

SET(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

TARGET_LINK_LIBRARIES(PyDealII deal_II.g)
TARGET_LINK_LIBRARIES(PyDealII ${Boost_PYTHON_LIBRARY})
TARGET_LINK_LIBRARIES(PyDealII ${PYTHON_LIBRARIES})

SET_TARGET_PROPERTIES(PyDealII PROPERTIES
  CXX_STANDARD 11
  )

INSTALL(TARGETS PyDealII DESTINATION ${PYTHON_INSTALL_PREFIX}/pydealii)

SET(PYTHON_SOURCES
  __init__.py
  dummy.py
  )                                                                                                                                                                                                               
FOREACH(PYTHON_SOURCE ${PYTHON_SOURCES})                                                                                                                                                                            
  ADD_CUSTOM_COMMAND(                                                                                                                                                                                             
    OUTPUT ${CMAKE_BINARY_DIR}/python/pydealii/${PYTHON_SOURCE}                                                                                                                                                    
    DEPENDS ${CMAKE_SOURCE_DIR}/source/python/${PYTHON_SOURCE}                                                                                                                                                  
    COMMAND ${CMAKE_COMMAND}                                                                                                                                                                                    
    ARGS -E copy ${CMAKE_SOURCE_DIR}/source/python/${PYTHON_SOURCE} ${CMAKE_BINARY_DIR}/python/pydealii/${PYTHON_SOURCE}                                                                                           
    COMMENT "Copying ${PYTHON_SOURCE}"                                                                                                                                                                          
    )                                                                                                                                                                                                           
  ADD_CUSTOM_TARGET(                                                                                                                                                                                              
    ${PYTHON_SOURCE} ALL
    DEPENDS ${CMAKE_BINARY_DIR}/python/pydealii/${PYTHON_SOURCE}                                                                                                                                                   
    )                                                                                                                                                                                                           
ENDFOREACH()                                                                                                                                                                                                        
INSTALL(                                                                                                                                                                                                            
  DIRECTORY ${CMAKE_BINARY_DIR}/python/pydealii
  DESTINATION ${PYTHON_INSTALL_PREFIX}                                                                                                                                                                            
  FILES_MATCHING PATTERN "*.py"                                                                                                                                                                                   
  ) 
