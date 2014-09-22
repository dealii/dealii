## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by Martin Steigemann and Wolfgang Bangerth
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#PROJECT(parameter_gui)

IF("${DEAL_II_EXECUTABLE_RELDIR}" STREQUAL "")
  SET(DEAL_II_EXECUTABLE_RELDIR "bin")
ENDIF()

FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtXml)
MARK_AS_ADVANCED(QT_QMAKE_EXECUTABLE)

INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

QT4_WRAP_CPP(SOURCE_MOC
  browse_lineedit.h
  info_message.h
  mainwindow.h
  parameter_delegate.h
  )

QT4_ADD_RESOURCES(SOURCE_RCC
  application.qrc
  )

ADD_EXECUTABLE(parameter_gui_exe
  browse_lineedit.cpp
  info_message.cpp
  main.cpp
  mainwindow.cpp
  parameter_delegate.cpp
  xml_parameter_reader.cpp
  xml_parameter_writer.cpp
  ${SOURCE_MOC}
  ${SOURCE_RCC}
  )
SET_TARGET_PROPERTIES(parameter_gui_exe
  PROPERTIES
  OUTPUT_NAME parameter_gui
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
  )
TARGET_LINK_LIBRARIES(parameter_gui_exe ${QT_LIBRARIES})

ADD_DEPENDENCIES(parameter_gui parameter_gui_exe)

INSTALL(TARGETS parameter_gui_exe
  RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT parameter_gui
  )
EXPORT(TARGETS parameter_gui_exe
  FILE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  APPEND
  )
