#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#PROJECT(parameter_gui)

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

ADD_EXECUTABLE(parameter_gui
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

TARGET_LINK_LIBRARIES(parameter_gui ${QT_LIBRARIES})

IF("${DEAL_II_EXECUTABLE_RELDIR}" STREQUAL "")
  SET(DEAL_II_EXECUTABLE_RELDIR "bin")
ENDIF()

INSTALL(TARGETS parameter_gui
  RUNTIME
  DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT contrib
  )

