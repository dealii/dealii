#####
##
## Copyright (C) 2013 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/algorithms/*.h
  )

SET(_src
  operator.cc
  timestep_control.cc
  )

SET(_inst
  operator.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_algorithms OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_algorithms "${_inst}")
