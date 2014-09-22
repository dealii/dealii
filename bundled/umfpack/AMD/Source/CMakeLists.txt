## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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

INCLUDE_DIRECTORIES(${UMFPACK_BUNDLED_INCLUDE_DIRS})

SET(src_amd
  amd_1.cc
  amd_2.cc
  amd_aat.cc
  amd_control.cc
  amd_defaults.cc
  amd_dump.cc
  amd_info.cc
  amd_order.cc
  amd_postorder.cc
  amd_post_tree.cc
  amd_preprocess.cc
  amd_valid.cc
  )

DEAL_II_ADD_LIBRARY(obj_amd_int OBJECT ${src_amd})
DEAL_II_ADD_DEFINITIONS(obj_amd_int "DINT")

DEAL_II_ADD_LIBRARY(obj_amd_long OBJECT ${src_amd})
DEAL_II_ADD_DEFINITIONS(obj_amd_long "DLONG")

DEAL_II_ADD_LIBRARY(obj_amd_global OBJECT "amd_global.cc")

