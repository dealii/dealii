#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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

