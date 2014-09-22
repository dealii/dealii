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

#
# Logic copied from the distributed UMFPACK Makefile.
# (... where everything is explained contrary to this... mess)
#

INCLUDE_DIRECTORIES(${UMFPACK_BUNDLED_INCLUDE_DIRS})

#
# *Shush!*
# (Remove locally in this subdirectory)
#
STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wall")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-sign-compare")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-write-strings")

SET(src_umfpack_UMF
  umf_2by2.cc
  umf_analyze.cc
  umf_apply_order.cc
  umf_assemble.cc
  umf_blas3_update.cc
  umf_build_tuples.cc
  umf_colamd.cc
  umf_create_element.cc
  umf_dump.cc
  umf_extend_front.cc
  umf_free.cc
  umf_fsize.cc
  umf_garbage_collection.cc
  umf_get_memory.cc
  umf_grow_front.cc
  umf_init_front.cc
  umf_is_permutation.cc
  umf_kernel.cc
  umf_kernel_init.cc
  umf_kernel_wrapup.cc
  umf_local_search.cc
  umf_lsolve.cc
  umf_ltsolve.cc
  umf_malloc.cc
  umf_mem_alloc_element.cc
  umf_mem_alloc_head_block.cc
  umf_mem_alloc_tail_block.cc
  umf_mem_free_tail_block.cc
  umf_mem_init_memoryspace.cc
  #
  # Don't ever think about including
  # umf_multicompile.cc, we do it by hand below
  #
  umf_realloc.cc
  umf_report_perm.cc
  umf_report_vector.cc
  umf_row_search.cc
  umf_scale.cc
  umf_scale_column.cc
  umf_set_stats.cc
  umf_singletons.cc
  umf_solve.cc
  umf_start_front.cc
  umf_store_lu.cc
  umf_symbolic_usage.cc
  umf_transpose.cc
  umf_triplet.cc
  umf_tuple_lengths.cc
  umf_usolve.cc
  umf_utsolve.cc
  umf_valid_numeric.cc
  umf_valid_symbolic.cc
  )

SET(src_umfpack_UMFPACK
  umfpack_col_to_triplet.cc
  umfpack_defaults.cc
  umfpack_free_numeric.cc
  umfpack_free_symbolic.cc
  umfpack_get_determinant.cc
  umfpack_get_lunz.cc
  umfpack_get_numeric.cc
  umfpack_get_symbolic.cc
  umfpack_load_numeric.cc
  umfpack_load_symbolic.cc
  umfpack_numeric.cc
  umfpack_qsymbolic.cc
  umfpack_report_control.cc
  umfpack_report_info.cc
  umfpack_report_matrix.cc
  umfpack_report_numeric.cc
  umfpack_report_perm.cc
  umfpack_report_status.cc
  umfpack_report_symbolic.cc
  umfpack_report_triplet.cc
  umfpack_report_vector.cc
  umfpack_save_numeric.cc
  umfpack_save_symbolic.cc
  umfpack_scale.cc
  umfpack_solve.cc
  umfpack_symbolic.cc
  umfpack_transpose.cc
  umfpack_triplet_to_col.cc
  )

SET(src_umfpack_GENERAL
  umfpack_global.cc
  umfpack_tictoc.cc
  umfpack_timer.cc
  )


#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

DEAL_II_ADD_LIBRARY(obj_umfpack_I_UMF OBJECT ${src_umfpack_UMF})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_I_UMF "DINT")

DEAL_II_ADD_LIBRARY(obj_umfpack_L_UMF OBJECT ${src_umfpack_UMF})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_L_UMF "DLONG")

DEAL_II_ADD_LIBRARY(obj_umfpack_I_UMFPACK OBJECT ${src_umfpack_UMFPACK})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_I_UMFPACK "DINT")

DEAL_II_ADD_LIBRARY(obj_umfpack_L_UMFPACK OBJECT ${src_umfpack_UMFPACK})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_L_UMFPACK "DLONG")

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_TSOLVE OBJECT umf_ltsolve.cc umf_utsolve.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TSOLVE "DINT" "CONJUGATE_SOLVE")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_TRIPLET_MAP_NOX OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_MAP_NOX "DINT" "DO_MAP")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_TRIPLET_MAP_X OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_MAP_X "DINT" "DO_MAP" "DO_VALUES")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_TRIPLET_NOMAP_X OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_NOMAP_X "DINT" "DO_VALUES")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_STORE OBJECT umf_store_lu.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_STORE "DINT" "DROP")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_ASSEMBLE OBJECT umf_assemble.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_ASSEMBLE "DINT" "FIXQ")

DEAL_II_ADD_LIBRARY(obj_umfpack_DI_SOLVE OBJECT umfpack_solve.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_SOLVE "DINT" "WSOLVE")

#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_TSOLVE OBJECT umf_ltsolve.cc umf_utsolve.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TSOLVE "DLONG" "CONJUGATE_SOLVE")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_TRIPLET_MAP_NOX OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_MAP_NOX "DLONG" "DO_MAP")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_TRIPLET_MAP_X OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_MAP_X "DLONG" "DO_MAP" "DO_VALUES")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_TRIPLET_NOMAP_X OBJECT umf_triplet.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_NOMAP_X "DLONG" "DO_VALUES")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_STORE OBJECT umf_store_lu.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_STORE "DLONG" "DROP")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_ASSEMBLE OBJECT umf_assemble.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_ASSEMBLE "DLONG" "FIXQ")

DEAL_II_ADD_LIBRARY(obj_umfpack_DL_SOLVE OBJECT umfpack_solve.cc)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_SOLVE "DLONG" "WSOLVE")

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

DEAL_II_ADD_LIBRARY(obj_umfpack_GENERIC OBJECT ${src_umfpack_GENERAL})

