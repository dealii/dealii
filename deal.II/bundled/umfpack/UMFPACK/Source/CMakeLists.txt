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

#
# Logic copied from the distributed UMFPACK Makefile.
# (... where everything is explained contrary to this... mess)
#

SET(src_umfpack_UMF
  umf_2by2.c
  umf_analyze.c
  umf_apply_order.c
  umf_assemble.c
  umf_blas3_update.c
  umf_build_tuples.c
  umf_colamd.c
  umf_create_element.c
  umf_dump.c
  umf_extend_front.c
  umf_free.c
  umf_fsize.c
  umf_garbage_collection.c
  umf_get_memory.c
  umf_grow_front.c
  umf_init_front.c
  umf_is_permutation.c
  umf_kernel.c
  umf_kernel_init.c
  umf_kernel_wrapup.c
  umf_local_search.c
  umf_lsolve.c
  umf_ltsolve.c
  umf_malloc.c
  umf_mem_alloc_element.c
  umf_mem_alloc_head_block.c
  umf_mem_alloc_tail_block.c
  umf_mem_free_tail_block.c
  umf_mem_init_memoryspace.c
  #
  # Don't ever think about including
  # umf_multicompile.c, we do it by hand below
  #
  umf_realloc.c
  umf_report_perm.c
  umf_report_vector.c
  umf_row_search.c
  umf_scale.c
  umf_scale_column.c
  umf_set_stats.c
  umf_singletons.c
  umf_solve.c
  umf_start_front.c
  umf_store_lu.c
  umf_symbolic_usage.c
  umf_transpose.c
  umf_triplet.c
  umf_tuple_lengths.c
  umf_usolve.c
  umf_utsolve.c
  umf_valid_numeric.c
  umf_valid_symbolic.c
  )

SET(src_umfpack_UMFPACK
  umfpack_col_to_triplet.c
  umfpack_defaults.c
  umfpack_free_numeric.c
  umfpack_free_symbolic.c
  umfpack_get_determinant.c
  umfpack_get_lunz.c
  umfpack_get_numeric.c
  umfpack_get_symbolic.c
  umfpack_load_numeric.c
  umfpack_load_symbolic.c
  umfpack_numeric.c
  umfpack_qsymbolic.c
  umfpack_report_control.c
  umfpack_report_info.c
  umfpack_report_matrix.c
  umfpack_report_numeric.c
  umfpack_report_perm.c
  umfpack_report_status.c
  umfpack_report_symbolic.c
  umfpack_report_triplet.c
  umfpack_report_vector.c
  umfpack_save_numeric.c
  umfpack_save_symbolic.c
  umfpack_scale.c
  umfpack_solve.c
  umfpack_symbolic.c
  umfpack_transpose.c
  umfpack_triplet_to_col.c
  )

SET(src_umfpack_GENERAL
  umfpack_global.c
  umfpack_tictoc.c
  umfpack_timer.c
  )


#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

DEAL_II_ADD_C_LIBRARY(obj_umfpack_I_UMF OBJECT ${src_umfpack_UMF})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_I_UMF "DINT")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_L_UMF OBJECT ${src_umfpack_UMF})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_L_UMF "DLONG")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_I_UMFPACK OBJECT ${src_umfpack_UMFPACK})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_I_UMFPACK "DINT")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_L_UMFPACK OBJECT ${src_umfpack_UMFPACK})
DEAL_II_ADD_DEFINITIONS(obj_umfpack_L_UMFPACK "DLONG")

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_TSOLVE OBJECT umf_ltsolve.c umf_utsolve.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TSOLVE "DINT" "CONJUGATE_SOLVE")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_TRIPLET_MAP_NOX OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_MAP_NOX "DINT" "DO_MAP")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_TRIPLET_MAP_X OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_MAP_X "DINT" "DO_MAP" "DO_VALUES")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_TRIPLET_NOMAP_X OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_TRIPLET_NOMAP_X "DINT" "DO_VALUES")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_STORE OBJECT umf_store_lu.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_STORE "DINT" "DROP")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_ASSEMBLE OBJECT umf_assemble.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_ASSEMBLE "DINT" "FIXQ")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DI_SOLVE OBJECT umfpack_solve.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DI_SOLVE "DINT" "WSOLVE")

#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_TSOLVE OBJECT umf_ltsolve.c umf_utsolve.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TSOLVE "DLONG" "CONJUGATE_SOLVE")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_TRIPLET_MAP_NOX OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_MAP_NOX "DLONG" "DO_MAP")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_TRIPLET_MAP_X OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_MAP_X "DLONG" "DO_MAP" "DO_VALUES")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_TRIPLET_NOMAP_X OBJECT umf_triplet.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_TRIPLET_NOMAP_X "DLONG" "DO_VALUES")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_STORE OBJECT umf_store_lu.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_STORE "DLONG" "DROP")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_ASSEMBLE OBJECT umf_assemble.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_ASSEMBLE "DLONG" "FIXQ")

DEAL_II_ADD_C_LIBRARY(obj_umfpack_DL_SOLVE OBJECT umfpack_solve.c)
DEAL_II_ADD_DEFINITIONS(obj_umfpack_DL_SOLVE "DLONG" "WSOLVE")

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

DEAL_II_ADD_C_LIBRARY(obj_umfpack_GENERIC OBJECT ${src_umfpack_GENERAL})

