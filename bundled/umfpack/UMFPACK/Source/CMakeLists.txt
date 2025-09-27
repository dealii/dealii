## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# The logic of this file is copied from the distributed UMFPACK Makefile.
# Rather than duplicate the documentation from that file, the following
# copies the GNUmakefile of UMFPACK 5.0.3 (which is part of SuiteSparse
# 2.4.0). The structure of the CMake implementation follows that of
# the GNUmakefile.
#
# ..................................................
#-#  #-------------------------------------------------------------------------------
#-#  # UMFPACK Makefile for compiling on Unix systems (for GNU Make)
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  default: ../Lib/libumfpack.a
#-#  
#-#  include ../../UFconfig/UFconfig.mk
#-#  
#-#  C = $(CC) $(CFLAGS) $(UMFPACK_CONFIG) \
#-#      -I../Include -I../../AMD/Include -I../../UFconfig
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # source files
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  # non-user-callable umf_*.[ch] files:
#-#  UMFCH = umf_assemble umf_blas3_update umf_build_tuples umf_create_element \
#-#  	umf_dump umf_extend_front umf_garbage_collection umf_get_memory \
#-#  	umf_init_front umf_kernel umf_kernel_init umf_kernel_wrapup \
#-#  	umf_local_search umf_lsolve umf_ltsolve umf_mem_alloc_element \
#-#  	umf_mem_alloc_head_block umf_mem_alloc_tail_block \
#-#  	umf_mem_free_tail_block umf_mem_init_memoryspace \
#-#  	umf_report_vector umf_row_search umf_scale_column \
#-#  	umf_set_stats umf_solve umf_symbolic_usage umf_transpose \
#-#  	umf_tuple_lengths umf_usolve umf_utsolve umf_valid_numeric \
#-#  	umf_valid_symbolic umf_grow_front umf_start_front umf_2by2 \
#-#  	umf_store_lu umf_scale
#-#  
#-#  # non-user-callable umf_*.[ch] files, int/UF_long versions only (no real/complex):
#-#  UMFINT = umf_analyze umf_apply_order umf_colamd umf_free umf_fsize \
#-#  	umf_is_permutation umf_malloc umf_realloc umf_report_perm \
#-#  	umf_singletons
#-#  
#-#  # non-user-callable, created from umf_ltsolve.c, umf_utsolve.c,
#-#  # umf_triplet.c, and umf_assemble.c , with int/UF_long and real/complex versions:
#-#  UMF_CREATED = umf_lhsolve umf_uhsolve umf_triplet_map_nox \
#-#  	umf_triplet_nomap_x umf_triplet_nomap_nox umf_triplet_map_x \
#-#  	umf_assemble_fixq umf_store_lu_drop
#-#  
#-#  # non-user-callable, int/UF_long and real/complex versions:
#-#  UMF = $(UMF_CREATED) $(UMFCH)
#-#  
#-#  # user-callable umfpack_*.[ch] files (int/UF_long and real/complex):
#-#  UMFPACK = umfpack_col_to_triplet umfpack_defaults umfpack_free_numeric \
#-#  	umfpack_free_symbolic umfpack_get_numeric umfpack_get_lunz \
#-#  	umfpack_get_symbolic umfpack_get_determinant umfpack_numeric \
#-#  	umfpack_qsymbolic umfpack_report_control umfpack_report_info \
#-#  	umfpack_report_matrix umfpack_report_numeric umfpack_report_perm \
#-#  	umfpack_report_status umfpack_report_symbolic umfpack_report_triplet \
#-#  	umfpack_report_vector umfpack_solve umfpack_symbolic \
#-#  	umfpack_transpose umfpack_triplet_to_col umfpack_scale \
#-#  	umfpack_load_numeric umfpack_save_numeric \
#-#  	umfpack_load_symbolic umfpack_save_symbolic
#-#  
#-#  # user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
#-#  # with int/UF_long and real/complex versions:
#-#  UMFPACKW = umfpack_wsolve
#-#  
#-#  USER = $(UMFPACKW) $(UMFPACK)
#-#  
#-#  # user-callable, only one version for int/UF_long, real/complex, *.[ch] files:
#-#  GENERIC = umfpack_timer umfpack_tictoc umfpack_global
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # include files:
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  INC = ../Include/umfpack.h ../../UFconfig/UFconfig.h \
#-#  	umf_config.h umf_version.h umf_internal.h umf_triplet.h \
#-#  	$(addsuffix .h,$(UMFCH)) \
#-#  	$(addsuffix .h,$(UMFINT)) \
#-#  	$(addprefix ../Include/, $(addsuffix .h,$(USER))) \
#-#  	$(addprefix ../Include/, $(addsuffix .h,$(GENERIC))) \
#-#  	../../AMD/Include/amd_internal.h ../../AMD/Include/amd.h
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # object files for each version
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  DI = $(addsuffix .o, $(subst umf_,umf_di_,$(UMF)) $(subst umfpack_,umfpack_di_,$(USER)))
#-#  DL = $(addsuffix .o, $(subst umf_,umf_dl_,$(UMF)) $(subst umfpack_,umfpack_dl_,$(USER)))
#-#  ZI = $(addsuffix .o, $(subst umf_,umf_zi_,$(UMF)) $(subst umfpack_,umfpack_zi_,$(USER)))
#-#  ZL = $(addsuffix .o, $(subst umf_,umf_zl_,$(UMF)) $(subst umfpack_,umfpack_zl_,$(USER)))
#-#  II = $(addsuffix .o, $(subst umf_,umf_i_,$(UMFINT)))
#-#  LL = $(addsuffix .o, $(subst umf_,umf_l_,$(UMFINT)))
#-#  GN = $(addsuffix .o, $(subst umfpack_,umfpack_gn_,$(GENERIC)))
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # compile each int and UF_long routine (with no real/complex version)
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umf_i_%.o: umf_%.c $(INC)
#-#  	$(C) -DDINT -c $< -o $@
#-#  
#-#  umf_l_%.o: umf_%.c $(INC)
#-#  	$(C) -DDLONG -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # compile each routine in the DI version
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umf_di_%.o: umf_%.c $(INC)
#-#  	$(C) -DDINT -c $< -o $@
#-#  
#-#  umf_di_%hsolve.o: umf_%tsolve.c $(INC)
#-#  	$(C) -DDINT -DCONJUGATE_SOLVE -c $< -o $@
#-#  
#-#  umf_di_triplet_map_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DDINT -DDO_MAP -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_di_triplet_map_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DDINT -DDO_MAP -c $< -o $@
#-#  
#-#  umf_di_triplet_nomap_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DDINT -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_di_triplet_nomap_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DDINT -c $< -o $@
#-#  
#-#  umf_di_assemble_fixq.o: umf_assemble.c $(INC)
#-#  	$(C) -DDINT -DFIXQ -c $< -o $@
#-#  
#-#  umf_di_store_lu_drop.o: umf_store_lu.c $(INC)
#-#  	$(C) -DDINT -DDROP -c $< -o $@
#-#  
#-#  umfpack_di_wsolve.o: umfpack_solve.c $(INC)
#-#  	$(C) -DDINT -DWSOLVE -c $< -o $@
#-#  
#-#  umfpack_di_%.o: umfpack_%.c $(INC)
#-#  	$(C) -DDINT -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # compile each routine in the DL version
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umf_dl_%.o: umf_%.c $(INC)
#-#  	$(C) -DDLONG -c $< -o $@
#-#  
#-#  umf_dl_%hsolve.o: umf_%tsolve.c $(INC)
#-#  	$(C) -DDLONG -DCONJUGATE_SOLVE -c $< -o $@
#-#  
#-#  umf_dl_triplet_map_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DDLONG -DDO_MAP -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_dl_triplet_map_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DDLONG -DDO_MAP -c $< -o $@
#-#  
#-#  umf_dl_triplet_nomap_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DDLONG -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_dl_triplet_nomap_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DDLONG -c $< -o $@
#-#  
#-#  umf_dl_assemble_fixq.o: umf_assemble.c $(INC)
#-#  	$(C) -DDLONG -DFIXQ -c $< -o $@
#-#  
#-#  umf_dl_store_lu_drop.o: umf_store_lu.c $(INC)
#-#  	$(C) -DDLONG -DDROP -c $< -o $@
#-#  
#-#  umfpack_dl_wsolve.o: umfpack_solve.c $(INC)
#-#  	$(C) -DDLONG -DWSOLVE -c $< -o $@
#-#  
#-#  umfpack_dl_%.o: umfpack_%.c $(INC)
#-#  	$(C) -DDLONG -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # compile each routine in the ZI version
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umf_zi_%.o: umf_%.c $(INC)
#-#  	$(C) -DZINT -c $< -o $@
#-#  
#-#  umf_zi_%hsolve.o: umf_%tsolve.c $(INC)
#-#  	$(C) -DZINT -DCONJUGATE_SOLVE -c $< -o $@
#-#  
#-#  umf_zi_triplet_map_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DZINT -DDO_MAP -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_zi_triplet_map_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DZINT -DDO_MAP -c $< -o $@
#-#  
#-#  umf_zi_triplet_nomap_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DZINT -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_zi_triplet_nomap_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DZINT -c $< -o $@
#-#  
#-#  umf_zi_assemble_fixq.o: umf_assemble.c $(INC)
#-#  	$(C) -DZINT -DFIXQ -c $< -o $@
#-#  
#-#  umf_zi_store_lu_drop.o: umf_store_lu.c $(INC)
#-#  	$(C) -DZINT -DDROP -c $< -o $@
#-#  
#-#  umfpack_zi_wsolve.o: umfpack_solve.c $(INC)
#-#  	$(C) -DZINT -DWSOLVE -c $< -o $@
#-#  
#-#  umfpack_zi_%.o: umfpack_%.c $(INC)
#-#  	$(C) -DZINT -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # compile each routine in the ZL version
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umf_zl_%.o: umf_%.c $(INC)
#-#  	$(C) -DZLONG -c $< -o $@
#-#  
#-#  umf_zl_%hsolve.o: umf_%tsolve.c $(INC)
#-#  	$(C) -DZLONG -DCONJUGATE_SOLVE -c $< -o $@
#-#  
#-#  umf_zl_triplet_map_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DZLONG -DDO_MAP -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_zl_triplet_map_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DZLONG -DDO_MAP -c $< -o $@
#-#  
#-#  umf_zl_triplet_nomap_x.o: umf_triplet.c $(INC)
#-#  	$(C) -DZLONG -DDO_VALUES -c $< -o $@
#-#  
#-#  umf_zl_triplet_nomap_nox.o: umf_triplet.c $(INC)
#-#  	$(C) -DZLONG -c $< -o $@
#-#  
#-#  umf_zl_assemble_fixq.o: umf_assemble.c $(INC)
#-#  	$(C) -DZLONG -DFIXQ -c $< -o $@
#-#  
#-#  umf_zl_store_lu_drop.o: umf_store_lu.c $(INC)
#-#  	$(C) -DZLONG -DDROP -c $< -o $@
#-#  
#-#  umfpack_zl_wsolve.o: umfpack_solve.c $(INC)
#-#  	$(C) -DZLONG -DWSOLVE -c $< -o $@
#-#  
#-#  umfpack_zl_%.o: umfpack_%.c $(INC)
#-#  	$(C) -DZLONG -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # Create the generic routines (GN) using a generic rule
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  umfpack_gn_%.o: umfpack_%.c $(INC)
#-#  	$(C) -c $< -o $@
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # Create the libumfpack.a library
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  ../Lib/libumfpack.a: $(II) $(LL) $(GN) $(DI) $(DL) $(ZI) $(ZL)
#-#  	$(AR) ../Lib/libumfpack.a $^
#-#  	- $(RANLIB) ../Lib/libumfpack.a
#-#  
#-#  #-------------------------------------------------------------------------------
#-#  # Remove all but the files in the original distribution
#-#  #-------------------------------------------------------------------------------
#-#  
#-#  purge: clean
#-#  	- $(RM) ../Lib/libumfpack.a
#-#  
#-#  clean:
#-#  	- $(RM) $(CLEAN)
# ..................................................
#

#
# *Shush!*
# (Remove locally in this subdirectory)
#
strip_flag(DEAL_II_WARNING_FLAGS "-pedantic")
strip_flag(DEAL_II_WARNING_FLAGS "-Wall")
strip_flag(DEAL_II_WARNING_FLAGS "-Wextra")
strip_flag(DEAL_II_WARNING_FLAGS "-Wfloat-conversion")
enable_if_supported(DEAL_II_WARNING_FLAGS "-Wno-sign-compare")
enable_if_supported(DEAL_II_WARNING_FLAGS "-Wno-write-strings")


#
# We want to create real- and complex-valued versions of UMFPACK. Because we only
# ever call these functions with 'long int' arguments (not regular 'int' indices),
# that means that we need to compile all of the files that the GNUmakefile above
# lists as parts of the DL, ZL, LL, and GN target collections. These target
# collections depend on the following file lists:
# - DL: $(UMF) and $(USER); $(USER) is $(UMFPACK) + $(UMFPACKW), where the latter
#       only consists of a single file (umfpack_wsolve) that we don't need.
# - ZL: same
# - LL: $(UMFINT)
# - GN: $(GENERIC)
#
#
# Start by definining which files these collections include.
#
set(src_umfpack_DL_AND_ZL
  # Files that originate in $(UMF), which is really just $(UMFCH) for our
  # purposes:
  umf_assemble.cc umf_blas3_update.cc umf_build_tuples.cc umf_create_element.cc
  umf_dump.cc umf_extend_front.cc umf_garbage_collection.cc umf_get_memory.cc
  umf_init_front.cc umf_kernel.cc umf_kernel_init.cc umf_kernel_wrapup.cc
  umf_local_search.cc umf_lsolve.cc umf_ltsolve.cc umf_mem_alloc_element.cc
  umf_mem_alloc_head_block.cc umf_mem_alloc_tail_block.cc
  umf_mem_free_tail_block.cc umf_mem_init_memoryspace.cc
  umf_report_vector.cc umf_row_search.cc umf_scale_column.cc
  umf_set_stats.cc umf_solve.cc umf_symbolic_usage.cc umf_transpose.cc
  umf_tuple_lengths.cc umf_usolve.cc umf_utsolve.cc umf_valid_numeric.cc
  umf_valid_symbolic.cc umf_grow_front.cc umf_start_front.cc umf_2by2.cc
  umf_store_lu.cc umf_scale.cc
  # Files that originate in $(USER), which as explained above is really
  # just $(UMFPACK)
  umfpack_col_to_triplet.cc umfpack_defaults.cc umfpack_free_numeric.cc
  umfpack_free_symbolic.cc umfpack_get_numeric.cc umfpack_get_lunz.cc
  umfpack_get_symbolic.cc umfpack_get_determinant.cc umfpack_numeric.cc
  umfpack_qsymbolic.cc umfpack_report_control.cc umfpack_report_info.cc
  umfpack_report_matrix.cc umfpack_report_numeric.cc umfpack_report_perm.cc
  umfpack_report_status.cc umfpack_report_symbolic.cc umfpack_report_triplet.cc
  umfpack_report_vector.cc umfpack_solve.cc umfpack_symbolic.cc
  umfpack_transpose.cc umfpack_triplet_to_col.cc umfpack_scale.cc
  umfpack_load_numeric.cc umfpack_save_numeric.cc
  umfpack_load_symbolic.cc umfpack_save_symbolic.cc
  )

set(src_umfpack_LL
  umf_analyze.cc umf_apply_order.cc umf_colamd.cc umf_free.cc umf_fsize.cc
  umf_is_permutation.cc umf_malloc.cc umf_realloc.cc umf_report_perm.cc
  umf_singletons.cc
  )

set(src_umfpack_GENERIC
  umfpack_global.cc
  umfpack_tictoc.cc
  umfpack_timer.cc
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../Include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../AMD/Include
  )

#-------------------------------------------------------------------------------
# Compile each of the following files as long int routines; UMFPACK also allows
# compiling with int indices, but we don't use these and so let's not
# bother -- it only takes compile time unnecessarily.
#-------------------------------------------------------------------------------

# Do what the GNUmakefile calls the DI targets:
define_object_library(bundled_umfpack_L_UMF OBJECT ${src_umfpack_DL_AND_ZL})
deal_ii_add_definitions(bundled_umfpack_L_UMF "DLONG")

# Repeat these steps for complex-valued matrices. This corresponds to
# the ZL targets:
define_object_library(bundled_umfpack_Z_UMF OBJECT ${src_umfpack_DL_AND_ZL})
deal_ii_add_definitions(bundled_umfpack_Z_UMF "ZLONG")

# Do what the GNUmakefile calls the LL targets:
define_object_library(bundled_umfpack_L_UMFPACK OBJECT ${src_umfpack_LL})
deal_ii_add_definitions(bundled_umfpack_L_UMFPACK "DLONG")


#-------------------------------------------------------------------------------
# The following files are special. We again compile each routine in the 'long int'
# version.
#-------------------------------------------------------------------------------

define_object_library(bundled_umfpack_DL_TSOLVE OBJECT umf_ltsolve.cc umf_utsolve.cc)
deal_ii_add_definitions(bundled_umfpack_DL_TSOLVE "DLONG" "CONJUGATE_SOLVE")

define_object_library(bundled_umfpack_DL_TRIPLET_MAP_NOX OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_DL_TRIPLET_MAP_NOX "DLONG" "DO_MAP")

define_object_library(bundled_umfpack_DL_TRIPLET_MAP_X OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_DL_TRIPLET_MAP_X "DLONG" "DO_MAP" "DO_VALUES")

define_object_library(bundled_umfpack_DL_TRIPLET_NOMAP_X OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_DL_TRIPLET_NOMAP_X "DLONG" "DO_VALUES")

define_object_library(bundled_umfpack_DL_TRIPLET_NOMAP_NOX OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_DL_TRIPLET_NOMAP_NOX "DLONG")

define_object_library(bundled_umfpack_DL_STORE OBJECT umf_store_lu.cc)
deal_ii_add_definitions(bundled_umfpack_DL_STORE "DLONG" "DROP")

define_object_library(bundled_umfpack_DL_ASSEMBLE OBJECT umf_assemble.cc)
deal_ii_add_definitions(bundled_umfpack_DL_ASSEMBLE "DLONG" "FIXQ")

define_object_library(bundled_umfpack_DL_SOLVE OBJECT umfpack_solve.cc)
deal_ii_add_definitions(bundled_umfpack_DL_SOLVE "DLONG" "WSOLVE")

# Again with complex arguments:
define_object_library(bundled_umfpack_ZL_TSOLVE OBJECT umf_ltsolve.cc umf_utsolve.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_TSOLVE "ZLONG" "CONJUGATE_SOLVE")

define_object_library(bundled_umfpack_ZL_TRIPLET_MAP_NOX OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_TRIPLET_MAP_NOX "ZLONG" "DO_MAP")

define_object_library(bundled_umfpack_ZL_TRIPLET_MAP_X OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_TRIPLET_MAP_X "ZLONG" "DO_MAP" "DO_VALUES")

define_object_library(bundled_umfpack_ZL_TRIPLET_NOMAP_X OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_TRIPLET_NOMAP_X "ZLONG" "DO_VALUES")

define_object_library(bundled_umfpack_ZL_TRIPLET_NOMAP_NOX OBJECT umf_triplet.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_TRIPLET_NOMAP_NOX "ZLONG")

define_object_library(bundled_umfpack_ZL_STORE OBJECT umf_store_lu.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_STORE "ZLONG" "DROP")

define_object_library(bundled_umfpack_ZL_ASSEMBLE OBJECT umf_assemble.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_ASSEMBLE "ZLONG" "FIXQ")

define_object_library(bundled_umfpack_ZL_SOLVE OBJECT umfpack_solve.cc)
deal_ii_add_definitions(bundled_umfpack_ZL_SOLVE "ZLONG" "WSOLVE")

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

define_object_library(bundled_umfpack_GENERIC OBJECT ${src_umfpack_GENERIC})

