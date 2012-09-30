#-------------------------------------------------------------------------------
# UMFPACK Makefile for the UMFPACK MATLAB mexFunction (old "make" only)
#-------------------------------------------------------------------------------

# This is a very ugly Makefile, and is only provided for those who do not
# have GNU make.  Note that it is not used if you have GNU make.  It ignores
# dependency checking and just compiles everything.

default: umfpack luflop

include ../../UFconfig/UFconfig.mk

MX = $(MEX) -I../Include -I../Source -I../../AMD/Include -I../../UFconfig

umfpack:
	$(MX) -c ../Source/umfpack_global.c
	$(MV) -f umfpack_global.o umfpack_m_global.o
	$(MX) -DDLONG -c ../Source/umf_analyze.c
	$(MV) -f umf_analyze.o umf_m_analyze.o
	$(MX) -DDLONG -c ../Source/umf_apply_order.c
	$(MV) -f umf_apply_order.o umf_m_apply_order.o
	$(MX) -DDLONG -c ../Source/umf_colamd.c
	$(MV) -f umf_colamd.o umf_m_colamd.o
	$(MX) -DDLONG -c ../Source/umf_free.c
	$(MV) -f umf_free.o umf_m_free.o
	$(MX) -DDLONG -c ../Source/umf_fsize.c
	$(MV) -f umf_fsize.o umf_m_fsize.o
	$(MX) -DDLONG -c ../Source/umf_is_permutation.c
	$(MV) -f umf_is_permutation.o umf_m_is_permutation.o
	$(MX) -DDLONG -c ../Source/umf_malloc.c
	$(MV) -f umf_malloc.o umf_m_malloc.o
	$(MX) -DDLONG -c ../Source/umf_realloc.c
	$(MV) -f umf_realloc.o umf_m_realloc.o
	$(MX) -DDLONG -c ../Source/umf_report_perm.c
	$(MV) -f umf_report_perm.o umf_m_report_perm.o
	$(MX) -DDLONG -c ../Source/umf_singletons.c
	$(MV) -f umf_singletons.o umf_m_singletons.o
	$(MX) -DDLONG -DCONJUGATE_SOLVE -c ../Source/umf_ltsolve.c
	$(MV) -f umf_ltsolve.o umf_md_lhsolve.o
	$(MX) -DDLONG -DCONJUGATE_SOLVE -c ../Source/umf_utsolve.c
	$(MV) -f umf_utsolve.o umf_md_uhsolve.o
	$(MX) -DDLONG -DDO_MAP -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_md_triplet_map_nox.o
	$(MX) -DDLONG -DDO_VALUES -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_md_triplet_nomap_x.o
	$(MX) -DDLONG -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_md_triplet_nomap_nox.o
	$(MX) -DDLONG -DDO_MAP -DDO_VALUES -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_md_triplet_map_x.o
	$(MX) -DDLONG -DFIXQ -c ../Source/umf_assemble.c
	$(MV) -f umf_assemble.o umf_md_assemble_fixq.o
	$(MX) -DDLONG -DDROP -c ../Source/umf_store_lu.c
	$(MV) -f umf_store_lu.o umf_md_store_lu_drop.o
	$(MX) -DDLONG -c ../Source/umf_assemble.c
	$(MV) -f umf_assemble.o umf_md_assemble.o
	$(MX) -DDLONG -c ../Source/umf_blas3_update.c
	$(MV) -f umf_blas3_update.o umf_md_blas3_update.o
	$(MX) -DDLONG -c ../Source/umf_build_tuples.c
	$(MV) -f umf_build_tuples.o umf_md_build_tuples.o
	$(MX) -DDLONG -c ../Source/umf_create_element.c
	$(MV) -f umf_create_element.o umf_md_create_element.o
	$(MX) -DDLONG -c ../Source/umf_dump.c
	$(MV) -f umf_dump.o umf_md_dump.o
	$(MX) -DDLONG -c ../Source/umf_extend_front.c
	$(MV) -f umf_extend_front.o umf_md_extend_front.o
	$(MX) -DDLONG -c ../Source/umf_garbage_collection.c
	$(MV) -f umf_garbage_collection.o umf_md_garbage_collection.o
	$(MX) -DDLONG -c ../Source/umf_get_memory.c
	$(MV) -f umf_get_memory.o umf_md_get_memory.o
	$(MX) -DDLONG -c ../Source/umf_init_front.c
	$(MV) -f umf_init_front.o umf_md_init_front.o
	$(MX) -DDLONG -c ../Source/umf_kernel.c
	$(MV) -f umf_kernel.o umf_md_kernel.o
	$(MX) -DDLONG -c ../Source/umf_kernel_init.c
	$(MV) -f umf_kernel_init.o umf_md_kernel_init.o
	$(MX) -DDLONG -c ../Source/umf_kernel_wrapup.c
	$(MV) -f umf_kernel_wrapup.o umf_md_kernel_wrapup.o
	$(MX) -DDLONG -c ../Source/umf_local_search.c
	$(MV) -f umf_local_search.o umf_md_local_search.o
	$(MX) -DDLONG -c ../Source/umf_lsolve.c
	$(MV) -f umf_lsolve.o umf_md_lsolve.o
	$(MX) -DDLONG -c ../Source/umf_ltsolve.c
	$(MV) -f umf_ltsolve.o umf_md_ltsolve.o
	$(MX) -DDLONG -c ../Source/umf_mem_alloc_element.c
	$(MV) -f umf_mem_alloc_element.o umf_md_mem_alloc_element.o
	$(MX) -DDLONG -c ../Source/umf_mem_alloc_head_block.c
	$(MV) -f umf_mem_alloc_head_block.o umf_md_mem_alloc_head_block.o
	$(MX) -DDLONG -c ../Source/umf_mem_alloc_tail_block.c
	$(MV) -f umf_mem_alloc_tail_block.o umf_md_mem_alloc_tail_block.o
	$(MX) -DDLONG -c ../Source/umf_mem_free_tail_block.c
	$(MV) -f umf_mem_free_tail_block.o umf_md_mem_free_tail_block.o
	$(MX) -DDLONG -c ../Source/umf_mem_init_memoryspace.c
	$(MV) -f umf_mem_init_memoryspace.o umf_md_mem_init_memoryspace.o
	$(MX) -DDLONG -c ../Source/umf_report_vector.c
	$(MV) -f umf_report_vector.o umf_md_report_vector.o
	$(MX) -DDLONG -c ../Source/umf_row_search.c
	$(MV) -f umf_row_search.o umf_md_row_search.o
	$(MX) -DDLONG -c ../Source/umf_scale_column.c
	$(MV) -f umf_scale_column.o umf_md_scale_column.o
	$(MX) -DDLONG -c ../Source/umf_set_stats.c
	$(MV) -f umf_set_stats.o umf_md_set_stats.o
	$(MX) -DDLONG -c ../Source/umf_solve.c
	$(MV) -f umf_solve.o umf_md_solve.o
	$(MX) -DDLONG -c ../Source/umf_symbolic_usage.c
	$(MV) -f umf_symbolic_usage.o umf_md_symbolic_usage.o
	$(MX) -DDLONG -c ../Source/umf_transpose.c
	$(MV) -f umf_transpose.o umf_md_transpose.o
	$(MX) -DDLONG -c ../Source/umf_tuple_lengths.c
	$(MV) -f umf_tuple_lengths.o umf_md_tuple_lengths.o
	$(MX) -DDLONG -c ../Source/umf_usolve.c
	$(MV) -f umf_usolve.o umf_md_usolve.o
	$(MX) -DDLONG -c ../Source/umf_utsolve.c
	$(MV) -f umf_utsolve.o umf_md_utsolve.o
	$(MX) -DDLONG -c ../Source/umf_valid_numeric.c
	$(MV) -f umf_valid_numeric.o umf_md_valid_numeric.o
	$(MX) -DDLONG -c ../Source/umf_valid_symbolic.c
	$(MV) -f umf_valid_symbolic.o umf_md_valid_symbolic.o
	$(MX) -DDLONG -c ../Source/umf_grow_front.c
	$(MV) -f umf_grow_front.o umf_md_grow_front.o
	$(MX) -DDLONG -c ../Source/umf_start_front.c
	$(MV) -f umf_start_front.o umf_md_start_front.o
	$(MX) -DDLONG -c ../Source/umf_2by2.c
	$(MV) -f umf_2by2.o umf_md_2by2.o
	$(MX) -DDLONG -c ../Source/umf_store_lu.c
	$(MV) -f umf_store_lu.o umf_md_store_lu.o
	$(MX) -DDLONG -c ../Source/umf_scale.c
	$(MV) -f umf_scale.o umf_md_scale.o
	$(MX) -DDLONG -DWSOLVE -c ../Source/umfpack_solve.c
	$(MV) -f umfpack_solve.o umfpack_md_wsolve.o
	$(MX) -DDLONG -c ../Source/umfpack_col_to_triplet.c
	$(MV) -f umfpack_col_to_triplet.o umfpack_md_col_to_triplet.o
	$(MX) -DDLONG -c ../Source/umfpack_defaults.c
	$(MV) -f umfpack_defaults.o umfpack_md_defaults.o
	$(MX) -DDLONG -c ../Source/umfpack_free_numeric.c
	$(MV) -f umfpack_free_numeric.o umfpack_md_free_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_free_symbolic.c
	$(MV) -f umfpack_free_symbolic.o umfpack_md_free_symbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_get_numeric.c
	$(MV) -f umfpack_get_numeric.o umfpack_md_get_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_get_lunz.c
	$(MV) -f umfpack_get_lunz.o umfpack_md_get_lunz.o
	$(MX) -DDLONG -c ../Source/umfpack_get_symbolic.c
	$(MV) -f umfpack_get_symbolic.o umfpack_md_get_symbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_get_determinant.c
	$(MV) -f umfpack_get_determinant.o umfpack_md_get_determinant.o
	$(MX) -DDLONG -c ../Source/umfpack_numeric.c
	$(MV) -f umfpack_numeric.o umfpack_md_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_qsymbolic.c
	$(MV) -f umfpack_qsymbolic.o umfpack_md_qsymbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_report_control.c
	$(MV) -f umfpack_report_control.o umfpack_md_report_control.o
	$(MX) -DDLONG -c ../Source/umfpack_report_info.c
	$(MV) -f umfpack_report_info.o umfpack_md_report_info.o
	$(MX) -DDLONG -c ../Source/umfpack_report_matrix.c
	$(MV) -f umfpack_report_matrix.o umfpack_md_report_matrix.o
	$(MX) -DDLONG -c ../Source/umfpack_report_numeric.c
	$(MV) -f umfpack_report_numeric.o umfpack_md_report_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_report_perm.c
	$(MV) -f umfpack_report_perm.o umfpack_md_report_perm.o
	$(MX) -DDLONG -c ../Source/umfpack_report_status.c
	$(MV) -f umfpack_report_status.o umfpack_md_report_status.o
	$(MX) -DDLONG -c ../Source/umfpack_report_symbolic.c
	$(MV) -f umfpack_report_symbolic.o umfpack_md_report_symbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_report_triplet.c
	$(MV) -f umfpack_report_triplet.o umfpack_md_report_triplet.o
	$(MX) -DDLONG -c ../Source/umfpack_report_vector.c
	$(MV) -f umfpack_report_vector.o umfpack_md_report_vector.o
	$(MX) -DDLONG -c ../Source/umfpack_solve.c
	$(MV) -f umfpack_solve.o umfpack_md_solve.o
	$(MX) -DDLONG -c ../Source/umfpack_symbolic.c
	$(MV) -f umfpack_symbolic.o umfpack_md_symbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_transpose.c
	$(MV) -f umfpack_transpose.o umfpack_md_transpose.o
	$(MX) -DDLONG -c ../Source/umfpack_triplet_to_col.c
	$(MV) -f umfpack_triplet_to_col.o umfpack_md_triplet_to_col.o
	$(MX) -DDLONG -c ../Source/umfpack_scale.c
	$(MV) -f umfpack_scale.o umfpack_md_scale.o
	$(MX) -DDLONG -c ../Source/umfpack_load_numeric.c
	$(MV) -f umfpack_load_numeric.o umfpack_md_load_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_save_numeric.c
	$(MV) -f umfpack_save_numeric.o umfpack_md_save_numeric.o
	$(MX) -DDLONG -c ../Source/umfpack_load_symbolic.c
	$(MV) -f umfpack_load_symbolic.o umfpack_md_load_symbolic.o
	$(MX) -DDLONG -c ../Source/umfpack_save_symbolic.c
	$(MV) -f umfpack_save_symbolic.o umfpack_md_save_symbolic.o
	$(MX) -DZLONG -DCONJUGATE_SOLVE -c ../Source/umf_ltsolve.c
	$(MV) -f umf_ltsolve.o umf_mz_lhsolve.o
	$(MX) -DZLONG -DCONJUGATE_SOLVE -c ../Source/umf_utsolve.c
	$(MV) -f umf_utsolve.o umf_mz_uhsolve.o
	$(MX) -DZLONG -DDO_MAP -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_mz_triplet_map_nox.o
	$(MX) -DZLONG -DDO_VALUES -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_mz_triplet_nomap_x.o
	$(MX) -DZLONG -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_mz_triplet_nomap_nox.o
	$(MX) -DZLONG -DDO_MAP -DDO_VALUES -c ../Source/umf_triplet.c
	$(MV) -f umf_triplet.o umf_mz_triplet_map_x.o
	$(MX) -DZLONG -DFIXQ -c ../Source/umf_assemble.c
	$(MV) -f umf_assemble.o umf_mz_assemble_fixq.o
	$(MX) -DZLONG -DDROP -c ../Source/umf_store_lu.c
	$(MV) -f umf_store_lu.o umf_mz_store_lu_drop.o
	$(MX) -DZLONG -c ../Source/umf_assemble.c
	$(MV) -f umf_assemble.o umf_mz_assemble.o
	$(MX) -DZLONG -c ../Source/umf_blas3_update.c
	$(MV) -f umf_blas3_update.o umf_mz_blas3_update.o
	$(MX) -DZLONG -c ../Source/umf_build_tuples.c
	$(MV) -f umf_build_tuples.o umf_mz_build_tuples.o
	$(MX) -DZLONG -c ../Source/umf_create_element.c
	$(MV) -f umf_create_element.o umf_mz_create_element.o
	$(MX) -DZLONG -c ../Source/umf_dump.c
	$(MV) -f umf_dump.o umf_mz_dump.o
	$(MX) -DZLONG -c ../Source/umf_extend_front.c
	$(MV) -f umf_extend_front.o umf_mz_extend_front.o
	$(MX) -DZLONG -c ../Source/umf_garbage_collection.c
	$(MV) -f umf_garbage_collection.o umf_mz_garbage_collection.o
	$(MX) -DZLONG -c ../Source/umf_get_memory.c
	$(MV) -f umf_get_memory.o umf_mz_get_memory.o
	$(MX) -DZLONG -c ../Source/umf_init_front.c
	$(MV) -f umf_init_front.o umf_mz_init_front.o
	$(MX) -DZLONG -c ../Source/umf_kernel.c
	$(MV) -f umf_kernel.o umf_mz_kernel.o
	$(MX) -DZLONG -c ../Source/umf_kernel_init.c
	$(MV) -f umf_kernel_init.o umf_mz_kernel_init.o
	$(MX) -DZLONG -c ../Source/umf_kernel_wrapup.c
	$(MV) -f umf_kernel_wrapup.o umf_mz_kernel_wrapup.o
	$(MX) -DZLONG -c ../Source/umf_local_search.c
	$(MV) -f umf_local_search.o umf_mz_local_search.o
	$(MX) -DZLONG -c ../Source/umf_lsolve.c
	$(MV) -f umf_lsolve.o umf_mz_lsolve.o
	$(MX) -DZLONG -c ../Source/umf_ltsolve.c
	$(MV) -f umf_ltsolve.o umf_mz_ltsolve.o
	$(MX) -DZLONG -c ../Source/umf_mem_alloc_element.c
	$(MV) -f umf_mem_alloc_element.o umf_mz_mem_alloc_element.o
	$(MX) -DZLONG -c ../Source/umf_mem_alloc_head_block.c
	$(MV) -f umf_mem_alloc_head_block.o umf_mz_mem_alloc_head_block.o
	$(MX) -DZLONG -c ../Source/umf_mem_alloc_tail_block.c
	$(MV) -f umf_mem_alloc_tail_block.o umf_mz_mem_alloc_tail_block.o
	$(MX) -DZLONG -c ../Source/umf_mem_free_tail_block.c
	$(MV) -f umf_mem_free_tail_block.o umf_mz_mem_free_tail_block.o
	$(MX) -DZLONG -c ../Source/umf_mem_init_memoryspace.c
	$(MV) -f umf_mem_init_memoryspace.o umf_mz_mem_init_memoryspace.o
	$(MX) -DZLONG -c ../Source/umf_report_vector.c
	$(MV) -f umf_report_vector.o umf_mz_report_vector.o
	$(MX) -DZLONG -c ../Source/umf_row_search.c
	$(MV) -f umf_row_search.o umf_mz_row_search.o
	$(MX) -DZLONG -c ../Source/umf_scale_column.c
	$(MV) -f umf_scale_column.o umf_mz_scale_column.o
	$(MX) -DZLONG -c ../Source/umf_set_stats.c
	$(MV) -f umf_set_stats.o umf_mz_set_stats.o
	$(MX) -DZLONG -c ../Source/umf_solve.c
	$(MV) -f umf_solve.o umf_mz_solve.o
	$(MX) -DZLONG -c ../Source/umf_symbolic_usage.c
	$(MV) -f umf_symbolic_usage.o umf_mz_symbolic_usage.o
	$(MX) -DZLONG -c ../Source/umf_transpose.c
	$(MV) -f umf_transpose.o umf_mz_transpose.o
	$(MX) -DZLONG -c ../Source/umf_tuple_lengths.c
	$(MV) -f umf_tuple_lengths.o umf_mz_tuple_lengths.o
	$(MX) -DZLONG -c ../Source/umf_usolve.c
	$(MV) -f umf_usolve.o umf_mz_usolve.o
	$(MX) -DZLONG -c ../Source/umf_utsolve.c
	$(MV) -f umf_utsolve.o umf_mz_utsolve.o
	$(MX) -DZLONG -c ../Source/umf_valid_numeric.c
	$(MV) -f umf_valid_numeric.o umf_mz_valid_numeric.o
	$(MX) -DZLONG -c ../Source/umf_valid_symbolic.c
	$(MV) -f umf_valid_symbolic.o umf_mz_valid_symbolic.o
	$(MX) -DZLONG -c ../Source/umf_grow_front.c
	$(MV) -f umf_grow_front.o umf_mz_grow_front.o
	$(MX) -DZLONG -c ../Source/umf_start_front.c
	$(MV) -f umf_start_front.o umf_mz_start_front.o
	$(MX) -DZLONG -c ../Source/umf_2by2.c
	$(MV) -f umf_2by2.o umf_mz_2by2.o
	$(MX) -DZLONG -c ../Source/umf_store_lu.c
	$(MV) -f umf_store_lu.o umf_mz_store_lu.o
	$(MX) -DZLONG -c ../Source/umf_scale.c
	$(MV) -f umf_scale.o umf_mz_scale.o
	$(MX) -DZLONG -DWSOLVE -c ../Source/umfpack_solve.c
	$(MV) -f umfpack_solve.o umfpack_mz_wsolve.o
	$(MX) -DZLONG -c ../Source/umfpack_col_to_triplet.c
	$(MV) -f umfpack_col_to_triplet.o umfpack_mz_col_to_triplet.o
	$(MX) -DZLONG -c ../Source/umfpack_defaults.c
	$(MV) -f umfpack_defaults.o umfpack_mz_defaults.o
	$(MX) -DZLONG -c ../Source/umfpack_free_numeric.c
	$(MV) -f umfpack_free_numeric.o umfpack_mz_free_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_free_symbolic.c
	$(MV) -f umfpack_free_symbolic.o umfpack_mz_free_symbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_get_numeric.c
	$(MV) -f umfpack_get_numeric.o umfpack_mz_get_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_get_lunz.c
	$(MV) -f umfpack_get_lunz.o umfpack_mz_get_lunz.o
	$(MX) -DZLONG -c ../Source/umfpack_get_symbolic.c
	$(MV) -f umfpack_get_symbolic.o umfpack_mz_get_symbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_get_determinant.c
	$(MV) -f umfpack_get_determinant.o umfpack_mz_get_determinant.o
	$(MX) -DZLONG -c ../Source/umfpack_numeric.c
	$(MV) -f umfpack_numeric.o umfpack_mz_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_qsymbolic.c
	$(MV) -f umfpack_qsymbolic.o umfpack_mz_qsymbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_report_control.c
	$(MV) -f umfpack_report_control.o umfpack_mz_report_control.o
	$(MX) -DZLONG -c ../Source/umfpack_report_info.c
	$(MV) -f umfpack_report_info.o umfpack_mz_report_info.o
	$(MX) -DZLONG -c ../Source/umfpack_report_matrix.c
	$(MV) -f umfpack_report_matrix.o umfpack_mz_report_matrix.o
	$(MX) -DZLONG -c ../Source/umfpack_report_numeric.c
	$(MV) -f umfpack_report_numeric.o umfpack_mz_report_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_report_perm.c
	$(MV) -f umfpack_report_perm.o umfpack_mz_report_perm.o
	$(MX) -DZLONG -c ../Source/umfpack_report_status.c
	$(MV) -f umfpack_report_status.o umfpack_mz_report_status.o
	$(MX) -DZLONG -c ../Source/umfpack_report_symbolic.c
	$(MV) -f umfpack_report_symbolic.o umfpack_mz_report_symbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_report_triplet.c
	$(MV) -f umfpack_report_triplet.o umfpack_mz_report_triplet.o
	$(MX) -DZLONG -c ../Source/umfpack_report_vector.c
	$(MV) -f umfpack_report_vector.o umfpack_mz_report_vector.o
	$(MX) -DZLONG -c ../Source/umfpack_solve.c
	$(MV) -f umfpack_solve.o umfpack_mz_solve.o
	$(MX) -DZLONG -c ../Source/umfpack_symbolic.c
	$(MV) -f umfpack_symbolic.o umfpack_mz_symbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_transpose.c
	$(MV) -f umfpack_transpose.o umfpack_mz_transpose.o
	$(MX) -DZLONG -c ../Source/umfpack_triplet_to_col.c
	$(MV) -f umfpack_triplet_to_col.o umfpack_mz_triplet_to_col.o
	$(MX) -DZLONG -c ../Source/umfpack_scale.c
	$(MV) -f umfpack_scale.o umfpack_mz_scale.o
	$(MX) -DZLONG -c ../Source/umfpack_load_numeric.c
	$(MV) -f umfpack_load_numeric.o umfpack_mz_load_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_save_numeric.c
	$(MV) -f umfpack_save_numeric.o umfpack_mz_save_numeric.o
	$(MX) -DZLONG -c ../Source/umfpack_load_symbolic.c
	$(MV) -f umfpack_load_symbolic.o umfpack_mz_load_symbolic.o
	$(MX) -DZLONG -c ../Source/umfpack_save_symbolic.c
	$(MV) -f umfpack_save_symbolic.o umfpack_mz_save_symbolic.o
	$(MX) -c ../Source/umfpack_timer.c
	$(MV) -f umfpack_timer.o umfpack_m_timer.o
	$(MX) -c ../Source/umfpack_tictoc.c
	$(MV) -f umfpack_tictoc.o umfpack_m_tictoc.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_global.c
	$(MV) -f amd_global.o amd_m_global.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_aat.c
	$(MV) -f amd_aat.o amd_m_aat.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_1.c
	$(MV) -f amd_1.o amd_m_1.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_2.c
	$(MV) -f amd_2.o amd_m_2.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_dump.c
	$(MV) -f amd_dump.o amd_m_dump.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_postorder.c
	$(MV) -f amd_postorder.o amd_m_postorder.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_post_tree.c
	$(MV) -f amd_post_tree.o amd_m_post_tree.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_defaults.c
	$(MV) -f amd_defaults.o amd_m_defaults.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_order.c
	$(MV) -f amd_order.o amd_m_order.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_control.c
	$(MV) -f amd_control.o amd_m_control.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_info.c
	$(MV) -f amd_info.o amd_m_info.o
	$(MX) -DDLONG -c ../../AMD/Source/amd_valid.c
	$(MV) -f amd_valid.o amd_m_valid.o
	$(MX) -output umfpack2 umfpackmex.c \
	    umf_m_analyze.o umf_m_apply_order.o umf_m_colamd.o umf_m_free.o \
	    umf_m_fsize.o umf_m_is_permutation.o umf_m_malloc.o \
	    umf_m_realloc.o umf_m_report_perm.o umf_m_singletons.o \
	    umf_md_lhsolve.o umf_md_uhsolve.o umf_md_triplet_map_nox.o \
	    umf_md_triplet_nomap_x.o umf_md_triplet_nomap_nox.o \
	    umf_md_triplet_map_x.o umf_md_assemble_fixq.o \
	    umf_md_store_lu_drop.o umf_md_assemble.o umf_md_blas3_update.o \
	    umf_md_build_tuples.o umf_md_create_element.o umf_md_dump.o \
	    umf_md_extend_front.o umf_md_garbage_collection.o \
	    umf_md_get_memory.o umf_md_init_front.o umf_md_kernel.o \
	    umf_md_kernel_init.o umf_md_kernel_wrapup.o umf_md_local_search.o \
	    umf_md_lsolve.o umf_md_ltsolve.o umf_md_mem_alloc_element.o \
	    umf_md_mem_alloc_head_block.o umf_md_mem_alloc_tail_block.o \
	    umf_md_mem_free_tail_block.o umf_md_mem_init_memoryspace.o \
	    umf_md_report_vector.o umf_md_row_search.o umf_md_scale_column.o \
	    umf_md_set_stats.o umf_md_solve.o umf_md_symbolic_usage.o \
	    umf_md_transpose.o umf_md_tuple_lengths.o umf_md_usolve.o \
	    umf_md_utsolve.o umf_md_valid_numeric.o umf_md_valid_symbolic.o \
	    umf_md_grow_front.o umf_md_start_front.o umf_md_2by2.o \
	    umf_md_store_lu.o umf_md_scale.o umfpack_md_wsolve.o \
	    umfpack_md_col_to_triplet.o umfpack_md_defaults.o \
	    umfpack_md_free_numeric.o umfpack_md_free_symbolic.o \
	    umfpack_md_get_numeric.o umfpack_md_get_lunz.o \
	    umfpack_md_get_symbolic.o umfpack_md_get_determinant.o \
	    umfpack_md_numeric.o \
	    umfpack_md_qsymbolic.o umfpack_md_report_control.o \
	    umfpack_md_report_info.o umfpack_md_report_matrix.o \
	    umfpack_md_report_numeric.o umfpack_md_report_perm.o \
	    umfpack_md_report_status.o umfpack_md_report_symbolic.o \
	    umfpack_md_report_triplet.o umfpack_md_report_vector.o \
	    umfpack_md_solve.o umfpack_md_symbolic.o umfpack_md_transpose.o \
	    umfpack_md_triplet_to_col.o umfpack_md_scale.o \
	    umfpack_md_load_numeric.o umfpack_md_save_numeric.o \
	    umfpack_md_load_symbolic.o umfpack_md_save_symbolic.o \
	    umf_mz_lhsolve.o umf_mz_uhsolve.o umf_mz_triplet_map_nox.o \
	    umf_mz_triplet_nomap_x.o umf_mz_triplet_nomap_nox.o \
	    umf_mz_triplet_map_x.o umf_mz_assemble_fixq.o \
	    umf_mz_store_lu_drop.o umf_mz_assemble.o umf_mz_blas3_update.o \
	    umf_mz_build_tuples.o umf_mz_create_element.o umf_mz_dump.o \
	    umf_mz_extend_front.o umf_mz_garbage_collection.o \
	    umf_mz_get_memory.o umf_mz_init_front.o umf_mz_kernel.o \
	    umf_mz_kernel_init.o umf_mz_kernel_wrapup.o umf_mz_local_search.o \
	    umf_mz_lsolve.o umf_mz_ltsolve.o umf_mz_mem_alloc_element.o \
	    umf_mz_mem_alloc_head_block.o umf_mz_mem_alloc_tail_block.o \
	    umf_mz_mem_free_tail_block.o umf_mz_mem_init_memoryspace.o \
	    umf_mz_report_vector.o umf_mz_row_search.o umf_mz_scale_column.o \
	    umf_mz_set_stats.o umf_mz_solve.o umf_mz_symbolic_usage.o \
	    umf_mz_transpose.o umf_mz_tuple_lengths.o umf_mz_usolve.o \
	    umf_mz_utsolve.o umf_mz_valid_numeric.o umf_mz_valid_symbolic.o \
	    umf_mz_grow_front.o umf_mz_start_front.o umf_mz_2by2.o \
	    umf_mz_store_lu.o umf_mz_scale.o umfpack_mz_wsolve.o \
	    umfpack_mz_col_to_triplet.o umfpack_mz_defaults.o \
	    umfpack_mz_free_numeric.o umfpack_mz_free_symbolic.o \
	    umfpack_mz_get_numeric.o umfpack_mz_get_lunz.o \
	    umfpack_mz_get_symbolic.o umfpack_mz_get_determinant.o \
	    umfpack_mz_numeric.o \
	    umfpack_mz_qsymbolic.o umfpack_mz_report_control.o \
	    umfpack_mz_report_info.o umfpack_mz_report_matrix.o \
	    umfpack_mz_report_numeric.o umfpack_mz_report_perm.o \
	    umfpack_mz_report_status.o umfpack_mz_report_symbolic.o \
	    umfpack_mz_report_triplet.o umfpack_mz_report_vector.o \
	    umfpack_mz_solve.o umfpack_mz_symbolic.o umfpack_mz_transpose.o \
	    umfpack_mz_triplet_to_col.o umfpack_mz_scale.o \
	    umfpack_mz_load_numeric.o umfpack_mz_save_numeric.o \
	    umfpack_mz_load_symbolic.o umfpack_mz_save_symbolic.o \
	    umfpack_m_timer.o umfpack_m_tictoc.o umfpack_m_global.o \
	    amd_m_global.o \
	    amd_m_aat.o amd_m_1.o amd_m_2.o amd_m_dump.o \
	    amd_m_postorder.o amd_m_post_tree.o amd_m_defaults.o amd_m_order.o \
	    amd_m_control.o amd_m_info.o amd_m_valid.o

luflop: luflopmex.c
	$(MX) -output luflop luflopmex.c

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) *.mex*

clean:
	- $(RM) $(CLEAN)
