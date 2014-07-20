#-------------------------------------------------------------------------------
# UMFPACK GNUmakefile for the UMFPACK MATLAB mexFunction (GNU "make" only)
#-------------------------------------------------------------------------------

default: umfpack2 luflop

include ../../UFconfig/UFconfig.mk

I =  -I../Include -I../Source -I../../AMD/Include -I../../UFconfig

MX = $(MEX) $(I)
    

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

# non-user-callable umf_*.[ch] files:
UMFCH = umf_assemble umf_blas3_update \
	umf_build_tuples umf_create_element \
	umf_dump umf_extend_front umf_garbage_collection \
	umf_get_memory umf_init_front umf_kernel \
	umf_kernel_init umf_kernel_wrapup \
	umf_local_search umf_lsolve umf_ltsolve \
	umf_mem_alloc_element umf_mem_alloc_head_block \
	umf_mem_alloc_tail_block umf_mem_free_tail_block \
	umf_mem_init_memoryspace \
	umf_report_vector umf_row_search umf_scale_column \
	umf_set_stats umf_solve umf_symbolic_usage umf_transpose \
	umf_tuple_lengths umf_usolve umf_utsolve umf_valid_numeric \
	umf_valid_symbolic umf_grow_front umf_start_front umf_2by2 \
	umf_store_lu umf_scale

# non-user-callable umf_*.[ch] files, int/UF_long versions only (no real/complex):
UMFINT = umf_analyze umf_apply_order umf_colamd umf_free umf_fsize \
	umf_is_permutation umf_malloc umf_realloc umf_report_perm \
	umf_singletons

# non-user-callable and user-callable amd_*.[ch] files (int/UF_long versions only):
AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess amd_global

# non-user-callable, created from umf_ltsolve.c, umf_utsolve.c,
# umf_triplet.c, and umf_assemble.c , with int/UF_long and real/complex versions:
UMF_CREATED = umf_lhsolve umf_uhsolve umf_triplet_map_nox \
	umf_triplet_nomap_x umf_triplet_nomap_nox umf_triplet_map_x \
	umf_assemble_fixq umf_store_lu_drop

# non-user-callable, int/UF_long and real/complex versions:
UMF = $(UMF_CREATED) $(UMFCH)

# user-callable umfpack_*.[ch] files (int/UF_long and real/complex):
UMFPACK = umfpack_col_to_triplet umfpack_defaults umfpack_free_numeric \
	umfpack_free_symbolic umfpack_get_numeric umfpack_get_lunz \
	umfpack_get_symbolic umfpack_get_determinant umfpack_numeric \
	umfpack_qsymbolic umfpack_report_control umfpack_report_info \
	umfpack_report_matrix umfpack_report_numeric umfpack_report_perm \
	umfpack_report_status umfpack_report_symbolic umfpack_report_triplet \
	umfpack_report_vector umfpack_solve umfpack_symbolic \
	umfpack_transpose umfpack_triplet_to_col umfpack_scale \
	umfpack_load_numeric umfpack_save_numeric \
	umfpack_load_symbolic umfpack_save_symbolic

# user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
# with int/UF_long and real/complex versions:
UMFPACKW = umfpack_wsolve

USER = $(UMFPACKW) $(UMFPACK)

# user-callable, only one version for int/UF_long, real/complex, *.[ch] files:
GENERIC = umfpack_timer umfpack_tictoc umfpack_global

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

AMDH = ../../AMD/Include/amd_internal.h ../../AMD/Include/amd.h \
    ../../UFconfig/UFconfig.h

INC1 = umf_config.h umf_version.h umf_internal.h umf_triplet.h

INC = ../Include/umfpack.h \
	$(addprefix ../Source/, $(INC1)) \
	$(addprefix ../Source/, $(addsuffix .h,$(UMFCH))) \
	$(addprefix ../Source/, $(addsuffix .h,$(UMFINT))) \
	$(addprefix ../Include/, $(addsuffix .h,$(USER))) \
	$(addprefix ../Include/, $(addsuffix .h,$(GENERIC))) \
	$(AMDH)

#-------------------------------------------------------------------------------
# Create the umfpack and amd mexFunctions for MATLAB (int versions only)
#-------------------------------------------------------------------------------

MEXI   = $(addsuffix .o, $(subst umf_,umf_m_,$(UMFINT)))
MEXDI  = $(addsuffix .o, $(subst umf_,umf_md_,$(UMF)) $(subst umfpack_,umfpack_md_,$(USER)))
MEXZI  = $(addsuffix .o, $(subst umf_,umf_mz_,$(UMF)) $(subst umfpack_,umfpack_mz_,$(USER)) )
MEXAMD = $(addsuffix .o, $(subst amd_,amd_m_,$(AMD)))
MEXGN  = $(addsuffix .o, $(subst umfpack_,umfpack_m_,$(GENERIC)))

MEXUMFPACK = $(MEXI) $(MEXDI) $(MEXZI) $(MEXGN)

# Note that mex has no "-o" option, thus the need for $(MV) commands.
# If it did, then the rules would be much simpler:
#	$(MX) -DDLONG -c $< -o $@

#----------------------------------------
# integer-only routines (no real/complex):
#----------------------------------------

amd_m_%.o: ../../AMD/Source/amd_%.c $(AMDH)
	$(MX) -DDLONG -c $<
	- $(MV) amd_$*.o $@

umf_m_%.o: ../Source/umf_%.c $(INC)
	$(MX) -DDLONG -c $<
	- $(MV) umf_$*.o $@

#----------------------------------------
# Double precision, int version, for MATLAB
#----------------------------------------

umf_md_%.o: ../Source/umf_%.c $(INC)
	$(MX) -DDLONG -c $<
	- $(MV) umf_$*.o $@

umf_md_%hsolve.o: ../Source/umf_%tsolve.c $(INC)
	$(MX) -DDLONG -DCONJUGATE_SOLVE -c $<
	- $(MV) umf_$*tsolve.o $@

umf_md_triplet_map_x.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DDLONG -DDO_MAP -DDO_VALUES -c $<
	- $(MV) umf_triplet.o $@

umf_md_triplet_map_nox.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DDLONG -DDO_MAP -c $<
	- $(MV) umf_triplet.o $@

umf_md_triplet_nomap_x.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DDLONG -DDO_VALUES -c $<
	- $(MV) umf_triplet.o $@

umf_md_triplet_nomap_nox.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DDLONG -c $<
	- $(MV) umf_triplet.o $@

umf_md_assemble_fixq.o: ../Source/umf_assemble.c $(INC)
	$(MX) -DDLONG -DFIXQ -c $<
	- $(MV) umf_assemble.o $@

umf_md_store_lu_drop.o: ../Source/umf_store_lu.c $(INC)
	$(MX) -DDLONG -DDROP -c $<
	- $(MV) umf_store_lu.o $@

umfpack_md_wsolve.o: ../Source/umfpack_solve.c $(INC)
	$(MX) -DDLONG -DWSOLVE -c $<
	- $(MV) umfpack_solve.o $@

umfpack_md_%.o: ../Source/umfpack_%.c $(INC)
	$(MX) -DDLONG -c $<
	- $(MV) umfpack_$*.o $@

#----------------------------------------
# Complex double precision, int version, for MATLAB
#----------------------------------------

umf_mz_%.o: ../Source/umf_%.c $(INC)
	$(MX) -DZLONG -c $<
	- $(MV) umf_$*.o $@

umf_mz_%hsolve.o: ../Source/umf_%tsolve.c $(INC)
	$(MX) -DZLONG -DCONJUGATE_SOLVE -c $<
	- $(MV) umf_$*tsolve.o $@

umf_mz_triplet_map_x.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DZLONG -DDO_MAP -DDO_VALUES -c $<
	- $(MV) umf_triplet.o $@

umf_mz_triplet_map_nox.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DZLONG -DDO_MAP -c $<
	- $(MV) umf_triplet.o $@

umf_mz_triplet_nomap_x.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DZLONG -DDO_VALUES -c $<
	- $(MV) umf_triplet.o $@

umf_mz_triplet_nomap_nox.o: ../Source/umf_triplet.c $(INC)
	$(MX) -DZLONG -c $<
	- $(MV) umf_triplet.o $@

umf_mz_assemble_fixq.o: ../Source/umf_assemble.c $(INC)
	$(MX) -DZLONG -DFIXQ -c $<
	- $(MV) umf_assemble.o $@

umf_mz_store_lu_drop.o: ../Source/umf_store_lu.c $(INC)
	$(MX) -DZLONG -DDROP -c $<
	- $(MV) umf_store_lu.o $@

umfpack_mz_wsolve.o: ../Source/umfpack_solve.c $(INC)
	$(MX) -DZLONG -DWSOLVE -c $<
	- $(MV) umfpack_solve.o $@

umfpack_mz_%.o: ../Source/umfpack_%.c $(INC)
	$(MX) -DZLONG -c $<
	- $(MV) umfpack_$*.o $@

#----------------------------------------
# Generic routines for MATLAB
#----------------------------------------

umfpack_m_timer.o: ../Source/umfpack_timer.c $(INC)
	$(MX) -c $<
	- $(MV) umfpack_timer.o $@

umfpack_m_tictoc.o: ../Source/umfpack_tictoc.c $(INC)
	$(MX) -c $<
	- $(MV) umfpack_tictoc.o $@

umfpack_m_global.o: ../Source/umfpack_global.c $(INC)
	$(MX) -c $<
	- $(MV) umfpack_global.o $@

#----------------------------------------
# umfpack mexFunction
#----------------------------------------

umfpack2: umfpackmex.c $(MEXUMFPACK) $(MEXAMD)
	$(MX) -output umfpack2 umfpackmex.c $(MEXUMFPACK) $(MEXAMD)

luflop: luflopmex.c
	$(MX) -output luflop luflopmex.c

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) *.mex*

clean:
	- $(RM) $(CLEAN)
