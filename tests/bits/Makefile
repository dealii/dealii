# Generated automatically from Makefile.in by configure.
############################################################
# Makefile,v 1.15 2002/06/13 12:51:13 hartmann Exp
# Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-deal2-1d.g) \
            $(lib-deal2-2d.g) \
            $(lib-deal2-3d.g) \
            $(lib-lac.g)  \
            $(lib-base.g)

default: run-tests

############################################################

tests_x = anna_? \
	  gerold_1 \
	  geometry_info_1 \
	  point_inside_? \
	  full_matrix_1 \
	  data_out_0* \
          data_out_faces_0* \
	  data_out_rotation_0* \
	  data_out_stack_0* \
	  dof_tools_[0-9]* \
	  fe_tools_[0-9]* \
	  roy_1 \
          denis_1 \
	  unit_support_points \
	  parameter_handler_* \
	  sparse_lu_decomposition_1 \
	  block_sparse_matrix_* \
	  block_vector_* \
	  hyper_ball_3d cylinder coarsening_3d \
	  mesh_3d_* \
	  normals_* \
	  q_point_sum_* \
	  volume_* \
	  mapping_cartesian_1 \
	  mapping_q4_3d \
	  q_points \
          find_cell_* \
          sparsity_pattern_* \
          sparse_matrix_* \
	  full_matrix_vector_* \
	  solver_* \
	  deal_solver_* \
	  vector_* \
	  subdomain_ids_* \
	  dof_constraints_* \
	  apply_boundary_values_* \
	  error_estimator_* \
	  cylinder_shell_* \
	  grid_in_* \
	  compressed_sparsity_pattern_* \
	  point_difference_* \
	  umfpack_*

# tests for the hp branch:
#	  fe_collection_*

ifeq ($(USE_CONTRIB_PETSC),yes)
  tests_x += petsc_*
endif

ifeq ($(USE_CONTRIB_METIS),yes)
  tests_x += metis_*
endif

ifeq ($(USE_CONTRIB_HSL),yes)
  tests_x += hsl_*
endif



# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests

.PHONY: default