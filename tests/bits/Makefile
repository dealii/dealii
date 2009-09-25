############################################################
# Makefile,v 1.15 2002/06/13 12:51:13 hartmann Exp
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
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

tests_x = grid_generator_?? \
	  grid_in_* \
	  point_inside_? \
	  full_matrix_* \
	  complex_real_full_matrix_* \
	  complex_complex_full_matrix_* \
	  block_matrix_array_* \
	  data_out_0* \
          data_out_faces_0* \
	  data_out_rotation_0* \
	  data_out_stack_0* \
	  data_out_curved_cells \
	  dof_tools_[0-9]* \
	  fe_tools_[0-9]* \
	  fe_tools_cpfqpm* \
	  fe_face_to_cell_index \
	  fe_field_function \
          anna_? \
	  gerold_? \
	  roy_1 \
          denis_1 \
	  unit_support_points \
	  parameter_handler_* \
	  parsed_function \
	  sparse_lu_decomposition_1 \
	  block_sparse_matrix_* \
	  block_vector_* \
	  hyper_ball_3d \
	  cylinder* \
	  cone* \
	  coarsening_3d \
	  mesh_3d_* \
	  normals_* \
	  q_point_sum_* \
	  volume_* \
	  mapping_cartesian_1 \
	  mapping_q4_3d \
	  q_points \
	  find_closest_vertex_* \
	  find_cells_adjacent_to_vertex_* \
          find_cell_* \
          sparsity_pattern_* \
          sparse_matrix_* \
          chunk_sparse_matrix_* \
	  full_matrix_vector_* \
	  solver_* \
	  deal_solver_* \
	  vector_* \
	  complex_vector_* \
	  subdomain_ids_* \
	  dof_constraints_* \
	  apply_boundary_values_* \
	  error_estimator_* \
	  cylinder_shell_* \
	  compressed_sparsity_pattern_* \
	  compressed_set_sparsity_pattern_* \
	  point_difference_* \
	  point_value_* \
	  fe_q_3d* \
          rt_* \
	  nedelec_* \
	  get_fe_from_name \
	  christian_* \
	  oliver_* \
	  accessor_equality* \
	  get_dof_indices_* \
	  n_dofs_per_object \
	  step-* \
	  hp_constraints* \
	  face_interpolation \
	  subface_interpolation \
	  face_orientation_crash \
	  refine_and_coarsen_* \
	  tria_crash_* \
	  face_orientation_* \
	  max_n_cells_* \
	  joa_* \
	  neighboring_cells_at_two_faces \
	  solution_transfer \
          make_boundary_constraints_* \
	  count_dofs_per_block* \
	  distorted_cells_*

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
