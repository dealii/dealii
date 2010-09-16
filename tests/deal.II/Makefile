############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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

tests_x = block_info block_matrices \
	user_data_* \
	constraint* \
	inhomogeneous_constraints \
	inhomogeneous_constraints_nonsymmetric \
	inhomogeneous_constraints_block \
	inhomogeneous_constraints_vector \
	data_out \
	derivative_* \
	derivatives \
	dof_test \
	dof_renumbering \
	dof_renumbering_* \
	error_estimator \
	grid_invert \
	grid_generator* \
	grid_out \
	grid_out_* \
	grid_output_input \
	grid_test \
	grid_transform \
	grid_transform_3d \
	grid_hyper_* \
	measure_et_al \
	extent_in_direction \
	normals_at_vertices \
	intergrid_constraints \
	intergrid_map \
	matrices \
	vectors_* \
	support_point_map \
	filtered_matrix \
	boundaries \
	sparsity_pattern* \
	grid_tools \
	subdomain_* \
	filtered_iterator* \
	grid_in* \
	grid_in_3d \
	face_orientations_3d \
	mg_dof_handler \
	subcelldata \
	interpolate_* \
	have_same_coarse_mesh_* \
	get_finest_common_cells_* \
	project_*[0-9] \
	project_boundary_* \
	project_to_surface_* \
	minimal_cell_diameter \
	maximal_cell_diameter \
	union_triangulation \
	iterators* \
	mesh_smoothing_* \
	create_* \
	mesh_worker_* \
	line_coarsening_3d \
	no_flux_* \
	fe_values_view_* \
	fe_values_extractor_* \
	rt_hessian \
	rt_covariant \
	coarsening_* \
	number_cache_* \
	kelly_crash_* \
	recursively_set_material_id \
	point_value_history_* \
	newton_*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
