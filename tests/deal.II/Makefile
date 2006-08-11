############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-lac.g)  \
            $(lib-base.g) \
            $(lib-deal2-1d.g) \
            $(lib-deal2-2d.g) \
            $(lib-deal2-3d.g)

default: run-tests 

############################################################

tests_x = block_matrices \
	constraints \
	data_out \
	derivative_* \
	derivatives \
	dof_test \
	dof_renumbering \
	error_estimator \
	grid_out \
	grid_test \
	grid_transform \
	grid_transform_3d \
	measure_et_al \
	normals_at_vertices \
	intergrid_constraints \
	intergrid_map \
	matrices \
	wave-test-3 \
	support_point_map \
	filtered_matrix \
	boundaries \
	sparsity_pattern \
	grid_tools \
	subdomain_ids \
	filtered_iterator \
	grid_in* \
	grid_in_3d \
	face_orientations_3d \
	mg_dof_handler \
	subcelldata \
	interpolate_* \
	have_same_coarse_mesh_* \
	get_finest_common_cells_* \
	project_*[0-9]

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
