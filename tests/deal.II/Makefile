# Generated automatically from Makefile.in by configure.
############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-lac.g)  \
            $(lib-base.g)
lib-1d    = $(lib-deal2-1d.g)
lib-2d    = $(lib-deal2-2d.g)
lib-3d    = $(lib-deal2-3d.g)

default: run-tests 

############################################################

block_matrices.exe           : block_matrices.g.$(OBJEXT)              $(lib-2d)           $(libraries)
constraints.exe              : constraints.g.$(OBJEXT)                 $(lib-2d) $(lib-3d) $(libraries) 
data_out.exe                 : data_out.g.$(OBJEXT)                    $(lib-2d) $(lib-3d) $(libraries)
derivative_approximation.exe : derivative_approximation.g.$(OBJEXT) $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
derivatives.exe              : derivatives.g.$(OBJEXT)       $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
dof_test.exe                 : dof_test.g.$(OBJEXT)                    $(lib-2d) $(lib-3d) $(libraries)
dof_renumbering.exe          : dof_renumbering.g.$(OBJEXT)   $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
error_estimator.exe          : error_estimator.g.$(OBJEXT)   $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
gradients.exe                : gradients.g.$(OBJEXT)                   $(lib-2d)           $(libraries)
grid_out.exe                 : grid_out.g.$(OBJEXT)          $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
grid_test.exe                : grid_test.g.$(OBJEXT)                   $(lib-2d) $(lib-3d) $(libraries)
grid_transform.exe           : grid_transform.g.$(OBJEXT)              $(lib-2d) $(lib-3d) $(libraries)
intergrid_constraints.exe    : intergrid_constraints.g.$(OBJEXT)    $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
intergrid_map.exe            : intergrid_map.g.$(OBJEXT)     $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
matrices.exe                 : matrices.g.$(OBJEXT)          $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
mg.exe                       : mg.g.$(OBJEXT)                          $(lib-2d)           $(libraries)
mglocal.exe                  : mglocal.g.$(OBJEXT)                     $(lib-2d)           $(libraries)
second_derivatives.exe       : second_derivatives.g.$(OBJEXT)          $(lib-2d)           $(libraries)
wave-test-3.exe              : wave-test-3.g.$(OBJEXT)                 $(lib-2d)           $(libraries)
support_point_map.exe        : support_point_map.g.$(OBJEXT) $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
filtered_matrix.exe          : filtered_matrix.g.$(OBJEXT)   $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
boundaries.exe               : boundaries.g.$(OBJEXT)        $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
sparsity_pattern.exe         : sparsity_pattern.g.$(OBJEXT)  $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
grid_tools.exe               : grid_tools.g.$(OBJEXT)        $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
subdomain_ids.exe            : subdomain_ids.g.$(OBJEXT)     $(lib-1d) $(lib-2d) $(lib-3d) $(libraries)
filtered_iterator.exe        : filtered_iterator.g.$(OBJEXT)           $(lib-2d)           $(libraries)
grid_in.exe                  : grid_in.g.$(OBJEXT)                     $(lib-2d)           $(libraries)
grid_in_3d.exe               : grid_in_3d.g.$(OBJEXT)                            $(lib-3d) $(libraries)


tests = $(sort \
          grid_test grid_transform dof_test data_out derivatives gradients \
	  constraints block_matrices second_derivatives derivative_approximation \
          matrices error_estimator intergrid_constraints intergrid_map \
          wave-test-3 dof_renumbering support_point_map filtered_matrix \
	  boundaries sparsity_pattern grid_tools subdomain_ids \
	  filtered_iterator grid_out grid_in grid_in_3d)

############################################################


include ../Makefile.rules
include Makefile.depend
