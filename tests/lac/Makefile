############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-lac.g)  \
            $(lib-base.g)
libraries_o = $(lib-lac.o)  \
            $(lib-base.o)

default: run-tests 

############################################################

block_matrices.exe        : block_matrices.g.$(OBJEXT)                  $(libraries)
block_vector.exe          : block_vector.g.$(OBJEXT) abort.g.$(OBJEXT)  $(libraries)
block_vector_iterator.exe : block_vector_iterator.g.$(OBJEXT)           $(libraries)
full_matrix.exe           : full_matrix.g.$(OBJEXT)                     $(libraries)
matrix_out.exe            : matrix_out.g.$(OBJEXT)                      $(libraries)
solver.exe                : solver.g.$(OBJEXT)   testmatrix.g.$(OBJEXT) $(libraries)
sparse_matrices.exe       : sparse_matrices.g.$(OBJEXT) testmatrix.g.$(OBJEXT) $(libraries)
sparse_matrices.opt.exe   : sparse_matrices.$(OBJEXT) testmatrix.$(OBJEXT) $(libraries_o)
sparsity_pattern.exe      : sparsity_pattern.g.$(OBJEXT) testmatrix.g.$(OBJEXT) $(libraries)
sparse_ilu.exe            : sparse_ilu.g.$(OBJEXT) testmatrix.g.$(OBJEXT) $(libraries)
vector-vector.exe         : vector-vector.g.$(OBJEXT)                   $(libraries)
#mgbase.exe               : mgbase.g.$(OBJEXT)                          $(libraries)
#mg.exe                   : mg.g.$(OBJEXT) testmatrix.g.$(OBJEXT)       $(libraries)


tests = $(sort \
          sparse_matrices block_matrices block_vector block_vector_iterator \
          full_matrix \
          matrix_out solver sparsity_pattern sparse_ilu vector-vector)

############################################################


include ../Makefile.rules
include Makefile.depend
