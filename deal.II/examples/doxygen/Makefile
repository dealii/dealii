###########################################################################
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 1998 - 2005 by the deal.II authors
#
#    This file is subject to QPL and may not be  distributed
#    without copyright and license information. Please refer
#    to the file deal.II/doc/license.html for the  text  and
#    further information on this license.
#
###########################################################################
include ../../common/Make.global_options

all: product_matrix.exe block_matrix_array.exe

######################################################################
# Compilation of source code
######################################################################

%.o : %.cc
	@echo =====debug======$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(CXXFLAGS) -c $< -o $@

######################################################################
# Specific targets
######################################################################

block_matrix_array.exe: block_matrix_array.o $(lib-base.g) $(lib-lac.g)
	@echo =====linking====$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(LDFLAGS) -o $@ $^ $(LIBS)

product_matrix.exe: product_matrix.o $(lib-base.g) $(lib-lac.g)
	@echo =====linking====$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(LDFLAGS) -o $@ $^ $(LIBS)

######################################################################
# Pseudo target for cleaning directory
######################################################################

clean:
	-rm *~ *.o *.exe
