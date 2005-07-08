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

all: product_matrix$(EXEEXT) block_matrix_array$(EXEEXT)

######################################################################
# Compilation of source code
######################################################################

%.$(OBJEXT) : %.cc
	@echo =====debug======$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(CXXFLAGS) -c $< -o $@

######################################################################
# Specific targets
######################################################################

block_matrix_array$(EXEEXT): block_matrix_array.$(OBJEXT) $(lib-lac.g) $(lib-base.g)
	@echo =====linking====$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(LDFLAGS) -o $@ $^ $(LIBS)

product_matrix$(EXEEXT): product_matrix.$(OBJEXT) $(lib-lac.g) $(lib-base.g)
	@echo =====linking====$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(LDFLAGS) -o $@ $^ $(LIBS)

######################################################################
# Pseudo target for cleaning directory
######################################################################

clean:
	-rm -f *~ *.$(OBJEXT) $(all)
