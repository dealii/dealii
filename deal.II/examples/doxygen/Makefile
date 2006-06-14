###########################################################################
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 1998 - 2006 by the deal.II authors
#
#    This file is subject to QPL and may not be  distributed
#    without copyright and license information. Please refer
#    to the file deal.II/doc/license.html for the  text  and
#    further information on this license.
#
###########################################################################
include ../../common/Make.global_options

all: product_matrix$(EXEEXT) block_matrix_array$(EXEEXT) \
     compressed_block_sparsity_pattern$(EXEEXT)

######################################################################
# Compilation of source code
######################################################################

%.$(OBJEXT) :
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

compressed_block_sparsity_pattern$(EXEEXT): compressed_block_sparsity_pattern.$(OBJEXT) $(lib-deal2-2d.g) $(lib-lac.g) $(lib-base.g)
	@echo =====linking====$(MT)== $<
	@$(CXX) $(CXXFLAGS.g) $(LDFLAGS) -o $@ $^ $(LIBS)

######################################################################
# Pseudo target for cleaning directory
######################################################################


clean:
	-rm -f *~ *.$(OBJEXT) $(all)


.PHONY: all clean

# Finally there is a rule which you normally need not care much about:
# since the executable depends on some include files from the library,
# besides the C++ application file of course, it is necessary to
# re-generate the executable when one of the files it depends on has
# changed. The following rule to created a dependency file
# `Makefile.dep', which `make' uses to determine when to regenerate
# the executable. This file is automagically remade whenever needed,
# i.e. whenever one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom of this file.
#
# If the creation of Makefile.dep fails, blow it away and fail
Makefile.dep: *.cc Makefile \
              $(shell echo $D/*/include/*/*.h)
	@echo ============================ Remaking $@
	@$D/common/scripts/make_dependencies  $(INCLUDE) -B. *.cc \
		> $@ \
	  || (rm -f $@ ; false)
	@if test -s $@ ; then : else rm $@ ; fi



# To make the dependencies known to `make', we finally have to include
# them:
include Makefile.dep
