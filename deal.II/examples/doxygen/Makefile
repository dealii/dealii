###########################################################################
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 1998 - 2007, 2010 by the deal.II authors
#
#    This file is subject to QPL and may not be  distributed
#    without copyright and license information. Please refer
#    to the file deal.II/doc/license.html for the  text  and
#    further information on this license.
#
###########################################################################
include ../../common/Make.global_options

SHELL=/bin/bash

tests = $(wildcard *.cc)
exe =  $(tests:%.cc=%$(EXEEXT))
ok = $(tests:%.cc=%.ok)

all: $(exe)

run: $(ok)
	for prog in $(exe) ; do ./$$prog ; done

######################################################################
# Compilation of source code
######################################################################

./%.g.$(OBJEXT) :
	@echo "==============debug========= $(<F)  ->  $@"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
./%.$(OBJEXT) :
	@echo ==============optimized===== $(<F)
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@

######################################################################
# Running the programs
######################################################################
%.ok: %$(EXEEXT)
	./$<

######################################################################
# Specific targets
######################################################################

block_matrix_array$(EXEEXT): block_matrix_array.g.$(OBJEXT) $(lib-deal2.g)
	@echo ============================ Linking $@
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

product_matrix$(EXEEXT): product_matrix.g.$(OBJEXT) $(lib-deal2.g)
	@echo ============================ Linking $@
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

compressed_block_sparsity_pattern$(EXEEXT): compressed_block_sparsity_pattern.g.$(OBJEXT) $(lib-deal2.g)
	@echo ============================ Linking $@
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

theta_timestepping$(EXEEXT): theta_timestepping.g.$(OBJEXT) $(lib-deal2.g)
	@echo ============================ Linking $@
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

######################################################################
# Pseudo target for cleaning directory
######################################################################


clean:
	-rm -f *~ *.$(OBJEXT) $(exe)


.PHONY: all clean

# Finally there is a rule which you normally need not care much about:
# since the executable depends on some include files from the library,
# besides the C++ application file of course, it is necessary to
# re-generate the executable when one of the files it depends on has
# changed. The following rule creates a dependency file
# `Makefile.dep', which `make' uses to determine when to regenerate
# the executable. This file is automagically remade whenever needed,
# i.e. whenever one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom of this file.
#
# If the creation of Makefile.dep fails, blow it away and fail
Makefile.dep: *.cc Makefile \
              $(shell echo $D/include/deal.II/*/*.h)
	@echo ============================ Remaking $@
	@$D/common/scripts/make_dependencies  $(INCLUDE) -B. *.cc \
		> $@ \
	  || (rm -f $@ ; false)
	@if test -s $@ ; then : else rm $@ ; fi



# To make the dependencies known to `make', we finally have to include
# them:
include Makefile.dep
