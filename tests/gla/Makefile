############################################################
# $Id: Makefile 17098 2008-10-03 19:27:13Z bangerth $
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options

default: run-tests

############################################################

# all .cc-files are tests by default.

ifeq ($(USE_CONTRIB_TRILINOS),no)
  tests_x_1 = $(shell for i in *.cc ; do \
			if grep -q TrilinosWrappers $$i ; then \
				 : ; \
			else \
				 echo $$i ; \
			fi ; done)
else
  tests_x_1 = $(wildcard *.cc)
endif

tests   = $(shell echo $(addsuffix /ncpu_*, $(basename $(tests_x_1))))


############################################################


OVERRIDE_OUTPUT_RULE = yes
include ../Makefile.rules
include Makefile.depend


Makefile.mpi.tests: Makefile $(shell echo *.cc) $(shell echo */*/cmp/generic)
	@echo =====Targets======= $@
	@(for i in $(basename $(wildcard *.cc)) ; do \
		echo "$$i/exe : $$i/obj.g.\$$(OBJEXT) \$$(libraries)"; \
	 done ; \
	 for i in $(basename $(wildcard *.cc)) ; do \
		for j in $$i/ncpu_*; do \
			echo "$$j/output : $$i/exe" ; \
		done ; \
	 done) \
		> $@

include Makefile.mpi.tests

# override the rule from ../Makefile.rules by running things with
# mpirun -np X
get_ncpus = $(shell echo $(1) | $(PERL) -pi -e 's/.*ncpu_(\d+).*/\1/g;')
get_exe = $(shell echo $(1) | $(PERL) -pi -e 's/ncpu_.*/exe/g;')

%/output:
	@echo =====Running==MPI== $@
	@echo Running > $(dir $@)/status
	@$(ULIMIT) -t 2400 ; \
		mpirun -np $(call get_ncpus, $@) ./$(call get_exe, $@) ; \
		if test ! $$? = 0 ; then rm $@ ; false ; fi
	@perl -pi $(normalize)  $@

