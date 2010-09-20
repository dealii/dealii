############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options


default: run-tests

############################################################
# rules to generate the threading tests. the .cc files are
# automatically generated; note that for the output generation, we
# have to sort the output since otherwise in MT mode things would come
# out in a random order. the same holds for a few other tests
threads_01.cc : threads.pl
	@echo =====regenerating== $@
	@$(PERL) $^ > $@

threads_02.cc : threads_01.cc
	@echo =====regenerating== $@
	@cat $^ \
	   | $(PERL) -pi -e 's/spawn \((.*), (.*)\) \(/new_thread (\2, \1, /g;' \
	   | $(PERL) -pi -e 's/(new_thread.*), \)/\1)/g;' \
	   | $(PERL) -pi -e 's/threads_01/threads_02/g;' \
	   > $@

threads_%/output : threads_%/exe
	@echo =====Running======= $<
	@echo Running > $(dir $@)/status
	@$(ULIMIT) -t 2400 ; ./$< ; \
	  if test ! $$? = 0 ; then rm $@ ; false ; fi
	@perl -pi $(normalize)  $@
	@sort $@ -o $@

thread_validity_07/output : thread_validity_07/exe
	@echo =====Running======= $<
	@echo Running > $(dir $@)/status
	@$(ULIMIT) -t 2400 ; ./$< ; \
	  if test ! $$? = 0 ; then rm $@ ; false ; fi
	@perl -pi $(normalize)  $@
	@sort $@ -o $@

task_%/output : task_%/exe
	@echo =====Running======= $<
	@echo Running > $(dir $@)/status
	@$(ULIMIT) -t 2400 ; ./$< ; \
	  if test ! $$? = 0 ; then rm $@ ; false ; fi
	@perl -pi $(normalize)  $@
	@sort $@ -o $@


############################################################

# all .cc-files are tests, though the test bdm.cc appears to be dysfunctional
# and we don't want the data_out_base_tecplot_bin test if we haven't found the
# corresponding libraries:
ifneq (tecio, $(findstring tecio, $LIBS))
  tests = $(filter-out bdm data_out_base_tecplot_bin, $(basename $(wildcard *)))
else
  tests = $(filter-out bdm, $(basename $(wildcard *)))
endif

# add threading tests. note that we have
# to list them individually, without wildcards, because the .cc files are
# generated and don't exist yet (so wildcard expansion will fail)
tests += threads_01 threads_02


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
