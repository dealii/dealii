############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries += $(lib-lac.g)  \
             $(lib-base.g)

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
tests_x =  geometry_info_* \
	   logtest \
	   logstream_end \
	   log_* \
	   reference \
	   quadrature_* qprojector \
	   path_search \
	   named_data \
	   slice_vector \
	   table \
	   tensor \
	   tensor_* \
	   symmetric_tensor_* \
	   full_tensor_* \
           timer \
	   polynomial* \
	   auto_derivative_function \
	   anisotropic_* \
           hierarchical \
	   data_out_base \
	   data_out_base_dx \
	   data_out_base_eps \
	   data_out_base_gmv \
	   data_out_base_gnuplot \
	   data_out_base_povray \
	   data_out_base_tecplot \
	   data_out_base_vtk \
	   function_* \
           functions_* \
	   utilities_* \
	   work_stream_* \
	   is_finite* \
	   number_traits_* \
	   thread_validity_* \
	   task_* \
	   parallel_* \
	   mutex_* \
	   scalar_*

# add threading tests. note that we have
# to list them individually, without wildcards, because the .cc files are
# generated and don't exist yet (so wildcard expansion will fail)
tests_x += threads_01 threads_02


ifeq (tecio, $(findstring tecio, $LIBS))
  tests_x += data_out_base_tecplot_bin
endif

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
