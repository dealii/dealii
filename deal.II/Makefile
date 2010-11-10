# $Id$
#
# Makefile for the whole library

SHELL=/bin/bash

help:
	@echo "========================================================================"
	@echo "=              Global Makefile for the deal.II libraries               ="
	@echo "========================================================================"
	@echo "=                                                                      ="
	@echo "= The following targets exist:                                         ="
	@echo "=    all        : debug and optimized libraries                        ="
	@echo "=    debug      : the debug version of the deal.II library             ="
	@echo "=    optimized  : the optimized version of the deal.II library         ="
	@echo "=                                                                      ="
	@echo "=    contrib    : additional libraries in contrib, if there are any    ="
	@echo "=                                                                      ="
	@echo "=    online-doc : generate the documentation in HTML format            ="
	@echo "=    TODO       : create a \"TODO\" file from the source files           ="
	@echo "=    TAGS       : create a TAGS file from include and source files     ="
	@echo "=                                                                      ="
	@echo "=    clean      : removes all object files in subdirs                  ="
	@echo "=    distclean  : removes all object files, libraries, etc in subdirs  ="
	@echo "========================================================================"

deps:
	cd common/scripts && $(MAKE) $(MAKEOPTIONS)
	@cd source && if $(MAKE) $(MAKEOPTIONS) Makefile.dep;\
	  then : ;\
	  else rm Makefile.dep && $(MAKE) $(MAKEOPTIONS) Makefile.dep; fi
	cd lib && $(MAKE) $(MAKEOPTIONS) external-links


all: debug optimized

debug: deps contrib
	@cd source && $(MAKE) debug

optimized: deps contrib
	@cd source && $(MAKE) optimized

contrib: contrib-functionparser
	cd contrib && $(MAKE)

contrib-functionparser:
	@cd contrib && $(MAKE) functionparser

online-doc doc:
	cd doc && $(MAKE)
	@echo
	@echo
	@echo =======================================================
	@echo "The online documentation can now be accessed through"
	@echo "  doc/index.html"
	@echo =======================================================
	@echo
	@echo

TODO:
	@egrep -i '// *todo' include/deal.II/*/*.h source/*/*.cc \
	| perl common/scripts/make_todo.pl > $@


TAGS:
	@cd common ; etags --language=c++ \
	../include/deal.II/*/*.h ../source/*/*.cc \
	../examples/*/*.cc ; \
	perl -pi -e 's/ \* .*//g; s&[\t]*/\*.*&&g; s&[\t]*\*/.*&&g;' TAGS


clean: clean-contrib clean-source \
       clean-doc clean-examples clean-lib clean-tests
	-rm -f TODO common/TAGS

clean-contrib:
	cd contrib && $(MAKE) clean

clean-source:
	cd source && $(MAKE) CLEAN=yes clean

clean-doc:
	cd doc && $(MAKE) CLEAN=yes clean

clean-examples:
	cd examples && $(MAKE) CLEAN=yes clean

clean-lib:
	cd lib && $(MAKE) CLEAN=yes clean

# if that directory exists, go into 'tests' and clean up there as well
clean-tests:
	-cd tests && $(MAKE) clean

clean-common-scripts:
	cd common/scripts && $(MAKE) clean

distclean: clean
	cd lib && $(MAKE) distclean
	-cd tests && $(MAKE) distclean
	-rm source/*/*.inst
	$(MAKE)	clean-common-scripts



# #############################################################################
# targets for build tests. these targets are supposed to be called for regular
# build tests that send their results back home to our build test server. note
# that this target is supposed to be called right after checking out a version
# from our CVS server, so at this time, common/Make.global_options doesn't
# exist yet (and isn't included at the top of this file, due to a
# conditional), and we first have to set up all this stuff.
#
# the way this target works is that it first generates a bunch of lines for
# the output file that will ultimately be sent home. it then calls ./configure
# and stashes away both the output and the return value. if ./configure
# succeeds, we extract a little more information concluding the
# header. irrespective of this, we then append the output of ./configure to
# the log file. if the call to ./configure failed, we send the result back
# home and exit.
#
# if ./configure succeeded, we then go on to the build-test-do target, which
# actually compiles something, and when we come back from there we send the
# whole thing over to the server
build-test:
	echo AUTOMATED DEAL.II BUILD TEST          | tee    build-test-log
	-(./configure $(BUILDTESTFLAGS) 2>&1) > build-test-config
	echo BEGIN HEADER `date -u '+%Y-%m-%d %T'`        | tee -a build-test-log
	(cd common ; $(MAKE) BUILDTEST=yes -f Make.global_options print-summary) | tee -a build-test-log
	echo END HEADER `date -u '+%Y-%m-%d %T'`             | tee -a build-test-log
	echo                                  | tee -a build-test-log
	echo BEGIN CONFIGURE OUTPUT `date -u '+%Y-%m-%d %T'` | tee -a build-test-log
	cat build-test-config                 | tee -a build-test-log
	echo END CONFIGURE OUTPUT `date -u '+%Y-%m-%d %T'`   | tee -a build-test-log
	echo                                  | tee -a build-test-log
	-$(MAKE) BUILDTEST= build-test-do 2>&1 | tee -a build-test-log

# target to do the actual compilation tests for an automated build test. build
# the library, the example programs, and the doxygen example programs. if the
# testsuite has been checked out, run that, too, and send the results to the
# testsuite server
build-test-do:
	@echo BEGIN MAKE LIB OUTPUT `date -u '+%Y-%m-%d %T'`
	$(MAKE) all
	@echo END MAKE LIB OUTPUT `date -u '+%Y-%m-%d %T'`
	@echo
	@echo BEGIN MAKE EXAMPLES OUTPUT `date -u '+%Y-%m-%d %T'`
	cd examples ; $(MAKE)
	@echo END MAKE EXAMPLES OUTPUT `date -u '+%Y-%m-%d %T'`
	@echo
	@echo BEGIN MAKE RUN EXAMPLES OUTPUT `date -u '+%Y-%m-%d %T'`
	cd examples/doxygen ; $(MAKE) run
	@echo END MAKE RUN EXAMPLES OUTPUT `date -u '+%Y-%m-%d %T'`
	@echo
	@echo BUILD TEST SUCCESSFUL `date -u '+%Y-%m-%d %T'`


.PHONY: all \
        online-doc doc printable-doc tex-doc contrib \
	contrib-functionparser \
        clean clean-contrib clean-base clean-lac clean-dealII \
	clean-doc clean-examples clean-lib clean-tests clean-common-scripts \
	distclean \
	build-test build-test-do \
	TAGS TODO

# mark the following two as phony since they mark two actual files but
# we don't want to list the full dependencies here. let the sub-makefiles
# do that
.PHONY: common/scripts/make_dependencies \
        common/scripts/expand_instantiations
