# $Id$
#
# Makefile for the whole library

# include global options and pathes. there's one exception where we don't want
# that, namely for build tests: in that case, the makefile is called before we
# have even run ./configure, so we can't include this file, but the target
# then called doesn't need the variables defined in there doesn't need them
# anyway and calls ./configure first
ifneq ($(BUILDTEST),yes)
  include common/Make.global_options
endif


help:
	@echo "========================================================================"
	@echo "=              Global Makefile for the deal.II libraries               ="
	@echo "========================================================================"
	@echo "=                                                                      ="
	@echo "= The following targets exist:                                         ="
	@echo "=    baseg lacg : debug libraries \"base\" and \"lac\"                     ="
	@echo "=    baseo laco : optimized libraries \"base\" and \"lac\"                 ="
	@echo "=    base lac   : debug and optimized libraries \"base\" and \"lac\"       ="
	@echo "=    1dg 2dg 3dg: deal.II debug version for specified dimension        ="
	@echo "=    1do 2do 3do: deal.II optimized version for specified dimension    ="
	@echo "=    1d 2d 3d   : deal.II debug and optimized for specified dimension  ="
	@echo "=                                                                      ="
	@echo "=    contrib    : additional libraries in contrib, if there are any    ="
	@echo "=                                                                      ="
	@echo "=    all        : base, lac, 1d, 2d and 3d, contrib                    ="
	@echo "=    debug      : baseg, lacg, 1dg, 2dg and 3dg, contrib               ="
	@echo "=    optimized  : baseo, laco, 1do, 2do and 3do, contrib               ="
	@echo "=                                                                      ="
	@echo "=    online-doc : generate the documentation in HTML format            ="
	@echo "=    TODO       : create a \"TODO\" file from the source files           ="
	@echo "=    TAGS       : create a TAGS file from include and source files     ="
	@echo "=                                                                      ="
	@echo "=    clean      : removes all object files in subdirs                  ="
	@echo "=    distclean  : removes all object files, libraries, etc in subdirs  ="
	@echo "========================================================================"



baseg: contrib-functionparser
	cd $D/base && $(MAKE) $(MAKEOPTIONS) libg

baseo: contrib-functionparser
	cd $D/base && $(MAKE) $(MAKEOPTIONS) libo

lacg: baseg contrib
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) libg

laco: baseo contrib
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) libo

1dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 1dg

2dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 2dg

3dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 3dg

1do:baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 1d

2do:baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 2d

3do: baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 3d

deps: common/scripts/make_dependencies
	cd $D/base && $(MAKE) $(MAKEOPTIONS) Makefile.dep
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) Makefile.dep
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) Makefile.dep
	cd $D/lib && $(MAKE) $(MAKEOPTIONS) external-links

common/scripts/make_dependencies:
	cd $D/common/scripts && $(MAKE) $(MAKEOPTIONS)

base: baseg baseo
lac: lacg laco
1d: 1dg 1do
2d: 2dg 2do
3d: 3dg 3do

$(LIBDIR):
	@mkdir $@

baseg baseo base lacg laco lac 1dg 1do 1d 2dg 2do 2d 3dg 3do 3d all debug contrib -functionparser : $(LIBDIR) deps

all: deps
	$(MAKE) $(MAKEOPTIONS) contrib base lac 2d 1d 3d

debug: deps
	$(MAKE) $(MAKEOPTIONS) contrib baseg lacg 2dg 1dg 3dg

optimized: deps
	$(MAKE) $(MAKEOPTIONS) contrib baseo laco 2do 1do 3do

# when HSL functions are used, we also compile the detached_ma27
# program which itself uses the optimized base library; so add this to
# the dependencies if necessary
#
# make sure the function parser is also built (it isn't in the default
# build list of the subdirectory, since it needs to be built before
# libbase
contrib: $(subst no,,$(subst yes,base,$(USE_CONTRIB_HSL))) \
         contrib-functionparser
	cd $D/contrib && $(MAKE)

contrib-functionparser:
	@cd $D/contrib && $(MAKE) functionparser

online-doc:
	cd $D/doc && $(MAKE)
	@echo
	@echo
	@echo =======================================================
	@echo "The online documentation can now be accessed through"
	@echo "  $D/doc/index.html"
	@echo =======================================================
	@echo
	@echo

TODO:
	@egrep -i '// *todo' $D/base/include/*/*.h $D/base/source/*.cc \
                     $D/lac/include/*/*.h $D/lac/source/*.cc   \
                     $D/deal.II/include/*/*.h $D/deal.II/source/*/*.cc \
	| perl -pi -e 's#$D/?#\n#;' \
	| perl -pi -e 's#(.*)(TODO:?)(\[[^\]]+\])#\3 \1\2\3#i;' \
	| perl -pi -e 's#\s*//\s*TODO:?\s*(\[[^\]]+\]):?\s*#\n\1     #i;' \
             > $@


TAGS:
	@cd $D/common ; etags $D/base/source/*.cc \
	$D/lac/source/*.cc $D/deal.II/source/*/*.cc $D/deal.II/source/*/*.h \
	$D/tests/*/*.cc $D/examples/*/*.cc \
	--language=c++ \
	$D/base/include/*/*.h $D/lac/include/*/*.h \
	$D/deal.II/include/*/*.h $D/tests/*/*.h; \
	perl -pi -e 's/ \* .*//g; s&[\t]*/\*.*&&g; s&[\t]*\*/.*&&g;' TAGS


clean: common/scripts/make_dependencies \
       clean-contrib clean-base clean-lac clean-dealII \
       clean-doc clean-examples clean-lib clean-tests
	-rm -f TODO common/TAGS

clean-contrib:
	cd $D/contrib && $(MAKE) clean

clean-base:
	cd $D/base && $(MAKE) clean

clean-lac:
	cd $D/lac && $(MAKE) clean

clean-dealII:
	cd $D/deal.II && $(MAKE) clean

clean-doc:
	cd $D/doc && $(MAKE) clean

clean-examples:
	cd $D/examples && $(MAKE) clean

clean-lib:
	cd $D/lib && $(MAKE) clean

# if that directory exists, go into 'tests' and clean up there as well
clean-tests:
	-cd tests && $(MAKE) clean

clean-common-scripts:
	cd $D/common/scripts && $(MAKE) clean

distclean: clean \
           clean-common-scripts
	cd lib && $(MAKE) distclean
	-cd tests && $(MAKE) distclean



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
	echo BEGIN HEADER `date -u '+%Y-%m-%d %T'`                | tee -a build-test-log
	echo USER=$(USER)                          | tee -a build-test-log
	echo HOST=`hostname`.`dnsdomainname`       | tee -a build-test-log
	echo BUILDTESTFLAGS=$(BUILDTESTFLAGS)      | tee -a build-test-log
	echo BUILDTESTNOTE=$(BUILDTESTNOTE)        | tee -a build-test-log
	(./configure $(BUILDTESTFLAGS) 2>&1) > build-test-config
	if test -f common/Make.global_options ; then \
	   grep '^enable-shared' common/Make.global_options | grep = ; \
	   grep enable-multithreading common/Make.global_options | grep = ; \
	   grep with-multithreading common/Make.global_options | grep = ; \
	   grep USE_CONTRIB_PETSC common/Make.global_options | grep = ; \
	   grep USE_CONTRIB_METIS common/Make.global_options | grep = ; \
	   grep USE_CONTRIB_HSL common/Make.global_options | grep = ; \
	   grep USE_CONTRIB_UMFPACK common/Make.global_options | grep = ; \
	 fi    | tee -a build-test-log
	echo END HEADER `date -u '+%Y-%m-%d %T'`             | tee -a build-test-log
	echo                                  | tee -a build-test-log
	echo BEGIN CONFIGURE OUTPUT `date -u '+%Y-%m-%d %T'` | tee -a build-test-log
	cat build-test-config                 | tee -a build-test-log
	echo END CONFIGURE OUTPUT `date -u '+%Y-%m-%d %T'`   | tee -a build-test-log
	echo                                  | tee -a build-test-log
	$(MAKE) BUILDTEST= build-test-do 2>&1 | tee -a build-test-log

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
	@echo BEGIN MAKE DOXYGENTUTORIALS OUTPUT `date -u '+%Y-%m-%d %T'`
	cd examples/doxygen ; $(MAKE)
	@echo END MAKE DOXYGENTUTORIALS OUTPUT `date -u '+%Y-%m-%d %T'`
	@echo 
#	if test -d tests ; then \
#	   @echo BEGIN MAKE TESTS OUTPUT `date -u '+%Y-%m-%d %T'` ; \
#	   cd tests ; $(MAKE) report+mail ; \
#	   @echo END MAKE TESTS OUTPUT `date -u '+%Y-%m-%d %T'` ; \
#	   @echo ; \
#	 fi
	@echo BUILD TEST SUCCESSFUL `date -u '+%Y-%m-%d %T'`


.PHONY: base baseg lac lacg 1d 2d 3d 1dg 2dg 3dg all \
        online-doc printable-doc tex-doc contrib \
        clean clean-contrib clean-base clean-lac clean-dealII \
	clean-doc clean-examples clean-lib clean-tests clean-common-scripts \
	distclean \
	build-test build-test-do \
	TAGS TODO
