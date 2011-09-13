# collect all header files, and for each try to do the following:
# generate a file tmp.cc that contains nothing but "#include <xxx>"
# with this particular header file. then try to compile it. this
# sometimes fails since this header does not include all other headers
# it actually needs. we don't usually see this since we may be
# including these missing headers into our files before, so that nobody
# notices the problem. expose it this way.

include ../Makefile.paths
include $D/common/Make.global_options


############################################################
# First how to create executables, including all necessary
# flags:
############################################################


%.OK :
	@echo "=================== $@"
	@if test -n "`echo $< | grep 2d`" ; then \
	  (($D/contrib/mesh_conversion/convert_mesh 2 $< $(@:OK=ucd) 2>&1) > $(@:OK=log)) ; \
	else \
	  (($D/contrib/mesh_conversion/convert_mesh 3 $< $(@:OK=ucd) 2>&1) > $(@:OK=log)) ; \
	fi
	@perl -pi -e 's#$D/contrib/mesh_conversion/##g;' $(@:OK=ucd)
	@cmd="diff $(<:%.inp=%.ucd) $(@:OK=ucd) > /dev/null" ; \
	 if eval $$cmd ; then \
			echo '=====OK============' $@ ; \
			echo 'diff ok' > $(dir $@)/status ; \
			echo $$i > $@ ; \
	 else \
			echo '+++++Error+++++++++' $@ ; \
	 fi



default: all

# collect all headers
MESH_INPUTS = $(shell echo $D/contrib/mesh_conversion/mesh/?d/*.inp)
TESTS_X   = $(shell for i in $(MESH_INPUTS) ; do echo $$i | $(PERL) -pi -e 's/.inp$$/.OK/g;' ; done)
TESTS     = $(shell echo $(TESTS_X:$D/%=%) | $(PERL) -pi -e 's/\//-/g;')

ID = `id -un`"@"`hostname`


# how to make the executable if that hasn't happened yet
$D/contrib/mesh_conversion/convert_mesh:
	@cd $(@D) && $(MAKE)

# next the target that creates the target rules. First, we list all the header
# files as the head of the rule, followed by a colon. Then the actual rule
# that includes writing the .cc file that only includes a single header file,
# compiles it, and then deletes everything again
Makefile.dep: $D/contrib/mesh_conversion/convert_mesh $(MESH_INPUTS)
	for i in $(MESH_INPUTS:$D/%=%) ; do \
		echo `echo $$i | $(PERL) -pi -e 's/\.inp/.OK/g; s/\//-/g;'` : $D/$$i ; \
	 done > Makefile.dep

include Makefile.dep


# here's what's to be done
all: $(TESTS)


# create a report, where we get exactly one line per test.
report:
	@-$(MAKE) -k
	@-svn info $D | grep '^Revision'
	@echo 'Date: ' `date -u +"%Y %j  %F  %U-%w"`
	@echo 'Id: ' $(ID)
	@for test in $(TESTS) ; do \
		testname=mesh-converter/`echo $$test | perl -pi -e 's/.*include.//g; s/\//-/g;'` ; \
		if test -f $$test ; then \
			echo "`date -u +"%Y-%m-%d %H:%M"`   +  $$testname" ; \
		else \
			echo "`date -u +"%Y-%m-%d %H:%M"`  0   $$testname" ; \
		fi ; \
	done

report+mail: all
	@$(MAKE) report | tee testresults
	@cat testresults | mail regression-tests@dealii.org
	@rm testresults


clean:
	-rm -f Makefile.dep *.$(OBJEXT) \
               *.g.$(OBJEXT) *.output *~ *.OK *.cc *.ucd *.log


.PHONY: report clean distclean all default report report+mail
