# $Id$


##########################################################################
# Change this section only, if you want to add a documentation library   #
##########################################################################
#
# specify which are the include files that a specific documentation libary
# depends upon
kdoc.base = $(wildcard ../../base/include/base/*.h)
kdoc.lac  = $(wildcard ../../lac/include/lac/*.h)
kdoc.basic= $(wildcard ../../deal.II/include/basic/*.h)
kdoc.grid = $(wildcard ../../deal.II/include/grid/tria*.h)\
            ../../deal.II/include/grid/geometry_info.h
kdoc.dof  = $(wildcard ../../deal.II/include/grid/dof*.h) \
            $(wildcard ../../deal.II/include/grid/mg_*.h) \
            $(wildcard ../../deal.II/include/fe/*.h)
kdoc.num  = $(wildcard ../../deal.II/include/numerics/*.h)

# user projects. the projects should be independant of each other, so
# we need not take care of cross-dependencies each time 
kdoc.transport     = $(wildcard ../../projects/transport/trans_dg/*.h)
kdoc.wave-equation = $(wildcard ../../projects/wave-equation/include/*.h)


# specify the names of the resulting doc libraries
kdoc.library-files = kdoc/base.kdoc         \
                     kdoc/lac.kdoc          \
                     kdoc/basic.kdoc        \
                     kdoc/grid.kdoc         \
                     kdoc/dof.kdoc          \
                     kdoc/numerics.kdoc     \
                     kdoc/transport.kdoc    \
                     kdoc/wave-equation.kdoc
##########################################################################
# End of changeable section                                              #
##########################################################################


# make rule names to rerun the generation of a library
kdoc.library-files.rerun = $(kdoc.library-files:%.kdoc=%.kdoc.rerun)


KDOCFLAGS   = -I../scripts/kdoc ../scripts/kdoc/kdoc -a -p


# by default: make docu libraries and generate an index of all names
default: kdoc kdoc/names.html

# alternatively: do the above and also generate a backlog of CVS
# changes in the last 100 days
all: default cvslog

cvslog:
	cd ../.. ; perl doc/auto/scripts/cvs2html -o doc/auto/cvs-backlog/newdeal -a -k -D 100


# make the .kdoc files twice: once without and a second time with
# crossreferencing the other libraries
kdoc: $(kdoc.library-files) $(kdoc.library-files.rerun)


# rule how to generate the documentation the first time, i.e. without
# crossreferencing the other doc-libraries. we do not try to make it
# dependent on anything, since it needs to be generated anyway
#
# since the rule is a rather complex one, involving the computation of
# variable names, here is the stripped down version for the base library:
# kdoc/base.kdoc:
#	@cd kdoc ; perl $(KDOCFLAGS) -dbase \
#		base $(kdoc.base:../..%=../../..%)
$(kdoc.library-files):
	@cd kdoc ;                                                 \
	 perl   $(KDOCFLAGS)                                       \
		-d$(subst .kdoc,,$(subst kdoc/,,$@))               \
		$(subst .kdoc,,$(subst kdoc/,,$@))                 \
		$($(subst /,.,$(subst .kdoc,,$@)):../..%=../../..%)

# rule how to generate the documentation the second time, i.e. with
# crossreferencing the other doc-libraries this time. the rule is
# similar to the above one, but includes the other libraries
$(kdoc.library-files.rerun):
	@cd kdoc ;                                                       \
	 perl   $(KDOCFLAGS) $(kdoc.library-files:kdoc/%.kdoc=-l%)       \
		-d$(subst .kdoc.rerun,,$(subst kdoc/,,$@))               \
		$(subst .kdoc.rerun,,$(subst kdoc/,,$@))                 \
		$($(subst /,.,$(subst .kdoc.rerun,,$@)):../..%=../../..%)


kdoc/names.html: $(kdoc.libray-files)
	@cd kdoc ; perl ../scripts/index.pl *.kdoc > names.html


clean:
	-rm -f kdoc/names.html kdoc/*.kdoc kdoc/*/*.html cvs-backlog/*.html


.PHONY: default all kdoc $(kdoc.library-files) $(kdoc.library-files.rerun) clean