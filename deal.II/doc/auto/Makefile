# Makefile,v 1.8 1999/04/15 17:36:14 wolf Exp


##########################################################################
# Change this section only, if you want to add a documentation library   #
##########################################################################
#
# specify which are the include files that a specific documentation libary
# depends upon
kdoc.base     = $(wildcard ../../base/include/base/*.h)
kdoc.lac      = $(wildcard ../../lac/include/lac/*.h)
kdoc.basic    = $(wildcard ../../deal.II/include/basic/*.h)
kdoc.grid     = $(wildcard ../../deal.II/include/grid/tria*.h)\
                $(wildcard ../../deal.II/include/grid/grid*.h)\
                ../../deal.II/include/grid/geometry_info.h
kdoc.dof      = $(wildcard ../../deal.II/include/grid/dof*.h) \
                $(wildcard ../../deal.II/include/grid/mg_*.h) \
                $(wildcard ../../deal.II/include/fe/*.h)
kdoc.numerics = $(wildcard ../../deal.II/include/numerics/*.h)

# for latex-output, we join the different libraries of the
# deal.II library into just one file
kdoc.dealII = $(kdoc.basic) $(kdoc.grid) $(kdoc.dof) $(kdoc.numerics)

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
ifeq ($(shell uname),SunOS)
DOCPP = /home/people/kanschat/bin/doc++
else
DOCPP = doc++
endif


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


# generate the input file for doc++
%.doc++: Makefile
	@echo $(kdoc.$(subst .doc++,,$(subst latex/,,$@)))  | \
	 perl -pi -e 's# #\n#g;  s#\.\./\.\./#//\@Include: ../../../#g;' > \
	 $@


# generate one of the latex-output files; this needs the respective
# include file for DOC++
#
# delete the .doc++ files immediately again, since they are no more 
# used and since this forces make to re-make the .tex-file each time
# we ask for it; the right way to do so would be to use dependencies,
# but it is hard to compute them from the '%'-sign used in the rule
%.tex: %.doc++
	$(DOCPP) -p -t -u -o $@ -ep a4wide $<
	rm $^


# generate temporary files from the doc++-latex-output which have
# headers and footers stripped and label names replaced by something
# file dependent
%.tex1: %.tex
	cat $<                           | \
	 perl scripts/process-tex        | \
	 perl -pi -e 's/cxx\./$(basename $(notdir $@)).cxx./g;' > $@

latex/deal_II_technical_reference.tex: latex/base.tex1 latex/lac.tex1 latex/dealII.tex1
	@cat latex/header.tex                              > $@
	@echo "\\\\chapter{Base library}"                 >> $@
	@echo "\\\\def\\\\presentlib{base}"               >> $@
	@cat latex/base.tex1                              >> $@
	@echo "\\\\chapter{Linear algebra library}"       >> $@
	@echo "\\\\def\\\\presentlib{lac}"                >> $@
	@cat latex/lac.tex1                               >> $@
	@echo "\\\\chapter{deal.II library}"              >> $@
	@echo "\\\\def\\\\presentlib{dealII}"             >> $@
	@cat latex/dealII.tex1                            >> $@
	@cat latex/footer.tex                             >> $@
	rm $^



# generate the latex-processed files from these
%.dvi: %.tex
	cd latex ; latex $(notdir $<)
	cd latex ; latex $(notdir $<)
	-rm $(subst .tex,.log,$<)
	-rm $(subst .tex,.aux,$<)


# generate postscript files from a latex output file
%.ps: %.dvi
	dvips $< -o $@


clean:
	-rm -f kdoc/names.html kdoc/*.kdoc kdoc/*/*.html cvs-backlog/*.html \
		latex/base.* latex/lac.* latex/dealII.* latex/deal_II_technical_reference.*


.PHONY: default all kdoc $(kdoc.library-files) $(kdoc.library-files.rerun) clean
