# Makefile,v 1.8 1999/04/15 17:36:14 wolf Exp


##########################################################################
# Change this section only, if you want to add a documentation library   #
##########################################################################
#
# specify which are the include files that a specific documentation libary
# depends upon
kdoc.base     = $(shell echo ../../base/include/base/*.h)
kdoc.lac      = $(shell echo ../../lac/include/lac/*.h)
kdoc.grid     = $(shell echo ../../deal.II/include/grid/*.h)
kdoc.dof      = $(shell echo ../../deal.II/include/dofs/*.h) \
                $(shell echo ../../deal.II/include/fe/*.h)
kdoc.numerics = $(shell echo ../../deal.II/include/numerics/*.h)

# for latex-output, we join the different libraries of the
# deal.II library into just one file
kdoc.dealII = $(kdoc.grid) $(kdoc.dof) $(kdoc.numerics)

# specify the names of the resulting doc libraries
kdoc.library-files = kdoc/base.kdoc         \
                     kdoc/lac.kdoc          \
                     kdoc/basic.kdoc        \
                     kdoc/grid.kdoc         \
                     kdoc/dof.kdoc          \
                     kdoc/numerics.kdoc

##########################################################################
# End of changeable section                                              #
##########################################################################


# include global paths and variables
include ../../common/Make.global_options

# make rule names to rerun the generation of a library
kdoc.library-files.rerun = $(kdoc.library-files:%.kdoc=%.kdoc.rerun)


KDOCFLAGS   = -I../scripts/kdoc ../scripts/kdoc/kdoc -a -p
#KDOCFLAGS   = /home/people/wolf/var/kdoc/kdoc2/bin/kdoc -a -p \
#              --html-logo=../../../pictures/deal.II_3.0_2.150.jpg
ifeq ($(shell uname),SunOS)
DOCPP = /home/people/wolf/Config/doc++/bin/doc++
#DOCPP = /home/people/kanschat/bin/doc++
else
DOCPP = doc++
endif


# by default: make library doc and generate an index of all names
default: kdoc kdoc/names.html

# alternatively: do the above and also generate a backlog of CVS
# changes in the last 100 days
all: default cvslog


#################################################################################
##  targets for cvslog generated HTML output    
#################################################################################

cvslog:
	cd ../.. ; $(PERL) doc/auto/scripts/cvs2html -o doc/auto/cvs-backlog/newdeal -a -k -D 100



#################################################################################
##  targets for kdoc generated HTML output    
#################################################################################

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
#	@cd kdoc ; $(PERL) $(KDOCFLAGS) -dbase \
#		base $(kdoc.base:../..%=../../..%)
$(kdoc.library-files):
	@cd kdoc ;                                                 \
	 $(PERL)   $(KDOCFLAGS)                                    \
		-d$(subst .kdoc,,$(subst kdoc/,,$@))               \
		$(subst .kdoc,,$(subst kdoc/,,$@))                 \
		$($(subst /,.,$(subst .kdoc,,$@)):../..%=../../..%)

# rule how to generate the documentation the second time, i.e. with
# crossreferencing the other doc-libraries this time. the rule is
# similar to the above one, but includes the other libraries
$(kdoc.library-files.rerun):
	@cd kdoc ;                                                       \
	 $(PERL)   $(KDOCFLAGS) $(kdoc.library-files:kdoc/%.kdoc=-l%)    \
		-d$(subst .kdoc.rerun,,$(subst kdoc/,,$@))               \
		$(subst .kdoc.rerun,,$(subst kdoc/,,$@))                 \
		$($(subst /,.,$(subst .kdoc.rerun,,$@)):../..%=../../..%)


kdoc/names.html: $(kdoc.libray-files)
	@cd kdoc ; $(PERL) ../scripts/index.pl *.kdoc > names.html




#################################################################################
##  targets for DOC++ generated latex output    
#################################################################################

# generate the input file for doc++
%.doc++: Makefile
	@echo $(kdoc.$(subst .doc++,,$(subst latex/,,$@)))  | \
	 $(PERL) -pi -e 's# #\n#g;  s#\.\./\.\./#//\@Include: ../../../#g;' > \
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
	 $(PERL) scripts/process-tex        | \
	 $(PERL) -pi -e 's/cxx\./$(basename $(notdir $@)).cxx./g;' > $@

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
	-rm $(subst .tex,.toc,$<)


# generate postscript files from a latex output file
%.ps: %.dvi
	cd latex ; dvips $(notdir $<) -o $(notdir $@)




#################################################################################
##  targets for DOC++ generated HTML output    
#################################################################################
doc++: doc++/base/index.html doc++/lac/index.html doc++/dealII/index.html

# generate one of the html-output files; this needs the respective
# include file for DOC++
#
# delete the .doc++ files immediately again, since they are no more 
# used and since this forces make to re-make the .html-file each time
# we ask for it; the right way to do so would be to use dependencies,
# but it is hard to compute them from the '%'-sign used in the rule
doc++/%/index.html: %.doc++
	cp $^ $(patsubst %/index.html,%,$@)
	cd $(patsubst %/index.html,%,$@) ; $(PERL) -pi -e 's!: ../!: ../../!g;' $^
	cd $(patsubst %/index.html,%,$@) ; $(DOCPP) -p -H -u -G -S -a $^
	cd $(patsubst %/index.html,%,$@) ; rm $^
	rm $^




#################################################################################
##  administrativa
#################################################################################
clean:
	-rm -f kdoc/names.html kdoc/*.kdoc kdoc/*/*.html cvs-backlog/*.html \
		latex/base.* latex/lac.* latex/dealII.* \
                latex/deal_II_technical_reference.* \
		doc++/base/* doc++/lac/* doc++/dealII/* \
		*.doc++


.PHONY: default all kdoc $(kdoc.library-files) $(kdoc.library-files.rerun) doc++ clean
