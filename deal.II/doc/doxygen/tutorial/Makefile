# $Id$
# Create proper HTML output from the different example programs

# include paths and global variables
include ../../../common/Make.global_options

# plain names of those files that we shall treat. simply take all those
# directories in examples/ that have a /doc subdirectory
example-names = $(subst /doc,,    \
                        $(shell cd $D/examples && echo */doc))



# prefix names for generated intermediate files
gen-example-names = $(addprefix generated/, $(example-names))

# names of commented program files, plain files, and tocs for the commented
# and commented and plain versions of the program as well as a combined one
example-toc-prog   = $(addsuffix _toc.prog,$(gen-example-names))
example-toc-plain  = $(addsuffix _toc.plain,$(gen-example-names))
example-toc-intro  = $(addsuffix _toc.intro,$(gen-example-names))
example-toc-results= $(addsuffix _toc.results,$(gen-example-names))
example-toc        = $(addsuffix _toc.combined,$(gen-example-names))

# same for intermediate files for doxygen output
example-dox-prog   = $(addsuffix _prog.dox,$(gen-example-names))
example-dox-plain  = $(addsuffix _plain.dox,$(gen-example-names))
example-plain      = $(addprefix plain/,$(addsuffix .cc,$(example-names)))

# finally a target for combined Doxygen files
example-doxygen     = $(addprefix doxygen/,$(addsuffix .h,$(example-names)))


# a makefile command that extracts the base name of a program
get_basename = $(shell echo $(1) | $(PERL) -pi -e 's/(generated|doxygen).//g; s/\..*//g; s/_(prog|plain|toc.*)//g;')



all: $(example-doxygen) toc.html

# Generation of the map of the example programs. Note, that the
# following command only merges the map into toc.html. The map itself
# is not generated automatically, but requires to call
#
# make steps.png
#
#  by hand. This way, we avoid testing for neato during
#  configuration. The map only needs to be generated after changes in
#  the structure of the tutorials.
toc.html: steps.cmapx toc.html.in
	@echo ================== Making $@
	@perl make_toc.pl > toc.html

# generate dot file
steps.dot: steps.pl $D/examples/*/doc/tooltip
	@echo ================== Making $@
	@perl $< > $@

steps.png: steps.dot
	@echo ================== Making $@
	@neato -Tpng -Tcmapx -O $<
	@mv steps.dot.png steps.png
	@mv steps.dot.cmapx steps.cmapx


$(example-toc-prog):
	@echo ================== Making $@
	@cat $D/examples/$(call get_basename, $@)/*.cc \
	  | $(PERL) program2toc \
	  > $@


$(example-toc-plain):
	@echo ================== Making $@
	@cat generated/$(call get_basename, $@)_toc.prog \
	  | $(PERL) -pi -e 's/href=\"\#/href=\"\#plain-/g;' \
	  > $@


$(example-toc-intro):
	@echo ================== Making $@
	@cat $D/examples/$(call get_basename, $@)/doc/intro.dox \
	  | egrep -i "<h[3456]>.*</h[3456]>" \
	  | $(PERL) ./intro2toc \
	  > $@


$(example-toc-results):
	@echo ================== Making $@
	@cat $D/examples/$(call get_basename, $@)/doc/results.dox \
	  | egrep -i "<h[3456]>.*</h[3456]>" \
	  | $(PERL) ./intro2toc \
	  > $@


$(example-toc):
	@echo ================== Making $@
	@echo "<table class=\"tutorial\" width=\"50%\">" > $@
	@echo "<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>" >> $@
	@echo "<tr><td width=\"50%\" valign=\"top\">" >> $@
	@echo "<ol>" >> $@
	@echo "  <li> <a href=\"#Intro\" class=bold>Introduction</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.intro >> $@
	@echo "  <li> <a href=\"#CommProg\" class=bold>The commented program</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.prog >> $@
	@echo "</ol></td><td width=\"50%\" valign=\"top\"><ol>" >> $@
	@echo "  <li value=\"3\"> <a href=\"#Results\" class=bold>Results</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.results >> $@
	@echo "</ol> </td> </tr> </table>" >> $@

# To revert, enter those lines before the last
#	@echo "  <li> <a href=\"#PlainProg\" class=bold>The plain program</a>" >> $@
#	@cat generated/$(call get_basename, $@)_toc.plain >> $@


$(example-dox-prog):
	@echo ================== Making $@
	@cat $D/examples/$(call get_basename, $@)/*.cc \
	  | $(PERL) program2doxygen \
	  > $@

$(example-plain):
	@echo ================== Making $@
	@cat $D/examples/$(call get_basename, $(@F))/$(@F) | $(PERL) program2plain >> $@

$(example-dox-plain):
	@echo ================== Making $@
	@echo " * <a name=\"PlainProg\"></a>" > $@
	@echo " * <h1> The plain program</h1>" >> $@
	@echo " * " >> $@
	@echo " * (If you are looking at a locally installed deal.II version, then the" >> $@
	@echo " * program can be found at <i>" >> $@
	@echo " * " $D/examples/$(call get_basename,$@)/*.cc | $(PERL) -pi -e 's!/! /!g;' >> $@
	@echo " * </i>. Otherwise, this is only" >> $@
	@echo " * the path on some remote server.)" >> $@
	@echo " @code" >> $@
	@cat $D/examples/$(call get_basename, $@)/*.cc | $(PERL) program2doxyplain >> $@
	@echo " @endcode" >> $@


$(example-doxygen):
	@echo ================== Assembling $@
	@echo "/**" > $@
	@echo " * @page $(shell echo $(call get_basename,$@) | $(PERL) -pi -e 's/-/_/g;') \
	                The $(call get_basename,$@) tutorial program" >> $@
	@echo "@htmlonly" >> $@
	@cat generated/$(call get_basename, $@)_toc.combined >> $@
	@echo "@endhtmlonly" >> $@
	@cat $D/examples/$(call get_basename, $@)/doc/intro.dox \
		| $(PERL) create_anchors >> $@
	@cat generated/$(call get_basename, $@)_prog.dox >> $@
	@cat $D/examples/$(call get_basename, $@)/doc/results.dox \
		| $(PERL) create_anchors >> $@
	@echo " <a name=\"PlainProg\"></a>" >> $@
	@echo " <h1> The plain program</h1>" >> $@
	@echo " @include \"$(call get_basename, $@).cc\"" >> $@
	@echo " */" >> $@

# Previous version of the lines after the last 'create_anchor'
#	@cat generated/$(call get_basename, $@)_plain.dox >> $@

validate-xrefs:
	@$(PERL) $D/common/scripts/validate-xrefs.pl \
                $(filter-out head.html foot.html, \
                             $(shell echo *.html))


clean:
	-rm -f  plain/*.cc $(example-toc) \
		$(example-toc-prog) $(example-toc-plain) \
		$(example-toc-intro) $(example-toc-results) \
		$(example-dox-prog) $(example-dox-plain) $(example-doxygen)
	-rm -f Makefile.dep


.PHONY: validate-xrefs clean


Makefile.dep: $(shell echo $D/examples/*/*.cc $D/examples/*/doc/*.dox) Makefile
	@echo ================== Generating $@
	@for i in $(example-names) ; do \
		echo generated/$${i}_prog.html generated/$${i}_plain.html \
		     generated/$${i}_toc.prog \
		     generated/$${i}_prog.dox generated/$${i}_plain.dox : \
			`echo $D/examples/$${i}/*.cc` ; \
		echo generated/$${i}_toc.plain : generated/$${i}_toc.prog ; \
		echo generated/$${i}_toc.combined : \
			generated/$${i}_toc.intro \
			generated/$${i}_toc.results \
			generated/$${i}_toc.prog \
			generated/$${i}_toc.plain ; \
		echo plain/$$i.cc : `echo $D/examples/$${i}/*.cc` ; \
		echo doxygen/$$i.h : generated/$${i}_prog.dox generated/$${i}_plain.dox \
				generated/$${i}_toc.combined \
				$D/examples/$$i/doc/intro.dox \
				$D/examples/$$i/doc/results.dox \
				plain/$$i.cc ; \
		echo generated/$${i}_toc.intro : \
				$D/examples/$$i/doc/intro.dox \
			 ; \
		echo generated/$${i}_toc.results : \
				$D/examples/$$i/doc/results.dox \
			 ; \
	 done \
	> $@

ifneq ($(CLEAN),yes)
include Makefile.dep
endif
