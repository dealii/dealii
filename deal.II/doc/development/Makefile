# $Id$


D = ../..
include $D/common/Make.global_options


default: makefiles.html makefile.small.html makefile.large.html


makefiles.html: Makefile makefiles.1.html $D/common/Make.global_options makefiles.2.html
	@cat makefiles.1.html > $@
	@echo '<li><code>CXX=$(CXX)</code>' >> $@
	@echo '<li><code>CC=$(CC)</code>' >> $@
	@echo '<li><code>PERL=$(PERL)</code>' >> $@
	@echo '<li><code>lib-path-base=$(lib-path-base)</code>' >> $@
	@echo '<li><code>lib-path-lac=$(lib-path-lac)</code>' >> $@
	@echo '<li><code>lib-path-deal2=$(lib-path-deal2)</code>' >> $@
	@echo '<li><code>LIBPATH=$(LIBPATH)</code>' >> $@
	@echo '<li><code>lib-base.o=$(lib-base.o)</code>' >> $@
	@echo '<li><code>lib-base.g=$(lib-base.g)</code>' >> $@
	@echo '<li><code>lib-lac.o=$(lib-lac.o)</code>' >> $@
	@echo '<li><code>lib-lac.g=$(lib-lac.g)</code>' >> $@
	@echo '<li><code>lib-deal2-1d.o=$(lib-deal2-1d.o)</code>' >> $@
	@echo '<li><code>lib-deal2-1d.g=$(lib-deal2-1d.g)</code>' >> $@
	@echo '<li><code>lib-deal2-2d.o=$(lib-deal2-2d.o)</code>' >> $@
	@echo '<li><code>lib-deal2-2d.g=$(lib-deal2-2d.g)</code>' >> $@
	@echo '<li><code>lib-deal2-3d.o=$(lib-deal2-3d.o)</code>' >> $@
	@echo '<li><code>lib-deal2-3d.g=$(lib-deal2-3d.g)</code>' >> $@
	@echo '<li><code>include-path-base=$(include-path-base)</code>' >> $@
	@echo '<li><code>include-path-lac=$(include-path-lac)</code>' >> $@
	@echo '<li><code>include-path-deal2=$(include-path-deal2)</code>' >> $@
	@echo '<li><code>INCLUDE=$(INCLUDE)</code>' >> $@
	@echo '<li><code>CXXFLAGS.g=$(CXXFLAGS.g)</code>' >> $@
	@echo '<li><code>CXXFLAGS.o=$(CXXFLAGS.o)</code>' >> $@
	@echo '<li><code>ACE_ROOT=$(ACE_ROOT)</code>' >> $@
	@echo '<li><code>lib-ACE=$(lib-ACE)</code>' >> $@
	@echo '<li><code>with-multithreading=$(with-multithreading)</code>' >> $@
	@cat makefiles.2.html >> $@


makefile.small.html: Makefile.small
	@$(PERL) makefile2html < $< > $@

makefile.large.html: Makefile.large
	@$(PERL) makefile2html < $< > $@


validate-xrefs:
	$(PERL) $D/doc/auto/scripts/validate-xrefs.pl \
		$(filter-out makefiles.1.html makefiles.2.html, \
		  $(shell echo *.html))



clean:
	-rm -f makefiles.html makefile.small.html makefile.large.html

.PHONY: default validate-xrefs clean
