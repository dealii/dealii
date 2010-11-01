# $Id$


D = ../..
include $D/common/Make.global_options


default: makefiles.html makefile.small.html makefile.large.html makefile.rules.html


# Have a list of variables defined in common/Make.global_options and expand them into
# something that looks decent in HTML. The special text "--break--" draws a horizontal
# line to separate different blocks of variables
VARIABLES =     D DEAL_II_VERSION DEAL_II_MAJOR DEAL_II_MINOR \
	--break-- \
		CXX GXX-VERSION GXX-VERSION-DETAILED \
		CC CC-VERSION \
		F77 F77-VERSION \
		SHLIBLD \
		AR \
		RANLIB \
		PERL \
		TARGET \
	--break-- \
		LIBPATH \
		LIBDIR \
		LIBS \
		OBJEXT EXEEXT static-lib-suffix shared-lib-suffix \
		lib-suffix \
	--break-- \
		lib-deal2.o \
		lib-deal2.g \
	--break-- \
		INCLUDE \
		CXXFLAGS.g \
		CXXFLAGS.o \
		CFLAGS.g \
		CFLAGS.o \
		F77FLAGS.g \
		F77FLAGS.o \
		SHLIBFLAGS \
		LDFLAGS \
		F77LIBS \
	--break-- \
		enable-shared \
		enable-parser \
		DEAL_II_USE_MPI \
	--break-- \
		enable-threads \
		lib-contrib-tbb.g \
		lib-contrib-tbb.o \
	--break-- \
		USE_CONTRIB_UMFPACK \
	--break-- \
		USE_CONTRIB_HSL \
		lib-contrib-hsl \
	--break-- \
		USE_CONTRIB_PETSC \
		DEAL_II_PETSC_DIR \
		DEAL_II_PETSC_ARCH \
		DEAL_II_PETSC_VERSION_MAJOR    \
		DEAL_II_PETSC_VERSION_MINOR    \
		DEAL_II_PETSC_VERSION_SUBMINOR \
		DEAL_II_PETSC_MPIUNI_LIB \
		lib-contrib-petsc.g \
		lib-contrib-petsc.o \
		include-path-petsc \
		include-path-petsc-bmake \
	--break-- \
		USE_CONTRIB_SLEPC \
		DEAL_II_SLEPC_DIR \
		DEAL_II_SLEPC_VERSION_MAJOR    \
		DEAL_II_SLEPC_VERSION_MINOR    \
		DEAL_II_SLEPC_VERSION_SUBMINOR \
		lib-contrib-slepc \
		include-path-slepc \
	--break-- \
		USE_CONTRIB_TRILINOS \
		DEAL_II_TRILINOS_INCDIR \
		DEAL_II_TRILINOS_LIBDIR \
		DEAL_II_TRILINOS_SHARED \
		DEAL_II_TRILINOS_STATIC \
		DEAL_II_TRILINOS_VERSION_MAJOR    \
		DEAL_II_TRILINOS_VERSION_MINOR    \
		DEAL_II_TRILINOS_VERSION_SUBMINOR \
		lib-contrib-trilinos \
		include-path-trilinos \
	--break-- \
		USE_CONTRIB_MUMPS  \
		DEAL_II_MUMPS_DIR  \
		DEAL_II_SCALAPACK_DIR \
		DEAL_II_BLACS_DIR  \
		DEAL_II_BLACS_ARCH \
		lib-contrib-mumps \
		include-path-mumps \
	--break-- \
		USE_CONTRIB_ARPACK  \
		DEAL_II_ARPACK_DIR  \
		DEAL_II_ARPACK_ARCH \
		lib-contrib-arpack \
	--break-- \
		USE_CONTRIB_METIS \
		DEAL_II_METIS_DIR \
		lib-contrib-metis \
		include-path-contrib-metis \
	--break-- \
		USE_CONTRIB_P4EST \
		DEAL_II_P4EST_DIR \
	--break-- \
		TECIO_INCLUDE \
		TECIO_LIBRARY

expand_variable = $(if $(1:--break--=), \
			'<tr><td><code>$$($1)</td> <td>$($1)</code></td></tr>',\
			'<tr><th colspan="2"><hr></th></tr>')

VARIABLE_VALUES = $(foreach v,$(VARIABLES),$(call expand_variable,$v))

makefiles.html: Makefile makefiles.1.html $D/common/Make.global_options makefiles.2.html
	@cat makefiles.1.html > $@
	@echo "<p align=\"center\"><table class=\"tutorial\" width=\"60%\">" >> $@
	@echo $(VARIABLE_VALUES) >> $@
	@echo "</table></p>" >> $@
	@cat makefiles.2.html >> $@


makefile.small.html: Makefile.small
	@$(PERL) makefile2html < $< > $@

makefile.large.html: Makefile.large
	@$(PERL) makefile2html < $< > $@

makefile.rules.html: $D/common/Makefile.template
	@$(PERL) makefile2html < $< > $@

validate-xrefs:
	@$(PERL) $D/common/scripts/validate-xrefs.pl \
		$(filter-out makefiles.1.html makefiles.2.html, \
		  $(shell echo *.html))



clean:
	-rm -f makefiles.html makefile.small.html makefile.large.html

.PHONY: default validate-xrefs clean
