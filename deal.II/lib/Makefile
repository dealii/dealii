# $Id$

D = ..

# get options which hold for all files of the project
include $D/common/Make.global_options


default:
	@echo There is no default target in this directory


######################### EXTERNAL LINK TARGETS ###################
# link in external libs, if any
external-links: external-links-petsc

ifeq ($(USE_CONTRIB_PETSC),yes)

  # PETSc has the nasty habit of putting debug and optimized builds into
  # different directories, but naming the respective libraries the same.
  # Since one doesn't want to always set and reset $LD_LIBRARY_PATH when
  # one switches from debug to optimized mode, they pass -Wl,-rpath to
  # the linker line. This is not really portable, and it has the additional
  # drawback that one has to pass different linker flags for debug and
  # optimized mode. The latter is a problem with our present build system,
  # which allows to use different libraries, but doesn't take into account
  # that one may want to pass different $LDFLAGS.
  #
  # To avoid such hassle, we build some kind of convenience library: a shared
  # library that has no content, but contains references to all the different
  # PETSc libraries, and we build one for debug and one for optimized mode.
  # We can't really get around the use of -Wl,-rpath, but at least we
  # have now centralized this part to one single place, and are back to where
  # we want to be: have different libraries for debug and optimized mode,
  # but in the same directory. And, voila, our build system works again :-)
  #
  # There is one gotcha: PETSc has a habit of making incompatible changes by
  # simply changing the names of libraries between different versions. That
  # really is too bad, so we have to conditionalize things.
ifeq ($(DEAL_II_PETSC_VERSION_MINOR),1)
  xlib-ksp=libpetscsles
else
  xlib-ksp=libpetscksp
endif
  # This is how it is done. Starting with PETSc versions 2.3.*:
ifeq ($(DEAL_II_PETSC_VERSION_MAJOR),2)
  xlib.g = $(lib-contrib-petsc-path.g)/libpetscksp$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetscdm$(lib-suffix)  \
           $(lib-contrib-petsc-path.g)/libpetscmat$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetscvec$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetsc$(lib-suffix)
  xlib.o = $(lib-contrib-petsc-path.o)/libpetscksp$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetscdm$(lib-suffix)  \
           $(lib-contrib-petsc-path.o)/libpetscmat$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetscvec$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetsc$(lib-suffix)
else
  # which is the same for PETSc 3.0.0:
  ifeq ($(DEAL_II_PETSC_VERSION_MAJOR)$(DEAL_II_PETSC_VERSION_MINOR),30)
    xlib.g = $(lib-contrib-petsc-path.g)/libpetscksp$(lib-suffix) \
             $(lib-contrib-petsc-path.g)/libpetscdm$(lib-suffix)  \
             $(lib-contrib-petsc-path.g)/libpetscmat$(lib-suffix) \
             $(lib-contrib-petsc-path.g)/libpetscvec$(lib-suffix) \
             $(lib-contrib-petsc-path.g)/libpetsc$(lib-suffix)
    xlib.o = $(lib-contrib-petsc-path.o)/libpetscksp$(lib-suffix) \
             $(lib-contrib-petsc-path.o)/libpetscdm$(lib-suffix)  \
             $(lib-contrib-petsc-path.o)/libpetscmat$(lib-suffix) \
             $(lib-contrib-petsc-path.o)/libpetscvec$(lib-suffix) \
             $(lib-contrib-petsc-path.o)/libpetsc$(lib-suffix)
  else
    # but after that (petsc-3.1++), we can use the simpler PETSc
    # default "--with-single-library=1" like this:
    xlib.g = $(lib-contrib-petsc-path.o)/libpetsc$(lib-suffix)
    xlib.o = $(lib-contrib-petsc-path.o)/libpetsc$(lib-suffix)
  endif
endif

ifeq ($(enable-shared),yes)

  external-links-petsc: libpetscall.g$(shared-lib-suffix) libpetscall$(shared-lib-suffix)

  libpetscall.g$(shared-lib-suffix): $(xlib.g)
	@echo "======================debug============= Linking library:   libpetscall.g$(shared-lib-suffix)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,petscall.g) $(call DEAL_II_ADD_SONAME,petscall.g) $(xlib.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,petscall.g) $@

  libpetscall$(shared-lib-suffix): $(xlib.o)
	@echo "======================optimized========= Linking library:   libpetscall$(shared-lib-suffix)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,petscall) $(call DEAL_II_ADD_SONAME,petscall)  $(xlib.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,petscall) $@

else

  # make static libraries. ensure that they are built one after the other because
  # both create temporary files with the same name
  external-links-petsc:
	@$(MAKE) libpetscall.g$(static-lib-suffix)
	@$(MAKE) libpetscall$(static-lib-suffix)

  # As I don't know how to merge archives, here's a workaround: Extract each of
  # the archives and append the object files to the libpetsc archive. Do the
  # same for debug and optimized mode.
  libpetscall.g$(static-lib-suffix): $(xlib.g)
	@echo "======================debug============= Linking library:   libpetscall.g$(static-lib-suffix)"
	@rm -f $@
	@rm -rf .ar_tmp
	@mkdir .ar_tmp; cd .ar_tmp; \
	 for arfile in $(xlib.g) ; do \
	   $(AR) x $$arfile; \
           $(AR) q ../libpetscall.g$(static-lib-suffix) *; \
           rm *; \
	 done
	@cd ..
	@rm -r .ar_tmp
	@$(RANLIB) libpetscall.g$(static-lib-suffix)

  libpetscall$(static-lib-suffix): $(xlib.o)
	@echo "======================optimized========= Linking library:   libpetscall$(static-lib-suffix)"
	@rm -f $@
	@rm -rf .ar_tmp
	@mkdir .ar_tmp; cd .ar_tmp; \
	 for arfile in $(xlib.o) ; do \
	   $(AR) x $$arfile; \
           $(AR) q ../libpetscall$(static-lib-suffix) *; \
           rm *; \
	 done
	@cd ..
	@rm -r .ar_tmp
	@$(RANLIB) libpetscall$(static-lib-suffix)
endif


else
  external-links-petsc:
endif



######################### CLEAN TARGETS ###################
# clean everything
clean: clean-objects

# only clean object files, not libraries and executables
clean-objects:  clean-debug clean-optimized \
		clean-contrib


# clean targets for the individual subdirectories. remove
# object files and template instantiation files.
# extract dir name from target
clean-debug clean-optimized:
	-cd $(@:clean-%=%) ; rm -f *.$(OBJEXT) *.ti *.ii

clean-contrib:
	-rm contrib/*/*.$(OBJEXT)
	-rm -r contrib/tbb/*debug
	-rm -r contrib/tbb/*release

distclean: clean
	-rm -f lib* bin/*
	-rm bin/detached_ma27
	-rm bin/dealii_parameter_gui
	-rm libtbb*

.PHONY: clean clean-objects
.PHONY: clean-base clean-lac clean-1d clean-2d clean-3d
.PHONY: clean-contrib
.PHONY: external-links external-links-petsc
