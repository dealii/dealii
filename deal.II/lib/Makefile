# $Id$
# Copyright W. Bangerth, University of Heidelberg, 2002

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
  lib-contrib-petsc-path.g = $(DEAL_II_PETSC_DIR)/lib/libg_c++/$(DEAL_II_PETSC_ARCH)
  lib-contrib-petsc-path.o = $(DEAL_II_PETSC_DIR)/lib/libO_c++/$(DEAL_II_PETSC_ARCH)
ifeq ($(DEAL_II_PETSC_VERSION_MINOR),1)
  xlib-ksp=libpetscsles
else
  xlib-ksp=libpetscksp
endif
  xlib.g = $(lib-contrib-petsc-path.g)/$(xlib-ksp)$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetscdm$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetscmat$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetscvec$(lib-suffix) \
           $(lib-contrib-petsc-path.g)/libpetsc$(lib-suffix)
  xlib.o = $(lib-contrib-petsc-path.o)/$(xlib-ksp)$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetscdm$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetscmat$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetscvec$(lib-suffix) \
           $(lib-contrib-petsc-path.o)/libpetsc$(lib-suffix)
  external-links-petsc:
	@echo Setting up libraries for PETSc
	@echo "======================debug============= Linking library:   libpetsc.g.so"
	@$(CXX) -shared -o libpetsc.g.so \
		-Wl,-rpath,$(lib-contrib-petsc-path.g) $(xlib.g)
	@echo "======================optimized========= Linking library:   libpetsc.so"
	@$(CXX) -shared -o libpetsc.so \
		-Wl,-rpath,$(lib-contrib-petsc-path.o) $(xlib.o)

else
  external-links-petsc:
endif



######################### CLEAN TARGETS ###################
# clean everything
clean: clean-objects clean-bin

# only clean object files, not libraries and executables
clean-objects:  clean-base clean-lac \
		clean-1d clean-2d clean-3d \
		clean-contrib


# clean targets for the individual subdirectories. remove
# object files and template instantiation files.
# extract dir name from target
clean-base clean-lac clean-1d clean-2d clean-3d:
	-cd $(@:clean-%=%) ; rm *.$(OBJEXT) *.ti *.ii

clean-contrib:
	-rm contrib/*/*.$(OBJEXT)

clean-bin:
	-rm bin/detached_ma27

distclean: clean
	-rm -f lib* bin/*

.PHONY: clean clean-objects 
.PHONY: clean-base clean-lac clean-1d clean-2d clean-3d
.PHONY: clean-contrib
.PHONY: external-links external-links-petsc
