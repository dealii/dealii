# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

all: base lac fe deal.II multigrid bits hp fail
	@-if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		cd umfpack ; $(MAKE) ; cd .. ; \
	fi
	@-if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) ; cd .. ; \
	fi
	@-if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		cd petsc ; $(MAKE) ; cd .. ; \
	fi

base:
	cd base ; $(MAKE)

lac:
	cd lac ; $(MAKE)

fe:
	cd fe ; $(MAKE)

bits:
	cd bits ; $(MAKE)

deal.II:
	cd deal.II ; $(MAKE)

multigrid:
	cd multigrid ; $(MAKE)

hp:
	cd hp ; $(MAKE)

fail:
	cd fail ; $(MAKE)

report:
	@for i in base lac fe deal.II multigrid bits hp fail all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	 done
	@if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		echo =======Report: umfpack ======= ; \
		cd umfpack ; $(MAKE) report ; cd .. ; \
	fi
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		echo =======Report: lapack ======= ; \
		cd lapack ; $(MAKE) report ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) report ; cd .. ; \
	fi

report+mail:
	@for i in base lac fe deal.II multigrid bits hp fail all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report+mail ; cd .. ; \
	 done
	@if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		echo =======Report: umfpack ======= ; \
		cd umfpack ; $(MAKE) report+mail ; cd .. ; \
	fi
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		echo =======Report: lapack ======= ; \
		cd lapack ; $(MAKE) report+mail ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) report+mail ; cd .. ; \
	fi

# compilation of tests only, no execution of tests

build: build-base build-lac build-fe build-deal.II build-multigrid build-bits build-hp build-fail
	@if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		cd umfpack ; $(MAKE) build ; cd .. ; \
	fi
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) build ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		cd petsc ; $(MAKE) build ; cd .. ; \
	fi

build-base:
	cd base ; $(MAKE) build

build-lac:
	cd lac ; $(MAKE) build

build-fe:
	cd fe ; $(MAKE) build

build-bits:
	cd bits ; $(MAKE) build

build-deal.II:
	cd deal.II ; $(MAKE) build

build-multigrid:
	cd multigrid ; $(MAKE) build

build-hp:
	cd hp ; $(MAKE) build

build-fail:
	cd fail ; $(MAKE) build


# compilation and execution of tests only, no comparison with
# reference values

output: output-base output-lac output-fe output-deal.II output-bits output-hp output-fail
	@if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		cd umfpack ; $(MAKE) output ; cd .. ; \
	fi
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) output ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) output ; cd .. ; \
	fi

output-base:
	cd base ; $(MAKE) output

output-lac:
	cd lac ; $(MAKE) output

output-fe:
	cd fe ; $(MAKE) output

output-bits:
	cd bits ; $(MAKE) output

output-deal.II:
	cd deal.II ; $(MAKE) output

output-multigrid:
	cd multigrid ; $(MAKE) output

output-hp:
	cd hp ; $(MAKE) output

output-fail:
	cd fail ; $(MAKE) output



# update/generate references

refs: refs-base refs-lac refs-fe refs-deal.II refs-bits refs-hp
	@if grep -q 'define HAVE_LIBUMFPACK' $D/base/include/base/config.h ; then \
		cd umfpack ; $(MAKE) refs ; cd .. ; \
	fi
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) refs ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) refs ; cd .. ; \
	fi

refs-base:
	cd base ; $(MAKE) refs

refs-lac:
	cd lac ; $(MAKE) refs

refs-fe:
	cd fe ; $(MAKE) refs

refs-bits:
	cd bits ; $(MAKE) refs

refs-deal.II:
	cd deal.II ; $(MAKE) refs

refs-multigrid:
	cd multigrid ; $(MAKE) refs

refs-hp:
	cd hp ; $(MAKE) refs



resultname:
	@echo $(TARGET)+$(GXX-VERSION)

clean: clean-base clean-lac clean-deal.II clean-fe clean-bits clean-hp clean-fail \
	clean-umfpack clean-lapack clean-petsc clean-multigrid

distclean: distclean-base distclean-lac distclean-deal.II distclean-fe distclean-bits distclean-hp distclean-fail \
	distclean-umfpack distclean-lapack distclean-petsc distclean-multigrid

clean-base:
	cd base ; $(MAKE) clean

distclean-base:
	cd base ; $(MAKE) distclean

clean-lac:
	cd lac ; $(MAKE) clean

distclean-lac:
	cd lac ; $(MAKE) distclean

clean-fe:
	cd fe ; $(MAKE) clean

distclean-fe:
	cd fe ; $(MAKE) distclean

clean-deal.II:
	cd deal.II ; $(MAKE) clean

distclean-deal.II:
	cd deal.II ; $(MAKE) distclean

clean-multigrid:
	cd multigrid ; $(MAKE) clean

distclean-multigrid:
	cd multigrid ; $(MAKE) distclean

clean-bits:
	cd bits ; $(MAKE) clean

distclean-bits:
	cd bits ; $(MAKE) distclean

clean-hp:
	cd hp ; $(MAKE) clean

distclean-hp:
	cd hp ; $(MAKE) distclean

clean-fail:
	cd fail ; $(MAKE) clean

distclean-fail:
	cd fail ; $(MAKE) distclean

clean-umfpack:
	cd umfpack ; $(MAKE) clean

distclean-umfpack:
	cd umfpack ; $(MAKE) distclean

clean-lapack:
	cd lapack ; $(MAKE) clean

distclean-lapack:
	cd lapack ; $(MAKE) distclean

clean-petsc:
	cd petsc ; $(MAKE) clean

distclean-petsc:
	cd petsc ; $(MAKE) distclean

.PHONY : all base lac fe deal.II multigrid bits hp fail\
	clean clean-base clean-lac clean-fe \
	clean-deal.II clean-multigrid clean-bits clean-hp clean-fail \
	distclean distclean-base distclean-lac distclean-fe \
	distclean-deal.II distclean-multigrid distclean-bits distclean-hp distclean-fail

include Makefile.paths
include $D/common/Make.global_options

