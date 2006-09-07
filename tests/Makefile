# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

all: base lac fe deal.II multigrid bits hp
	@if grep 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
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

report:
	@for i in base lac fe deal.II multigrid bits hp all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	 done
	@if grep 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		echo =======Report: lapack ======= ; \
		cd lapack ; $(MAKE) report ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) report ; cd .. ; \
	fi

report+mail:
	@for i in base lac fe deal.II multigrid bits hp all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report+mail ; cd .. ; \
	 done
	@if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		echo =======Report: lapack ======= ; \
		cd lapack ; $(MAKE) report+mail ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
		cd petsc ; $(MAKE) report+mail ; cd .. ; \
	fi

# compilation of tests only, no execution of tests

build: build-base build-lac build-fe build-deal.II build-multigrid build-bits build-hp
	@if grep 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
		cd lapack ; $(MAKE) build ; cd .. ; \
	fi
	@if grep -q 'define DEAL_II_USE_PETSC' $D/base/include/base/config.h ; then \
		echo =======Report: petsc ======= ; \
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


# compilation and execution of tests only, no comparison with
# reference values

output: output-base output-lac output-fe output-deal.II output-bits output-hp
	@if grep 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
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



# update/generate references

refs: refs-base refs-lac refs-fe refs-deal.II refs-bits refs-hp
	@if grep 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ; then \
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

clean: clean-base clean-lac clean-deal.II clean-fe clean-bits clean-hp clean-lapack clean-petsc

distclean: distclean-base distclean-lac distclean-deal.II distclean-fe distclean-bits distclean-hp distclean-lapack distclean-petsc

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

clean-lapack:
	cd lapack ; $(MAKE) clean

distclean-lapack:
	cd lapack ; $(MAKE) distclean

clean-petsc:
	cd petsc ; $(MAKE) clean

distclean-petsc:
	cd petsc ; $(MAKE) distclean

.PHONY : all base lac fe deal.II multigrid bits \
	clean clean-base clean-lac clean-fe \
	clean-deal.II clean-multigrid clean-bits clean-hp \
	distclean distclean-base distclean-lac distclean-fe \
	distclean-deal.II distclean-multigrid distclean-bits distclean-hp

include Makefile.paths
include $D/common/Make.global_options

