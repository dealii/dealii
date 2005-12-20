# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

all: base lac fe deal.II multigrid bits

base: abort.o
	cd base ; $(MAKE)

lac: abort.o
	cd lac ; $(MAKE)

fe: abort.o
	cd fe ; $(MAKE)

bits: abort.o
	cd bits ; $(MAKE)

deal.II: abort.o
	cd deal.II ; $(MAKE)

multigrid: abort.o
	cd multigrid ; $(MAKE)

hp: abort.o
	cd hp ; $(MAKE)

report: abort.o
	@for i in base lac fe deal.II multigrid bits all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	 done

report+mail: abort.o
	@for i in base lac fe deal.II multigrid bits all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report+mail ; cd .. ; \
	 done


# compilation of tests only, no execution of tests

build: build-base build-lac build-fe build-deal.II build-multigrid build-bits

build-base: abort.o
	cd base ; $(MAKE) build

build-lac: abort.o
	cd lac ; $(MAKE) build

build-fe: abort.o
	cd fe ; $(MAKE) build

build-bits: abort.o
	cd bits ; $(MAKE) build

build-deal.II: abort.o
	cd deal.II ; $(MAKE) build

build-multigrid: abort.o
	cd multigrid ; $(MAKE) build


# compilation and execution of tests only, no comparison with
# reference values

output: output-base output-lac output-fe output-deal.II output-bits

output-base: abort.o
	cd base ; $(MAKE) output

output-lac: abort.o
	cd lac ; $(MAKE) output

output-fe: abort.o
	cd fe ; $(MAKE) output

output-bits: abort.o
	cd bits ; $(MAKE) output

output-deal.II: abort.o
	cd deal.II ; $(MAKE) output

output-multigrid: abort.o
	cd multigrid ; $(MAKE) output



# update/generate references

refs: refs-base refs-lac refs-fe refs-deal.II refs-bits

refs-base: abort.o
	cd base ; $(MAKE) refs

refs-lac: abort.o
	cd lac ; $(MAKE) refs

refs-fe: abort.o
	cd fe ; $(MAKE) refs

refs-bits: abort.o
	cd bits ; $(MAKE) refs

refs-deal.II: abort.o
	cd deal.II ; $(MAKE) refs

refs-multigrid: abort.o
	cd multigrid ; $(MAKE) refs



resultname:
	@echo $(TARGET)+$(GXX-VERSION)

clean: clean-base clean-lac clean-deal.II clean-fe clean-bits
	-rm -f abort.o

distclean: distclean-base distclean-lac distclean-deal.II distclean-fe distclean-bits
	-rm -f abort.o

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

abort.o: abort.cc
	 @$(CXX) $(CXXFLAGS.g) -o $@ -c $<

.PHONY : all base lac fe deal.II multigrid bits \
	clean clean-base clean-lac clean-fe clean-deal.II clean-multigrid clean-bits \
	distclean distclean-base distclean-lac distclean-fe \
	distclean-deal.II distclean-multigrid distclean-bits

include Makefile.paths
include $D/common/Make.global_options

