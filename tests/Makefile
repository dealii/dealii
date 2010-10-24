# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

include Makefile.paths
include $D/common/Make.global_options

# first define the directories we should work on
nofail-dirs = a-framework       \
	      base              \
              lac               \
              fe                \
              deal.II           \
              multigrid         \
              bits              \
              hp                \
              aniso     \
              codim_one \
              serialization \
	      $(shell if test "$(DEAL_II_USE_MPI)" = yes ; then echo mpi ; fi) \
	      $(shell if test "$(USE_CONTRIB_P4EST)" = yes ; then echo distributed_grids ; fi) \
	      $(shell if test "$(USE_CONTRIB_PETSC)" = yes ; then echo petsc ; fi) \
	      $(shell if test "$(USE_CONTRIB_TRILINOS)" = yes ; then echo trilinos ; fi) \
	      $(shell if test "$(USE_CONTRIB_UMFPACK)" = yes ; then echo umfpack ; fi) \
	      $(shell if test "$(USE_CONTRIB_METIS)" = yes ; then echo metis ; fi) \
	      $(shell if test "$(USE_CONTRIB_HSL)" = yes ; then echo hsl ; fi) \
	      $(shell if grep -q 'define HAVE_LIBLAPACK' $D/base/include/base/config.h ;\
                then echo lapack ; fi)

# now define the targets that actually do something
all: nofail do-fail

nofail: $(addprefix do-,$(nofail-dirs))

do-%:
	cd $(@:do-%=%) && $(MAKE)


report:
	-@$(MAKE) -k nofail
	@for i in $(nofail-dirs) fail all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	 done

report+mail:
	-@$(MAKE) -k nofail
	@for i in $(nofail-dirs) fail all-headers ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report+mail ; cd .. ; \
	 done

# compilation of tests only, no execution of tests
build: $(addprefix build-,$(nofail-dirs) fail)
build-%:
	cd $(@:build-%=%) && $(MAKE) build


# compilation and execution of tests only, no comparison with
# reference values
output: $(addprefix output-,$(nofail-dirs) fail)
output-%:
	cd $(@:output-%=%) && $(MAKE) output


# update/generate references
refs: $(addprefix output-,$(nofail-dirs))
refs-%:
	cd $(@:refs-%=%) && $(MAKE) refs


# clean up output files
clean:     $(addprefix clean-,$(nofail-dirs) fail)
distclean: $(addprefix distclean-,$(nofail-dirs) fail)
clean-%:
	cd $(@:clean-%=%) && $(MAKE) clean
distclean-%:
	cd $(@:distclean-%=%) && $(MAKE) distclean



resultname:
	@echo $(TARGET)+$(GXX-VERSION)


.PHONY : all do-% report report+mail build build-% refs refs-% \
         clean clean-% \
	 distclean distclean-% \
	 resultname


