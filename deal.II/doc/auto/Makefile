# $Id$

KDOCFLAGS   = -I../scripts/kdoc ../scripts/kdoc/kdoc -a -p
KDOCINCLUDE = -lbase -llac -lbasic -lgrid -ldof -lnumerics -ltransport

kdoc.base = $(wildcard ../../base/include/base/*.h)
kdoc.lac  = $(wildcard ../../lac/include/lac/*.h)
kdoc.basic= $(wildcard ../../deal.II/include/basic/*.h)
kdoc.grid = $(wildcard ../../deal.II/include/grid/tria*.h) ../../deal.II/include/grid/geometry_info.h
kdoc.dof  = $(wildcard ../../deal.II/include/grid/dof*.h) \
            $(wildcard ../../deal.II/include/grid/mg_*.h) \
            $(wildcard ../../deal.II/include/fe/*.h)
kdoc.num  = $(wildcard ../../deal.II/include/numerics/*.h)

#user projects
kdoc.transport = $(wildcard ../../projects/transport/trans_dg/*.h)




kdoc.all  = $(kdoc.base) $(kdoc.lac) $(kdoc.basic) $(kdoc.grid) $(kdoc.dof) $(kdoc.num) \
            $(kdoc.ntrans)

all-include-files = $(wildcard ../../deal.II/include/*/*.h \
                               ../../lac/include/lac/*.h   \
                               ../../base/include/base/*.h \
                               ../../projects/transport/trans_dg/*.h \
                     )



default: kdoc kdoc/names.html
all: default cvslog

cvslog:
	cd ../.. ; perl doc/auto/scripts/cvs2html -o doc/auto/cvs-backlog/newdeal -a -k -D 100


# make sure that the *.kdoc files exist by
# using the dependancies and the following rules.
kdoc: kdoc/base.kdoc kdoc/lac.kdoc kdoc/basic.kdoc kdoc/grid.kdoc kdoc/dof.kdoc kdoc/numerics.kdoc $(kdoc.all) kdoc-rerun

kdoc-rerun:
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dbase \
		 base $(kdoc.base:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dlac \
		 lac $(kdoc.lac:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dbasic \
		 basic $(kdoc.basic:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dgrid  \
		grid  $(kdoc.grid:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -ddof   \
		dof   $(kdoc.dof:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dnumerics   \
		numerics $(kdoc.num:../..%=../../..%)
	@cd kdoc ; perl $(KDOCFLAGS) $(KDOCINCLUDE) -dtransport   \
		transport $(kdoc.transport:../..%=../../..%)
# output the names of the files which were not processed since they should
# be added to one of the above lists.
	@echo -------------------------------------------------------------
	@echo Unprocessed include files:
	@echo "  " $(filter-out $(kdoc.all), $(all-include-files))


kdoc/base.kdoc: $(kdoc.base)
	@cd kdoc ; perl $(KDOCFLAGS) -dbase \
		base $(kdoc.base:../..%=../../..%)

kdoc/lac.kdoc: $(kdoc.lac)
	@cd kdoc ; perl $(KDOCFLAGS) -dlac  \
		lac   $(kdoc.lac:../..%=../../..%)

kdoc/basic.kdoc: $(kdoc.basic)
	@cd kdoc ; perl $(KDOCFLAGS) -dbasic \
		basic $(kdoc.basic:../..%=../../..%)

kdoc/grid.kdoc: $(kdoc.grid)
	@cd kdoc ; perl $(KDOCFLAGS) -dgrid  \
		grid   $(kdoc.grid:../..%=../../..%)

kdoc/dof.kdoc: $(kdoc.tria)
	@cd kdoc ; perl $(KDOCFLAGS) -ddof   \
		dof     $(kdoc.dof:../..%=../../..%)

kdoc/numerics.kdoc: $(kdoc.num)
	@cd kdoc ; perl $(KDOCFLAGS) -dnumerics   \
		numerics $(kdoc.num:../..%=../../..%)

kdoc/transport.kdoc: $(kdoc.transport)
	@cd kdoc ; perl $(KDOCFLAGS) -dtransport   \
		transport $(kdoc.transport:../..%=../../..%)


kdoc/names.html: kdoc/base.kdoc kdoc/lac.kdoc kdoc/basic.kdoc kdoc/grid.kdoc kdoc/dof.kdoc kdoc/numerics.kdoc kdoc/transport.kdoc
	@cd kdoc ; perl ../scripts/index.pl *.kdoc > names.html


clean:
	-rm -f kdoc/names.html kdoc/*.kdoc kdoc/*/*.html cvs-backlog/*.html



.PHONY: default all kdoc-rerun clean