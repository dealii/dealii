# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998


# list the directories we want to visit
subdirs = grid/ dof/ poisson/ convergence/ error-estimation/ multigrid/

# define lists of targets: for each directory we produce a target name
# for compilation, running and cleaning by appending the action to
# the directory name (replacing the slash by ".action")
compile = $(subdirs:/=.compile)
run     = $(subdirs:/=.run)
clean   = $(subdirs:/=.clean)

# define global targets which are to be excuted in every subdirectory
compile: $(compile)
run    : $(run)
# for cleaning up: do this also for the present directory
clean  : $(clean)
	-rm -f *~


# define the action of the targets for the specific subdirectories
$(compile) :
	cd $(@:.compile=) ; $(MAKE) D=$D

$(run) :
	cd $(@:.run=) ; $(MAKE) run D=$D

$(clean) :
	-cd $(@:.clean=) ; $(MAKE) clean D=$D


