//----------------------------  mg_base.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_base.cc  ---------------------------


#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>
#include <iostream>
#include <cmath>


MGBase::~MGBase () 
{}


MGBase::MGBase(const MGTransferBase &transfer,
	       const unsigned        minlevel,
	       const unsigned        maxlevel)
		:
		maxlevel(maxlevel),
		minlevel(minlevel),
		defect(minlevel,maxlevel),
		solution(minlevel,maxlevel),
		transfer(&transfer)
{
  Assert(minlevel <= maxlevel,
	 ExcSwitchedLevels(minlevel, maxlevel));
}


void
MGBase::vcycle(const MGSmootherBase     &pre_smooth,
	       const MGSmootherBase     &post_smooth,
	       const MGCoarseGridSolver &coarse_grid_solver)
{
  level_mgstep (maxlevel, pre_smooth, post_smooth, coarse_grid_solver);
  abort ();
}


void
MGBase::level_mgstep(const unsigned int        level,
		     const MGSmootherBase     &pre_smooth,
		     const MGSmootherBase     &post_smooth,
		     const MGCoarseGridSolver &coarse_grid_solver)
{
  char *name = new char[100];

  sprintf(name, "MG%d-defect",level);
  print_vector(level, defect[level], name);
  
  solution[level] = 0.;
  
  if (level == minlevel)
    {
      coarse_grid_solver(level, solution[level], defect[level]);
      return;
    }

			   // smoothing of the residual by modifying s
  pre_smooth.smooth(level, solution[level], defect[level]);
				   // t = d-As

  sprintf(name, "MG%d-pre",level);
  print_vector(level, solution[level], name);
  
  t.reinit(solution[level].size());
  level_vmult(level, t, solution[level], defect[level]);
  
				   // make t rhs of lower level
//TODO: this function adds the restricted t to defect[level-1].
//TODO: why don't we have to clear it before?  
  transfer->restrict_and_add (level, defect[level-1], t);
  
				   // do recursion
  level_mgstep(level-1, pre_smooth, post_smooth, coarse_grid_solver);

				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_mgstep directly above
  t.reinit(solution[level].size());

				   // do coarse grid correction

  transfer->prolongate(level, t, solution[level-1]);

  sprintf(name, "MG%d-cgc",level);
  print_vector(level, t, name);

  solution[level] += t;
  
				   // smoothing (modify solution again)
//TODO: what happens here? smooth overwrites the solution[level],
//TODO: so the previous two statements should have no effect. No?  
  post_smooth.smooth(level, solution[level], defect[level]);

  sprintf(name, "MG%d-post",level);
  print_vector(level, solution[level], name);
}


//////////////////////////////////////////////////////////////////////

MGCoarseGridSolver::~MGCoarseGridSolver()
{};


//////////////////////////////////////////////////////////////////////

MGTransferBase::~MGTransferBase()
{};


