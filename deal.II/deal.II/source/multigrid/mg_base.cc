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


//TODO: this function is only for debugging purposes and should be removed sometimes
template<class VECTOR>
static
void print_vector(ostream& s, const VECTOR& v)
{
  const unsigned int n = (unsigned int)(sqrt(v.size())+.3);
  unsigned int k=0;

				   // write the vector entries in a
				   // kind of square
  for (unsigned int i=0;i<n;++i)
    {
      for (unsigned int j=0;j<n;++j)
	s << '\t' << v(k++);
      s << endl;
    }
  s << endl;
}


//////////////////////////////////////////////////////////////////////


MGBase::~MGBase () 
{};


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
};


void
MGBase::vcycle(const MGSmootherBase     &pre_smooth,
	       const MGSmootherBase     &post_smooth,
	       const MGCoarseGridSolver &coarse_grid_solver)
{
  level_mgstep(maxlevel, pre_smooth, post_smooth, coarse_grid_solver);
};


void
MGBase::level_mgstep(const unsigned int        level,
		     const MGSmootherBase     &pre_smooth,
		     const MGSmootherBase     &post_smooth,
		     const MGCoarseGridSolver &coarse_grid_solver)
{
  solution[level] = 0.;
  
  if (level == minlevel)
    {
      coarse_grid_solver(level, solution[level], defect[level]);
      return;
    }

			   // smoothing of the residual by modifying s
  pre_smooth.smooth(level, solution[level], defect[level]);
				   // t = d-As
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
  solution[level] += t;
  
				   // smoothing (modify solution again)
//TODO: what happens here? smooth overwrites the solution[level],
//TODO: so the previous two statements should have no effect. No?  
  post_smooth.smooth(level, solution[level], defect[level]);
};


//////////////////////////////////////////////////////////////////////

MGCoarseGridSolver::~MGCoarseGridSolver()
{};


//////////////////////////////////////////////////////////////////////

MGTransferBase::~MGTransferBase()
{};


