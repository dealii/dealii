// $Id$

#include <lac/mgbase.h>

MGBase::~MGBase () 
{};



MGBase::MGBase(const MGTransferBase& transfer,
			     unsigned maxlevel, unsigned minlevel)
		:
		d(maxlevel-minlevel),
		s(maxlevel-minlevel),
		transfer(&transfer),
		maxlevel(maxlevel), minlevel(minlevel)
{}

void
MGBase::level_mgstep(unsigned level,
		     const MGSmootherBase& pre_smooth,
		     const MGSmootherBase& post_smooth,
		     const MGCoarseGridSolver& coarse_grid_solver)
{
  if (level == minlevel)
    {
      coarse_grid_solver(level, s[level], d[level]);
      return;
    }
  
				   // smoothing of the residual by modifying s
  pre_smooth.smooth(level, s[level], d[level]);

				   // t = d-As
  level_residual(level, t, s[level], d[level]);
  
				   // make t rhs of lower level
  transfer->restrict(level, t, d[level-1]);
				   // do recursion
  level_mgstep(level-1, pre_smooth, post_smooth, coarse_grid_solver);
				   // do coarse grid correction
  transfer->prolongate(level, s[level], s[level-1]);

				   // smoothing (modify s again)
  post_smooth.smooth(level, s[level], d[level]);
}

MGTransferBase::~MGTransferBase()
{}

MGSmootherBase::~MGSmootherBase()
{}
