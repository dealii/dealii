// $Id$
// Copyright Guido Kanschat, Universitaet Heidelberg, 1999

#include <numerics/multigrid.h>
#include <lac/vector.h>

MultiGridBase::MultiGridBase(MGTransferBase& transfer,
			     unsigned maxlevel, unsigned minlevel,
			     unsigned pre_smooth, unsigned post_smooth)
		:
		maxlevel(maxlevel), minlevel(minlevel),
		transfer(&transfer),
		n_pre_smooth(pre_smooth), n_post_smooth(post_smooth)
{}

void
MultiGridBase::level_mgstep(unsigned level)
{
  if (level == minlevel)
    {
      coarse_grid_solution(level, d[level], s[level]);
      return;
    }
  
				   // smoothing of the residual by modifying s
  smooth(level, d[level], s[level], n_pre_smooth);

				   // t = d-As
  level_residual(level, t, s[level], d[level]);
  
				   // make t rhs of lower level
  transfer->restrict(level, t, d[level-1]);
				   // do recursion
  level_mgstep(level-1);
				   // do coarse grid correction
  transfer->prolongate(level, s[level], s[level-1]);

				   // smoothing (modify s again)
  post_smooth(level, d[level], s[level], n_post_smooth);
}


void
MultiGridBase::post_smooth(unsigned level,
			   Vector<float>& dst, const Vector<float>& src,
			   unsigned steps)
{
  smooth(level, dst, src, steps);
}

void
MultiGridBase::coarse_grid_solution(unsigned level,
			  Vector<float>& dst,
			  const Vector<float>& src)
{
  smooth(level, dst, src, 10 * (n_pre_smooth + n_post_smooth));
}

// ab hier Wolfgang

