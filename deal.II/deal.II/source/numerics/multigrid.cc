// $Id$
// Copyright Guido Kanschat, Universitaet Heidelberg, 1999

#include <numerics/multigrid.h>
#include <lac/vector.h>

void
MultiGridBase::level_mgstep(unsigned level)
{
  if (level == minlevel)
    {
      coarse_grid_solution(level, d[level], s[level]);
      return;
    }
  
  smooth(level, d[level], s[level], n_pre_smooth);

  post_smooth(level, d[level], s[level], n_post_smooth);
}


void MultiGridBase::post_smooth(unsigned level,
				Vector<float>& dst, const Vector<float>& src,
				unsigned steps)
{
  smooth(level, dst, src, steps);
}

// ab hier Wolfgang

