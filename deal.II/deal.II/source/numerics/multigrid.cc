// $Id$
// Copyright Guido Kanschat, Universitaet Heidelberg, 1999

#include <numerics/multigrid.h>
#include <lac/dvector.h>

void
MultiGrid::vmult(dVector& dst, const dVector& src) const
{
  dst = 0.;
  
  copy_to_mg(s,src);
  
  for (unsigned l=0;l<maxlevel;++l)
  {
    level_active_vmult(l,d[l],s[l]);
  }
  copy_from_mg(dst,d);
}


void
MultiGrid::precondition(dVector& dst, const dVector& src) const
{
  copy_to_mg(s,src);
  copy_to_mg(d,dst);
  level_mgstep(maxlevel);
}


void
MultiGrid::level_mgstep(unsigned level) const
{
  if (level == minlevel)
  {
    coarse_grid_solution(level, d[level], s[level]);
    return;
  }
  
  smooth(level, d[level], s[level], n_pre_smooth);

  post_smooth(level, d[level], s[level], n_post_smooth);
}


void MultiGrid::post_smooth(unsigned level,
			    dVector& dst, const dVector& src,
			    unsigned steps) const
{
  smooth(level, dst, src, steps);
}
