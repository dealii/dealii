// $Id$
// Copyright Guido Kanschat, Universitaet Heidelberg, 1999

#include <numerics/multigrid.h>

void
MGMatrix::vmult(dVector& dst, const dVector& src) const
{
  dst = 0.;
  
  copy_to_mg(s,src);
  
  for (int l=0;l<maxlevel;++l)
  {
    level_active_vmult(l,d,s);
  }
  copy_from_mg(dst,d);
}


void
MGMatrix::precondition(dVector& dst, const dVector& src) const
{
  copy_to_mg(s,src);
  d = 0.;
  
    
}

void
MGMatrix::level_mgstep(int level)
{
  for (int i=0; i< n_presmooth; ++i)
  {
    smooth(level, d, s);
  }
}

