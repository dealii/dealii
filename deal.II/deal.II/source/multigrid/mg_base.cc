// $Id$

#include <lac/mgbase.h>
#include <iostream>
#include <cmath>

template<class VECTOR>
void print_vector(ostream& s, const VECTOR& v)
{
  const unsigned int n = (unsigned int)(sqrt(v.size())+.3);
  unsigned int k=0;
  
  for (unsigned int i=0;i<n;++i)
    {
      for (unsigned int j=0;j<n;++j)
//	s << ((float)i)/n << '\t' << ((float)j)/n << '\t' << v(k++) << endl;
	s << '\t' << v(k++);
      s << endl;
    }
  s << endl;
}

MGBase::~MGBase () 
{};



MGBase::MGBase(const MGTransferBase& transfer,
			     unsigned minlevel, unsigned maxlevel)
		:
		maxlevel(maxlevel), minlevel(minlevel),
		d(minlevel,maxlevel),
		s(minlevel,maxlevel),
		transfer(&transfer)
{
  Assert(minlevel <= maxlevel, ExcSwitchedLevels(minlevel, maxlevel));
}

void
MGBase::vcycle(const MGSmootherBase& pre_smooth,
	       const MGSmootherBase& post_smooth,
	       const MGCoarseGridSolver& coarse_grid_solver)
{
  level_mgstep(maxlevel, pre_smooth, post_smooth, coarse_grid_solver);
}


void
MGBase::level_mgstep(unsigned level,
		     const MGSmootherBase& pre_smooth,
		     const MGSmootherBase& post_smooth,
		     const MGCoarseGridSolver& coarse_grid_solver)
{
  s[level] = 0.;
  
  if (level == minlevel)
    {
      coarse_grid_solver(level, s[level], d[level]);
      return;
    }

			   // smoothing of the residual by modifying s
  pre_smooth.smooth(level, s[level], d[level]);
				   // t = d-As
  t.reinit(s[level].size());
  level_vmult(level, t, s[level], d[level]);

				   // make t rhs of lower level
  transfer->restrict(level, d[level-1], t);
  
				   // do recursion
  level_mgstep(level-1, pre_smooth, post_smooth, coarse_grid_solver);
				   // do coarse grid correction
  t.reinit(s[level].size());
  transfer->prolongate(level, t, s[level-1]);
  s[level].add(t);
  
				   // smoothing (modify s again)
  post_smooth.smooth(level, s[level], d[level]);
}


//////////////////////////////////////////////////////////////////////

MGCoarseGridSolver::~MGCoarseGridSolver()
{}


//////////////////////////////////////////////////////////////////////

MGTransferBase::~MGTransferBase()
{}

MGSmootherBase::~MGSmootherBase()
{}


//////////////////////////////////////////////////////////////////////


void
MGSmootherIdentity::smooth (const unsigned int,
			    Vector<double>       &,
			    const Vector<double> &) const
{}


