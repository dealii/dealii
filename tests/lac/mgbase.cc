// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <lac/mgbase.h>

class TransferTest
  :
  public MGTransferBase
{
    virtual void prolongate(unsigned l,
			    Vector<float>& dst,
			    const Vector<float>& src) const;
    virtual void restrict(unsigned l,
			  Vector<float>& dst,
			  const Vector<float>& src) const;
    friend class MGBase;
};

class SmoothTest
  :
  public MGSmootherBase
{
  public:
        virtual void smooth (const unsigned int  level,
			     Vector<float>      &u,
			     const Vector<float>& rhs) const;
};

class MGCoarseGridTest
  :
  public MGCoarseGridSolver
{
  public:
    virtual void operator () (unsigned int level,
			      Vector<float>& dst,
			      const Vector<float>& src) const ;
};


class MGTest
  :
  public MGBase
{
    MGSmootherBase& smoother;
    MGCoarseGridSolver& solver;

  public:
    MGTest(TransferTest& t, SmoothTest& sm, MGCoarseGridTest& cs)
		    : MGBase(t, 9, 3),
		      smoother(sm),
		      solver(cs)
      {
	for (unsigned int l = minlevel; l <= maxlevel; ++l)
	  {
	    d[l].reinit(1);
	    s[l].reinit(1);
	  }
      }
    
    virtual void level_residual(unsigned level,
				Vector<float>& dst,
				const Vector<float>& src,
				const Vector<float>& rhs);
    
    void doit()
      {
	level_mgstep(9, smoother, smoother, solver);
      }
};

main()
{
  TransferTest tr;
  SmoothTest s;
  MGCoarseGridTest c;
  
  MGTest mg(tr, s, c);
  mg.doit();
  
}

void
TransferTest::prolongate(unsigned l,
			 Vector<float>&,
			 const Vector<float>&) const
{
  cout << "Prolongating " << l-1 << " to " << l << endl;
}

void
TransferTest::restrict(unsigned l,
		       Vector<float>&,
		       const Vector<float>&) const
{
  cout << "Restricting  " << l << " to " << l-1 << endl;
}

void
SmoothTest::smooth (const unsigned int  level,
		    Vector<float>      &,
		    const Vector<float>&) const
{
  cout << "Smoothing on " << level << endl;
}

void
MGTest::level_residual(unsigned l,
		       Vector<float>&,
		       const Vector<float>&,
		       const Vector<float>&)
{
  cout << "Residual on  " << l << endl;
}

void
MGCoarseGridTest::operator() (unsigned int level,
			      Vector<float>&,
			      const Vector<float>&) const
{
  cout << "Solving on   " << level << endl;
}
