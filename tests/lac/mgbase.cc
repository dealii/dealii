// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <lac/mgbase.h>
#include <base/logstream.h>

#include <fstream>

class TransferTest
  :
  public MGTransferBase
{
    virtual void prolongate(unsigned l,
			    Vector<double>& dst,
			    const Vector<double>& src) const;
    virtual void restrict(unsigned l,
			  Vector<double>& dst,
			  const Vector<double>& src) const;
    friend class MGBase;
};

class SmoothTest
  :
  public MGSmootherBase
{
  public:
        virtual void smooth (const unsigned int  level,
			     Vector<double>      &u,
			     const Vector<double>& rhs) const;
};

class MGCoarseGridTest
  :
  public MGCoarseGridSolver
{
  public:
    virtual void operator () (unsigned int level,
			      Vector<double>& dst,
			      const Vector<double>& src) const ;
};


class MGTest
  :
  public MGBase
{
    MGSmootherBase& smoother;
    MGCoarseGridSolver& solver;

  public:
    MGTest(TransferTest& t, SmoothTest& sm, MGCoarseGridTest& cs)
		    : MGBase(t, 3, 9),
		      smoother(sm),
		      solver(cs)
      {
	for (unsigned int l = minlevel; l <= maxlevel; ++l)
	  {
	    d[l].reinit(1);
	    s[l].reinit(1);
	  }
      }
    
    virtual void level_vmult(unsigned level,
				Vector<double>& dst,
				const Vector<double>& src,
				const Vector<double>& rhs);
    
    void doit()
      {
	level_mgstep(9, smoother, smoother, solver);
      }
};

int main()
{
  ofstream logfile("mgbase.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  TransferTest tr;
  SmoothTest s;
  MGCoarseGridTest c;
  
  MGTest mg(tr, s, c);
  mg.doit();
  
}

void
TransferTest::prolongate(unsigned l,
			 Vector<double>&,
			 const Vector<double>&) const
{
  deallog << "Prolongating " << l-1 << " to " << l << endl;
}

void
TransferTest::restrict(unsigned l,
		       Vector<double>&,
		       const Vector<double>&) const
{
  deallog << "Restricting  " << l << " to " << l-1 << endl;
}

void
SmoothTest::smooth (const unsigned int  level,
		    Vector<double>      &,
		    const Vector<double>&) const
{
  deallog << "Smoothing on " << level << endl;
}

void
MGTest::level_vmult(unsigned l,
		       Vector<double>&,
		       const Vector<double>&,
		       const Vector<double>&)
{
  deallog << "Residual on  " << l << endl;
}

void
MGCoarseGridTest::operator() (unsigned int level,
			      Vector<double>&,
			      const Vector<double>&) const
{
  deallog << "Solving on   " << level << endl;
}
