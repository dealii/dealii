// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <numerics/multigrid.h>
#include <numerics/mg_smoother.h>

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
    friend class MultiGridBase;
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


class MGTest
  :
  public MultiGridBase
{
    MGSmootherBase& smoother;
    
  public:
    MGTest(TransferTest& t, SmoothTest& s)
		    : MultiGridBase(t, 9, 3),
		      smoother(s)
      {}
    
    virtual void level_residual(unsigned level,
				Vector<float>& dst,
				const Vector<float>& src,
				const Vector<float>& rhs);
    
    void doit()
      {
	level_mgstep(9, smoother, smoother);
      }
};

main()
{
  TransferTest tr;
  SmoothTest s;
  MGTest mg(tr, s);
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
  cout << "Smoothing on" << level << endl;
}

void
MGTest::level_residual(unsigned l,
		       Vector<float>&,
		       const Vector<float>&,
		       const Vector<float>&)
{
  cout << "Residual on  " << l << endl;
}
