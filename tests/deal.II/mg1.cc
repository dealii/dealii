// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <numerics/multigrid.h>

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

class MGTest
  :
  public MultiGridBase
{
  public:
    MGTest(TransferTest& t)
		    : MultiGridBase(t, 9, 3, 2, 4)
      {}
    
    virtual void level_residual(unsigned level,
				Vector<float>& dst,
				const Vector<float>& src,
				const Vector<float>& rhs);
    
    virtual void smooth(unsigned level,
			Vector<float>& x,
			const Vector<float>& b,
			unsigned steps);
    void doit()
      {
	level_mgstep(9);
      }
};

main()
{
  TransferTest tr;
  MGTest mg(tr);
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
MGTest::level_residual(unsigned l,
		       Vector<float>&,
		       const Vector<float>&,
		       const Vector<float>&)
{
  cout << "Residual on  " << l << endl;
}

void
MGTest::smooth(unsigned l,
	       Vector<float>&,
	       const Vector<float>&,
	       unsigned steps)
{
  cout << "Smooth on    " << l << " : " << steps << " steps" << endl;
}
