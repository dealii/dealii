//----------------------------  mgbase.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mgbase.cc  ---------------------------


#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>
#include <base/logstream.h>

#include <fstream>
#include <iostream>

class TransferTest
  :
  public MGTransferBase
{
    virtual void prolongate(unsigned l,
			    Vector<double>& dst,
			    const Vector<double>& src) const;
    virtual void restrict_and_add (unsigned l,
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
	    defect[l].reinit(1);
	    solution[l].reinit(1);
	  }
      }
    
    virtual void level_vmult(unsigned level,
			     Vector<double>& dst,
			     const Vector<double>& src,
			     const Vector<double>& rhs);

    virtual void print_vector(unsigned int,
			      const Vector<double> &,
			      const char *) const;
    
    void doit()
      {
	level_mgstep(9, smoother, smoother, solver);
      }
};

int main()
{
  std::ofstream logfile("mgbase.output");
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
  deallog << "Prolongating " << l-1 << " to " << l << std::endl;
}

void
TransferTest::restrict_and_add (unsigned l,
				Vector<double>&,
				const Vector<double>&) const
{
  deallog << "Restricting  " << l << " to " << l-1 << std::endl;
}

void
SmoothTest::smooth (const unsigned int  level,
		    Vector<double>      &,
		    const Vector<double>&) const
{
  deallog << "Smoothing on " << level << std::endl;
}

void
MGTest::level_vmult(unsigned l,
		       Vector<double>&,
		       const Vector<double>&,
		       const Vector<double>&)
{
  deallog << "Residual on  " << l << std::endl;
}

void
MGTest::print_vector(unsigned int,
		     const Vector<double> &,
		     const char *) const
{}



void
MGCoarseGridTest::operator() (unsigned int level,
			      Vector<double>&,
			      const Vector<double>&) const
{
  deallog << "Solving on   " << level << std::endl;
}
