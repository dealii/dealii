// $Id$

// Test program for multigrid with test matrices (fd)

#include <cmath>
#include <fstream>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparsematrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/precondition.h>
#include <lac/mgbase.h>

template<class VECTOR>
void print_vector(ostream& s, const VECTOR& v)
{
  const unsigned int n = (unsigned int)(sqrt(v.size())+.3);
  unsigned int k=0;
  
  for (unsigned int i=0;i<n;++i)
    {
      for (unsigned int j=0;j<n;++j)
//	s << ((FLOAT)i)/n << '\t' << ((FLOAT)j)/n << '\t' << v(k++) << endl;
	s << '\t' << v(k++);
      s << endl;
    }
  s << endl;
}

class FDMG
  :
  public MGBase
{
  private:
				     /**
				      * Pointer to the level matrices.
				      */
    SmartPointer<MGMatrix<SparseMatrix<FLOAT> > >matrices;
  public:
    FDMG(unsigned int maxlevel, MGMatrix<SparseMatrix<FLOAT> >& matrices,
	 FDMGTransfer& transfer);
    
    
    virtual void level_residual (unsigned int level,
				 Vector<FLOAT>& dst,
				 const Vector<FLOAT>& src,
				 const Vector<FLOAT>& rhs);

    void copy_to_mg(const Vector<double>& rhs);
    void copy_from_mg(Vector<double>& lsg);
};

class MGSmootherLAC
  :
  public MGSmootherBase
{
  private:
    SmartPointer<MGMatrix<SparseMatrix<FLOAT> > >matrices;
  public:
    MGSmootherLAC(MGMatrix<SparseMatrix<FLOAT> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<FLOAT> &u,
			 const Vector<FLOAT> &rhs) const;
    
};

typedef MGCoarseGridLACIteration<SolverCG<SparseMatrix<FLOAT> , Vector<FLOAT>  >,
SparseMatrix<FLOAT>, /*PreconditionRelaxation<SparseMatrix<FLOAT> ,*/
  PreconditionIdentity<Vector<FLOAT> > >
Coarse;

class MGPrecondition
{
  private:
    FDMG& mg;
    MGSmootherBase& smooth;
    Coarse& coarse;
  public:
    MGPrecondition(FDMG& m, MGSmootherBase& s, Coarse& c)
		    :
		    mg(m), smooth(s), coarse(c)
      {}
    
    
    void operator () (Vector<double>& dst, const Vector<double>& src) const
      {
//	print_vector(cout, src);
	
	mg.copy_to_mg(src);
	mg.vcycle(smooth,smooth,coarse);
	mg.copy_from_mg(dst);
      }
    
    
};


main()
{
  ofstream logfile("mg.output");
  deallog.attach(logfile);
  
  PrimitiveVectorMemory<Vector<double>  > mem;
  SolverControl control(100, 1.e-14, true);

  const unsigned int base = 3;
  const unsigned int maxlevel = 8;
  
  const unsigned int maxsize = base * (1 << maxlevel);
  
				   // grid transfer
  FDMGTransfer transfer(maxsize, maxsize, maxlevel);

				   // coarse grid solver
  PrimitiveVectorMemory<Vector<FLOAT> > cgmem;
  ReductionControl cgcontrol(100, 1.e-30, 1.e-2, false, false);
  SolverCG<SparseMatrix<FLOAT> , Vector<FLOAT> > cgcg(cgcontrol,cgmem);

  
  for (unsigned int level = 0; level <= maxlevel; ++level)
    {
      const unsigned int minlevel = 0;

      const unsigned int size = base * (1 << level);

      const unsigned int dim = (size-1)*(size-1);
      deallog << "Level " << level << " size " << size << endl;

				       // Make matrix
      vector<SparseMatrixStruct >  structure(maxlevel+1);
      MGMatrix<SparseMatrix<FLOAT> > A(minlevel,maxlevel);

      FDMatrix testproblem(size, size);

      for(unsigned ll = minlevel; ll <= level; ++ll)
	{
	  const unsigned int size = base * (1 << ll);
	  const unsigned int dim = (size-1)*(size-1);
	  FDMatrix testproblem(size, size);
	  structure[ll].reinit(dim, dim, 5);
	  testproblem.build_structure(structure[ll]);
	  structure[ll].compress();
	  A[ll].reinit(structure[ll]);
	  testproblem.laplacian(A[ll]);
	}
      
      FDMG multigrid(level, A, transfer);

//      PreconditionRelaxation<SparseMatrix<FLOAT> , Vector<FLOAT> >
      PreconditionIdentity<Vector<FLOAT> >
	cgprec;//(A[minlevel], &SparseMatrix<FLOAT> ::template precondition_SSOR<FLOAT>, 1.2);
      
      Coarse coarsegrid(cgcg, A[minlevel], cgprec);
      
      MGSmootherLAC smoother(A);
//      MGSmootherIdentity smoother;

      MGPrecondition precondition(multigrid, smoother, coarsegrid);

//      SolverRichardson<SparseMatrix<FLOAT> , Vector<double> > solver(control, mem);
      SolverCG<SparseMatrix<FLOAT> , Vector<double> > solver(control, mem);

      Vector<double> u(dim);
      Vector<double> f(dim);
      u = 0.;
      f = 1.;//(size*size);

      solver.solve(A[level], u, f, precondition);
      ofstream out("T");
      testproblem.gnuplot_print(out,u);
    }
}

FDMG::FDMG(unsigned int maxlevel, MGMatrix<SparseMatrix<FLOAT> >& matrices,
	   FDMGTransfer& transfer)
		:
		MGBase(transfer, maxlevel, 0),
		matrices(&matrices)
{
  for (unsigned int level = minlevel; level<=maxlevel ; ++level)
    {
      s[level].reinit(matrices[level].m());
      d[level].reinit(matrices[level].m());
    }
}


void
FDMG::level_residual (unsigned int level,
		      Vector<FLOAT>& dst,
		      const Vector<FLOAT>& src,
		      const Vector<FLOAT>& rhs)
{
  (*matrices)[level].residual(dst, src, rhs);
}

void
FDMG::copy_to_mg(const Vector<double>& v)
{
  d[maxlevel] = v;
}

void
FDMG::copy_from_mg(Vector<double>& v)
{
  v = s[maxlevel];
}


MGSmootherLAC::MGSmootherLAC(MGMatrix<SparseMatrix<FLOAT> >& matrix)
		:
		matrices(&matrix)
{}


void
MGSmootherLAC::smooth (const unsigned int level,
		       Vector<FLOAT> &u,
		       const Vector<FLOAT> &rhs) const
{
  SolverControl control(1,1.e-300,false,false);
  PrimitiveVectorMemory<Vector<FLOAT> > mem;
  SolverRichardson<SparseMatrix<FLOAT> , Vector<FLOAT>  > rich(control, mem);
  PreconditionRelaxation<SparseMatrix<FLOAT> , Vector<FLOAT> >
    prec((*matrices)[level], &SparseMatrix<FLOAT> ::template precondition_SSOR<FLOAT>, 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
