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

class FDMG
  :
  public MGBase
{
  private:
    SmartPointer<MGMatrix<SparseMatrix<float> > >matrices;
  public:
    FDMG(unsigned int maxlevel, MGMatrix<SparseMatrix<float> >& matrices,
	 FDMGTransfer& transfer)
		    :
		    MGBase(transfer, maxlevel, 0),
		    matrices(&matrices)
      {}
    
    
    virtual void level_residual (unsigned int level,
				 Vector<float>& dst,
				 const Vector<float>& src,
				 const Vector<float>& rhs);
};

class MGSmootherLAC
  :
  MGSmootherBase
{
  private:
    SmartPointer<MGMatrix<SparseMatrix<float> > >matrices;
  public:
    MGSmootherLAC(MGMatrix<SparseMatrix<float> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<float> &u,
			 const Vector<float> &rhs) const;
    
};

typedef MGCoarseGridLACIteration<SolverCG<SparseMatrix<float> , Vector<float>  >,
SparseMatrix<float>, PreconditionRelaxation<SparseMatrix<float> , Vector<float> > >
Coarse;

class MGPrecondition
{
  private:
    FDMG& mg;
    MGSmootherLAC& smooth;
    Coarse& coarse;
  public:
    MGPrecondition(FDMG& m, MGSmootherLAC& s, Coarse& c)
		    :
		    mg(m), smooth(s), coarse(c)
      {}
    
    
    void operator () (Vector<double>& dst, const Vector<double>& src) const
      {
	
      }
    
    
};


main()
{
  ofstream logfile("mg.output");
  deallog.attach(logfile);
  
  PrimitiveVectorMemory<Vector<double>  > mem;
  SolverControl control(100, 1.e-5, true);

  const unsigned int base = 3;
  const unsigned int maxlevel = 6;
  
  const unsigned int maxsize = base * (1 << maxlevel);
  
				   // grid transfer
  FDMGTransfer transfer(maxsize, maxsize, maxlevel);

				   // coarse grid solver
  PrimitiveVectorMemory<Vector<float> > cgmem;
  SolverControl cgcontrol(100, 1.e-12);
  SolverCG<SparseMatrix<float> , Vector<float> > cgcg(cgcontrol,cgmem);

  vector<SparseMatrixStruct >  structure(maxlevel+1);
  MGMatrix<SparseMatrix<float> > A(maxlevel+1);
  
  for (unsigned int level = 0; level <= maxlevel; ++level)
    {
      const unsigned int size = base * (1 << level);

      const unsigned int dim = (size+1)*(size+1);

				       // Make matrix
      FDMatrix testproblem(size, size);
      structure[level].reinit(dim, dim, 5);
      testproblem.build_structure(structure[level]);
      structure[level].compress();
      A[level].reinit(structure[level]);
      testproblem.laplacian(A[level]);

      FDMG multigrid(level, A, transfer);

      PreconditionRelaxation<SparseMatrix<float> , Vector<float> >
	cgprec(A[0], &SparseMatrix<float> ::template precondition_SSOR<float>, 1.2);
      
      Coarse coarsegrid(cgcg, A[0], cgprec);
      
      MGSmootherLAC smoother(A);
      
//      MGPrecondition prec(multigrid);
    }
}

void
FDMG::level_residual (unsigned int level,
		      Vector<float>& dst,
		      const Vector<float>& src,
		      const Vector<float>& rhs)
{
  (*matrices)[level].residual(dst, src, rhs);
}

MGSmootherLAC::MGSmootherLAC(MGMatrix<SparseMatrix<float> >& matrix)
		:
		matrices(&matrix)
{}


void
MGSmootherLAC::smooth (const unsigned int level,
		       Vector<float> &u,
		       const Vector<float> &rhs) const
{
  SolverControl control(0.,1);
  PrimitiveVectorMemory<Vector<float> > mem;
  SolverRichardson<SparseMatrix<float> , Vector<float>  > rich(control, mem);
  PreconditionRelaxation<SparseMatrix<float> , Vector<float> >
    prec((*matrices)[level], &SparseMatrix<float> ::template precondition_SSOR<float>, 1.2);

  rich.solve((*matrices)[level], u, rhs, prec);
}
