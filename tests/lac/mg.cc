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

#define TYPE  double
#define ACCURACY 1.e-5

// template<class VECTOR>
// void print_vector(ostream& s, const VECTOR& v)
// {
//   const unsigned int n = (unsigned int)(sqrt(v.size())+.3);
//   unsigned int k=0;
  
//   for (unsigned int i=0;i<n;++i)
//     {
//       for (unsigned int j=0;j<n;++j)
// //	s << ((float)i)/n << '\t' << ((float)j)/n << '\t' << v(k++) << endl;
// 	s << '\t' << v(k++);
//       s << endl;
//     }
//   s << endl;
// }

class FDMG
  :
  public MGBase
{
  private:
				     /**
				      * Pointer to the level matrices.
				      */
    SmartPointer<MGMatrix<SparseMatrix<double> > >matrices;
  public:
    FDMG(unsigned int maxlevel, MGMatrix<SparseMatrix<double> >& matrices,
	 FDMGTransfer& transfer);
    
    
    virtual void level_vmult (unsigned int level,
				 Vector<double>& dst,
				 const Vector<double>& src,
				 const Vector<double>& rhs);

    void copy_to_mg(const Vector<TYPE>& rhs);
    void copy_from_mg(Vector<TYPE>& lsg);
};

class MGSmootherLAC
  :
  public MGSmootherBase
{
  private:
    SmartPointer<MGMatrix<SparseMatrix<double> > >matrices;
  public:
    MGSmootherLAC(MGMatrix<SparseMatrix<double> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<double> &u,
			 const Vector<double> &rhs) const;
    
};

typedef MGCoarseGridLACIteration<SolverCG<>,
                                 SparseMatrix<double>,
                                 PreconditionIdentity > Coarse;


int main()
{
  ofstream logfile("mg.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  PrimitiveVectorMemory<Vector<TYPE>  > mem;
  SolverControl control(10000, ACCURACY);

  const unsigned int base = 3;
  const unsigned int maxlevel = 5;
  
  const unsigned int maxsize = base * (1 << maxlevel);
  
				   // grid transfer
  FDMGTransfer transfer(maxsize, maxsize, maxlevel);

				   // coarse grid solver
  PrimitiveVectorMemory<> cgmem;
  ReductionControl cgcontrol(100, 1.e-30, 1.e-2, false, false);
  SolverCG<> cgcg(cgcontrol,cgmem);

  
  for (unsigned int level = 0; level <= maxlevel; ++level)
    {
      const unsigned int minlevel = 0;

      const unsigned int size = base * (1 << level);

      const unsigned int dim = (size-1)*(size-1);
      deallog << "Level " << level << " size " << size << endl;

				       // Make matrix
      vector<SparseMatrixStruct >  structure(maxlevel+1);
      MGMatrix<SparseMatrix<double> > A(minlevel,maxlevel);

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

//      PreconditionRelaxation<>
      PreconditionIdentity
	cgprec;//(A[minlevel], &SparseMatrix<double> ::template precondition_SSOR<double>, 1.2);
      
      Coarse coarsegrid(cgcg, A[minlevel], cgprec);
      
      MGSmootherLAC smoother(A);
//      MGSmootherIdentity smoother;

      PreconditionMG<FDMG, Vector<TYPE> >
	precondition(multigrid, smoother, smoother, coarsegrid);

//      SolverRichardson<SparseMatrix<double> , Vector<TYPE> > solver(control, mem);
      SolverCG<SparseMatrix<double> , Vector<TYPE> > solver(control, mem);

      Vector<TYPE> u(dim);
      Vector<TYPE> f(dim);
      u = 0.;
      f = 1.;//(size*size);

      solver.solve(A[level], u, f, precondition);
      ofstream out("T");
      testproblem.gnuplot_print(out,u);
    }
}

FDMG::FDMG(unsigned int maxlevel, MGMatrix<SparseMatrix<double> >& matrices,
	   FDMGTransfer& transfer)
		:
		MGBase(transfer, 0, maxlevel),
		matrices(&matrices)
{
  for (unsigned int level = minlevel; level<=maxlevel ; ++level)
    {
      solution[level].reinit(matrices[level].m());
      defect[level].reinit(matrices[level].m());
    }
}


void
FDMG::level_vmult (unsigned int level,
		      Vector<double>& dst,
		      const Vector<double>& src,
		      const Vector<double>&)
{
  (*matrices)[level].vmult(dst, src);
}

void
FDMG::copy_to_mg(const Vector<TYPE>& v)
{
  defect[maxlevel] = v;
}

void
FDMG::copy_from_mg(Vector<TYPE>& v)
{
  v = solution[maxlevel];
}


MGSmootherLAC::MGSmootherLAC(MGMatrix<SparseMatrix<double> >& matrix)
		:
		matrices(&matrix)
{}


void
MGSmootherLAC::smooth (const unsigned int level,
		       Vector<double> &u,
		       const Vector<double> &rhs) const
{
  SolverControl control(1,1.e-300,false,false);
  PrimitiveVectorMemory<> mem;
  SolverRichardson<> rich(control, mem);
  PreconditionRelaxation<>
    prec((*matrices)[level], &SparseMatrix<double> ::template precondition_SSOR<double>, 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
