//----------------------------  mg.cc  ---------------------------
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
//----------------------------  mg.cc  ---------------------------


#include <cmath>
#include <fstream>
#include <iostream>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/precondition.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>


#define TYPE  double
#define ACCURACY 1.e-5

// template<class VECTOR>
// void print_vector(std::ostream& s, const VECTOR& v)
// {
//   const unsigned int n = (unsigned int)(sqrt(v.size())+.3);
//   unsigned int k=0;
  
//   for (unsigned int i=0;i<n;++i)
//     {
//       for (unsigned int j=0;j<n;++j)
// //	s << ((float)i)/n << '\t' << ((float)j)/n << '\t' << v(k++) << std::endl;
// 	s << '\t' << v(k++);
//       s << std::endl;
//     }
//   s << std::endl;
// }

class FDMG
  :
  public MGBase
{
  private:
				     /**
				      * Pointer to the level matrices.
				      */
    SmartPointer<MGLevelObject<SparseMatrix<double> > >matrices;
  public:
    FDMG(unsigned int maxlevel, MGLevelObject<SparseMatrix<double> >& matrices,
	 FDMGTransfer& transfer);


virtual void level_vmult (unsigned int level,
				 Vector<double>& dst,
				 const Vector<double>& src,
				 const Vector<double>& rhs);

    void copy_to_mg(const Vector<TYPE>& rhs);
    void copy_from_mg(Vector<TYPE>& lsg);
    virtual void print_vector(unsigned int,
			      const Vector<double> &,
			      const char *) const;

};

class MGSmootherLAC
  :
  public MGSmootherBase
{
  private:
    SmartPointer<MGLevelObject<SparseMatrix<double> > >matrices;
  public:
    MGSmootherLAC(MGLevelObject<SparseMatrix<double> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<double> &u,
			 const Vector<double> &rhs) const;
    
};

typedef MGCoarseGridLACIteration<SolverCG<>,
                                 SparseMatrix<double>,
                                 PreconditionIdentity > Coarse;


int main()
{
  std::ofstream logfile("mg.output");
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
      deallog << "Level " << level << " size " << size << std::endl;

				       // Make matrix
      std::vector<SparsityPattern>  structure(maxlevel+1);
      MGLevelObject<SparseMatrix<double> > A(minlevel,maxlevel);

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

//      SolverRichardson<Vector<TYPE> > solver(control, mem);
      SolverCG<Vector<TYPE> > solver(control, mem);

      Vector<TYPE> u(dim);
      Vector<TYPE> f(dim);
      u = 0.;
      f = 1.;//(size*size);

      solver.solve(A[level], u, f, precondition);
      std::ofstream out("T");
      testproblem.gnuplot_print(out,u);
    }
}

FDMG::FDMG(unsigned int maxlevel, MGLevelObject<SparseMatrix<double> >& matrices,
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

void
FDMG::print_vector(unsigned int,
		   const Vector<double> &,
		   const char *) const
{}



MGSmootherLAC::MGSmootherLAC(MGLevelObject<SparseMatrix<double> >& matrix)
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
  PreconditionSSOR<> prec;
  prec.initialize((*matrices)[level], 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
