// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/sparsematrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_richardson.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/mg_dof.h>
#include <grid/dof_accessor.h>
#include <grid/mg_dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <numerics/multigrid.h>
#include <numerics/multigrid.templates.h>
#include <numerics/mg_smoother.h>

#include "helmholtz.h"

template<int dim>
class RHSFunction
  :
  public Function<dim>
{
  public:
    virtual double operator() (const Point<dim>&) const;
//    virtual Tensor<1,dim> gradient (const Point<dim> &p) const;
};

class MGSmootherLAC
  :
  public MGSmootherBase
{
  private:
    SmartPointer<MGMatrix<SparseMatrix<float> > >matrices;
  public:
    MGSmootherLAC(MGMatrix<SparseMatrix<float> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<float> &u,
			 const Vector<float> &rhs) const;
    
};

main(int argc, const char** argv)
{
  unsigned int refine;
  if (argc<=1)
    {
      cerr << "Usage:" << argv[0] << " <refinement>" << endl;
      refine = 0;
    }
  else
    refine = atoi(argv[1]);
  

  Helmholtz equation;
  RHSFunction<2> rhs;
  QGauss3<2> quadrature;
  
  FELinear<2> fe;
  Triangulation<2> tr;
  MGDoFHandler<2> mgdof(&tr);
  DoFHandler<2>& dof(mgdof);
  
  GridGenerator::hyper_cube(tr,-1.,1.);
  tr.refine_global(refine);
  dof.distribute_dofs(fe);

  const unsigned int size = dof.n_dofs();
  
  SparseMatrixStruct structure(size, dof.max_couplings_between_dofs());
  dof.make_sparsity_pattern(structure);
  structure.compress();
  
  SparseMatrix<double> A(structure);
  Vector<double> f(size);

  equation.build_all(A,f,dof, quadrature, rhs);

  Vector<double> u;
  u.reinit(f);
  PrimitiveVectorMemory<Vector<double> > mem;
  SolverControl control(100, 1.e-15);
  SolverCG<SparseMatrix<double>, Vector<double> >solver(control, mem);
  
  PreconditionIdentity<Vector<double> > precondition;
  solver.solve(A,u,f,precondition);
  
  u = 0.;
  vector<SparseMatrixStruct> mgstruct(tr.n_levels());
  MGMatrix<SparseMatrix<float> > mgA(0,tr.n_levels()-1);
  for (unsigned int i=0;i<tr.n_levels();++i)
    {
      mgstruct[i].reinit(mgdof.n_dofs(i), mgdof.n_dofs(i),
			 mgdof.max_couplings_between_dofs());
      mgdof.make_sparsity_pattern(i, mgstruct[i]);
      mgstruct[i].compress();

      mgA[i].reinit(mgstruct[i]);
    }
  equation.build_mgmatrix(mgA, mgdof, quadrature);
  
  SolverControl cgcontrol(10,0., false, false);
  
  PrimitiveVectorMemory<Vector<float> > cgmem;
  SolverCG<SparseMatrix<float>, Vector<float> > cgsolver(cgcontrol, cgmem);
  PreconditionIdentity<Vector<float> > cgprec;
  MGCoarseGridLACIteration<SolverCG<SparseMatrix<float>, Vector<float> >,
    SparseMatrix<float>, PreconditionIdentity<Vector<float> > >
    coarse(cgsolver, mgA[0], cgprec);

  MGSmootherLAC smoother(mgA);
  MGTransferPrebuilt transfer;
  transfer.build_matrices(mgdof);
  
  deallog << "Levels: " << tr.n_levels() << endl;
  
  MG<2> multigrid(mgdof, mgA, transfer);
  PreconditionMG<MG<2>, Vector<double> >
    mgprecondition(multigrid, smoother, smoother, coarse);
  
  solver.solve(A, u, f, mgprecondition);
  
  cerr << "OK" << endl;
}

template<int dim>
double
RHSFunction<dim>::operator() (const Point<dim>&) const
{
  return 1.;
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
  SolverControl control(1,1.e-300,false,false);
  PrimitiveVectorMemory<Vector<float> > mem;
  SolverRichardson<SparseMatrix<float> , Vector<float>  > rich(control, mem);
  PreconditionRelaxation<SparseMatrix<float> , Vector<float> >
    prec((*matrices)[level], &SparseMatrix<float> ::template precondition_SSOR<float>, 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
