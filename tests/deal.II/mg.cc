// $Id$
/* (c) Guido Kanschat, 1999 */

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
#include <grid/dof_constraints.h>
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
    SmartPointer<MGMatrix<SparseMatrix<double> > >matrices;
  public:
    MGSmootherLAC(MGMatrix<SparseMatrix<double> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<double> &u,
			 const Vector<double> &rhs) const;
    
};

main()
{
  Helmholtz equation;
  RHSFunction<2> rhs;
  QGauss5<2> quadrature;
  
  FEQ1<2> fe1;
  FEQ2<2> fe2;
  FEQ3<2> fe3;
  FEQ4<2> fe4;
  for (unsigned int degree=1;degree<=3;degree++)
    {
      Triangulation<2> tr;
      MGDoFHandler<2> mgdof(&tr);
      DoFHandler<2>& dof(mgdof);
  
      GridGenerator::hyper_cube(tr,-1.,1.);
  
      FiniteElement<2>* fe;
      switch(degree)
	{
	  case 1: fe = &fe1; deallog.push("Q1"); break;
	  case 2: fe = &fe2; deallog.push("Q2"); break;
	  case 3: fe = &fe3; deallog.push("Q3"); break;
	  case 4: fe = &fe4; deallog.push("Q4"); break;
	}
      
      for (unsigned int step=0;step < 3; ++step)
	{
	  tr.refine_global(1);
	  dof.distribute_dofs(*fe);

	  ConstraintMatrix hanging_nodes;
	  dof.make_hanging_node_constraints(hanging_nodes);
	  hanging_nodes.close();

	  const unsigned int size = dof.n_dofs();
	  deallog << "DoFs " << size << endl;
	  deallog << "Levels: " << tr.n_levels() << endl;
	  
	  SparseMatrixStruct structure(size, dof.max_couplings_between_dofs());
	  dof.make_sparsity_pattern(structure);
	  structure.compress();
	  
	  SparseMatrix<double> A(structure);
	  Vector<double> f(size);
	  
	  equation.build_all(A,f,dof, quadrature, rhs);
	  
	  Vector<double> u;
	  u.reinit(f);
	  PrimitiveVectorMemory<Vector<double> > mem;
	  SolverControl control(100, 1.e-12);
	  SolverCG<SparseMatrix<double>, Vector<double> >solver(control, mem);
	  
	  PreconditionIdentity<Vector<double> > precondition;
	  solver.solve(A,u,f,precondition);
	  
	  u = 0.;
	  vector<SparseMatrixStruct> mgstruct(tr.n_levels());
	  MGMatrix<SparseMatrix<double> > mgA(0,tr.n_levels()-1);
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
	  PrimitiveVectorMemory<Vector<double> > cgmem;
	  SolverCG<SparseMatrix<double>, Vector<double> > cgsolver(cgcontrol, cgmem);
	  PreconditionIdentity<Vector<double> > cgprec;
	  MGCoarseGridLACIteration<SolverCG<SparseMatrix<double>, Vector<double> >,
	    SparseMatrix<double>, PreconditionIdentity<Vector<double> > >
	    coarse(cgsolver, mgA[0], cgprec);
	  
	  MGSmootherLAC smoother(mgA);
	  MGTransferPrebuilt transfer;
	  transfer.build_matrices(mgdof);
	  
	  
	  MG<2> multigrid(mgdof, hanging_nodes, mgA, transfer);
	  PreconditionMG<MG<2>, Vector<double> >
	    mgprecondition(multigrid, smoother, smoother, coarse);
	  
	  solver.solve(A, u, f, mgprecondition);
	}
      deallog.pop();
    }
}

template<int dim>
double
RHSFunction<dim>::operator() (const Point<dim>&) const
{
  return 1.;
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
  SolverControl control(2,1.e-300,false,false);
  PrimitiveVectorMemory<Vector<double> > mem;
  SolverRichardson<SparseMatrix<double> , Vector<double>  > rich(control, mem);
  PreconditionRelaxation<SparseMatrix<double> , Vector<double> >
    prec((*matrices)[level], &SparseMatrix<double> ::template precondition_SSOR<double>, 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
