//----------------------------  mg.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg.cc  ---------------------------


#include <base/function.h>
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_richardson.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <dofs/dof_tools.h>

#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_tools.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/multigrid.h>

#include "helmholtz.h"

#include <fstream>


template<int dim>
class RHSFunction : public Function<dim>
{
  public:
    virtual double value (const Point<dim>&,
			  const unsigned int) const;
};


class MGSmootherLAC : public MGSmootherBase
{
  private:
    SmartPointer<MGLevelObject<SparseMatrix<double> > >matrices;
  public:
    MGSmootherLAC(MGLevelObject<SparseMatrix<double> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<double> &u,
			 const Vector<double> &rhs) const;
    
};

int main()
{
  ofstream logfile("mg.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

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
      MGDoFHandler<2> mgdof(tr);
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
	  DoFTools::make_hanging_node_constraints(dof, hanging_nodes);
	  hanging_nodes.close();

	  const unsigned int size = dof.n_dofs();
	  deallog << "DoFs " << size << endl;
	  deallog << "Levels: " << tr.n_levels() << endl;
	  
	  SparsityPattern structure(size, dof.max_couplings_between_dofs());
	  DoFTools::make_sparsity_pattern(dof, structure);
	  structure.compress();
	  
	  SparseMatrix<double> A(structure);
	  Vector<double> f(size);
	  
	  equation.build_all(A,f,dof, quadrature, rhs);
	  
	  Vector<double> u;
	  u.reinit(f);
	  PrimitiveVectorMemory<> mem;
	  SolverControl control(100, 1.e-12);
	  SolverCG<>    solver(control, mem);
	  
	  u = 0.;
	  MGLevelObject<SparsityPattern> mgstruct(0, tr.n_levels()-1);
	  MGLevelObject<SparseMatrix<double> > mgA(0,tr.n_levels()-1);
	  for (unsigned int i=0;i<tr.n_levels();++i)
	    {
	      mgstruct[i].reinit(mgdof.n_dofs(i), mgdof.n_dofs(i),
				 mgdof.max_couplings_between_dofs());
	      MGDoFTools::make_sparsity_pattern(mgdof, mgstruct[i], i);
	      mgstruct[i].compress();
	      
	      mgA[i].reinit(mgstruct[i]);
	    }
	  equation.build_mgmatrix(mgA, mgdof, quadrature);
	  
	  SolverControl cgcontrol(10,0., false, false);
	  PrimitiveVectorMemory<> cgmem;
	  SolverCG<> cgsolver(cgcontrol, cgmem);
	  PreconditionIdentity cgprec;
	  MGCoarseGridLACIteration<SolverCG<>, SparseMatrix<double>, PreconditionIdentity>
	    coarse(cgsolver, mgA[0], cgprec);
	  
	  MGSmootherLAC smoother(mgA);
	  MGTransferPrebuilt transfer;
	  transfer.build_matrices(mgdof);


Multigrid<2> multigrid(mgdof, hanging_nodes, mgstruct, mgA, transfer);
	  PreconditionMG<Multigrid<2> >
	    mgprecondition(multigrid, smoother, smoother, coarse);
	  
	  solver.solve(A, u, f, mgprecondition);
	}
      deallog.pop();
    }
}

template<int dim>
double
RHSFunction<dim>::value (const Point<dim>&,
			 const unsigned int) const
{
  return 1.;
}

MGSmootherLAC::MGSmootherLAC(MGLevelObject<SparseMatrix<double> >& matrix)
		:
		matrices(&matrix)
{}

void
MGSmootherLAC::smooth (const unsigned int level,
		       Vector<double> &u,
		       const Vector<double> &rhs) const
{
  SolverControl control(2,1.e-300,false,false);
  PrimitiveVectorMemory<> mem;
  SolverRichardson<> rich(control, mem);
  PreconditionRelaxation<>
    prec((*matrices)[level], &SparseMatrix<double> ::template precondition_SSOR<double>, 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}
