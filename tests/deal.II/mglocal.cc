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
#include <basic/data_io.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <numerics/multigrid.h>
#include <numerics/mg_smoother.h>
MGSmoother* smoother_object;
#include <numerics/multigrid.templates.h>

#include <fstream>

// Does anybody know why cmath does not do this?

#ifndef M_PI_2
#define	M_PI_2		1.57079632679489661923
#endif



int step_counter = 0;

#include "helmholtz.h"

template<int dim>
class RHSFunction
  :
  public Function<dim>
{
  public:
    virtual double value (const Point<dim>&,
			  const unsigned int ) const;
};

extern void write_gnuplot (const MGDoFHandler<2>& dofs,
			   const Vector<double>& v,
			   unsigned int level,
			   ostream &out);

int main()
{
  Helmholtz equation;
  RHSFunction<2> rhs;
  QGauss5<2> quadrature;
  
  FEQ1<2> fe1;
  FEQ2<2> fe2;
  FEQ3<2> fe3;
  FEQ4<2> fe4;
  for (unsigned int degree=1;degree<=1;degree++)
    {
      Triangulation<2> tr;
      MGDoFHandler<2> mgdof(&tr);
      DoFHandler<2>& dof(mgdof);
  
      GridGenerator::hyper_cube(tr,-M_PI_2,M_PI_2);
  
      FiniteElement<2>* fe;
      switch(degree)
	{
	  case 1: fe = &fe1; deallog.push("Q1"); break;
	  case 2: fe = &fe2; deallog.push("Q2"); break;
	  case 3: fe = &fe3; deallog.push("Q3"); break;
	  case 4: fe = &fe4; deallog.push("Q4"); break;
	}

      tr.refine_global(1);
      Triangulation<2>::active_cell_iterator cell = tr.begin_active();
      cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      tr.refine_global(2);
      dof.distribute_dofs(*fe);
      const unsigned int size = dof.n_dofs();
      deallog << "DoFs " << size << endl;
      deallog << "Levels: " << tr.n_levels() << endl;
      for (unsigned int step=1;step < 5; ++step)
	{
	  deallog << "smoothing-steps" << step << endl;
	  SparseMatrixStruct structure(size, dof.max_couplings_between_dofs());
	  dof.make_sparsity_pattern(structure);

	  ConstraintMatrix hanging_nodes;
	  dof.make_hanging_node_constraints(hanging_nodes);
	  hanging_nodes.close();
	  hanging_nodes.condense(structure);

	  structure.compress();
	  
	  SparseMatrix<double> A(structure);
	  Vector<double> f(size);
	  
	  equation.build_all(A,f,dof, quadrature, rhs);
	  
	  hanging_nodes.condense(A);
	  hanging_nodes.condense(f);

	  Vector<double> u;
	  u.reinit(f);
	  PrimitiveVectorMemory<Vector<double> > mem;
	  SolverControl control(20, 1.e-12, true);
	  SolverRichardson<SparseMatrix<double>, Vector<double> >solver(control, mem);
	  
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
	  
	  SolverControl cgcontrol(20,0., false, false);
	  PrimitiveVectorMemory<Vector<double> > cgmem;
	  SolverCG<SparseMatrix<double>, Vector<double> > cgsolver(cgcontrol, cgmem);
	  PreconditionIdentity cgprec;
	  MGCoarseGridLACIteration<SolverCG<SparseMatrix<double>, Vector<double> >,
	    SparseMatrix<double>, PreconditionIdentity>
	    coarse(cgsolver, mgA[tr.n_levels()-2], cgprec);
	  
	  MGSmootherRelaxation<double>
	    smoother(mgdof, mgA, &SparseMatrix<double>::template precondition_SSOR<double>,
		     step, 1.);
	  smoother_object = &smoother;
	  
	  MGTransferPrebuilt transfer;
	  transfer.build_matrices(mgdof);
	  
	  
	  MG<2> multigrid(mgdof, hanging_nodes, mgA, transfer, tr.n_levels()-2);
	  PreconditionMG<MG<2>, Vector<double> >
	    mgprecondition(multigrid, smoother, smoother, coarse);

	   u = 0.;
	   
	  solver.solve(A, u, f, mgprecondition);
	  hanging_nodes.distribute(u);
	}
      deallog.pop();
    }
}

template<int dim>
double
RHSFunction<dim>::value (const Point<dim>&p,
			 const unsigned int) const
{
//  return 1.;
  
  return p(0)*p(0)+p(1)*p(1);
  
  return (2.1)*(sin(p(0))* sin(p(1)));
}


void write_gnuplot (const MGDoFHandler<2>& dofs, const Vector<double>& v, unsigned int level,
		    ostream &out)
{
  MGDoFHandler<2>::cell_iterator cell;
  MGDoFHandler<2>::cell_iterator endc = dofs.end(level);

  Vector<double> values(dofs.get_fe().total_dofs);
  
  unsigned int cell_index=0;
  for (cell=dofs.begin(level); cell!=endc; ++cell, ++cell_index) 
    {
      cell->get_mg_dof_values(v, values);
      
      out << cell->vertex(0) << "  " << values(0) << endl;
      out << cell->vertex(1) << "  " << values(1) << endl;
      out << endl;
      out << cell->vertex(3) << "  " << values(3) << endl;
      out << cell->vertex(2) << "  " << values(2) << endl;
      out << endl;
      out << endl;
    }
  out << endl;
}
