//----------------------------  mglocal.cc  ---------------------------
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
//----------------------------  mglocal.cc  ---------------------------


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
#include <multigrid/mg_dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/grid_generator.h>
#include <numerics/data_out.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <multigrid/multigrid.h>
#include <multigrid/mg_smoother.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_tools.h>

#include <fstream>

MGSmoother* smoother_object;

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
  ofstream logfile("mglocal.output");
//  logfile.setf(ios::fixed);
//  logfile.precision (3);
  deallog.attach(logfile);
//  deallog.depth_console(0);

  Helmholtz equation;
  RHSFunction<2> rhs;
  QGauss5<2> quadrature;
  
  FEQ1<2> fe1;
  FEQ2<2> fe2;
  FEQ3<2> fe3;
  FEQ4<2> fe4;
  for (unsigned int degree=1;degree<=4;degree++)
    {
      Triangulation<2> tr;
      MGDoFHandler<2> mgdof(tr);
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
      for (unsigned int step=14;step < 15; ++step)
	{
	  deallog << "smoothing-steps" << step << endl;
	  SparsityPattern structure(size, dof.max_couplings_between_dofs());
	  DoFTools::make_sparsity_pattern(dof, structure);

	  ConstraintMatrix hanging_nodes;
	  DoFTools::make_hanging_node_constraints (dof, hanging_nodes);
	  hanging_nodes.close();
	  hanging_nodes.condense(structure);

	  structure.compress();
	  
	  SparseMatrix<double> A(structure);
	  Vector<double> f(size);
	  
	  equation.build_all(A,f,dof, quadrature, rhs);
	  
	  hanging_nodes.condense(A);
	  hanging_nodes.condense(f);

	  if (true)
	    {
	      ofstream out_file("MGf");
	      DataOut<2> out;
	      out.attach_dof_handler(dof);
	      out.add_data_vector(f, "v");
	      out.build_patches(5);
	      out.write_gnuplot(out_file);
	    }
	  
	  Vector<double> u;
	  u.reinit(f);
	  PrimitiveVectorMemory<> mem;
	  SolverControl control(20, 1.e-8, true);
	  SolverCG<> solver(control, mem);
	  
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
	  
	  SolverControl cgcontrol(2000,1.e-14, false, false);
	  PrimitiveVectorMemory<> cgmem;
	  SolverCG<> cgsolver(cgcontrol, cgmem);
	  PreconditionIdentity cgprec;
	  MGCoarseGridLACIteration<SolverCG<>, SparseMatrix<double>, PreconditionIdentity>
	    coarse(cgsolver, mgA[tr.n_levels()-2], cgprec);
	  
	  MGSmootherRelaxation<double>
	    smoother(mgdof, mgA, &SparseMatrix<double>::template precondition_SSOR<double>,
		     step, 1.);
	  smoother_object = &smoother;
	  
	  MGTransferPrebuilt transfer;
	  transfer.build_matrices(mgdof);


Multigrid<2> multigrid(mgdof, hanging_nodes, mgstruct, mgA, transfer, tr.n_levels()-2);
	  PreconditionMG<Multigrid<2> >
	    mgprecondition(multigrid, smoother, smoother, coarse);

	   u = 0.;
	   
	  solver.solve(A, u, f, mgprecondition);
	  hanging_nodes.distribute(u);

	  DataOut<2> out;
	  char* name = new char[100];

	  sprintf(name, "MG-Q%d-%d", degree, step);
	  
	  ofstream ofile(name);
	  out.attach_dof_handler(dof);
	  out.add_data_vector(u,"u");
	  out.add_data_vector(f,"f");
	  out.build_patches(5);
	  out.write_gnuplot(ofile);

	  delete[] name;
	}
      deallog.pop();
    }
}

template<int dim>
double
RHSFunction<dim>::value (const Point<dim>&p,
			 const unsigned int) const
{
  return 1.;
  
  return p(0)*p(0)+p(1)*p(1);
  
  return (2.1)*(sin(p(0))* sin(p(1)));
}
