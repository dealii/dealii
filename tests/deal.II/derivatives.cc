//----------------------------------------------------------------------
// $Id$
// Version: $Name$
//
// (c) 2001 the deal.II authors
//
// purpose of this test:
//
// check 1st and 2nd order derivatives of continuous elements.
//
//
// Method:
//
// solve  (u,v) = (f,v)
//
// Compute L2, H1 and H2 norm of the errors. Does not use VectorTools.
//
// The logfile checked in for comparison does not contain the errors.
// There, we compute the L^2-norm of the derivatives. Change the variable
// errors below to true to compute errors.
//
//----------------------------------------------------------------------

const bool errors = false;


#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/mapping_cartesian.h>
#include <fe/mapping_q1.h>
#include <fe/mapping_q.h>
#include <fe/fe_values.h>
#include <numerics/data_out.h>
#include <vector>
#include <fstream>
#include <string>


template <int dim>
void
check (const unsigned int level,
       const Mapping<dim>& mapping,
       const FiniteElement<dim>& element,
       const Quadrature<dim>& quadrature)
{
  Triangulation<dim> tr;

  CosineFunction<dim> cosine;
  
  DoFHandler<dim> dof(tr);
  if (dim==2)
    GridGenerator::hyper_ball(tr);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  
  tr.refine_global (level);
  
  dof.distribute_dofs(element);

  FEValues<dim> fe (mapping, element, quadrature,
		    update_values
		    | update_q_points | update_JxW_values);

  vector <unsigned int> global_dofs (element.dofs_per_cell);
  vector <double> function (quadrature.n_quadrature_points);

  Vector<double> u (dof.n_dofs ());
  Vector<double> f (dof.n_dofs ());

  SparsityPattern A_pattern (dof.n_dofs (), dof.n_dofs (),
			     dof.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof, A_pattern);
  A_pattern.compress ();
  
  SparseMatrix<double> A(A_pattern);
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  const DoFHandler<dim>::cell_iterator end = dof.end();

  for (; cell != end;++cell)
    {
      fe.reinit(cell);
      cell->get_dof_indices (global_dofs);
      cosine.value_list (fe.get_quadrature_points(), function);

      for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
	{
	  double dx = fe.JxW (k);

	  for (unsigned int i=0;i<element.dofs_per_cell;++i)
	    {
	      const double v = fe.shape_value (i,k);
	      double rhs = dx * v * (function[k]);

	      f(global_dofs[i]) += rhs;
	      for (unsigned int j=0;j<element.dofs_per_cell;++j)
		{
		  const double u = fe.shape_value (j,k);
		  double el = dx * (u*v/* + (Du*Dv) */);
		  A.add (global_dofs[i], global_dofs[j], el);
		}
	    }
	}
    }

  SolverControl control (1000, 1.e-10, false, false);
  PrimitiveVectorMemory<Vector<double> > mem;
  SolverCG<Vector<double> > solver (control, mem);
  PreconditionIdentity prec;
  
  solver.solve (A, u, f, prec);

  FEValues<dim> fe2 (mapping, element, quadrature,
		     update_values | update_gradients
		     | update_second_derivatives
		     | update_q_points | update_JxW_values);

  double l2 = 0.;
  double h1 = 0.;
  double h2 = 0.;

  vector<double> u_local (quadrature.n_quadrature_points);
  vector<Tensor<1,dim> > Du (quadrature.n_quadrature_points);
  vector<Tensor<1,dim> > Df (quadrature.n_quadrature_points);
  vector<Tensor<2,dim> > DDu (quadrature.n_quadrature_points);
  vector<Tensor<2,dim> > DDf (quadrature.n_quadrature_points);
  
  for (cell = dof.begin_active(); cell != end;++cell)
    {
      fe2.reinit (cell);
      
      cosine.value_list (fe2.get_quadrature_points(), function);
      cosine.gradient_list (fe2.get_quadrature_points(), Df);
      cosine.hessian_list (fe2.get_quadrature_points(), DDf);
      fe2.get_function_values (u, u_local);
      fe2.get_function_grads (u, Du);
      fe2.get_function_2nd_derivatives (u,DDu);
      
      for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
	{
	  const double dx = fe.JxW (k);
	  double e =  u_local[k];
	  if (errors) e -= function[k];
	  l2 += dx*e*e;
	  for (unsigned int i=0;i<dim;++i)
	    {
	      e = Du[k][i];
	      if (errors) e -= Df[k][i];
	      h1 += dx*e*e;
	      for (unsigned int j=0;j<dim;++j)
		{
		  e = DDu[k][i][j];
		  if (errors) e -= DDf[k][i][j];
		  h2 += dx*e*e;
		}
	    }
	}
    }
      deallog << "L2: " << l2 << endl;
      deallog << "H1: " << h1 << endl;
      deallog << "H2: " << h2 << endl;  
}

template <int dim>
void loop ()
{
  QGauss<dim> gauss((dim<3) ? 5 : 3);

  vector<Mapping<dim>*> maps;
//  maps.push_back (new MappingCartesian<dim>);
  maps.push_back (new MappingQ1<dim>);
  maps.push_back (new MappingQ<dim>(2));

  vector<FiniteElement<dim>*> elements;
  elements.push_back (new FE_Q<dim> (1));
  elements.push_back (new FE_Q<dim> (2));
  if (dim<3)
    {
      elements.push_back (new FE_Q<dim> (3));
      elements.push_back (new FE_Q<dim> (4));
    }
  
  elements.push_back (new FE_DGQ<dim> (1));
  elements.push_back (new FE_DGQ<dim> (2));
  if (dim<3)
    {
      elements.push_back (new FE_DGQ<dim> (3));
      elements.push_back (new FE_DGQ<dim> (4));
    }
  
  for (unsigned int m=0;m<maps.size();++m)
    for (unsigned int e=0;e<elements.size();++e)
    {
      check (1, *maps[m], *elements[e], gauss);
//      check (2, *maps[m], *elements[e], gauss);
    }
}


int main ()
{
  ofstream logfile ("derivatives.output");
  logfile.precision (2);
  logfile.setf(ios::fixed);  
  deallog.attach(logfile);
  if (!errors) deallog.depth_console(0);

  deallog.push ("1d");
  loop<1> ();
  deallog.pop ();
  deallog.push ("2d");
  loop<2> ();
  deallog.pop ();
  deallog.push ("3d");
  loop<3> ();
  deallog.pop ();
}
