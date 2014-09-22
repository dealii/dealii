// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Author: Guido Kanschat

#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/full_matrix.templates.h>
#include <lac/solver_richardson.h>
#include <lac/vector_memory.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_system.h>
#include <fe/mapping_cartesian.h>
#include <fe/fe_values.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <string>

/*
 * Enter name of the finite element family here.
 */

#define ELNAME "FE_DGP"
#define PUSH(k) FE_DGP<dim> q ## k (k);\
 elements.push_back (ElPair(&q ## k, "dgp" # k))



template <int dim>
void compute_embedding (unsigned int degree,
			Triangulation<dim>& tr_coarse,
			Triangulation<dim>& tr_fine,
			const FiniteElement<dim>& fe_coarse,
			const FiniteElement<dim>& fe_fine,
			const char* name)
{
  deallog.push(name);
  
  DoFHandler<dim> dof_coarse(tr_coarse);
  dof_coarse.distribute_dofs(fe_coarse);
  DoFHandler<dim> dof_fine(tr_fine);
  dof_fine.distribute_dofs(fe_fine);
  
  FullMatrix<long double> A(dof_fine.n_dofs());
  Vector<long double> f(dof_fine.n_dofs());
  Vector<long double> u(dof_fine.n_dofs());
  vector<vector<vector<double> > >
    result(2<<dim,
	   vector<vector<double> >(fe_coarse.dofs_per_cell,
				   vector<double>(fe_fine.dofs_per_cell, 0.)));

  DoFHandler<dim>::active_cell_iterator coarse_cell
    = dof_coarse.begin_active();
  vector<unsigned int> indices(fe_fine.dofs_per_cell);

  MappingCartesian<dim> mapping;
  QGauss<dim> q_fine(degree+1);
  
  FEValues<dim> fine (mapping, fe_fine, q_fine,
		      update_q_points
		      | update_JxW_values
		      | update_values);


  DoFHandler<dim>::active_cell_iterator c;

  for (unsigned int coarse_no=0;coarse_no<fe_coarse.dofs_per_cell;
       ++coarse_no)
    {
      A.clear();
      f.clear();
      u.clear();

      for (c = dof_fine.begin_active();
	   c != dof_fine.end();
	   ++c)
	{
	  c->get_dof_indices(indices);
	  fine.reinit(c);
					   // Build mass matrix and RHS
	  
	  Quadrature<dim> q_coarse (fine.get_quadrature_points(),
				    fine.get_JxW_values());
	  FEValues<dim> coarse (mapping, fe_coarse, q_coarse,
				update_values);
	  coarse.reinit(coarse_cell);
	  
	  for (unsigned int k=0;k<fine.n_quadrature_points;++k)
	    {
	      double dx = fine.JxW(k);
	      for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
		{
		  double v = fine.shape_value(i,k);
		  f(indices[i]) += dx *
		    v * coarse.shape_value(coarse_no, k);
		  for (unsigned int j=0;j<fe_fine.dofs_per_cell;++j)
		    A(indices[i],indices[j]) += dx*v*fine.shape_value(j,k);
		}
	    }
	}
//      A.print_formatted(cout, 2, false, 4, "~", 9);
      FullMatrix<double> P(A.n());
      P.invert(A);
      SolverControl control (100, 1.e-24, true, false);
      PrimitiveVectorMemory<Vector<long double> > mem;
      SolverRichardson<Vector<long double> > solver(control, mem);
      solver.solve(A,u,f,P);
      
      unsigned int cell_i=0;
      for (c = dof_fine.begin_active();
	   c != dof_fine.end();
	   ++c, ++cell_i)
	{
	  c->get_dof_indices(indices);
	  for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
	    result[cell_i][coarse_no][i] = (fabs(u(indices[i])) > 1.e-16)
	      ? (27.*u(indices[i])) : 0.;
	}
    }

  for (unsigned int cell=0;cell<tr_fine.n_active_cells();++cell)
    {
      cout << "static const double "
	   << name << "_"
	   << cell << "[] =" << endl << '{' << endl;
      for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
	{
	  for (unsigned int j=0;j<fe_coarse.dofs_per_cell;++j)
	    {
	      double s = result[cell][j][i];
	      if (fabs(s) < 1.e-13)
		cout << " 0,";
	      else
		cout << ' ' << setprecision(10) << s << "/27.,";
	    }
	  cout << endl;
	}
      cout << "};" << endl << endl;
    }
  deallog.pop();
}


template <int dim>
void loop ()
{
  char prefix[3];
  sprintf(prefix, "%dd", dim);
  deallog.push(prefix);
  
  cout << "namespace " << ELNAME << "_" << dim << "d\n{";

  Triangulation<dim> tr_coarse;
  Triangulation<dim> tr_fine;
  GridGenerator::hyper_cube (tr_coarse, 0, 1);
  GridGenerator::hyper_cube (tr_fine, 0, 1);
  tr_fine.refine_global(1);

  typedef pair<const FiniteElement<dim>*, const char*> ElPair ;
  vector <ElPair> elements(0);

				   /*
				    * List element degrees for
				    * computation here.
				    */
  PUSH(0);
  PUSH(1);
  PUSH(2);
  PUSH(3);
  PUSH(4);
  PUSH(5);
  PUSH(6);

  char* name = new char[100];

				   // Embed all lower spaces into
				   // higher or just the same degree
				   // on different grids.
  bool same_degree_only = true;
  
  unsigned int n = elements.size();
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=((same_degree_only) ? i : 0);j<=i;++j)
      {
	ostrstream os (name, 99);
	os << elements[i].second << "_into_"
	   << elements[j].second << "_refined" << ends;
	
	compute_embedding (i, tr_coarse, tr_fine,
			   *elements[i].first,
			   *elements[j].first,
			   name);
      }

  delete [] name;
  cout << "};\n";
  deallog.pop();
  
}

int main ()
{
  loop<1> ();
  loop<2> ();
  loop<3> ();
}
