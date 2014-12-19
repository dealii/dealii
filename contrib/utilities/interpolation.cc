// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <string>


template <int dim>
inline bool in_cube (const Point<dim>& p)
{
  if ((p(0) < -1e-10) || (p(0) > 1.000000001))
    return false;
  if ((dim>1) && ((p(1) < -1e-10) || (p(1) > 1.000000001)))
    return false;
  if ((dim>2) && ((p(2) < -1e-10) || (p(2) > 1.000000001)))
    return false;
  return true;
}

template <int dim>
inline void
move_points (vector<Point<dim> >& v, unsigned int i)
{
  Point<dim> p;
  switch (dim)
    {
    case 1:
      if (i) p(0) = -1.;
      break;
    case 2:
      if ((i==1) || (i==2))
	p(0) = -1.;
      if ((i==2) || (i==3))
	p(1) = -1.;
      break;
    case 3:
      if ((i==1) || (i==2) || (i==5) || (i==6))
	p(0) = -1.;
      if ((i==4) || (i==5) || (i==6) || (i==7))
	p(1) = -1.;
      if ((i==2) || (i==3) || (i==6) || (i==7))
	p(2) = -1.;
      break;
    default:
      Assert(false, ExcNotImplemented());
    }
  for (unsigned int j=0;j<v.size();++j)
    v[j] += p;
}


template <int dim>
void compute_interpolation (Triangulation<dim>& tr_coarse,
			    Triangulation<dim>& tr_fine,
			    const FiniteElement<dim>& fe_coarse,
			    const FiniteElement<dim>& fe_fine,
			    const char* name)
{
  DoFHandler<dim> dof_coarse(tr_coarse);
  dof_coarse.distribute_dofs(fe_coarse);
  DoFHandler<dim> dof_fine(tr_fine);
  dof_fine.distribute_dofs(fe_fine);
  
  vector<vector<vector<double> > >
    result(2<<dim,
	   vector<vector<double> >(fe_coarse.dofs_per_cell,
				   vector<double>(fe_fine.dofs_per_cell, 0.)));

  vector<unsigned int> indices(fe_fine.dofs_per_cell);

  MappingQ1<dim> mapping;

  vector<Point<dim> > points;
  fe_coarse.get_unit_support_points (points);
  vector<Point<dim> > q_points(points.size());
  vector<double> dummy_weights(points.size(), 0.);
  
				   // Coarse level cell of fine grid
  DoFHandler<dim>::cell_iterator father = dof_fine.begin (0);

  for (unsigned int child_no = 0; child_no < GeometryInfo<dim>::children_per_cell;
       ++child_no)
    {
      DoFHandler<dim>::active_cell_iterator c = father->child(child_no);
      
      c->get_dof_indices(indices);
				       // Evaluate at support points
				       // of father cell.
      for (unsigned int i=0;i<points.size();++i)
	q_points[i] = 2.*points[i];
      move_points (q_points, child_no);
      
      Quadrature<dim> q_fine (q_points, dummy_weights);
      
      FEValues<dim> fine (mapping, fe_fine, q_fine,
			  update_values);
      fine.reinit(c);
					   // Build mass matrix and RHS
	  
	  for (unsigned int k=0;k<fine.n_quadrature_points;++k)
	    {
//	      cerr << k << '\t' << q_points[k];
	      if (in_cube(q_points[k]))
		{
		  for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
		    {
		      double v = fine.shape_value(i,k);
//		      cerr << '[' << i << ' ' << v << ']';
		      result[child_no][k][i] = v;
		    }
		}
//	      cerr << endl;      
	    }
    }
  
  for (unsigned int cell=0;cell<tr_fine.n_active_cells();++cell)
    {
      cout << "static double[] "
	   << name << "_"
	   << dim << "d_"
	   << cell << " =" << endl << '{' << endl;
      for (unsigned int i=0;i<fe_coarse.dofs_per_cell;++i)
	{
	  for (unsigned int j=0;j<fe_fine.dofs_per_cell;++j)
	    {
	      double a = result[cell][i][j];
	      if ((a<1.e-14) && (a>-1.e-14))
		a = 0.;
//	      cout << ' ' << setprecision(10) << a << ",";
	      if (a != 0.)
		cout << "restriction[" << cell << "]("
		     << i << ',' << j << ") = "
		     << a << ';' << endl;
	    }
	  cout << endl;
	}
      cout << "};" << endl << endl;
    }
}


template <int dim>
void compute_projection (Triangulation<dim>& tr_coarse,
			 Triangulation<dim>& tr_fine,
			 const FiniteElement<dim>& fe_coarse,
			 const FiniteElement<dim>& fe_fine,
			 const char* name)
{
  DoFHandler<dim> dof_coarse(tr_coarse);
  dof_coarse.distribute_dofs(fe_coarse);
  DoFHandler<dim> dof_fine(tr_fine);
  dof_fine.distribute_dofs(fe_fine);
  DoFHandler<dim>::cell_iterator cell_coarse = dof_coarse.begin_active();

  vector<vector<vector<double> > >
    result(2<<dim,
	   vector<vector<double> >(fe_coarse.dofs_per_cell,
				   vector<double>(fe_fine.dofs_per_cell, 0.)));

  vector<unsigned int> indices(fe_fine.dofs_per_cell);

  MappingQ1<dim> mapping;
  QGauss<dim> quadrature (9);
  FullMatrix<double> mass (fe_coarse.dofs_per_cell);
  
  FEValues<dim> coarse (mapping, fe_coarse, quadrature,
			update_values
			| update_JxW_values);
  coarse.reinit (cell_coarse);

				   // Build coarse level mass matrix
  for (unsigned int k=0;k<coarse.n_quadrature_points;++k)
    {
      for (unsigned int i=0;i<fe_coarse.dofs_per_cell;++i)
	for (unsigned int j=0;j<fe_coarse.dofs_per_cell;++j)
	  mass(i,j) += coarse.JxW(k)
	    * coarse.shape_value(i,k) * coarse.shape_value(j,k);
    }
  
  FullMatrix<double> inverse (fe_coarse.dofs_per_cell);
  inverse.invert(mass);
  Vector<double> rhs (fe_coarse.dofs_per_cell);
  Vector<double> u (fe_coarse.dofs_per_cell);

				   // Coarse level cell of fine grid
  DoFHandler<dim>::cell_iterator father = dof_fine.begin (0);

  FEValues<dim> fine (mapping, fe_fine, quadrature,
		      update_values
		      | update_JxW_values
		      | update_q_points);

				   // loop over all fine shape functions
  for (unsigned int n=0;n<fe_fine.dofs_per_cell;++n)
    {
      rhs = 0.;
      DoFHandler<dim>::active_cell_iterator c = father->child(0);
      
      c->get_dof_indices(indices);

      fine.reinit(c);

				       // Build a quadrature formula
				       // for the coarse cell
      Quadrature<dim> q_coarse (fine.get_quadrature_points(),
				fine.get_JxW_values());
      FEValues<dim> coarse (mapping, fe_coarse, q_coarse,
			    update_values);
      coarse.reinit(cell_coarse);
					   // Build RHS
	  
      for (unsigned int k=0;k<fine.n_quadrature_points;++k)
	{
	  for (unsigned int i=0;i<fe_coarse.dofs_per_cell;++i)
	    {
	      double v = coarse.shape_value(i,k);
	      double f = fine.shape_value(n,k);
	      rhs(i) += fine.JxW(k) * v*f;
	    }
	}
      inverse.vmult(u,rhs);
      for (unsigned int i=0;i<fe_coarse.dofs_per_cell;++i)
	result[0][i][n] = u(i);
    }
  
  for (unsigned int cell=0;cell<1;++cell)
    {
      cout << "static double "
	   << name << "_"
	   << dim << "d[] =" << endl << '{' << endl;
      for (unsigned int i=0;i<fe_coarse.dofs_per_cell;++i)
	{
	  for (unsigned int j=0;j<fe_fine.dofs_per_cell;++j)
	    {
	      double a = result[cell][i][j];
	      if ((a<1.e-14) && (a>-1.e-14)) a = 0.;
	      cout << ' ' << setprecision(8) << a << ",";
	    }
	  cout << endl;
	}
      cout << "};" << endl << endl;
    }
}

#define PUSH_Q(k) FE_Q<dim> q ## k(k);\
 elements.push_back (ElPair(&q ## k, "dgq" # k))
#define PUSH_DGQ(k) FE_DGQ<dim> dgq ## k(k);\
 dgelements.push_back (ElPair(&dgq ## k, "dgq" # k))

template <int dim>
void loop ()
{
  Triangulation<dim> tr_coarse;
  Triangulation<dim> tr_fine;
  GridGenerator::hyper_cube (tr_coarse, 0, 1);
  GridGenerator::hyper_cube (tr_fine, 0, 1);
  tr_fine.refine_global(1);

  typedef pair<const FiniteElement<dim>*, const char*> ElPair ;
  vector <ElPair> elements(0);
  vector <ElPair> dgelements(0);

  PUSH_Q(1);
  PUSH_Q(2);
  PUSH_Q(3);
  PUSH_Q(4);
  PUSH_DGQ(0);
  PUSH_DGQ(1);
  PUSH_DGQ(2);
  PUSH_DGQ(3);
  PUSH_DGQ(4);
  PUSH_DGQ(5);
  PUSH_DGQ(6);
  PUSH_DGQ(7);

  char* name = new char[100];

  unsigned int n = elements.size();
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=i;j<=i;++j)
      {
//  	ostrstream os (name, 99);
//  	os << "interpolate " << elements[j].second << " refined onto "
//  	   << elements[i].second << ends;
	
//  	compute_interpolation (tr_coarse, tr_fine,
//  			       *elements[i].first,
//  			       *elements[j].first,
//  			       name);
      }

  n = dgelements.size();
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=i;j<=i;++j)
      {
	ostrstream os (name, 99);
	os << "project_" << dgelements[j].second << "_refined_onto_"
	   << dgelements[i].second << ends;
	
	compute_projection (tr_coarse, tr_fine,
			    *dgelements[i].first,
			    *dgelements[j].first,
			    name);
      }

  delete [] name;
}

int main ()
{
  loop<1> ();
  loop<2> ();
  loop<3> ();
}
