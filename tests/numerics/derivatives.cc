// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2017 by the deal.II authors
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


const bool errors = false;


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <vector>
#include <string>


template <int dim>
void
check (const unsigned int level,
       const Mapping<dim> &mapping,
       const FiniteElement<dim> &element,
       const Quadrature<dim> &quadrature)
{
  Triangulation<dim> tr;

  Functions::CosineFunction<dim> cosine;

  DoFHandler<dim> dof(tr);
  if (dim==2)
    {
      GridGenerator::hyper_ball(tr, Point<dim>(), 1);
      tr.reset_manifold(0);
    }
  else
    GridGenerator::hyper_cube(tr, -1,1);

  tr.refine_global (level);

  dof.distribute_dofs(element);

  FEValues<dim> fe (mapping, element, quadrature,
                    update_values
                    | update_quadrature_points | update_JxW_values);

  std::vector <types::global_dof_index> global_dofs (element.dofs_per_cell);
  std::vector <double> function (quadrature.size());

  Vector<double> u (dof.n_dofs ());
  Vector<double> f (dof.n_dofs ());

  SparsityPattern A_pattern (dof.n_dofs (), dof.n_dofs (),
                             dof.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof, A_pattern);
  A_pattern.compress ();

  SparseMatrix<double> A(A_pattern);

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  const typename DoFHandler<dim>::cell_iterator end = dof.end();

  for (; cell != end; ++cell)
    {
      fe.reinit(cell);
      cell->get_dof_indices (global_dofs);
      cosine.value_list (fe.get_quadrature_points(), function);

      for (unsigned int k=0; k<quadrature.size(); ++k)
        {
          double dx = fe.JxW (k);

          for (unsigned int i=0; i<element.dofs_per_cell; ++i)
            {
              const double v = fe.shape_value (i,k);
              double rhs = dx * v * (function[k]);

              f(global_dofs[i]) += rhs;
              for (unsigned int j=0; j<element.dofs_per_cell; ++j)
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
                     | update_hessians
                     | update_quadrature_points | update_JxW_values);

  double l2 = 0.;
  double h1 = 0.;
  double h2 = 0.;

  std::vector<double> u_local (quadrature.size());
  std::vector<Tensor<1,dim> > Du (quadrature.size());
  std::vector<Tensor<1,dim> > Df (quadrature.size());
  std::vector<Tensor<2,dim> > DDu (quadrature.size());
  std::vector<SymmetricTensor<2,dim> > DDf (quadrature.size());

  for (cell = dof.begin_active(); cell != end; ++cell)
    {
      fe2.reinit (cell);

      cosine.value_list (fe2.get_quadrature_points(), function);
      cosine.gradient_list (fe2.get_quadrature_points(), Df);
      cosine.hessian_list (fe2.get_quadrature_points(), DDf);
      fe2.get_function_values (u, u_local);
      fe2.get_function_gradients (u, Du);
      fe2.get_function_hessians (u,DDu);

      for (unsigned int k=0; k<quadrature.size(); ++k)
        {
          const double dx = fe.JxW (k);
          double e =  u_local[k];
          if (errors) e -= function[k];
          l2 += dx*e*e;
          for (unsigned int i=0; i<dim; ++i)
            {
              e = Du[k][i];
              if (errors) e -= Df[k][i];
              h1 += dx*e*e;
              for (unsigned int j=0; j<dim; ++j)
                {
                  e = DDu[k][i][j];
                  if (errors) e -= DDf[k][i][j];
                  h2 += dx*e*e;
                }
            }
        }
    }
  deallog << "L2: " << l2 << std::endl;
  deallog << "H1: " << h1 << std::endl;
  deallog << "H2: " << h2 << std::endl;
}

template <int dim>
void loop ()
{
  QGauss<dim> gauss((dim<3) ? 5 : 3);

  std::vector<Mapping<dim>*> maps;
//  maps.push_back (new MappingCartesian<dim>);
  maps.push_back (new MappingQGeneric<dim>(1));
  maps.push_back (new MappingQ<dim>(2));

  std::vector<FiniteElement<dim>*> elements;
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

  for (unsigned int m=0; m<maps.size(); ++m)
    for (unsigned int e=0; e<elements.size(); ++e)
      {
        check (1, *maps[m], *elements[e], gauss);
//      check (2, *maps[m], *elements[e], gauss);
      }

  for (unsigned int m=0; m<maps.size(); ++m)
    delete maps[m];

  for (unsigned int e=0; e<elements.size(); ++e)
    delete elements[e];
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);

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
