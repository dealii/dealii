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




#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <vector>
#include <string>




void
solve_filtered (std::map<types::global_dof_index,double> &bv,
                SparseMatrix<double>          &A,
                Vector<double>                &u,
                Vector<double>                &f)
{
  FilteredMatrix<Vector<double> > A1 (A);
  A1.add_constraints (bv);

  SolverControl control (1000, 1.e-10, false, false);
  PrimitiveVectorMemory<Vector<double> > mem;
  SolverCG<Vector<double> > solver (control, mem);
  PreconditionJacobi<SparseMatrix<double> > prec;
  FilteredMatrix<Vector<double> > fprec;
  prec.initialize (A, 1.2);
  fprec.initialize(prec);

  Vector<double> f1 (f.size());
  f1 = f;
  A1.apply_constraints (f1, true);

  solver.solve (A1, u, f1, fprec);

  for (std::map<types::global_dof_index,double>::const_iterator i=bv.begin();
       i!=bv.end(); ++i)
    AssertThrow (std::fabs(u(i->first) - i->second) < 1e-8,
                 ExcInternalError());
}



template <int dim>
void
solve_eliminated (std::map<types::global_dof_index,double> &bv,
                  SparseMatrix<double>          &A,
                  Vector<double>                &u,
                  Vector<double>                &f)
{
  MatrixTools::apply_boundary_values (bv, A, u, f);

  SolverControl control (1000, 1.e-10, false, false);
  PrimitiveVectorMemory<Vector<double> > mem;
  SolverCG<Vector<double> > solver (control, mem);
  PreconditionJacobi<> prec;
  prec.initialize (A, 1.2);

  solver.solve (A, u, f, prec);
}



template <int dim>
void
check ()
{
  Triangulation<dim> tr;

  Functions::CosineFunction<dim> cosine;

  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.reset_manifold(0);

  tr.refine_global (5-dim);

  MappingQ<dim> mapping(2);
  FE_Q<dim> element(1);
  QGauss<dim> quadrature(4);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  FEValues<dim> fe (mapping, element, quadrature,
                    update_values | update_gradients
                    | update_quadrature_points | update_JxW_values);

  std::vector <types::global_dof_index> global_dofs (element.dofs_per_cell);
  std::vector <double> function (quadrature.size());

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
              const Tensor<1,dim> grad_v = fe.shape_grad(i,k);

              double rhs = dx * v * (function[k]);

              f(global_dofs[i]) += rhs;
              for (unsigned int j=0; j<element.dofs_per_cell; ++j)
                {
                  const Tensor<1,dim> grad_u = fe.shape_grad (j,k);
                  double el = dx * (grad_u*grad_v);
                  A.add (global_dofs[i], global_dofs[j], el);
                }
            }
        }
    }

  // interpolate boundary values
  std::map<types::global_dof_index,double> bv;
  VectorTools::interpolate_boundary_values (mapping, dof, 0, cosine, bv, std::vector<bool>());
  // the cosine has too many zero
  // values on the boundary of the
  // domain, so reset the elements to
  // some other value
  for (typename std::map<types::global_dof_index,double>::iterator i=bv.begin();
       i!=bv.end(); ++i)
    i->second = std::sin(i->second+0.5)+1.0;

  // first solve filtered. this does
  // not change the matrix
  Vector<double> u_filtered (dof.n_dofs ());
  solve_filtered (bv, A, u_filtered, f);

  // then solve by eliminating in the
  // matrix. since this changes the
  // matrix, this call must come
  // second
  Vector<double> u_eliminated (dof.n_dofs ());
  solve_eliminated<dim> (bv, A, u_eliminated, f);

  // output and check
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      deallog << u_filtered(i) << std::endl;
      Assert (std::fabs(u_filtered(i) - u_eliminated(i)) < 1e-8,
              ExcInternalError());
    };
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
