// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



// Same as create_laplace_matrix_0[1234] but passing an additional constraint
// matrix and comparing results without constraints



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  // create a system element composed
  // of one Q1 and one Q2 element
  FE_Q<dim> element(2);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping (3);
  QGauss<dim> quadrature(6);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  VectorTools::interpolate_boundary_values(dof, 0,
                                           ConstantFunction<dim>(1,1.),
                                           constraints);
  constraints.close ();

  // create sparsity pattern. note
  // that different components should
  // not couple, so use pattern
  SparsityPattern sparsity;
  {
    CompressedSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (dof, csp, constraints);
    sparsity.copy_from(csp);
  }

  SparseMatrix<double> matrix, matrix_ref;
  matrix.reinit (sparsity);
  matrix_ref.reinit (sparsity);

  Functions::ExpFunction<dim> rhs_function;

  Vector<double> rhs (dof.n_dofs()), rhs_ref(dof.n_dofs());

  MatrixTools::
  create_laplace_matrix (mapping, dof,
                         quadrature, matrix_ref,
                         rhs_function, rhs_ref,
                         &rhs_function);
  constraints.condense(matrix_ref, rhs_ref);

  MatrixTools::create_laplace_matrix (mapping, dof, quadrature, matrix,
                                      rhs_function, rhs, &rhs_function,
                                      constraints);

  // compute reference: need to cancel constrained entries as these will in
  // general get different values
  matrix.add(-1., matrix_ref);
  for (unsigned int i=0; i<matrix.m(); ++i)
    if (constraints.is_constrained(i)==true)
      matrix.diag_element(i) = 0;
  deallog << "Matrix error Frobenius: " << matrix.frobenius_norm() << std::endl;

  rhs -= rhs_ref;
  deallog << "RHS vector error l2: " << rhs.l2_norm() << std::endl;
}



int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);
  deallog.attach(logfile);
  deallog.threshold_double(1e-10);
  deallog.depth_console (0);

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
