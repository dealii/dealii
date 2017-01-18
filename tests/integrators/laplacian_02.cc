// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
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


// Test the tangential functions in integrators/laplace.h
// Output matrices and assert consistency of residuals

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/integrators/laplace.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

using namespace LocalIntegrators::Laplace;



template <int dim>
void test_boundary(const FiniteElement<dim> &fe, bool diff=false)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_zero_boundary_constraints(dof, constraints);

  QGauss<dim-1> quadrature(fe.tensor_degree()+1);
  FEFaceValues<dim> fev(fe, quadrature, update_values | update_gradients | update_normal_vectors | update_JxW_values);

  FullMatrix<double> M(fe.dofs_per_cell);
  FullMatrix<double> Mglobal(dof.n_dofs());
  std::vector<unsigned int> indices(fe.dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  cell->get_dof_indices(indices);
  for (unsigned i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      fev.reinit(cell, i);
      nitsche_tangential_matrix(M, fev, 10.);
      if (diff)
        nitsche_matrix(M, fev, 10., -1.);
      constraints.distribute_local_to_global(M, indices, indices, Mglobal);
    }
  deallog << fe.get_name() << ": bdry norm " << Mglobal.frobenius_norm() << std::endl;
}



template <int dim>
void test_face(const FiniteElement<dim> &fe, bool diff=false)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global();

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_zero_boundary_constraints(dof, constraints);

  QGauss<dim-1> quadrature(fe.tensor_degree()+1);
  FEFaceValues<dim> fev1(fe, quadrature, update_values | update_gradients | update_normal_vectors | update_JxW_values);
  FEFaceValues<dim> fev2(fe, quadrature, update_values | update_gradients | update_normal_vectors | update_JxW_values);

  FullMatrix<double> M11(fe.dofs_per_cell);
  FullMatrix<double> M12(fe.dofs_per_cell);
  FullMatrix<double> M21(fe.dofs_per_cell);
  FullMatrix<double> M22(fe.dofs_per_cell);
  FullMatrix<double> Mglobal(dof.n_dofs());
  std::vector<unsigned int> indices1(fe.dofs_per_cell);
  std::vector<unsigned int> indices2(fe.dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell1 = dof.begin_active();
  typename DoFHandler<dim>::active_cell_iterator cell2 = ++dof.begin_active();

  cell1->get_dof_indices(indices1);
  cell2->get_dof_indices(indices2);

  fev1.reinit(cell1, 1);
  fev2.reinit(cell2, 0);
  ip_tangential_matrix(M11, M12, M21, M22, fev1, fev2, 10);

  if (diff)
    ip_matrix(M11, M12, M21, M22, fev1, fev2, 10, -1.);

  constraints.distribute_local_to_global(M11, indices1, indices1, Mglobal);
  constraints.distribute_local_to_global(M21, indices2, indices1, Mglobal);
  constraints.distribute_local_to_global(M12, indices1, indices2, Mglobal);
  constraints.distribute_local_to_global(M22, indices2, indices2, Mglobal);
  deallog << fe.get_name() << ": face norm " << Mglobal.frobenius_norm() << std::endl;
}



template <int dim>
void
test()
{
  // In each dimension, the first four outputs should be zero: the
  // tangential jump of Nedelec elements is zero and for
  // Raviart-Thomas, the tangential and the full jump are the same.
  FE_Nedelec<dim> n1(2);
  test_boundary(n1);
  test_face(n1);

  FE_RaviartThomas<dim> r1(2);
  test_boundary(r1, true);
  test_face(r1,true);
  test_boundary(r1);
  test_face(r1);

  FE_DGQ<dim> q1(1);
  FESystem<dim> sys1(q1,dim);
  test_boundary(sys1);
  test_boundary(sys1, true);
  test_face(sys1);
  test_face(sys1, true);
  deallog << std::endl;
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}
