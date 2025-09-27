// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the functions in integrators/divergence.h
// Output matrices and assert consistency of residuals
#include <deal.II/base/vector_slice.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/integrators/divergence.h>

#include "../tests.h"

#include "../test_grids.h"

using namespace LocalIntegrators::Divergence;


template <int dim>
void
test_cell(const FEValuesBase<dim> &fev, const FEValuesBase<dim> &fes)
{
  AssertDimension(fev.get_fe().n_components(), dim);
  AssertDimension(fes.get_fe().n_components(), 1);

  const unsigned int nv = fev.dofs_per_cell;
  const unsigned int ns = fes.dofs_per_cell;

  FullMatrix<double> Md(ns, nv);
  cell_matrix(Md, fev, fes);
  Md.print(deallog, 16, 8);

  deallog << "gradient_matrix" << std::endl;
  FullMatrix<double> Mg(nv, ns);
  gradient_matrix(Mg, fes, fev);
  Mg.print(deallog, 16, 8);


  Vector<double>                           u(nv), v(ns), w(ns);
  std::vector<std::vector<Tensor<1, dim>>> ugrad(
    dim, std::vector<Tensor<1, dim>>(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(nv);
  for (unsigned int i = 0; i < nv; ++i)
    indices[i] = i;

  deallog << "Divergence-Residuals";
  for (unsigned int i = 0; i < nv; ++i)
    {
      u    = 0.;
      u(i) = 1.;
      w    = 0.;
      fev.get_function_gradients(u, indices, make_array_view(ugrad), true);
      cell_residual<dim>(w, fes, ugrad);
      Md.vmult(v, u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
    }
  deallog << std::endl;

  v.reinit(nv);
  indices.resize(ns);
  for (unsigned int i = 0; i < ns; ++i)
    indices[i] = i;

  deallog << "Gradient-Residuals";
  for (unsigned int i = 0; i < ns; ++i)
    {
      w    = 0.;
      w(i) = 1.;
      u    = 0.;
      fes.get_function_gradients(w, indices, ugrad[0]);
      gradient_residual(u, fev, ugrad[0]);
      Mg.vmult(v, w);
      u.add(-1., v);
      deallog << ' ' << u.l2_norm();
    }
  deallog << std::endl;
}


template <int dim>
void
test_boundary(const FEValuesBase<dim> &fev, const FEValuesBase<dim> &fes)
{
  AssertDimension(fev.get_fe().n_components(), dim);
  AssertDimension(fes.get_fe().n_components(), 1);

  const unsigned int nv = fev.dofs_per_cell;
  const unsigned int ns = fes.dofs_per_cell;

  FullMatrix<double> M(ns, nv);
  u_dot_n_matrix(M, fev, fes);
  M.print(deallog, 16, 8);

  Vector<double>                   u(nv), v(ns), w(ns);
  std::vector<std::vector<double>> uval(
    dim, std::vector<double>(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(nv);
  for (unsigned int i = 0; i < nv; ++i)
    indices[i] = i;

  deallog << "Residuals u dot n";
  for (unsigned int i = 0; i < nv; ++i)
    {
      u    = 0.;
      u(i) = 1.;
      w    = 0.;
      fev.get_function_values(u, indices, make_array_view(uval), true);
      u_dot_n_residual(w, fev, fes, uval);
      M.vmult(v, u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
    }
  deallog << std::endl;
}


template <int dim>
void
test_face(const FEValuesBase<dim> &fev1,
          const FEValuesBase<dim> &fev2,
          const FEValuesBase<dim> &fes1,
          const FEValuesBase<dim> &fes2)
{
  AssertDimension(fev1.get_fe().n_components(), dim);
  AssertDimension(fes1.get_fe().n_components(), 1);
  AssertDimension(fev2.get_fe().n_components(), dim);
  AssertDimension(fes2.get_fe().n_components(), 1);

  const unsigned int nv1 = fev1.dofs_per_cell;
  const unsigned int nv2 = fev1.dofs_per_cell;
  const unsigned int ns1 = fes1.dofs_per_cell;
  const unsigned int ns2 = fes1.dofs_per_cell;

  FullMatrix<double> M11(ns1, nv1);
  FullMatrix<double> M12(ns1, nv2);
  FullMatrix<double> M21(ns2, nv1);
  FullMatrix<double> M22(ns2, nv2);

  u_dot_n_matrix(M11, M12, M21, M22, fev1, fev2, fes1, fes2);
  deallog << "M11" << std::endl;
  M11.print(deallog, 16, 8);
  deallog << "M22" << std::endl;
  M22.print(deallog, 16, 8);
  deallog << "M12" << std::endl;
  M12.print(deallog, 16, 8);
  deallog << "M21" << std::endl;
  M21.print(deallog, 16, 8);
}


template <int dim>
void
test_fe(Triangulation<dim> &tr, FiniteElement<dim> &fv, FiniteElement<dim> &fs)
{
  deallog << fv.get_name() << " x " << fs.get_name() << std::endl
          << "cell matrix" << std::endl;
  QGauss<dim>   quadrature(fv.tensor_degree() + 1);
  FEValues<dim> fev(fv,
                    quadrature,
                    update_values | update_gradients | update_JxW_values);
  FEValues<dim> fes(fs,
                    quadrature,
                    update_values | update_gradients | update_JxW_values);

  typename Triangulation<dim>::cell_iterator cell1 = tr.begin(1);
  fev.reinit(cell1);
  fes.reinit(cell1);
  test_cell(fev, fes);

  QGauss<dim - 1>   face_quadrature(fv.tensor_degree() + 1);
  FEFaceValues<dim> fev1(fv,
                         face_quadrature,
                         update_values | update_gradients |
                           update_normal_vectors | update_JxW_values);
  FEFaceValues<dim> fes1(fs,
                         face_quadrature,
                         update_values | update_gradients |
                           update_normal_vectors | update_JxW_values);
  for (const unsigned int i : GeometryInfo<dim>::face_indices())
    {
      deallog << "boundary_matrix " << i << std::endl;
      fev1.reinit(cell1, i);
      fes1.reinit(cell1, i);
      test_boundary(fev1, fes1);
    }

  FEFaceValues<dim>                          fev2(fv,
                         face_quadrature,
                         update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim>                          fes2(fs,
                         face_quadrature,
                         update_values | update_gradients | update_JxW_values);
  typename Triangulation<dim>::cell_iterator cell2 = cell1->neighbor(1);

  deallog << "face_matrix " << 0 << std::endl;
  cell1 = tr.begin(1);
  fev1.reinit(cell1, 1);
  fev2.reinit(cell2, 0);
  fes1.reinit(cell1, 1);
  fes2.reinit(cell2, 0);
  test_face(fev1, fev2, fes1, fes2);
}


template <int dim>
void
test(Triangulation<dim> &tr)
{
  FE_DGQ<dim>           q1(1);
  FE_Nedelec<dim>       n1(1);
  FE_RaviartThomas<dim> r1(1);
  FESystem<dim>         s1(q1, dim);

  test_fe(tr, r1, q1);
  test_fe(tr, n1, q1);
  test_fe(tr, s1, q1);
}


int
main()
{
  const std::string logname = "output";
  std::ofstream     logfile(logname);
  deallog.attach(logfile);

  Triangulation<2> tr2;
  TestGrids::hypercube(tr2, 1);
  test(tr2);
}
