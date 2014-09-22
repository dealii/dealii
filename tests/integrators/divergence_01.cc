// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


// Test the functions in integrators/divergence.h
// Output matrices and assert consistency of residuals

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>
#include <deal.II/integrators/divergence.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

using namespace LocalIntegrators::Divergence;

template <int dim>
void test_cell(const FEValuesBase<dim> &fev, const FEValuesBase<dim> &fes)
{
  AssertDimension(fev.get_fe().n_components(), dim);
  AssertDimension(fes.get_fe().n_components(), 1);

  const unsigned int nv = fev.dofs_per_cell;
  const unsigned int ns = fes.dofs_per_cell;

  FullMatrix<double> Md(ns,nv);
  cell_matrix(Md,fev,fes);
  Md.print(deallog,8);

  deallog << "gradient_matrix" << std::endl;
  FullMatrix<double> Mg(nv,ns);
  gradient_matrix(Mg,fes,fev);
  Mg.print(deallog,8);


  Vector<double> u(nv), v(ns), w(ns);
  std::vector<std::vector<Tensor<1,dim> > >
  ugrad(dim,std::vector<Tensor<1,dim> >(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(nv);
  for (unsigned int i=0; i<nv; ++i)
    indices[i] = i;

  deallog << "Divergence-Residuals";
  for (unsigned int i=0; i<nv; ++i)
    {
      u = 0.;
      u(i) = 1.;
      w = 0.;
      fev.get_function_gradients(u, indices, ugrad, true);
      cell_residual(w, fes, make_slice(ugrad));
      Md.vmult(v,u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
    }
  deallog << std::endl;

  v.reinit(nv);
  indices.resize(ns);
  for (unsigned int i=0; i<ns; ++i)
    indices[i] = i;

  deallog << "Gradient-Residuals";
  for (unsigned int i=0; i<ns; ++i)
    {
      w = 0.;
      w(i) = 1.;
      u = 0.;
      fes.get_function_gradients(w, indices, ugrad[0]);
      gradient_residual(u, fev, ugrad[0]);
      Mg.vmult(v,w);
      u.add(-1., v);
      deallog << ' ' << u.l2_norm();
    }
  deallog << std::endl;
}


template <int dim>
void test_boundary(const FEValuesBase<dim> &fev, const FEValuesBase<dim> &fes)
{
  AssertDimension(fev.get_fe().n_components(), dim);
  AssertDimension(fes.get_fe().n_components(), 1);

  const unsigned int nv = fev.dofs_per_cell;
  const unsigned int ns = fes.dofs_per_cell;

  FullMatrix<double> M(ns,nv);
  u_dot_n_matrix(M, fev, fes);
  M.print(deallog,8);

  Vector<double> u(nv), v(ns), w(ns);
  std::vector<std::vector<double> >
  uval    (dim,std::vector<double>(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(nv);
  for (unsigned int i=0; i<nv; ++i)
    indices[i] = i;

  deallog << "Residuals u dot n";
  for (unsigned int i=0; i<nv; ++i)
    {
      u = 0.;
      u(i) = 1.;
      w = 0.;
      fev.get_function_values(u, indices, uval, true);
      u_dot_n_residual(w, fev, fes, make_slice(uval));
      M.vmult(v,u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
    }
  deallog << std::endl;
}


template <int dim>
void test_face(const FEValuesBase<dim> &fev1,
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

  FullMatrix<double> M11(ns1,nv1);
  FullMatrix<double> M12(ns1,nv2);
  FullMatrix<double> M21(ns2,nv1);
  FullMatrix<double> M22(ns2,nv2);

  u_dot_n_matrix(M11, M12, M21, M22, fev1, fev2, fes1, fes2);
  deallog << "M11" << std::endl;
  M11.print(deallog,8);
  deallog << "M22" << std::endl;
  M22.print(deallog,8);
  deallog << "M12" << std::endl;
  M12.print(deallog,8);
  deallog << "M21" << std::endl;
  M21.print(deallog,8);

  // Vector<double> u1(n1), v1(n1), w1(n1);
  // Vector<double> u2(n2), v2(n2), w2(n2);
  // std::vector<std::vector<double> >
  //   u1val (d,std::vector<double>(fev1.n_quadrature_points)),
  //   u2val (d,std::vector<double>(fev2.n_quadrature_points));
  // std::vector<std::vector<Tensor<1,dim> > >
  //   u1grad(d,std::vector<Tensor<1,dim> >(fev1.n_quadrature_points)),
  //   u2grad(d,std::vector<Tensor<1,dim> >(fev2.n_quadrature_points));

  // std::vector<unsigned int> indices1(n1), indices2(n2);
  // for (unsigned int i=0;i<n1;++i) indices1[i] = i;
  // for (unsigned int i=0;i<n2;++i) indices2[i] = i;

  // deallog << "Residuals";  for (unsigned int i=0;i<n1;++i) indices1[i] = i;

  // for (unsigned int i1=0;i1<n1;++i1)
  //   for (unsigned int i2=0;i2<n2;++i2)
  //     {
  //  u1 = 0.;
  //  u1(i1) = 1.;
  //  w1 = 0.;
  //  fev1.get_function_values(u1, indices1, u1val, true);
  //  fev1.get_function_gradients(u1, indices1, u1grad, true);
  //  u2 = 0.;
  //  u2(i2) = 1.;
  //  w2 = 0.;
  //  fev2.get_function_values(u2, indices2, u2val, true);
  //  fev2.get_function_gradients(u2, indices2, u2grad, true);
  //  ip_residual(w1, w2, fev1, fev2,
  //        make_slice(u1val), make_slice(u1grad),
  //        make_slice(u2val), make_slice(u2grad),
  //        17);
  //  M11.vmult(v1,u1); w1.add(-1., v1);
  //  M12.vmult(v1,u2); w1.add(-1., v1);
  //  M21.vmult(v2,u1); w2.add(-1., v2);
  //  M22.vmult(v2,u2); w2.add(-1., v2);
  //  deallog << ' ' << w1.l2_norm() + w2.l2_norm();
  //  if (d==1)
  //    {
  //      ip_residual(w1, w2, fev1, fev2, u1val[0], u1grad[0], u2val[0], u2grad[0], 17);
  //      M11.vmult(v1,u1); w1.add(-1., v1);
  //      M12.vmult(v1,u2); w1.add(-1., v1);
  //      M21.vmult(v2,u1); w2.add(-1., v2);
  //      M22.vmult(v2,u2); w2.add(-1., v2);
  //      deallog << " e" << w1.l2_norm() + w2.l2_norm();
  //    }
  //     }
  // deallog << std::endl;
}


template <int dim>
void
test_fe(Triangulation<dim> &tr, FiniteElement<dim> &fv, FiniteElement<dim> &fs)
{
  deallog << fv.get_name() << " x " << fs.get_name() << std::endl << "cell matrix" << std::endl;
  QGauss<dim> quadrature(fv.tensor_degree()+1);
  FEValues<dim> fev(fv, quadrature, update_values | update_gradients);
  FEValues<dim> fes(fs, quadrature, update_values | update_gradients);

  typename Triangulation<dim>::cell_iterator cell1 = tr.begin(1);
  fev.reinit(cell1);
  fes.reinit(cell1);
  test_cell(fev, fes);

  QGauss<dim-1> face_quadrature(fv.tensor_degree()+1);
  FEFaceValues<dim> fev1(fv, face_quadrature, update_values | update_gradients | update_normal_vectors);
  FEFaceValues<dim> fes1(fs, face_quadrature, update_values | update_gradients | update_normal_vectors);
  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      deallog << "boundary_matrix " << i << std::endl;
      fev1.reinit(cell1, i);
      fes1.reinit(cell1, i);
      test_boundary(fev1, fes1);
    }

  FEFaceValues<dim> fev2(fv, face_quadrature, update_values | update_gradients);
  FEFaceValues<dim> fes2(fs, face_quadrature, update_values | update_gradients);
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
  FE_DGQ<dim> q1(1);
  FE_Nedelec<dim> n1(1);
  FE_RaviartThomas<dim> r1(1);
  FESystem<dim> s1(q1,dim);

  test_fe(tr, r1, q1);
  test_fe(tr, n1, q1);
  test_fe(tr, s1, q1);

}


int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tr2;
  TestGrids::hypercube(tr2, 1);
  test(tr2);

//  Triangulation<3> tr3;
//  TestGrids::hypercube(tr3, 1);
//  test(tr3);

}
