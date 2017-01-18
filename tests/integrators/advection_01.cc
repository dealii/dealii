// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2016 by the deal.II authors
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


// Test the functions in integrators/laplace.h
// Output matrices and assert consistency of residuals

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>
#include <deal.II/integrators/advection.h>
#include <deal.II/integrators/laplace.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>

using namespace LocalIntegrators::Advection;

template <int dim>
void test_cell(const FEValuesBase<dim> &fev)
{
  const unsigned int n = fev.dofs_per_cell;
  unsigned int d=fev.get_fe().n_components();
  FullMatrix<double> M(n,n);
  std::vector<std::vector<double> > vel(dim,std::vector<double>(1));
  vel[0][0]=1.;
  vel[1][0]=1.;
  if (dim > 2)
    vel[2][0]=1.;
  cell_matrix(M,fev,fev,vel);
  {
    LogStream::Prefix pre("cell");
    M.print(deallog,8);
  }

  Vector<double> u(n), v(n), w(n);
  std::vector<std::vector<double > >
  uval    (d,std::vector<double>(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(n);
  for (unsigned int i=0; i<n; ++i)
    indices[i] = i;

  {
    LogStream::Prefix pre("Residuals");
    for (unsigned int i=0; i<n; ++i)
      {
        u = 0.;
        u(i) = 1.;
        w = 0.;
        fev.get_function_values(u, indices, VectorSlice<std::vector<std::vector<double> > >(uval), true);
        cell_residual(w, fev, make_slice(uval),vel);
        M.vmult(v,u);
        w.add(-1., v);
        deallog << ' ' << w.l2_norm();
        if (d==1)
          {
            cell_residual(w, fev, uval[0],vel);
            M.vmult(v,u);
            w.add(-1., v);
            deallog << " e" << w.l2_norm();
          }
      }
    deallog << std::endl;
  }
}


template <int dim>
void test_boundary(const FEValuesBase<dim> &fev)
{
  const unsigned int n = fev.dofs_per_cell;
  unsigned int d=fev.get_fe().n_components();
  FullMatrix<double> M(n,n);
  std::vector<std::vector<double> > vel(dim ,std::vector<double>(1));
  vel[0][0]=1.;
  vel[1][0]=1.;
  if (dim > 2)
    vel[2][0]=1.;
  upwind_value_matrix(M, fev, fev, vel);
  {
    LogStream::Prefix pre("bdry");
    M.print(deallog,8);
  }

  Vector<double> u(n), v(n), w(n);
  std::vector<std::vector<double> >
  uval(d, std::vector<double>(fev.n_quadrature_points));
  std::vector<std::vector<double> >
  null_val(d, std::vector<double>(fev.n_quadrature_points, 0.));
  std::vector<std::vector<Tensor<1,dim> > >
  ugrad   (d,std::vector<Tensor<1,dim> >(fev.n_quadrature_points));

  std::vector<types::global_dof_index> indices(n);
  for (unsigned int i=0; i<n; ++i)
    indices[i] = i;

  {
    LogStream::Prefix pre("Residuals");
    for (unsigned int i=0; i<n; ++i)
      {
        u = 0.;
        u(i) = 1.;
        w = 0.;
        fev.get_function_values(u, indices, VectorSlice<std::vector<std::vector<double> > >(uval), true);
        upwind_value_residual(w, fev, make_slice(uval), make_slice(null_val), vel);
        M.vmult(v,u);
        w.add(-1., v);
        deallog << ' ' << w.l2_norm();
        if (d==1)
          {
            upwind_value_residual(w, fev, uval[0], null_val[0], vel);
            M.vmult(v,u);
            w.add(-1., v);
            deallog << " e" << w.l2_norm();
          }
      }
    deallog << std::endl;
  }
}

template <int dim>
void test_face(const FEValuesBase<dim> &fev1,
               const FEValuesBase<dim> &fev2)
{
  const unsigned int n1 = fev1.dofs_per_cell;
  const unsigned int n2 = fev2.dofs_per_cell;
  unsigned int d=fev1.get_fe().n_components();
  FullMatrix<double> M11(n1,n1);
  FullMatrix<double> M12(n1,n2);
  FullMatrix<double> M21(n2,n1);
  FullMatrix<double> M22(n2,n2);
  std::vector<std::vector<double> > vel(dim ,std::vector<double>(1));
  vel[0][0]=1.;
  vel[1][0]=1.;
  if (dim > 2)
    vel[2][0]=1.;

  upwind_value_matrix(M11, M12, M21, M22, fev1, fev2, fev1, fev2, vel);

  {
    LogStream::Prefix pre("M11");
    M11.print(deallog,8);
  }
  {
    LogStream::Prefix pre("M12");
    M12.print(deallog,8);
  }
  {
    LogStream::Prefix pre("M21");
    M21.print(deallog,8);
  }
  {
    LogStream::Prefix pre("M22");
    M22.print(deallog,8);
  }

  Vector<double> u1(n1), v1(n1), w1(n1);
  Vector<double> u2(n2), v2(n2), w2(n2);
  std::vector<std::vector<double> >
  u1val (d,std::vector<double>(fev1.n_quadrature_points)),
        nullval (d,std::vector<double>(fev2.n_quadrature_points, 0.));

  std::vector<types::global_dof_index> indices1(n1), indices2(n2);
  for (unsigned int i=0; i<n1; ++i) indices1[i] = i;
  for (unsigned int i=0; i<n2; ++i) indices2[i] = i;

  {
    LogStream::Prefix pre("Residuals");
    for (unsigned int i1=0; i1<n1; ++i1)
      {
        u1 = 0.;
        u1(i1) = 1.;
        w1 = 0.;
        w2 = 0.;
        fev1.get_function_values(u1, indices1, VectorSlice<std::vector<std::vector<double> > >(u1val), true);
        upwind_face_residual(w1, w2, fev1, fev2,
                             make_slice(u1val),
                             make_slice(nullval),
                             vel);
        M11.vmult(v1,u1);
        w1.add(-1., v1);
        M21.vmult(v2,u1);
        w2.add(-1., v2);
        deallog << "sys  a " << w1.l2_norm() << ' ' << w2.l2_norm() << std::endl;
        if (d==1)
          {
            w1 = 0.;
            w2 = 0.;
            upwind_face_residual(w1, w2, fev1, fev2, u1val[0], nullval[0], vel);
            M11.vmult(v1,u1);
            w1.add(-1., v1);
            M21.vmult(v2,u1);
            w2.add(-1., v2);
            deallog << "sing a " << w1.l2_norm() << ' ' << w2.l2_norm() << std::endl;
          }
        w1 = 0.;
        w2 = 0.;
        fev2.get_function_values(u1, indices2, VectorSlice<std::vector<std::vector<double> > >(u1val), true);
        upwind_face_residual(w1, w2, fev1, fev2,
                             make_slice(nullval),
                             make_slice(u1val),
                             vel);
        M12.vmult(v1,u1);
        w1.add(-1., v1);
        M22.vmult(v2,u1);
        w2.add(-1., v2);
        deallog << "sys  b " << w1.l2_norm() << ' ' << w2.l2_norm() << std::endl;
        if (d==1)
          {
            w1 = 0.;
            w2 = 0.;
            upwind_face_residual(w1, w2, fev1, fev2, nullval[0], u1val[0], vel);
            M12.vmult(v1,u1);
            w1.add(-1., v1);
            M22.vmult(v2,u1);
            w2.add(-1., v2);
            deallog << "sing b " << w1.l2_norm() << ' ' << w2.l2_norm() << std::endl;
          }
      }
    deallog << std::endl;
  }
}


template <int dim>
void
test_fe(Triangulation<dim> &tr, FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl << "cell matrix" << std::endl;
  QGauss<dim> quadrature(fe.tensor_degree());
  FEValues<dim> fev(fe, quadrature, update_values | update_gradients | update_JxW_values);

  typename Triangulation<dim>::cell_iterator cell1 = tr.begin(1);
  fev.reinit(cell1);
  test_cell(fev);

  QGauss<dim-1> face_quadrature(fe.tensor_degree()+1);
  FEFaceValues<dim> fef1(fe, face_quadrature, update_values | update_gradients | update_normal_vectors | update_JxW_values);
  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      deallog << "boundary_matrix " << i << std::endl;
      fef1.reinit(cell1, i);
      test_boundary(fef1);
    }

  FEFaceValues<dim> fef2(fe, face_quadrature, update_values | update_gradients | update_normal_vectors | update_JxW_values);
  typename Triangulation<dim>::cell_iterator cell2 = cell1->neighbor(1);

  deallog << "face_matrix " << 0 << std::endl;
  cell1 = tr.begin(1);
  fef1.reinit(cell1, 1);
  fef2.reinit(cell2, 0);
  test_face(fef1, fef2);
  test_face(fef2, fef1);
}


template <int dim>
void
test(Triangulation<dim> &tr)
{
  FE_DGQ<dim> q1(1);
  test_fe(tr, q1);

  FE_DGQ<dim> q2(2);
  test_fe(tr, q2);

  FE_Nedelec<dim> n1(1);
  test_fe(tr, n1);
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  Triangulation<2> tr2;
  TestGrids::hypercube(tr2, 1);
  test(tr2);

  Triangulation<3> tr3;
  TestGrids::hypercube(tr3, 1);
  test(tr3);
}
