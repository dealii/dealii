// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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


// Test nitsche_tangential in Elasticity
// Output matrices and assert consistency of residuals

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>
#include <deal.II/integrators/elasticity.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>

using namespace LocalIntegrators::Elasticity;

template <int dim>
void test_boundary(const FEValuesBase<dim> &fev)
{
  const unsigned int n = fev.dofs_per_cell;
  unsigned int d=fev.get_fe().n_components();
  FullMatrix<double> M(n,n);
  nitsche_tangential_matrix(M, fev, 17);
  {
    LogStream::Prefix pre("bdry");
    M.print(deallog,12,8);
  }

  Vector<double> u(n), v(n), w(n);
  std::vector<std::vector<double> >
  uval    (d,std::vector<double>(fev.n_quadrature_points)),
          null_val(d,std::vector<double>(fev.n_quadrature_points, 0.));
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
        fev.get_function_gradients(u, indices, VectorSlice<std::vector<std::vector<Tensor<1,dim> > > >(ugrad), true);
        nitsche_tangential_residual(w, fev, make_slice(uval), make_slice(ugrad), make_slice(null_val), 17);
        M.vmult(v,u);
        w.add(-1., v);
        deallog << ' ' << w.l2_norm();
      }
    deallog << std::endl;
  }
}


template <int dim>
void
test_fe(Triangulation<dim> &tr, FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl << "cell matrix" << std::endl;
  typename Triangulation<dim>::cell_iterator cell1 = tr.begin(1);

  QGauss<dim-1> face_quadrature(fe.tensor_degree()+1);
  FEFaceValues<dim> fef1(fe, face_quadrature, update_values | update_gradients | update_normal_vectors);
  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      deallog << "boundary_matrix " << i << std::endl;
      fef1.reinit(cell1, i);
      test_boundary(fef1);
    }
}



template <int dim>
void
test(Triangulation<dim> &tr)
{
  FE_DGQ<dim> q1(1);
  FESystem<dim> fe1(q1,dim);
  test_fe(tr, fe1);

  FE_DGQ<dim> q2(2);
  FESystem<dim> fe2(q2,dim);
  test_fe(tr, fe2);

  FE_Nedelec<dim> n1(1);
  test_fe(tr, n1);

  FE_RaviartThomas<dim> rt1(1);
  test_fe(tr, rt1);
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);
  deallog.precision(8);

  Triangulation<2> tr2;
  TestGrids::hypercube(tr2, 1);
  test(tr2);

  Triangulation<3> tr3;
  TestGrids::hypercube(tr3, 1);
  test(tr3);
}
