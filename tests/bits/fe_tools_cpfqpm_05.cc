// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <memory>
#include <string>

#include "../tests.h"

// check
//   FETools::compute_projection_from_quadrature_points_matrix
// we put this into the fe_tools_common framework for simplicity, but
// in fact we ignore the second FE it passes to the check_this() function
// and we can only test as well for scalar elements, since this is all
// the function presently supports.
//
// this test makes sure that projecting onto a finite element space
// sufficiently fine to hold the quadrature point data, then interpolating
// back to the quadrature points is an identity operation
//
// This is the same as fe_tools_cpfqpm_02 but testing dim<spacedim

// output some indicators for a given matrix. we don't write out the
// entire matrix since this would blow up our output files beyond
// reasonable limits
void
output_matrix(const FullMatrix<double> &m)
{
  if ((m.m() == 0) || (m.n() == 0))
    {
      deallog << "(Empty matrix)" << std::endl;
      return;
    }

  deallog << m.l1_norm() << ' ' << m.linfty_norm() << std::endl;
  if (m.m() == m.n())
    deallog << m.frobenius_norm() << std::endl;

  for (unsigned int i = 0; i < std::min(m.m(), m.n()); ++i)
    deallog << m(i, i) << ' ' << m(i, std::min(m.m(), m.n()) - i - 1) << ' ';
  deallog << std::endl;
}



template <int dim, int spacedim>
void
check_this(const FiniteElement<dim, spacedim> &fe,
           const FiniteElement<dim, spacedim> & /*fe2*/)
{
  // only check if both elements have
  // support points. otherwise,
  // interpolation doesn't really
  // work
  if (fe.n_components() != 1)
    return;

  // ignore this check if this fe has already
  // been treated
  static std::set<std::string> already_checked;
  if (already_checked.find(fe.get_name()) != already_checked.end())
    return;
  already_checked.insert(fe.get_name());


  // test with the same quadrature formulas
  // of a degree that is high enough to
  // exactly capture the data
  QGauss<dim> q_lhs(fe.degree + 1);
  QGauss<dim> q_rhs(fe.degree + 1);

  // this test can only succeed if there are
  // at least as many degrees of freedom in
  // the finite element as there are
  // quadrature points
  if (fe.dofs_per_cell < q_rhs.size())
    return;

  deallog << "dofs_per_cell=" << fe.dofs_per_cell
          << ", n_q_points=" << q_rhs.size() << std::endl;

  FullMatrix<double> X(fe.dofs_per_cell, q_rhs.size());

  FETools::compute_projection_from_quadrature_points_matrix(fe,
                                                            q_lhs,
                                                            q_rhs,
                                                            X);

  // then compute the matrix that
  // interpolates back to the quadrature
  // points
  FullMatrix<double> I_q(q_rhs.size(), fe.dofs_per_cell);
  FETools::compute_interpolation_to_quadrature_points_matrix(fe, q_rhs, I_q);

  FullMatrix<double> product(q_rhs.size(), q_rhs.size());
  I_q.mmult(product, X);

  // the product should be the identity
  // matrix now. make sure that this is
  // indeed the case
  for (unsigned int i = 0; i < product.m(); ++i)
    product(i, i) -= 1;

  output_matrix(product);
  AssertThrow(product.frobenius_norm() < 1e-10, ExcInternalError());
}


template <int dim, int spacedim>
void
check(const FiniteElement<dim, spacedim> &fe1,
      const FiniteElement<dim, spacedim> &fe2,
      const std::string &                 name)
{
  deallog << "Checking " << name << " in " << dim << "d and spacedim "
          << spacedim << std::endl;

  // call main function in .cc files
  check_this(fe1, fe2);
}



#define CHECK(EL1, deg1, EL2, deg2, dim, spacedim)      \
  {                                                     \
    FE_##EL1<dim, spacedim> fe1(deg1);                  \
    FE_##EL2<dim, spacedim> fe2(deg2);                  \
    check(fe1, fe2, #EL1 #deg1 " against " #EL2 #deg2); \
    check(fe2, fe1, #EL2 #deg2 " against " #EL1 #deg1); \
  }

#define CHECK_ALL(EL1, deg1, EL2, deg2) \
  {                                     \
    CHECK(EL1, deg1, EL2, deg2, 1, 2);  \
    CHECK(EL1, deg1, EL2, deg2, 1, 3);  \
    CHECK(EL1, deg1, EL2, deg2, 2, 3);  \
  }


int
main()
{
  try
    {
      initlog();
      deallog << std::setprecision(8);
      deallog.depth_console(0);

      CHECK_ALL(Q, 1, Q, 1);
      CHECK_ALL(Q, 1, Q, 2);
      CHECK_ALL(Q, 1, Q, 3);
      CHECK_ALL(Q, 2, Q, 2);
      CHECK_ALL(Q, 2, Q, 3);
      CHECK_ALL(Q, 3, Q, 3);

      CHECK_ALL(DGQ, 0, DGQ, 0);
      CHECK_ALL(DGQ, 0, DGQ, 1);
      CHECK_ALL(DGQ, 0, DGQ, 2);
      CHECK_ALL(DGQ, 0, DGQ, 4);
      CHECK_ALL(DGQ, 1, DGQ, 1);
      CHECK_ALL(DGQ, 1, DGQ, 3);
      CHECK_ALL(DGQ, 2, DGQ, 2);
      CHECK_ALL(DGQ, 2, DGQ, 2);
      CHECK_ALL(DGQ, 2, DGQ, 4);
      CHECK_ALL(DGQ, 3, DGQ, 3);

      CHECK_ALL(DGP, 0, DGP, 0);
      CHECK_ALL(DGP, 0, DGP, 1);
      CHECK_ALL(DGP, 0, DGP, 2);
      CHECK_ALL(DGP, 0, DGP, 4);
      CHECK_ALL(DGP, 1, DGP, 1);
      CHECK_ALL(DGP, 1, DGP, 3);
      CHECK_ALL(DGP, 2, DGP, 2);
      CHECK_ALL(DGP, 2, DGP, 2);
      CHECK_ALL(DGP, 2, DGP, 4);
      CHECK_ALL(DGP, 3, DGP, 3);

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
}
