// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FETools::compute_nodal_quadrature()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

template <int dim>
void
print_quadrature(const FE_Q_Base<dim> &fe, std::ostream &out)
{
  out << fe.get_name() << std::endl;
  const Quadrature<dim> nodal_quadrature =
    FETools::compute_nodal_quadrature(fe);

  out << "points:" << std::endl;
  for (const auto &point : nodal_quadrature.get_points())
    out << "  " << point << std::endl;
  out << "weights:" << std::endl;
  for (const auto &weight : nodal_quadrature.get_weights())
    out << "  " << weight << std::endl;

  // Check that we are close enough to QGaussLobatto
  const Quadrature<dim> lobatto = QGaussLobatto<dim>(fe.tensor_degree() + 1);

  for (unsigned int q = 0; q < lobatto.size(); ++q)
    {
      AssertThrow((lobatto.point(q) - nodal_quadrature.point(q)).norm() < 1e-6,
                  ExcInternalError());
      AssertThrow(std::abs(lobatto.weight(q) - nodal_quadrature.weight(q)) <
                    1e-6,
                  ExcInternalError());
    }
}

template <int dim>
void
print_quadrature(const FiniteElement<dim> &fe, std::ostream &out)
{
  out << fe.get_name() << std::endl;
  const Quadrature<dim> nodal_quadrature =
    FETools::compute_nodal_quadrature(fe);

  out << "points:" << std::endl;
  for (const auto &point : nodal_quadrature.get_points())
    out << "  " << point << std::endl;
  out << "weights:" << std::endl;
  for (const auto &weight : nodal_quadrature.get_weights())
    out << "  " << weight << std::endl;
}

int
main()
{
  initlog();

  print_quadrature(FE_SimplexP<1>(1), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;
  print_quadrature(FE_SimplexP<1>(2), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;
  print_quadrature(FE_SimplexP<2>(1), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;
  print_quadrature(FE_SimplexP<2>(2), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;

  print_quadrature(FE_WedgeP<3>(1), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;

  print_quadrature(FE_DGQ<2>(1), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;
  print_quadrature(FE_DGQ<2>(2), deallog.get_file_stream());
  deallog.get_file_stream() << std::endl;
  print_quadrature(FE_Q<3>(1), deallog.get_file_stream());
}
