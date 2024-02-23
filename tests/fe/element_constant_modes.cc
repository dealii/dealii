// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_system.h>

#include <string>

#include "../tests.h"


template <int dim>
void
print_constant_modes(const FiniteElement<dim> &fe)
{
  deallog << "Testing " << fe.get_name() << std::endl;

  Table<2, bool> constant_modes = fe.get_constant_modes().first;
  for (unsigned int r = 0; r < constant_modes.n_rows(); ++r)
    {
      for (unsigned int c = 0; c < constant_modes.n_cols(); ++c)
        deallog << constant_modes(r, c) << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
test()
{
  print_constant_modes(FE_Q<dim>(1));
  print_constant_modes(FE_Q<dim>(2));
  print_constant_modes(FE_DGQ<dim>(1));
  print_constant_modes(FE_DGQLegendre<dim>(2));
  print_constant_modes(FE_DGQHermite<dim>(2));
  print_constant_modes(FE_DGQHermite<dim>(3));
  print_constant_modes(FE_DGP<dim>(2));
  print_constant_modes(FE_Q_Hierarchical<dim>(1));
  print_constant_modes(FE_Q_Hierarchical<dim>(2));
  print_constant_modes(FE_FaceQ<dim>(1));
  print_constant_modes(FE_FaceP<dim>(1));
  print_constant_modes(FESystem<dim>(FE_Q<dim>(1), 2, FE_Q<dim>(2), 1));
  print_constant_modes(
    FESystem<dim>(FE_DGP<dim>(1), 1, FE_Q_iso_Q1<dim>(2), 1));
  print_constant_modes(FE_Q_DG0<dim>(1));
  print_constant_modes(FESystem<dim>(FE_Q_DG0<dim>(2), 1, FE_Q<dim>(1), 2));
  print_constant_modes(FESystem<dim>(FE_Q<dim>(1), 2, FE_Q_DG0<dim>(1), 2));
}

template <>
void
test<1>()
{
  print_constant_modes(FE_Q<1>(1));
  print_constant_modes(FESystem<1>(FE_Q<1>(1), 2, FE_Q<1>(2), 1));
  print_constant_modes(FESystem<1>(FE_DGP<1>(1), 1, FE_Q_iso_Q1<1>(2), 1));
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
