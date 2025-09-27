// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// FESystem::get_unit_support_points would return an empty array if
// one base element had no support points. but that's not necessary
// that way: it should still return something useful if the element
// has no support points but in fact also has no degrees of freedom on
// faces at all. in that case we would simply not care
//
// check this by comparing
//    FESystem(FE_Q, 1,  FE_DGQ, 1)
// with
//    FESystem(FE_Q, 1,  FE_DGP, 1)
// (the latter not having any support points for the second component).

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"


template <int dim>
void
check2(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  const std::vector<Point<dim - 1>> unit_f_s_p =
    fe.get_unit_face_support_points();
  for (unsigned int i = 0; i < unit_f_s_p.size(); ++i)
    deallog << i << ' ' << unit_f_s_p[i] << std::endl;
}


template <int dim>
void
check()
{
  check2(FESystem<dim>(FE_Q<dim>(2), 1, FE_DGQ<dim>(2), 1));
  check2(FESystem<dim>(FE_Q<dim>(2), 1, FE_DGP<dim>(2), 1));
}



int
main()
{
  initlog();

  check<2>();
  check<3>();
  return 0;
}
