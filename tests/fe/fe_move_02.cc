// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check that we can move a FESystem<dim> object in a reasonable way.

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  deallog << "components: " << fe.n_components() << std::endl;
  deallog << "blocks: " << fe.n_blocks() << std::endl;
  deallog << "conforms H1: " << fe.conforms(FiniteElementData<dim>::H1)
          << std::endl;
  deallog << "n_base_elements: " << fe.n_base_elements() << std::endl;
  deallog << std::endl;
}

template <int dim>
void
move()
{
  FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  check(fe);
  FESystem<dim> fe2(std::move(fe));
  check(fe);
  check(fe2);
}

int
main()
{
  initlog();

  move<1>();
  move<2>();
  move<3>();

  return 0;
}
