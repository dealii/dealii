// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this is the same as system_02 but using the variadic constructor:
// check what happens with an FE_System if we hand it 0 components

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(FESystem<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  deallog << "components: " << fe.n_components() << std::endl;
  deallog << "blocks: " << fe.n_blocks() << std::endl;
  deallog << "conforms H1: " << fe.conforms(FiniteElementData<dim>::H1)
          << std::endl;
  deallog << "n_base_elements: " << fe.n_base_elements() << std::endl;
}

int
main()
{
  initlog();

  {
    FESystem<2> fe = {FE_Q<2>(1) ^ 2, FE_Q<2>(1) ^ 1};
    check<2>(fe);
  }
  {
    FESystem<2> fe = {FE_Q<2>(1) ^ 2, FE_DGQ<2>(2) ^ 0, FE_Q<2>(1) ^ 1};
    check<2>(fe);
  }
  {
    FESystem<2> fe = {FESystem<2>(FE_Q<2>(1) ^ 2) ^ 1, FE_Q<2>(1) ^ 1};
    check<2>(fe);
  }
  {
    FESystem<2> fe = {FESystem<2>(FE_Q<2>(1) ^ 2) ^ 1,
                      FE_DGQ<2>(2) ^ 0,
                      FE_Q<2>(1) ^ 1};
    check<2>(fe);
  }

  return 0;
}
