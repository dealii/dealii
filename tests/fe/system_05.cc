// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Document wrong behavior with multiplicity 0 FESystem and incorrectly reported
// degree

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

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

  deallog << "degree: " << fe.degree << std::endl;
}

int
main()
{
  initlog();

  {
    FESystem<2> fe(FE_Q<2>(1) ^ 2, FE_Q<2>(1), FE_Q<2>(1), FE_Q<2>(2) ^ 0);
    check<2>(fe);
    Assert(fe.degree == 1, ExcInternalError());
    auto copy = fe.clone();
    check<2>(*copy.get());
    Assert(fe.degree == 1, ExcInternalError());
  }
  return 0;
}
