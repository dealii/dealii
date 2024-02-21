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

// Test FiniteElement::get_sub_fe(), example used in documentation

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
check()
{
  FESystem<2> fe_velocity(FE_Q<2>(2), 2);
  FE_Q<2>     fe_pressure(1);
  FE_DGP<2>   fe_dg(0);
  FE_BDM<2>   fe_nonprim(1);

  FESystem<2> fe(fe_velocity, 1, fe_pressure, 1, fe_dg, 2, fe_nonprim, 1);

  // same using component masks to copy over:
  auto run = [&](const unsigned int first,
                 const unsigned int n,
                 const std::string &desc) {
    const unsigned int n_components = fe.n_components();

    ComponentMask mask(n_components, false);
    for (unsigned int i = first; i < first + n; ++i)
      mask.set(i, true);

    deallog << "<tr>\n"
            << "<td><code>" << mask << "</code></td>\n"
            << "<td><code>" << fe.get_sub_fe(mask).get_name()
            << "</code></td>\n"
            << "<td>" << desc << "</td>\n"
            << "</tr>\n";

    // we should be able to use ComponentMask or indices:
    Assert(fe.get_sub_fe(mask) == fe.get_sub_fe(first, n), ExcInternalError());
  };

  deallog
    << "\n<table>\n"
    << "<tr>\n<th>ComponentMask</th>\n<th>Result</th>\n<th>Description</th>\n</tr>\n";
  run(0, 7, "@p fe itself, the whole @p FESystem");
  run(0, 2, "just the @p fe_velocity");
  run(0, 1, "The first component in @p fe_velocity");
  run(1, 1, "The second component in @p fe_velocity");
  run(2, 1, "@p fe_pressure");
  run(3, 1, "first copy of @p fe_dg");
  run(4, 1, "second copy of @p fe_dg");
  run(5, 2, "both components of @p fe_nonprim");
  deallog << "</table>" << std::endl;
}

int
main()
{
  initlog();

  check<2>();
}
