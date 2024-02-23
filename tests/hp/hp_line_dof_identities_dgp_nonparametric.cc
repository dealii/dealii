// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FE_DGPNonparametric::hp_line_dof_identities


#include <deal.II/fe/fe_dgp_nonparametric.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;
  for (unsigned int i = 1; i < 8 - dim; ++i)
    fe_collection.push_back(FE_DGPNonparametric<dim>(i));

  for (unsigned int i = 0; i < fe_collection.size(); ++i)
    for (unsigned int j = 0; j < fe_collection.size(); ++j)
      {
        const std::vector<std::pair<unsigned int, unsigned int>> identities =
          fe_collection[i].hp_line_dof_identities(fe_collection[j]);

        deallog << "Identities for " << fe_collection[i].get_name() << " and "
                << fe_collection[j].get_name() << ": " << identities.size()
                << std::endl;

        for (unsigned int k = 0; k < identities.size(); ++k)
          {
            Assert(identities[k].first < fe_collection[i].dofs_per_line,
                   ExcInternalError());
            Assert(identities[k].second < fe_collection[j].dofs_per_line,
                   ExcInternalError());

            deallog << identities[k].first << ' ' << identities[k].second
                    << std::endl;
          }
      }
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
