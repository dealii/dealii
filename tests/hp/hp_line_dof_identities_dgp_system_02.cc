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



// check FESystem(FE_DGP)::hp_line_dof_identities, but with a different
// arrangement of base elements and multiplicities than in the _01 test


#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;
  for (unsigned int i = 1; i < 8 - dim; ++i)
    {
      // add the system three times, with
      // different numbers of base elements
      // and multiplicities
      fe_collection.push_back(FESystem<dim>(FE_DGP<dim>(i), 3));
      fe_collection.push_back(
        FESystem<dim>(FE_DGP<dim>(i), 2, FE_DGP<dim>(i), 1));
      fe_collection.push_back(
        FESystem<dim>(FE_DGP<dim>(i), 1, FE_DGP<dim>(i), 2));
    }

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

        // make sure the identities are the
        // same whatever the arrangement of
        // the elements to their base
        // elements (the operation i/3*3
        // brings us back to the first of the
        // three identical elements in the
        // collection)
        Assert(identities == fe_collection[i / 3 * 3].hp_line_dof_identities(
                               fe_collection[j / 3 * 3]),
               ExcInternalError());
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
