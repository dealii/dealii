// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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



// check FESystem(FE_Q,FE_DGQ)::hp_line_dof_identities. i thought this was
// what would show what happened crash_14, but it turned out to be a wrong
// lead. doesn't matter, keep it anyway :-)


#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;
  for (unsigned int k = 1; k < 4; ++k)
    {
      for (unsigned int i = 1; i < 4; ++i)
        fe_collection.push_back(
          FESystem<dim>(FE_Q<dim>(k), 1, FE_DGQ<dim>(i), 1));

      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        for (unsigned int j = 0; j < fe_collection.size(); ++j)
          {
            const std::vector<std::pair<unsigned int, unsigned int>>
              identities =
                fe_collection[i].hp_line_dof_identities(fe_collection[j]);

            deallog << "Identities for " << fe_collection[i].get_name()
                    << " and " << fe_collection[j].get_name() << ": "
                    << identities.size() << std::endl;

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
}



int
main()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
