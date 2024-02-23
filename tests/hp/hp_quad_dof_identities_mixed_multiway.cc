// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FE_Collection::hp_quad_dof_identities with multiway
// identities for a combination of FE_Q element and FE_DGQ
// elements. Only the FE_Q elements have DoFs on vertices so only they
// will participate in the identities.


#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  hp::FECollection<dim> fe_collection;
  for (unsigned int i = 1; i <= 8; ++i)
    {
      fe_collection.push_back(FE_Q<dim>(i));
      fe_collection.push_back(FE_DGQ<dim>(i));
    }

  // construct the complete set of fe indices
  std::set<unsigned int> fe_indices;
  for (unsigned int i = 0; i < fe_collection.size(); ++i)
    fe_indices.emplace(i);

  const auto identities = fe_collection.hp_quad_dof_identities(fe_indices);

  for (unsigned int i = 0; i < identities.size(); ++i)
    {
      deallog << "Identity set #" << i << std::endl;
      for (const auto &p : identities[i])
        {
          deallog << "  " << fe_collection[p.first].get_name()
                  << ": quad dof index " << p.second << std::endl;
        }
    }

  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<3>();

  deallog << "OK" << std::endl;
}
