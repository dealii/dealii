// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FE_Collection::hp_line_dof_identities with multiway
// identities for a system of FE_Q elements. Each vector component of
// the system will form its own identity set.


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  hp::FECollection<dim> fe_collection;
  for (unsigned int i = 1; i <= 8; ++i)
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(i), 2));

  // construct the complete set of fe indices
  std::set<unsigned int> fe_indices;
  for (unsigned int i = 0; i < fe_collection.size(); ++i)
    fe_indices.emplace(i);

  const auto identities = fe_collection.hp_line_dof_identities(fe_indices);

  for (unsigned int i = 0; i < identities.size(); ++i)
    {
      deallog << "Identity set #" << i << std::endl;
      for (const auto &p : identities[i])
        {
          deallog << "  " << fe_collection[p.first].get_name()
                  << ": line dof index " << p.second << std::endl;
        }
    }

  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
