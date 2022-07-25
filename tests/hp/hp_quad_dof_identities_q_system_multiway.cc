// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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



// check FE_Collection::hp_quad_dof_identities with multiway
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
