// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Verify sequences for a custom hierarchy registered on a hp::FECollection.


#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;

  // register custom hierarchy with two sequences: odd and even indices
  fe_collection.set_hierarchy(
    [](const hp::FECollection<dim> &fe_collection,
       const unsigned int           fe_index) {
      return ((fe_index + 2) < fe_collection.size()) ? fe_index + 2 : fe_index;
    },
    [](const hp::FECollection<dim> &fe_collection,
       const unsigned int           fe_index) {
      return (fe_index > 1) ? fe_index - 2 : fe_index;
    });

  // add dummy FEs to collection
  while (fe_collection.size() < 6)
    fe_collection.push_back(FE_Q<dim>(1));
  deallog << "size: " << fe_collection.size() << std::endl;

  // verify sequence for each FE index
  for (unsigned int fe_index = 0; fe_index < fe_collection.size(); ++fe_index)
    {
      const auto sequence = fe_collection.get_hierarchy_sequence(fe_index);

      deallog << " idx: " << fe_index << ", sequence:";
      for (const auto index : sequence)
        deallog << ' ' << index;
      deallog << std::endl;
    }
}


int
main()
{
  initlog();

  test<1>();

  deallog << "OK" << std::endl;
}
