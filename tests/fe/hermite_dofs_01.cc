// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/fe/fe_hermite.h>

#include "../tests.h"



/**
 * Test case to check the function
 * FE_Hermite::get_dofs_corresponding_to_outward_normal_derivatives().
 */


template <int dim>
void
print_hermite_face_dof_tables()
{
  for (unsigned int degree = 1; degree < 10 - dim; degree += 2)
    {
      FE_Hermite<dim> herm(degree);
      deallog << "Testing " << herm.get_name() << std::endl;
      for (unsigned int r = 0; r < (degree + 1) / 2; ++r)
        {
          const Table<2, unsigned int> table =
            herm.get_dofs_corresponding_to_outward_normal_derivatives(r);
          deallog << "Checking derivative order " << r << std::endl;
          for (unsigned int i = 0; i < table.size(0); ++i)
            {
              deallog << "Face " << i << ": ";
              for (unsigned int j = 0; j < table.size(1); ++j)
                deallog << table(i, j) << " ";
              deallog << std::endl;
            }
        }
      deallog << std::endl;
    }
}



int
main()
{
  std::ofstream logfile("output");

  deallog.attach(logfile);

  print_hermite_face_dof_tables<1>();
  print_hermite_face_dof_tables<2>();
  print_hermite_face_dof_tables<3>();

  return 0;
}
