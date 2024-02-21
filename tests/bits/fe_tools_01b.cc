// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// common framework for the various fe_tools_*.cc tests

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"


// check
//   FETools::hierarchic_to_lexicographic_numbering


template <int dim>
void
check(const FE_Q<dim> &fe, const std::string &name)
{
  deallog << "Checking " << name << " in " << dim << "d:" << std::endl;

  const std::vector<unsigned int> n =
    FETools::hierarchic_to_lexicographic_numbering<dim>(fe.degree);
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    deallog << n[i] << ' ';
  deallog << std::endl;
}



#define CHECK(EL, deg, dim) \
  {                         \
    FE_##EL<dim> EL(deg);   \
    check(EL, #EL #deg);    \
  }

#define CHECK_ALL(EL, deg) \
  {                        \
    CHECK(EL, deg, 1);     \
    CHECK(EL, deg, 2);     \
    CHECK(EL, deg, 3);     \
  }


int
main()
{
  try
    {
      initlog();
      deallog << std::setprecision(2);

      CHECK_ALL(Q, 1);
      CHECK_ALL(Q, 2);
      CHECK_ALL(Q, 3);
      CHECK_ALL(Q, 4);

      return 0;
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
