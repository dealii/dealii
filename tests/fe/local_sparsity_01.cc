// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check FiniteElement::get_local_dof_sparsity_pattern() for FESystem.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include "deal.II/fe/fe_system.h"
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  const auto &pattern = fe.get_local_dof_sparsity_pattern();
  deallog << fe.get_name() << ":\n";

  if (!pattern.empty())
    {
      for (unsigned int i = 0; i < pattern.size(0); ++i)
        {
          for (unsigned int j = 0; j < pattern.size(1); ++j)
            deallog << (pattern(i, j) == true ? "X" : ".");
          deallog << "\n";
        }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  // The simplest case is a single iso Q1:
  test<1>(FE_Q_iso_Q1<1>(2));
  test<1>(FE_Q_iso_Q1<1>(3));

  // It also works with custom support points:
  std::vector<Point<1>> points = {Point<1>(0.0), Point<1>(0.1), Point<1>(1.0)};
  test<1>(FE_Q_iso_Q1<1>(points));

  std::vector<Point<1>> points2 = {Point<1>(0.0),
                                   Point<1>(0.1),
                                   Point<1>(0.2),
                                   Point<1>(1.0)};
  test<1>(FE_Q_iso_Q1<1>(points2));

  // The 2 iso Q1 elements couple using their pattern:
  test<1>(FESystem<1, 1>(FE_DGQ<1>(1), 1, FE_Q_iso_Q1<1>(2), 2));
  // The coupling between the first two to the third copy is a full coupling
  // currently, because we don't detect this yet:
  test<1>(FESystem<1, 1>(FE_Q_iso_Q1<1>(2), 2, FE_Q_iso_Q1<1>(2), 1));
  // Different iso_Q1 degrees always couple fully (off diagonal blocks):
  test<1>(FESystem<1, 1>(FE_Q_iso_Q1<1>(2), 1, FE_Q_iso_Q1<1>(3), 1));
}
