// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test lexicographic renumbering on a patch

#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
test(const unsigned int fe_degree, unsigned int n_components)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);

  deallog << "Testing dim=" << dim << ", degree=" << fe_degree
          << ", components=" << n_components << std::endl;

  FESystem<dim>   fe(FE_Q<dim>(fe_degree), n_components);
  DoFHandler<dim> dof_handler(tria);

  const unsigned int n_dof_1d = (2 * fe_degree + 1);

  dof_handler.distribute_dofs(fe);
  DoFRenumbering::lexicographic(dof_handler);
  DoFRenumbering::component_wise(dof_handler);

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  std::map<types::global_dof_index, Point<dim>> support_points;
  MappingQ1<dim>                                mapping;
  DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
  deallog << "Support points:" << std::endl;
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      deallog << "  " << i << " (" << support_points[i] << ")";
      if ((i + 1) % n_dof_1d == 0)
        deallog << std::endl;
    }

  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.precision(2);

  test<2>(1, 1);
  test<2>(2, 1);

  test<2>(1, 2);
  test<2>(2, 2);

  test<3>(2, 1);

  return 0;
}
