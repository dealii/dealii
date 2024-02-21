/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Test that when copy_triangulation is used also the information about
 * periodic faces is copied.
 */
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim           = 2;
  const unsigned int n_components  = 1;
  const unsigned int n_refinements = 1;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;

  GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);

  tria.add_periodicity(periodic_faces);
  tria.refine_global(n_refinements);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(1), n_components));

  deallog << "Tria:" << std::endl;
  for (const auto &elem : tria.get_periodic_face_map())
    {
      deallog << elem.first.first << " " << elem.first.second << " "
              << elem.second.first.first << " " << elem.second.first.second
              << std::endl;
    }

  Triangulation<dim> other_tria;
  other_tria.copy_triangulation(tria);

  deallog << "Copy of tria:" << std::endl;
  for (const auto &elem : other_tria.get_periodic_face_map())
    {
      deallog << elem.first.first << " " << elem.first.second << " "
              << elem.second.first.first << " " << elem.second.first.second
              << std::endl;
    }

  deallog << "OK" << std::endl;
}
