/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Verify that DoFs are placed in the same locations for
// FESystem(FE_SimplexP<dim>(3)) and FE_SimplexP<dim>(3). Previously this didn't
// work in 2d since FESystem only populated line orientations in 3d (but for
// simplices we need them in 2d).

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_p1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  FE_SimplexP<dim> fe_1(3);
  FESystem<dim>    fe_2(fe_1, 1);

  DoFHandler<dim> dof_1(tria), dof_2(tria);
  dof_1.distribute_dofs(fe_1);
  dof_2.distribute_dofs(fe_2);

  deallog << "checking " << fe_1.get_name() << " vs. " << fe_2.get_name()
          << std::endl;

  deallog << "check face_to_cell_index:" << std::endl;
  for (const unsigned int &face_no : fe_1.reference_cell().face_indices())
    for (types::geometric_orientation combined_orientation = 0;
         combined_orientation <
         fe_1.reference_cell().n_face_orientations(face_no);
         ++combined_orientation)
      for (unsigned int face_dof = 0; face_dof < fe_1.max_dofs_per_face();
           ++face_dof)
        {
          Assert(
            fe_1.face_to_cell_index(face_dof, face_no, combined_orientation) ==
              fe_2.face_to_cell_index(face_dof, face_no, combined_orientation),
            ExcInternalError());
          deallog
            << fe_1.face_to_cell_index(face_dof, face_no, combined_orientation)
            << " = "
            << fe_2.face_to_cell_index(face_dof, face_no, combined_orientation)
            << std::endl;
        }

  // lines only have two orientations
  deallog << "check adjust_line_dof_index_for_line_orientation:" << std::endl;
  for (types::geometric_orientation combined_orientation = 0;
       combined_orientation < 2;
       ++combined_orientation)
    for (unsigned int line_dof = 0; line_dof < fe_1.n_dofs_per_line();
         ++line_dof)
      {
        Assert(fe_1.adjust_line_dof_index_for_line_orientation(
                 line_dof, combined_orientation) ==
                 fe_2.adjust_line_dof_index_for_line_orientation(
                   line_dof, combined_orientation),
               ExcInternalError());
        deallog << fe_1.adjust_line_dof_index_for_line_orientation(
                     line_dof, combined_orientation)
                << " = "
                << fe_2.adjust_line_dof_index_for_line_orientation(
                     line_dof, combined_orientation)
                << std::endl;
      }

  MappingP1<dim> mapping;
  auto points_1 = DoFTools::map_dofs_to_support_points(mapping, dof_1);
  auto points_2 = DoFTools::map_dofs_to_support_points(mapping, dof_2);

  for (const auto &cell : tria.active_cell_iterators())
    {
      typename DoFHandler<dim>::active_cell_iterator cell_1(&tria,
                                                            cell->level(),
                                                            cell->index(),
                                                            &dof_1),
        cell_2(&tria, cell->level(), cell->index(), &dof_2);

      std::vector<types::global_dof_index> dofs_1(fe_1.n_dofs_per_cell()),
        dofs_2(fe_2.n_dofs_per_cell());

      cell_1->get_dof_indices(dofs_1);
      cell_2->get_dof_indices(dofs_2);

      deallog << "dofs 1:" << std::endl;
      for (const auto &dof : dofs_1)
        deallog << "  " << std::setw(2) << dof << " -> " << points_1[dof]
                << std::endl;

      deallog << "dofs 2:" << std::endl;
      for (const auto &dof : dofs_2)
        deallog << "  " << std::setw(2) << dof << " -> " << points_2[dof]
                << std::endl;
    }
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}
