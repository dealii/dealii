// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "./simplex_grids.h"

void
test(const ReferenceCell<3> ref_cell)
{
  constexpr int              dim = 3;
  dealii::Triangulation<dim> tria;
  dealii::GridGenerator::reference_cell(tria, ref_cell);

  {
    deallog << "Reference cells before refinement given by "
               "tria.get_reference_cells():"
            << std::endl;

    const auto &ref_cells = tria.get_reference_cells();
    for (const auto &r : ref_cells)
      deallog << r.to_string() << ", ";
    deallog << std::endl;

    deallog << "Reference cells before refinement from cells:" << std::endl;
    for (const auto &cell : tria.active_cell_iterators())
      deallog << cell->reference_cell().to_string() << ", ";
    deallog << std::endl;
    deallog << std::endl;
  }

  tria.refine_global(1);

  {
    deallog << "Reference cells after refinement given by "
               "tria.get_reference_cells():"
            << std::endl;

    const auto &ref_cells = tria.get_reference_cells();
    for (const auto &r : ref_cells)
      deallog << r.to_string() << ", ";
    deallog << std::endl;

    deallog << "Reference cells after refinement from cells:" << std::endl;
    for (const auto &cell : tria.active_cell_iterators())
      deallog << cell->reference_cell().to_string() << ", ";
    deallog << std::endl;
    deallog << std::endl;
  }
}


int
main()
{
  initlog();

  test(ReferenceCells::Tetrahedron);
  deallog << "======================================" << std::endl << std::endl;
  test(ReferenceCells::Pyramid);
  deallog << "======================================" << std::endl << std::endl;
  test(ReferenceCells::Wedge);
  deallog << "======================================" << std::endl << std::endl;
  test(ReferenceCells::Hexahedron);
}
