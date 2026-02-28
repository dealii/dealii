// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::face_indices_by_type()

#include <deal.II/base/geometry_info.h>

#include <deal.II/grid/reference_cell.h>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../tests.h"

template <unsigned int dim>
std::vector<ReferenceCell>
get_all_ref_cells()
{
  switch (dim)
    {
      case 0:
        return {ReferenceCells::Vertex};
      case 1:
        return {ReferenceCells::Line};
      case 2:
        return {ReferenceCells::Triangle, ReferenceCells::Quadrilateral};
      case 3:
        return {ReferenceCells::Tetrahedron,
                ReferenceCells::Pyramid,
                ReferenceCells::Wedge,
                ReferenceCells::Hexahedron};
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}

template <unsigned int dim>
void
test_generic()
{
  for (auto cell_type : get_all_ref_cells<dim>())
    {
      deallog << std::setw(7) << cell_type.to_string() << " ";
      for (auto ref_choice : RefinementCase<dim>::all_refinement_cases())
        deallog << static_cast<int>(ref_choice) << " "
                << cell_type.n_children(ref_choice) << "   ";
      deallog << std::endl;
    }
}

void
test_special_tet()
{
  deallog << std::setw(7) << "Tet"
          << " ";
  auto choices = {
    // let's not test isotropic_refinement because it really doesn't match the
    // code structure (see https://github.com/dealii/dealii/pull/19261)
    // IsotropicRefinementChoice::isotropic_refinement,
    IsotropicRefinementChoice::cut_tet_68,
    IsotropicRefinementChoice::cut_tet_57,
    IsotropicRefinementChoice::cut_tet_49,
  };
  for (auto ref_choice : choices)
    deallog << static_cast<int>(ref_choice) << " "
            << ReferenceCells::Tetrahedron.n_children(
                 RefinementCase<3>(static_cast<std::uint8_t>(ref_choice)))
            << "   ";
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << std::left;

  // test_generic<0>(); // no refinement for vertex
  test_generic<1>();
  test_generic<2>();
  test_generic<3>();

  deallog << std::endl
          << std::endl
          << "Special cases" << std::endl
          << std::endl;
  test_special_tet();
}
