// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Checks dynamic patch distributor by printing the mapping
// between cell dofs and patch dofs

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_patch_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/patch_distributors.h>
#include <deal.II/matrix_free/patch_storage.h>

#include "../tests.h"


template <int dim, int degree>
void
test()
{
  deallog << "\n Running test in " << dim << "D with degree " << degree << "..."
          << std::endl;
  const unsigned int n_cell_dofs =
    (degree + 1) * (degree + 1) * (dim == 3 ? (degree + 1) : 1);
  constexpr const unsigned int n_cells = 1 << dim;

  using DistInterior = PatchDistributors::Dynamic<dim, degree>;

  DistInterior distributor;

  // Print an empty line every (degree+1) printed lines.
  std::size_t        printed_lines   = 0;
  const unsigned int lines_per_block = degree + 1;
  auto               maybe_newline   = [&]() {
    ++printed_lines;
    if ((printed_lines % lines_per_block) == 0)
      deallog << std::endl;
    if (printed_lines % (lines_per_block * lines_per_block) == 0)
      deallog << "------------\n";
  };

  // Simple map dump: print cell first, then dof, then patch dof (index).
  const auto print_map =
    [&](unsigned patch_index, unsigned cell, unsigned dof) {
      deallog << "Regular   C " << cell << " D " << dof << " -> P "
              << patch_index << std::endl;
      maybe_newline();
    };

  const auto print_map_duplicate =
    [&](unsigned patch_index, unsigned cell, unsigned dof) {
      deallog << "Duplicate C " << cell << " D " << dof << " -> P "
              << patch_index << std::endl;
      maybe_newline();
    };

  const auto print_skipped = [&](unsigned cell, unsigned dof) {
    deallog << "Skipped   C " << cell << " D " << dof << std::endl;
    maybe_newline();
  };

  deallog
    << "Printing cell, dof, then patch (primary + duplicates), then skipped:\n"
    << std::endl;
  distributor.loop(print_map, print_map_duplicate, print_skipped);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  initlog();


  test<1, 1>();
  test<1, 2>();
  test<1, 3>();


  test<2, 1>();
  test<2, 2>();
  test<2, 3>();

  // test<2, 4>();
  // test<2, 5>();
  // test<2, 6>();
  // test<2, 7>();


  test<3, 1>();
  test<3, 2>();
  test<3, 3>();

  // test<3, 4>();
  // test<3, 5>();
  // test<3, 6>();
  // test<3, 7>();
  deallog << "Tests finished." << std::endl;

  return 0; // Indicate success
}
