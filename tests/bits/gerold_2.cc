// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// apply SparsityTools::reorder_Cuthill_McKee to the cell connectivity
// graph for the mesh used in gerold_2. apparently the mesh consists
// of two or more non-connected parts, and the reordering algorithm
// trips over that

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"


template <int dim>
class LaplaceProblem
{
public:
  void
  run();

private:
  Triangulation<dim> triangulation;
};


template <int dim>
void
LaplaceProblem<dim>::run()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);

  std::ifstream input_file(SOURCE_DIR "/gerold_1.inp");
  grid_in.read_ucd(input_file);

  DynamicSparsityPattern cell_connectivity;
  GridTools::get_face_connectivity_of_cells(triangulation, cell_connectivity);
  std::vector<types::global_dof_index> permutation(
    triangulation.n_active_cells());
  SparsityTools::reorder_Cuthill_McKee(cell_connectivity, permutation);

  for (unsigned int i = 0; i < permutation.size(); ++i)
    deallog << permutation[i] << std::endl;
}


int
main()
{
  initlog();

  try
    {
      LaplaceProblem<3> laplace_problem_3d;
      laplace_problem_3d.run();
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

  return 0;
}
