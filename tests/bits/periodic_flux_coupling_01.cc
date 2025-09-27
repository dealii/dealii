// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/* checks that make_flux_sparsity_pattern() with masks runs without
 * calling inactive cells
 * on an adapted mesh with periodic boundary conditions
 */

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class MakeFlux
{
public:
  MakeFlux();
  void
  run();

private:
  void
  make_grid();

  Triangulation<dim> triangulation;
  FE_DGQ<dim>        fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern sparsity_pattern;
};

template <int dim>
MakeFlux<dim>::MakeFlux()
  : fe(1)
  , dof_handler(triangulation)
{}


template <int dim>
void
MakeFlux<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  using CellIteratorTria = typename dealii::Triangulation<dim>::cell_iterator;
  std::vector<dealii::GridTools::PeriodicFacePair<CellIteratorTria>>
                     periodic_faces;
  const unsigned int b_id1     = 2;
  const unsigned int b_id2     = 3;
  const unsigned int direction = 1;

  dealii::GridTools::collect_periodic_faces(triangulation,
                                            b_id1,
                                            b_id2,
                                            direction,
                                            periodic_faces,
                                            dealii::Tensor<1, dim>());
  triangulation.add_periodicity(periodic_faces);
  triangulation.refine_global(1);
}


template <int dim>
void
MakeFlux<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 3; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        {
          /* refine the top-left cell to ensure we test this cell being looped
           * over in make_flux_sparsity_pattern after the bottom-left cell and
           * before the top-right cell, both of which will be less refined.
           */

          Point<dim> refn_point;
          refn_point = Point<dim>(0.005, 0.995);
          typename Triangulation<dim>::active_cell_iterator cell_it =
            triangulation.begin_active();
          for (; cell_it != triangulation.end(); ++cell_it)
            {
              if (cell_it->is_locally_owned() &&
                  cell_it->point_inside(refn_point))
                {
                  cell_it->set_refine_flag();
                  break;
                }
            }
          triangulation.execute_coarsening_and_refinement();
        }

      dof_handler.distribute_dofs(fe);

      DynamicSparsityPattern dsp(dof_handler.n_dofs());

      // set up full mask not doing anything
      const unsigned int n_components = dof_handler.get_fe().n_components();
      Table<2, DoFTools::Coupling> mask(n_components, n_components);
      for (unsigned int i = 0; i < n_components; ++i)
        for (unsigned int j = 0; j < n_components; ++j)
          mask(i, j) = DoFTools::always;
      DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, mask, mask);
      deallog << Utilities::int_to_string(dof_handler.n_dofs(), 2) << std::endl;
      deallog.pop();
    }
  deallog << "PASSED" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  try
    {
      MakeFlux<2> test;
      test.run();
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
