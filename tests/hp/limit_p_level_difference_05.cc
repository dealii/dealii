// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check results of limit_p_level_difference on cells that
// share edges but no faces in 3D.
//
// 3D: n_refinements: 1           n_refinements: 2
//     mesh:       future FE:     mesh:       future FE:
//     +---+-+-+   +---+++++      +-+-+++++   +-+-+++++
//     |   | | |   |   |2|2|      | | +++++   |0|1+222+
//     |   +-+-+   | 1 +-+-+      +-+-+++++   +-+-+222+
//     |   | | |   |   |2|2|      | | +++++   |0|1+222+
//     +---X-+-+   +---X-+-+      +-+-X++++   +-+-X++++
//     |   |       |   |          | | |       |0|0|
//     |   |       | 0 |          +-+-+       |-+-|
//     |   |       |   |          | | |       |0|0|
//     +---+       +---+          +---+       +---+

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

#include "../tests.h"


template <int dim>
void
write(const DoFHandler<dim> &dofh, const unsigned int n_refinements)
{
  const auto &tria = dofh.get_triangulation();

  Vector<float> active_fe_indices(tria.n_active_cells());
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      active_fe_indices(cell->active_cell_index()) = cell->active_fe_index();

  Vector<float> future_fe_indices(tria.n_active_cells());
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      future_fe_indices(cell->active_cell_index()) = cell->future_fe_index();

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dofh);
  data_out.add_data_vector(active_fe_indices, "active_fe_indices");
  data_out.add_data_vector(future_fe_indices, "future_fe_indices");
  data_out.build_patches();

  std::ofstream ofile("fe_indices-" + Utilities::to_string(dim) + "D-nref" +
                      Utilities::to_string(n_refinements) + ".vtk");
  data_out.write_vtk(ofile);
}


template <int dim>
void
test(const unsigned int n_refinements)
{
  hp::FECollection<dim> fes;
  for (unsigned int p = 1; p <= 3; ++p)
    fes.push_back(FE_Q<dim>(p));

  Triangulation<dim> tria;
  DoFHandler<dim>    dofh(tria);
  {
    // create L-shaped mesh
    std::vector<unsigned int> repetitions(dim);
    Point<dim>                bottom_left, top_right;
    for (unsigned int d = 0; d < dim; ++d)
      if (d < 2)
        {
          repetitions[d] = 2;
          bottom_left[d] = -1.;
          top_right[d]   = 1.;
        }
      else
        {
          repetitions[d] = 1;
          bottom_left[d] = 0.;
          top_right[d]   = 1.;
        }

    std::vector<int> cells_to_remove(dim, 1);
    cells_to_remove[0] = -1;

    GridGenerator::subdivided_hyper_L(
      tria, repetitions, bottom_left, top_right, cells_to_remove);

    // refine grid according to sketch above
    for (unsigned int i = 0; i < n_refinements; ++i)
      {
        for (const auto &cell : dofh.active_cell_iterators())
          if (cell->center()[0] > 0)
            {
              cell->set_refine_flag();
              cell->set_active_fe_index(2);
            }

        tria.execute_coarsening_and_refinement();
      }
  }

  dofh.distribute_dofs(fes);
  hp::Refinement::limit_p_level_difference(dofh,
                                           /*max_difference=*/1,
                                           /*contains=*/0);

  // output FE indices
  for (const auto &cell : dofh.active_cell_iterators())
    deallog << cell->id().to_string() << ": active:" << cell->active_fe_index()
            << " future:" << cell->future_fe_index() << std::endl;

#if false
  write(dofh, n_refinements);
#endif

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>(1);
  test<2>(2);
  deallog.pop();
  deallog.push("3d");
  test<3>(1);
  test<3>(2);
  deallog.pop();
}
