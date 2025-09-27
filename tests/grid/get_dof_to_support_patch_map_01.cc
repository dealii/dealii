

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test GridTools::map_dof_to_support_patch () using a quadratic FE_Q<dim>
// space.  We output the barycenter of each cell in the patch around the dof
// for a few dofs in the triangulation.  We have three layers of refinement in
// the triangulation which represents the various patches that could arise in
// practice.



#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation, 0, 1);

  triangulation.refine_global(2);
  {
    // choose barycenters of cells to refine or coarsen to end up with
    // a three level triangulation.
    std::vector<Point<dim>> refine_centers;
    std::vector<Point<dim>> coarsen_centers;

    if (dim == 1)
      {
        refine_centers.push_back(Point<dim>(1. / 8.));

        coarsen_centers.push_back(Point<dim>(5. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8.));
      }
    else if (dim == 2)
      {
        refine_centers.push_back(Point<dim>(1. / 8., 7. / 8.));

        coarsen_centers.push_back(Point<dim>(5. / 8., 5. / 8.));
        coarsen_centers.push_back(Point<dim>(5. / 8., 7. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8., 5. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8., 7. / 8.));
      }
    else if (dim == 3)
      {
        refine_centers.push_back(Point<dim>(1. / 8., 7. / 8., 1. / 8.));

        coarsen_centers.push_back(Point<dim>(7. / 8., 7. / 8., 7. / 8.));
        coarsen_centers.push_back(Point<dim>(5. / 8., 7. / 8., 7. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8., 5. / 8., 7. / 8.));
        coarsen_centers.push_back(Point<dim>(5. / 8., 5. / 8., 7. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8., 7. / 8., 5. / 8.));
        coarsen_centers.push_back(Point<dim>(5. / 8., 7. / 8., 5. / 8.));
        coarsen_centers.push_back(Point<dim>(7. / 8., 5. / 8., 5. / 8.));
        coarsen_centers.push_back(Point<dim>(5. / 8., 5. / 8., 5. / 8.));
      }
    else
      DEAL_II_NOT_IMPLEMENTED();

    unsigned int index = 0;
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell, ++index)

      {
        Point<dim> cell_bary = cell->barycenter();

        // refine cells
        for (unsigned int i = 0; i < refine_centers.size(); ++i)
          {
            double diff = 0;
            for (unsigned int d = 0; d < dim; ++d)
              diff += pow(cell_bary[d] - refine_centers[i][d], 2.);
            diff = std::sqrt(diff);

            if (diff < 1e-14)
              {
                cell->set_refine_flag();
                break;
              }
          }
        // coarsen cells
        for (unsigned int i = 0; i < coarsen_centers.size(); ++i)
          {
            double diff = 0;
            for (unsigned int d = 0; d < dim; ++d)
              diff += pow(cell_bary[d] - coarsen_centers[i][d], 2.);
            diff = std::sqrt(diff);

            if (diff < 1e-14)
              {
                cell->set_coarsen_flag();
                break;
              }
          }
      }
    triangulation.execute_coarsening_and_refinement();
  }


  DoFHandler<dim> dof_handler(triangulation);
  unsigned int    iFEDeg = 2;
  FE_Q<dim>       finite_element(iFEDeg);
  dof_handler.distribute_dofs(finite_element);

  std::map<types::global_dof_index,
           std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator>>
    dof_to_cell_map = GridTools::get_dof_to_support_patch_map(dof_handler);

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      // loop through and print out barycenter of cells in patch of certain dofs
      // in system
      if (i % (dim == 1 ? 1 : (dim == 2 ? 5 : 25)) == 0)
        {
          deallog << "Patch around dof " << i << ": ";
          typename std::vector<
            typename dealii::DoFHandler<dim>::active_cell_iterator>::iterator
            patch_iter     = dof_to_cell_map[i].begin(),
            patch_iter_end = dof_to_cell_map[i].end();
          for (; patch_iter != patch_iter_end; ++patch_iter)
            {
              typename dealii::DoFHandler<dim>::active_cell_iterator
                         patch_cell = *(patch_iter);
              Point<dim> cell_bary  = patch_cell->barycenter();
              deallog << '(';
              for (unsigned int d = 0; d < dim - 1; ++d)
                deallog << cell_bary[d] << ", ";
              deallog << cell_bary[dim - 1] << ") ";
            }
          deallog << std::endl;
        }
    }


  // clean up data
  dof_handler.clear();
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
