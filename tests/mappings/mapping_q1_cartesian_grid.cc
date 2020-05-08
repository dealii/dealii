// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <iostream>

#include "../tests.h"

// Test that MappingQ1 works on Cartesian grids. This is a regression test for
// a bug where, due to the quadratic equation solver in
// transform_unit_to_real_cell, we lost nearly all digits of precision due to
// roundoff.


class RefineTest
{
public:
  static constexpr int dim = 2;

public:
  RefineTest(const unsigned int refine_case)
    : triangulation()
    , dof_handler(triangulation)
    , refine_case(refine_case)
  {
    make_grid();
    test_mapping();
  }

private:
  void
  make_grid()
  {
    if (refine_case == 1)
      {
        std::vector<unsigned int> repetitions(dim);
        repetitions[0] = repetitions[1] =
          5; // Non-unity number of initial partitions
        const std::array<double, dim> L = {{1.0, 1.0}};
        GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                  repetitions,
                                                  Point<dim>(0.0, 0.0),
                                                  Point<dim>(L[0], L[1]),
                                                  true);
        triangulation.refine_global(1); // With refinement
      }
    else if (refine_case == 2)
      {
        std::vector<unsigned int> repetitions(dim);
        repetitions[0] = repetitions[1] =
          5; // Non-unity number of initial partitions
        const std::array<double, dim> L = {{1.0, 1.0}};
        GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                  repetitions,
                                                  Point<dim>(0.0, 0.0),
                                                  Point<dim>(L[0], L[1]),
                                                  true);
        triangulation.refine_global(0); // With no refinement
      }
    else if (refine_case == 3)
      {
        std::vector<unsigned int> repetitions(dim);
        repetitions[0] = repetitions[1] = 10; // Finer initial mesh
        const std::array<double, dim> L = {{1.0, 1.0}};
        GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                  repetitions,
                                                  Point<dim>(0.0, 0.0),
                                                  Point<dim>(L[0], L[1]),
                                                  true);
        triangulation.refine_global(0); // With no refinement
      }
    else if (refine_case == 4)
      {
        std::vector<unsigned int> repetitions(dim);
        repetitions[0] = repetitions[1] = 1; // Single cell
        const std::array<double, dim> L = {{1.0, 1.0}};
        GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                  repetitions,
                                                  Point<dim>(0.0, 0.0),
                                                  Point<dim>(L[0], L[1]),
                                                  true);
        triangulation.refine_global(4); // With pure global refinement
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }
  }

  void
  test_mapping()
  {
    const double               tol = 1e-8;
    const MappingQGeneric<dim> mapping(1);

    deallog << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;

    auto       cell = dof_handler.begin_active();
    const auto endc = dof_handler.end();
    for (unsigned int index = 0; cell != endc; ++cell, ++index)
      {
        const auto P1 = cell->vertex(0);
        const auto P2 = cell->vertex(1);
        const auto P3 = cell->vertex(2);

        // Find the mid-point of the cell
        Point<dim> test;
        test[0] = (P1[0] + P2[0]) / 2.0;
        test[1] = (P2[1] + P3[1]) / 2.0;

        // Map to and from the unit cell
        const Point<dim> dp = mapping.transform_real_to_unit_cell(cell, test);
        const Point<dim> dp_real =
          mapping.transform_unit_to_real_cell(cell, dp);

        // Now we start off at the mid-point of the unit cell
        Point<dim> test_unit_mid;
        test_unit_mid[0] = 0.5;
        test_unit_mid[1] = 0.5;

        // And map to and from the real cell
        const Point<dim> dp_real_mid =
          mapping.transform_unit_to_real_cell(cell, test_unit_mid);
        const Point<dim> dp_test_mid =
          mapping.transform_real_to_unit_cell(cell, dp_real_mid);

        if ((test - dp_real).norm() > tol)
          {
            deallog << " " << std::endl;
            deallog << "ERROR" << std::endl;
            deallog << "cell = " << index << std::endl;
            deallog << "cell vertex(0): " << cell->vertex(0) << std::endl;
            deallog << "cell vertex(1): " << cell->vertex(1) << std::endl;
            deallog << "cell vertex(2): " << cell->vertex(2) << std::endl;
            deallog << "cell vertex(3): " << cell->vertex(3) << std::endl;
            deallog << " " << std::endl;
            deallog << "test point =   " << test << "  mapped point =   " << dp
                    << "  back-mapping =   " << dp_real << std::endl;
            deallog << "test point (unit mid) =   " << test_unit_mid
                    << "  mapped point =   " << dp_real_mid
                    << "  forward-mapping =   " << dp_test_mid << std::endl;
            deallog
              << "======================================================================================="
              << std::endl;
          }
      }
  }

private:
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler;
  const unsigned int refine_case;
};

int main(/*int argc, char *argv[]*/)
{
  initlog();

  {
    deallog << "Refine case = 1" << std::endl;
    RefineTest refineglobal(1);
  }
  {
    deallog << "Refine case = 2" << std::endl;
    RefineTest refineglobal(2);
  }
  {
    deallog << "Refine case = 3" << std::endl;
    RefineTest refineglobal(3);
  }
  {
    deallog << "Refine case = 4" << std::endl;
    RefineTest refineglobal(4);
  }
  return 0;
}
