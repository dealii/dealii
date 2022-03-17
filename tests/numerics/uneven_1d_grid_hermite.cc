// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

#include <deal.II/base/config.h>

#define PRECISION 8


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>



/**
 * Test case for Hermite on an irregular 1D grid.
 * <code>FE_Hermite<1>(poly_degree)<\code> should be able to perfectly represent
 * any
 * polynomial function up to degree @p poly_degree. If all basis functions
 * are correctly scaled according to element size, then projecting a
 * polynomial of this form onto the Hermite FE space will produce negligible
 * pointwise errors.
 */
using namespace dealii;

// Define a function to project onto the domain
class TestPoly : public Function<1>
{
public:
  virtual double
  value(const Point<1> &p, unsigned int c = 0) const override
  {
    return p(c) * (1.0 + p(c) * (0.5 - p(c)));
  }



  std::string
  get_polynomial_string()
  {
    return "X + 0.5 X^2 - X^3";
  }
};



void
test_fe_on_domain(const unsigned int regularity)
{
  Triangulation<1> tr;
  DoFHandler<1>    dof(tr);

  double   left = -1.0, right = 1.0;
  Point<1> left_point(left), right_point(right);
  GridGenerator::hyper_cube(tr, left, right);

  // Refine the right-most cell three times to get the elements
  // [-1,0],[0,0.5],[0.5,0.75],[0.75,1]
  for (unsigned int i = 0; i < 3; ++i)
    {
      for (auto &cell : tr.active_cell_iterators())
        {
          const double distance = right_point.distance(cell->vertex(1));
          if (distance < 1e-6)
            {
              cell->set_refine_flag();
              // break;
            }
        }
      tr.execute_coarsening_and_refinement();
    }

  FE_Hermite<1> herm(2 * regularity + 1);
  dof.distribute_dofs(herm);

  MappingCartesian<1> mapping;

  QGauss<1>      quadr(2 * regularity + 2);
  Vector<double> solution(dof.n_dofs());
  TestPoly       rhs_func;

  AffineConstraints<double> constraints;
  constraints.close();

  FEValues<1> fe_herm(mapping,
                      herm,
                      quadr,
                      update_values | update_quadrature_points |
                        update_JxW_values);

  std::vector<types::global_dof_index> local_to_global(herm.n_dofs_per_cell());

  VectorTools::project_hermite(
    mapping, dof, constraints, quadr, rhs_func, solution, false);

  double err_sq = 0;

  for (auto &cell : dof.active_cell_iterators())
    {
      fe_herm.reinit(cell);
      cell->get_dof_indices(local_to_global);
      for (const unsigned int q : fe_herm.quadrature_point_indices())
        {
          double sol_at_point = 0;
          for (const unsigned int i : fe_herm.dof_indices())
            sol_at_point +=
              fe_herm.shape_value(i, q) * solution(local_to_global[i]);

          sol_at_point -= rhs_func.value(fe_herm.quadrature_point(q));
          err_sq += sol_at_point * sol_at_point * fe_herm.JxW(q);
        }
    }

  err_sq = std::sqrt(err_sq);

  char fname[50];
  sprintf(fname, "Cell-1d-Hermite-%d", regularity);
  deallog.push(fname);

  deallog << "Test polynomial:" << std::endl;
  deallog << rhs_func.get_polynomial_string() << std::endl;
  deallog << std::endl;

  deallog << "Grid cells:" << std::endl;
  for (const auto &cell : tr.active_cell_iterators())
    {
      deallog << "(\t" << cell->vertex(0) << ","
              << "\t" << cell->vertex(1) << "\t)" << std::endl;
    }
  deallog << std::endl;

  deallog << "Interpolation error:" << std::endl;
  deallog << err_sq << "\n\n" << std::endl;
  deallog.pop();
}



int
main()
{
  std::ofstream logfile("output");

  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  test_fe_on_domain(0);
  test_fe_on_domain(1);
  test_fe_on_domain(2);
  test_fe_on_domain(3);

  return 0;
}
