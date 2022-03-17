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
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
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
 * are correctly scaled according to element size, then solving the Laplace
 * equation for a polynomial solution of this form in the Hermite FE space
 * will produce negligible pointwise errors.
 */

using namespace dealii;

// Define a function to project onto the domain [-1,1]
class Solution : public Function<1>
{
public:
  virtual double
  value(const Point<1> &p, unsigned int c = 0) const override
  {
    return p(c) * (1.0 + p(c) * (0.5 - p(c)));
  }



  std::string
  get_function_string()
  {
    return "X + 0.5 X^2 - X^3";
  }
};



class RHSFunction : public Function<1>
{
public:
  virtual double
  value(const Point<1> &p, unsigned int c = 0) const override
  {
    return -1.0 + 6 * p(c);
  }
};



void
test_fe_on_domain(const unsigned int regularity)
{
  deallog << std::endl;
  char fname[50];
  sprintf(fname, "Cell-1d-Hermite-%d", regularity);
  deallog.push(fname);
  deallog.depth_file(2);

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
            }
        }
      tr.execute_coarsening_and_refinement();
    }

  FE_Hermite<1> herm(2 * regularity + 1);
  dof.distribute_dofs(herm);

  MappingCartesian<1> mapping;

  QGauss<1> quadr(2 * regularity + 2);

  Vector<double> sol(dof.n_dofs());
  Vector<double> rhs(dof.n_dofs());

  Solution    sol_object;
  RHSFunction rhs_object;

  AffineConstraints<double> constraints;
  constraints.close();

  DynamicSparsityPattern dsp(dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, dsp);
  SparsityPattern sp;
  sp.copy_from(dsp);

  SparseMatrix<double> stiffness_matrix;
  stiffness_matrix.reinit(sp);

  MatrixCreator::create_laplace_matrix(mapping, dof, quadr, stiffness_matrix);

  FEValues<1> fe_herm(mapping,
                      herm,
                      quadr,
                      update_values | update_gradients |
                        update_quadrature_points | update_JxW_values);

  std::vector<types::global_dof_index> local_to_global(herm.n_dofs_per_cell());

  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_herm.reinit(cell);
      cell->get_dof_indices(local_to_global);
      for (const unsigned int i : fe_herm.dof_indices())
        {
          double rhs_temp = 0;
          for (const unsigned int q : fe_herm.quadrature_point_indices())
            rhs_temp += fe_herm.shape_value(i, q) *
                        rhs_object.value(fe_herm.quadrature_point(q)) *
                        fe_herm.JxW(q);
          rhs(local_to_global[i]) += rhs_temp;
        }
    }

  std::map<types::global_dof_index, double> bound_vals;

  std::map<types::boundary_id, const Function<1, double> *> bound_map;
  bound_map.emplace(std::make_pair(0U, &sol_object));
  bound_map.emplace(std::make_pair(1U, &sol_object));

  VectorTools::project_hermite_boundary_values(
    mapping, dof, bound_map, QGauss<0>(1), 0, bound_vals);

  MatrixTools::apply_boundary_values(bound_vals, stiffness_matrix, sol, rhs);

  IterationNumberControl solver_control_values(50, 1e-11);
  SolverCG<>             solver(solver_control_values);

  solver.solve(stiffness_matrix, sol, rhs, PreconditionIdentity());

  double err_sq = 0;

  for (auto &cell : dof.active_cell_iterators())
    {
      fe_herm.reinit(cell);
      cell->get_dof_indices(local_to_global);
      for (const unsigned int q : fe_herm.quadrature_point_indices())
        {
          double sol_at_point = 0;
          for (const unsigned int i : fe_herm.dof_indices())
            sol_at_point += fe_herm.shape_value(i, q) * sol(local_to_global[i]);

          sol_at_point -= sol_object.value(fe_herm.quadrature_point(q));
          err_sq += sol_at_point * sol_at_point * fe_herm.JxW(q);
        }
    }

  err_sq = std::sqrt(err_sq);

  deallog << "Test polynomial:" << std::endl;
  deallog << sol_object.get_function_string() << std::endl;
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
