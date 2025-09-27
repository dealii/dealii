// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
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

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>

#include "../tests.h"

#define PRECISION 8

/**
 * Test case for Hermite interpolation elements on an regular square grid,
 * to check the assembly of the mass matrix works as intended and
 * representation of polynomials up to $1D$ degree @p poly_degree by
 * <code>FE_Hermite<dim>(poly_degree)<\code> is exact.
 */

using namespace dealii;

// Define a function to project onto the domain [-1,1]^2
class Solution : public Function<2>
{
public:
  virtual double
  value(const Point<2> &p, unsigned int c = 0) const override
  {
    (void)c;
    return p[0] * p[1];
  }


  std::string
  get_function_string()
  {
    return "XY";
  }
};


void
test_fe_on_domain(const unsigned int regularity)
{
  deallog << std::endl;
  char fname[50];
  sprintf(fname, "Cell-%dd-Hermite-%d", 2, 2 * regularity + 1);
  deallog.push(fname);
  deallog.depth_file(2);

  Triangulation<2> tr;
  DoFHandler<2>    dof(tr);

  double   left = -1.0, right = 1.0;
  Point<2> left_point, right_point;
  for (unsigned int i = 0; i < 2; ++i)
    left_point[i] = left, right_point[i] = right;
  GridGenerator::subdivided_hyper_cube(tr, 4, left, right);

  FE_Hermite<2> herm(2 * regularity + 1);
  dof.distribute_dofs(herm);

  MappingCartesian<2> mapping;

  QGauss<2> quadr(2 * regularity + 2);

  Vector<double> sol(dof.n_dofs());
  Vector<double> rhs(dof.n_dofs());

  Solution sol_object;

  AffineConstraints<double> constraints;
  constraints.close();

  DynamicSparsityPattern dsp(dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, dsp);
  SparsityPattern sp;
  sp.copy_from(dsp);

  SparseMatrix<double> mass_matrix;
  mass_matrix.reinit(sp);
  MatrixCreator::create_mass_matrix(mapping, dof, quadr, mass_matrix);

  FEValues<2> fe_herm(mapping,
                      herm,
                      quadr,
                      update_values | update_quadrature_points |
                        update_JxW_values);

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
                        sol_object.value(fe_herm.quadrature_point(q)) *
                        fe_herm.JxW(q);
          rhs(local_to_global[i]) += rhs_temp;
        }
    }

  IterationNumberControl solver_control_values(8000, 1e-11);
  SolverCG<>             solver(solver_control_values);

  solver.solve(mass_matrix, sol, rhs, PreconditionIdentity());

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
  MultithreadInfo::set_thread_limit(1);

  test_fe_on_domain(0);
  test_fe_on_domain(1);
  test_fe_on_domain(2);

  return 0;
}
