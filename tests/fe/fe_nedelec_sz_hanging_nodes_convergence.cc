// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//
// This convergence test verifies that hanging nodes on FE_NedelecSZ
// elements are handled correctly. The orientation of edges and faces is
// automatically adapted in the presence of hanging edges and hanging
// faces. Furthermore, the function make_hanging_node_constraints()
// considers the orientation of the underlying edges and faces to prevent
// sign conflicts.
//
// This test solves the real-valued Maxwell equation in 2D and 3D:
//
// curl(curl(E)) - E = 0,
//
// where we consider Dirichlet boundary data
// such that we know the solution E.
// In 2D:
// E = (sin(y), sin(x))
//
// In 3D:
// E = (sin(y), sin(z), 0)
//
// In the first step, we solve the Maxwell equation on a coarse grid.
// After that, we refine the center of the domain so the grid contains
// hanging edges (and faces in 3D). We expect that the L2 difference
// between the numerical solution and the exact solution stays the same
// or gets smaller in each step, i.e., with a finer and finer grid, we
// converge more and more to the exact solution.
// The geometry is chosen in such a way that it covers the most important
// cases and is as small as possible.
// To cover all aspects of make_hanging_node_constraints(), we have to
// consider at least cubic base functions (that corresponds to
// polynomial_degree = 2).

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <cmath>
#include <iostream>

#include "../tests.h"

namespace ConvergenceTest
{

  template <int dim>
  class SolutionValues : public Function<dim>
  {
  public:
    SolutionValues()
      : Function<dim>(dim)
    {}

    void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>>   &value_list) const override
    {
      Assert(value_list.size() == points.size(),
             ExcDimensionMismatch(value_list.size(), points.size()));

      for (unsigned int p = 0; p < points.size(); p++)
        {
          SolutionValues<dim>::vector_value(points[p], value_list[p]);
        }
    }
  };

  // here we are using artificial boundary values,
  // the benefit is, that we know the exact solution for the
  // electric field E, i.e. this is also the exact solution
  // to the partial differential equation we aim to solve
  template <>
  void
  SolutionValues<2>::vector_value(const Point<2> &p,
                                  Vector<double> &values) const
  {
    values(0) = sin(p[1]);
    values(1) = sin(p[0]);
  }

  template <>
  void
  SolutionValues<3>::vector_value(const Point<3> &p,
                                  Vector<double> &values) const
  {
    values(0) = sin(p[1]);
    values(1) = sin(p[2]);
    values(2) = 0.0;
  }



  template <int dim>
  class MaxwellProblem
  {
  public:
    MaxwellProblem(unsigned int refinements,
                   unsigned int poly_degree,
                   unsigned int n_iterations);

    void
    run();

  private:
    void
    make_grid();
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid(double radius);
    double
    output_error();

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FE_NedelecSZ<dim>  fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution, system_rhs;

    const unsigned int refinements;
    const unsigned int poly_degree;
    const unsigned int n_iterations;
  };

  template <int dim>
  MaxwellProblem<dim>::MaxwellProblem(const unsigned int refinements,
                                      const unsigned int poly_degree,
                                      const unsigned int n_iterations)
    : dof_handler(triangulation)
    , fe(poly_degree)
    , refinements(refinements)
    , poly_degree(poly_degree)
    , n_iterations(n_iterations)
  {}

  template <>
  void
  MaxwellProblem<2>::make_grid()
  {
    const unsigned int dim = 2;

    const Point<dim> center(0.0, 0.0);
    const double     radius = 1.00;

    GridGenerator::hyper_ball_balanced(triangulation, center, radius);

    triangulation.refine_global(refinements);
  }

  template <>
  void
  MaxwellProblem<3>::make_grid()
  {
    const double radius      = 1.00;
    const double x_thickness = 0.75;

    GridGenerator::cylinder(triangulation, radius, x_thickness / 2.0);

    for (auto &cell : triangulation.active_cell_iterators())
      for (auto &face : cell->face_iterators())
        if (face->at_boundary())
          face->set_boundary_id(0);

    triangulation.refine_global(refinements);
  }

  template <int dim>
  void
  MaxwellProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0 /* vector component*/,
      SolutionValues<dim>(),
      0 /* boundary id*/,
      constraints);
    constraints.close();

    constraints.condense(dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  void
  MaxwellProblem<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int curl_dim = (dim == 2) ? 1 : 3;

    // choose the quadrature formulas
    QGauss<dim> quadrature_formula(fe.degree + 2);

    // get the number of quadrature points and dofs
    const unsigned int n_q_points    = quadrature_formula.size(),
                       dofs_per_cell = fe.dofs_per_cell;

    // set update flags
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    // Extractors for the real part
    const FEValuesExtractors::Vector E_re(0);

    // create the local left hand side and right hand side
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // loop over all cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        // initialize values:
        cell_matrix = 0;
        cell_rhs    = 0;
        fe_values.reinit(cell);

        for (const unsigned int i : fe_values.dof_indices())
          {
            // only compute this once
            std::vector<Tensor<1, dim>>      phi_i(n_q_points);
            std::vector<Tensor<1, curl_dim>> curl_phi_i(n_q_points);
            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              {
                phi_i[q_point]      = fe_values[E_re].value(i, q_point);
                curl_phi_i[q_point] = fe_values[E_re].curl(i, q_point);
              }

            // we use here, that the problem is symmetrical
            for (unsigned int j = i; j < dofs_per_cell; j++)
              {
                double mass_part = 0;
                double curl_part = 0;

                for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
                  {
                    Tensor<1, dim> phi_j = fe_values[E_re].value(j, q_point);
                    Tensor<1, curl_dim> curl_phi_j =
                      fe_values[E_re].curl(j, q_point);

                    curl_part +=
                      curl_phi_i[q_point] * curl_phi_j * fe_values.JxW(q_point);

                    mass_part +=
                      phi_i[q_point] * phi_j * fe_values.JxW(q_point);
                  }

                double mass_term  = curl_part - mass_part;
                cell_matrix(i, j) = mass_term;
                cell_matrix(j, i) = mass_term;
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }


  template <int dim>
  void
  MaxwellProblem<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
    constraints.distribute(solution);
  }

  template <int dim>
  void
  MaxwellProblem<dim>::refine_grid(const double radius)
  {
    for (auto &cell : triangulation.cell_iterators())
      {
        if (!cell->is_active())
          continue;

        double distance;
        if (dim == 3)
          distance = std::sqrt(std::pow(cell->center()[1], 2) +
                               std::pow(cell->center()[2], 2));
        else if (dim == 2)
          distance = std::sqrt(std::pow(cell->center()[0], 2) +
                               std::pow(cell->center()[1], 2));

        if (distance < radius)
          cell->set_refine_flag();
      }

    triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  double
  MaxwellProblem<dim>::output_error()
  {
    SolutionValues<dim> exact_solution;
    Vector<double>      diff_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      exact_solution,
                                      diff_per_cell,
                                      QGauss<dim>(poly_degree + 2),
                                      VectorTools::L2_norm);
    return diff_per_cell.l2_norm();
  }

  template <int dim>
  void
  MaxwellProblem<dim>::run()
  {
    deallog << "Testing for dim = " << dim
            << ", polynomial_degree p = " << poly_degree << std::endl;
    std::vector<double> L2_error(n_iterations);
    bool                passed = false;

    for (unsigned int cycle = 0; cycle < n_iterations; ++cycle)
      {
        if (cycle == 0)
          make_grid();
        else
          refine_grid(0.25);

        setup_system();
        assemble_system();
        solve();
        L2_error[cycle] = output_error();

        // check that we are close to the analytical solution
        passed = (std::fabs(L2_error[cycle]) < 1e-3);

        if (!passed)
          {
            deallog << "FAILED" << std::endl;
            deallog << "L2 Error = " << L2_error << "(threshold = 1e-3)"
                    << std::endl;
          }


        // check for convergence
        if (cycle > 0)
          {
            passed = (L2_error[cycle - 1] >= L2_error[cycle]);

            if (!passed)
              {
                deallog << "FAILED" << std::endl;
                deallog << "Convergence error in step " << cycle << std::endl;
              }
          }
      }

    if (passed)
      deallog << "OK" << std::endl;
  }
} // namespace ConvergenceTest

int
main()
{
  initlog();

  using namespace ConvergenceTest;

  const unsigned int refinements  = 1;
  const unsigned int poly_degree  = 2;
  const unsigned int n_iterations = (running_in_debug_mode ? 2 : 3);

  MaxwellProblem<2> maxwell_2d(refinements, poly_degree, n_iterations);
  maxwell_2d.run();

  MaxwellProblem<3> maxwell_3d(refinements, poly_degree, n_iterations);
  maxwell_3d.run();

  return 0;
}
