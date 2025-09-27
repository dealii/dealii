// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-08 on a simplex mesh. Following incompatible modifications had to be
// made:
//  - Change the FE_Q to FE_SimplexP.
//  - Change QGauss to QGaussSimplex.
//  - Use MappingFE (Do not use default mapping).
//  - Convert triangulation to a triangulation based on simplices.
//  - Use refine_global() instead of execute_coarsening_and_refinement().
//  - Reduce number of cycles from 7 to 4 to reduce runtime.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
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
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

namespace Step8
{

  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid();
    void
    output_results(const unsigned int cycle) const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    MappingFE<dim>     mapping;

    FESystem<dim> fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };

  template <int dim>
  void
  right_hand_side(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>>   &values)
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

    Point<dim> point_1, point_2;
    point_1[0] = 0.5;
    point_2[0] = -0.5;

    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
        if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
            ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
          values[point_n][0] = 1.0;
        else
          values[point_n][0] = 0.0;

        if (points[point_n].norm_square() < 0.2 * 0.2)
          values[point_n][1] = 1.0;
        else
          values[point_n][1] = 0.0;
      }
  }


  template <int dim>
  ElasticProblem<dim>::ElasticProblem()
    : dof_handler(triangulation)
    , mapping(FE_SimplexP<dim>(1))
    , fe(FE_SimplexP<dim>(1), dim)
  {}

  template <int dim>
  void
  ElasticProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Functions::ZeroFunction<dim>(dim), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }

  template <int dim>
  void
  ElasticProblem<dim>::assemble_system()
  {
    QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);

    Functions::ConstantFunction<dim> lambda(1.), mu(1.);

    std::vector<Tensor<1, dim>> rhs_values(n_q_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j =
                  fe.system_to_component_index(j).first;

                for (const unsigned int q_point :
                     fe_values.quadrature_point_indices())
                  {
                    cell_matrix(i, j) +=
                      ((fe_values.shape_grad(i, q_point)[component_i] *
                        fe_values.shape_grad(j, q_point)[component_j] *
                        lambda_values[q_point]) +
                       (fe_values.shape_grad(i, q_point)[component_j] *
                        fe_values.shape_grad(j, q_point)[component_i] *
                        mu_values[q_point]) +
                       ((component_i == component_j) ?        //
                          (fe_values.shape_grad(i, q_point) * //
                           fe_values.shape_grad(j, q_point) * //
                           mu_values[q_point]) :              //
                          0)                                  //
                       ) *                                    //
                      fe_values.JxW(q_point);                 //
                  }
              }
          }

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int q_point :
                 fe_values.quadrature_point_indices())
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             rhs_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }

  template <int dim>
  void
  ElasticProblem<dim>::solve()
  {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }

  template <int dim>
  void
  ElasticProblem<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    QGaussSimplex<dim - 1> quadrature(fe.degree + 1);
    KellyErrorEstimator<dim>::estimate(
      mapping, dof_handler, quadrature, {}, solution, estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.refine_global();
  }

  template <int dim>
  void
  ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names;
    switch (dim)
      {
        case 1:
          solution_names.emplace_back("displacement");
          break;
        case 2:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          break;
        case 3:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          solution_names.emplace_back("z_displacement");
          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }

    data_out.add_data_vector(solution, solution_names);
    data_out.build_patches(mapping);

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
    data_out.write_vtk(output);
  }

  template <int dim>
  void
  ElasticProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 5; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            Triangulation<dim> temp;
            GridGenerator::hyper_cube(temp, -1, 1);
            GridGenerator::convert_hypercube_to_simplex_mesh(temp,
                                                             triangulation);
          }
        else
          refine_grid();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

        assemble_system();
        solve();
        output_results(cycle);
      }
  }
} // namespace Step8

int
main()
{
  try
    {
      Step8::ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
