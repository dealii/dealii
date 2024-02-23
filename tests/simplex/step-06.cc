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

// Step-06 on a simplex mesh.
//
// The following changes had to be made due to incompatibilities:
//  - Manifolds are disabled
//    All manifold ids from the triangulation object generated with
//    GridGenerator::hyper_ball() had to be removed
//  - Error estimation has to be skipped
//    KellyErrorEstimator requires QProjector::project_to_all_subfaces()
//    to work with triangles
//  - GridOut::write_gnuplot() has been skipped


#define USE_SIMPLEX
// #define OUTPUT_RESULTS


#ifdef USE_SIMPLEX
#  include <deal.II/base/quadrature_lib.h>
#  include <deal.II/base/types.h>

#  include <deal.II/fe/fe_pyramid_p.h>
#  include <deal.II/fe/fe_simplex_p.h>
#  include <deal.II/fe/fe_simplex_p_bubbles.h>
#  include <deal.II/fe/fe_wedge_p.h>
#else
#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/fe/fe_q.h>
#endif


#ifdef OUTPUT_RESULTS
#  include <deal.II/numerics/data_out.h>

#  include <fstream>
#endif


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class Step6
{
public:
  Step6();

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

#ifdef USE_SIMPLEX
  FE_SimplexP<dim> fe;
#else
  FE_Q<dim>         fe;
#endif

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler;

  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparsityPattern      sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;
};



template <int dim>
double
coefficient(const Point<dim> &p)
{
  if (p.square() < 0.5 * 0.5)
    return 20;
  else
    return 1;
}



template <int dim>
Step6<dim>::Step6()
  : fe(2)
  , dof_handler(triangulation)
{}



template <int dim>
void
Step6<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
}



template <int dim>
void
Step6<dim>::assemble_system()
{
#ifdef USE_SIMPLEX
  const QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
#else
  const QGauss<dim> quadrature_formula(fe.degree + 1);
#endif

  MappingFE<dim> mapping(fe);
  FEValues<dim>  fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const double current_coefficient =
            coefficient(fe_values.quadrature_point(q_index));
          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (current_coefficient *              // a(x_q)
                   fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                   fe_values.JxW(q_index));           // dx
              cell_rhs(i) += (1.0 *                   // f(x)
                              fe_values.shape_value(i, q_index) * // phi_i(x_q)
                              fe_values.JxW(q_index));            // dx
            }
        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}



template <int dim>
void
Step6<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  constraints.distribute(solution);
}



template <int dim>
void
Step6<dim>::refine_grid()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

#if false
#  ifdef USE_SIMPLEX
  KellyErrorEstimator<dim>::estimate(MappingFE<dim>(fe),
                                     dof_handler,
                                     QGaussSimplex<dim - 1>(fe.degree + 1),
                                     {},
                                     solution,
                                     estimated_error_per_cell);
#  else
  KellyErrorEstimator<dim>::estimate(MappingFE<dim>(fe),
                                     dof_handler,
                                     QGauss<dim - 1>(fe.degree + 1),
                                     {},
                                     solution,
                                     estimated_error_per_cell);
#  endif

  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.3,
                                                  0.03);
#else
  // incompatibility: KellyErrorEstimator does not yet work for simplices
  // instead, we will mark cells for refinement by hand
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int v = 0; v < cell->n_vertices(); ++v)
      {
        const double vertex_square = cell->vertex(v).square();
        if (vertex_square > 0.48 * 0.48 && vertex_square < 0.52 * 0.52)
          cell->set_refine_flag();
      }
#endif

  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
Step6<dim>::output_results(const unsigned int cycle) const
{
#ifdef OUTPUT_RESULTS
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches(MappingFE<dim>(fe));

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
    data_out.write_vtu(output);
  }
#else
  (void)cycle;
#endif
}



template <int dim>
void
Step6<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 5; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
#ifdef USE_SIMPLEX
          Triangulation<dim> tria_quad;
          GridGenerator::hyper_ball(tria_quad);

          // incompatibility: remove all manifold ids
          for (const auto &cell : tria_quad.active_cell_iterators())
            {
              cell->set_manifold_id(numbers::flat_manifold_id);
              for (const auto &face : cell->face_iterators())
                face->set_manifold_id(numbers::flat_manifold_id);
            }

          GridGenerator::convert_hypercube_to_simplex_mesh(tria_quad,
                                                           triangulation);
#else
          GridGenerator::hyper_ball(triangulation);
#endif

          triangulation.refine_global(1);
        }
      else
        refine_grid();

      deallog << "   Number of active cells:       "
              << triangulation.n_active_cells() << std::endl;

      setup_system();

      deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

      assemble_system();
      solve();
      output_results(cycle);
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      Step6<2> laplace_problem_2d;
      laplace_problem_2d.run();
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
