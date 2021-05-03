/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004
 *         Wolfgang Bangerth, Texas A&M University, 2016
 */

// Step 17 on simplex mesh.
// Modifications with respect to original step-17:
// - grid generation via GridGenerator::subdivided_hyper_rectangle and
//   GridGenerator::subdivided_hyper_rectangle_with_simplices (instead of
//   GridGenerator::hyper_cube);
// - removed local, adaptive refinement.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// simplex
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

// Flag to toggle between hexes and simplices.
// #define HEX

namespace Step17
{
  using namespace dealii;

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

    unsigned int
    solve();

    void
    output_results() const;

    MPI_Comm mpi_communicator;

    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

    ConditionalOStream pcout;

#ifdef HEX
    MappingQGeneric<dim, dim> mapping;
#else
    MappingFE<dim, dim> mapping;
#endif

    Triangulation<dim> triangulation;
    FESystem<dim>      fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> hanging_node_constraints;

    PETScWrappers::MPI::SparseMatrix system_matrix;

    PETScWrappers::MPI::Vector solution;
    PETScWrappers::MPI::Vector system_rhs;
  };


  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
      Assert(dim >= 2, ExcInternalError());

      Point<dim> point_1, point_2;
      point_1(0) = 0.5;
      point_2(0) = -0.5;

      if (((p - point_1).norm_square() < 0.2 * 0.2) ||
          ((p - point_2).norm_square() < 0.2 * 0.2))
        values(0) = 1;
      else
        values(0) = 0;

      if (p.square() < 0.2 * 0.2)
        values(1) = 1;
      else
        values(1) = 0;
    }

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override
    {
      const unsigned int n_points = points.size();

      Assert(value_list.size() == n_points,
             ExcDimensionMismatch(value_list.size(), n_points));

      for (unsigned int p = 0; p < n_points; ++p)
        RightHandSide<dim>::vector_value(points[p], value_list[p]);
    }
  };


  template <int dim>
  ElasticProblem<dim>::ElasticProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, (this_mpi_process == 0))
#ifdef HEX
    , mapping(1)
    , fe(FE_Q<dim>(1), dim)
#else
    , mapping(FE_SimplexP<dim>(1))
    , fe(FE_SimplexP<dim>(1), dim)
#endif
    , dof_handler(triangulation)
  {}


  template <int dim>
  void
  ElasticProblem<dim>::setup_system()
  {
    GridTools::partition_triangulation_zorder(n_mpi_processes, triangulation);

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::subdomain_wise(dof_handler);

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    hanging_node_constraints,
                                    false);

    const std::vector<IndexSet> locally_owned_dofs_per_proc =
      DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
    const IndexSet locally_owned_dofs =
      locally_owned_dofs_per_proc[this_mpi_process];

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    solution.reinit(locally_owned_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
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

    RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->subdomain_id() == this_mpi_process)
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);

          lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
          mu.value_list(fe_values.get_quadrature_points(), mu_values);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const unsigned int component_j =
                    fe.system_to_component_index(j).first;

                  for (unsigned int q_point = 0; q_point < n_q_points;
                       ++q_point)
                    {
                      cell_matrix(i, j) +=
                        ((fe_values.shape_grad(i, q_point)[component_i] *
                          fe_values.shape_grad(j, q_point)[component_j] *
                          lambda_values[q_point]) +
                         (fe_values.shape_grad(i, q_point)[component_j] *
                          fe_values.shape_grad(j, q_point)[component_i] *
                          mu_values[q_point]) +
                         ((component_i == component_j) ?
                            (fe_values.shape_grad(i, q_point) *
                             fe_values.shape_grad(j, q_point) *
                             mu_values[q_point]) :
                            0)) *
                        fe_values.JxW(q_point);
                    }
                }
            }

          right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                            rhs_values);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                cell_rhs(i) += fe_values.shape_value(i, q_point) *
                               rhs_values[q_point](component_i) *
                               fe_values.JxW(q_point);
            }

          cell->get_dof_indices(local_dof_indices);
          hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                              cell_rhs,
                                                              local_dof_indices,
                                                              system_matrix,
                                                              system_rhs);
        }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(mapping,
                                             dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values);
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, false);
  }

  template <int dim>
  unsigned int
  ElasticProblem<dim>::solve()
  {
    SolverControl solver_control(solution.size(), 1e-8 * system_rhs.l2_norm());
    PETScWrappers::SolverCG cg(solver_control, mpi_communicator);

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    Vector<double> localized_solution(solution);
    hanging_node_constraints.distribute(localized_solution);
    solution = localized_solution;

    return solver_control.last_step();
  }

  template <int dim>
  void
  ElasticProblem<dim>::output_results() const
  {
    const Vector<double> localized_solution(solution);

    if (this_mpi_process == 0)
      {
        std::ofstream output("solution.vtk");

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
              Assert(false, ExcInternalError());
          }

        data_out.add_data_vector(localized_solution, solution_names);

        std::vector<unsigned int> partition_int(triangulation.n_active_cells());
        GridTools::get_subdomain_association(triangulation, partition_int);

        const Vector<double> partitioning(partition_int.begin(),
                                          partition_int.end());

        data_out.add_data_vector(partitioning, "partitioning");

        data_out.build_patches(mapping);
        data_out.write_vtk(output);
      }
  }

  template <int dim>
  void
  ElasticProblem<dim>::run()
  {
    const unsigned int n_subdivisions = 10;
    Point<dim>         a, b;

    switch (dim)
      {
        case 1:
          a = Point<dim>(-1);
          b = Point<dim>(1);
          break;
        case 2:
          a = Point<dim>(-1, -1);
          b = Point<dim>(1, 1);
          break;
        case 3:
          a = Point<dim>(-1, -1, -1);
          b = Point<dim>(1, 1, 1);
          break;
        default:
          Assert(false, ExcInternalError());
      }

#ifdef HEX
    GridGenerator::subdivided_hyper_rectangle(
      triangulation, std::vector<unsigned int>(dim, n_subdivisions), a, b);
#else
    GridGenerator::subdivided_hyper_rectangle_with_simplices(
      triangulation, std::vector<unsigned int>(dim, n_subdivisions), a, b);
#endif

    pcout << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl;

    setup_system();

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << " (by partition:";
    for (unsigned int p = 0; p < n_mpi_processes; ++p)
      pcout << (p == 0 ? ' ' : '+')
            << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
                                                                p));
    pcout << ")" << std::endl;

    assemble_system();
    const unsigned int n_iterations = solve();

    pcout << "   Solver converged in " << n_iterations << " iterations."
          << std::endl;

    output_results();
  }
} // namespace Step17


int
main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step17;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      ElasticProblem<2> elastic_problem;
      elastic_problem.run();
    }
  catch (std::exception &exc)
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
