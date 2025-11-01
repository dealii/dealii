/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2018 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------

 *
 * a light modification of step-37 tutorial to test MappingFEField
 * in parallel with non-zero BC and zero body load.
 */


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


namespace Step37
{
  using namespace dealii;


  const unsigned int degree_finite_element = 2;
  const unsigned int dimension             = 3;


  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    void
    run();

    ~LaplaceProblem();

  private:
    void
    setup_system();
    void
    assemble_rhs();
    void
    solve();
    void
    output_results(const unsigned int cycle) const;

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;

    FE_Q<dim>       fe_euler;
    FESystem<dim>   fe_system;
    DoFHandler<dim> dof_euler;
    std::shared_ptr<
      MappingFEField<dim, dim, LinearAlgebra::distributed::Vector<double>>>
                                               mapping;
    AffineConstraints<double>                  constraints_euler;
    LinearAlgebra::distributed::Vector<double> euler_positions;

    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;
    AffineConstraints<double> non_homogeneous_constraints;
    using SystemMatrixType = MatrixFreeOperators::LaplaceOperator<
      dim,
      degree_finite_element,
      degree_finite_element + 1,
      1,
      LinearAlgebra::distributed::Vector<double>>;
    SystemMatrixType system_matrix;

    MGConstrainedDoFs mg_constrained_dofs;
    using LevelMatrixType = MatrixFreeOperators::LaplaceOperator<
      dim,
      degree_finite_element,
      degree_finite_element + 1,
      1,
      LinearAlgebra::distributed::Vector<float>>;
    MGLevelObject<LevelMatrixType> mg_matrices;

    LinearAlgebra::distributed::Vector<double> solution;
    LinearAlgebra::distributed::Vector<double> system_rhs;

    ConditionalOStream pcout;
  };



  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
    ,
#else
    triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    ,
#endif
    fe(degree_finite_element)
    , dof_handler(triangulation)
    , fe_euler(degree_finite_element)
    , fe_system(fe_euler, dim)
    , dof_euler(triangulation)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {}



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    mapping.reset();
    dof_euler.clear();
  }



  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    system_matrix.clear();
    mg_matrices.clear_elements();

    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    dof_euler.distribute_dofs(fe_system);
    {
      const IndexSet locally_relevant_euler =
        DoFTools::extract_locally_relevant_dofs(dof_euler);
      euler_positions.reinit(dof_euler.locally_owned_dofs(),
                             locally_relevant_euler,
                             MPI_COMM_WORLD);
    }

    // Set up a quadrature formula based on the support points of FE_Q and go
    // through the mesh with a MappingQ of your favorite degree; on each cell,
    // evaluate the quadrature points and write them into the vector to be
    // associated with MappingFEField
    VectorTools::get_position_vector(dof_euler, euler_positions);

    // Apply hanging node constraints on that vector
    constraints_euler.clear();
    constraints_euler.reinit(dof_euler.locally_owned_dofs(),
                             DoFTools::extract_locally_relevant_dofs(
                               dof_euler));
    DoFTools::make_hanging_node_constraints(dof_euler, constraints_euler);
    constraints_euler.close();
    constraints_euler.distribute(euler_positions);

    euler_positions.update_ghost_values();

    mapping = std::make_shared<
      MappingFEField<dim, dim, LinearAlgebra::distributed::Vector<double>>>(
      dof_euler, euler_positions);

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    constraints.clear();
    constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    {
      typename MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      additional_data.mapping_update_flags =
        (update_gradients | update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
        new MatrixFree<dim, double>());
      system_mf_storage->reinit(*mapping.get(),
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);
      system_matrix.initialize(system_mf_storage);
    }

    system_matrix.compute_diagonal();

    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(system_rhs);

    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);

    const std::set<types::boundary_id> dirichlet_boundary = {0};
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary);

    for (unsigned int level = 0; level < nlevels; ++level)
      {
        const IndexSet relevant_dofs =
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
        AffineConstraints<double> level_constraints;
        level_constraints.reinit(dof_handler.locally_owned_mg_dofs(level),
                                 relevant_dofs);
        level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim, float>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, float>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level = level;
        std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level(
          new MatrixFree<dim, float>());
        mg_mf_storage_level->reinit(MappingQ1<dim>{},
                                    dof_handler,
                                    level_constraints,
                                    QGauss<1>(fe.degree + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);
        mg_matrices[level].compute_diagonal();
      }
  }



  template <int dim>
  class PotentialBCFunction : public Function<dim>
  {
  public:
    PotentialBCFunction(const double &charge, const Point<dim> &dipole)
      : Function<dim>(1)
      , charge(charge)
      , x0(std::abs(charge) < 1e-10 ? Point<dim>() : dipole / charge)
    {}

    virtual ~PotentialBCFunction() = default;

    virtual double
    value(const Point<dim> &p, const unsigned int) const override
    {
      const double r = p.distance(x0);
      Assert(r > 0, ExcDivideByZero());
      return charge / r;
    }

  private:
    const double     charge;
    const Point<dim> x0;
  };



  template <int dim>
  void
  LaplaceProblem<dim>::assemble_rhs()
  {
    Timer time;

    non_homogeneous_constraints.clear();
    {
      AffineConstraints<double> hanging_nodes_laplace_constraints;
      hanging_nodes_laplace_constraints.reinit(dof_handler.locally_owned_dofs(),
                                               locally_relevant_dofs);
      non_homogeneous_constraints.reinit(dof_handler.locally_owned_dofs(),
                                         locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints(
        dof_handler, hanging_nodes_laplace_constraints);

      const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
      std::map<types::boundary_id, const Function<dim> *>
                               dirichlet_boundary_functions;
      PotentialBCFunction<dim> bc_func(240, Point<dim>());
      dirichlet_boundary_functions[0] = &bc_func;
      VectorTools::interpolate_boundary_values(*mapping.get(),
                                               dof_handler,
                                               dirichlet_boundary_functions,
                                               non_homogeneous_constraints);
      // make sure hanging nodes override Dirichlet
      non_homogeneous_constraints.merge(
        hanging_nodes_laplace_constraints,
        AffineConstraints<double>::MergeConflictBehavior::right_object_wins);
      non_homogeneous_constraints.close();
    }

    solution = 0;
    non_homogeneous_constraints.distribute(solution);
    system_rhs = 0;
    solution.update_ghost_values();
    FEEvaluation<dim, degree_finite_element> phi(
      *system_matrix.get_matrix_free());
    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches();
         ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values_plain(solution);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            phi.submit_gradient(-phi.get_gradient(q), q);
            // phi.submit_value(make_vectorized_array<double>(1.0), q);
          }
        // phi.integrate(EvaluationFlags::values|EvaluationFlags::gradients);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(system_rhs);
      }
    system_rhs.compress(VectorOperation::add);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);

    using SmootherType =
      PreconditionChebyshev<LevelMatrixType,
                            LinearAlgebra::distributed::Vector<float>>;
    mg::SmootherRelaxation<SmootherType,
                           LinearAlgebra::distributed::Vector<float>>
                                                         mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        if (level > 0)
          {
            smoother_data[level].smoothing_range     = 15.;
            smoother_data[level].degree              = 4;
            smoother_data[level].eig_cg_n_iterations = 10;
          }
        else
          {
            smoother_data[0].smoothing_range = 1e-3;
            smoother_data[0].degree          = numbers::invalid_unsigned_int;
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
          }
        mg_matrices[level].compute_diagonal();
        smoother_data[level].preconditioner =
          mg_matrices[level].get_matrix_diagonal_inverse();
      }
    mg_smoother.initialize(mg_matrices, smoother_data);

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
      mg_coarse;
    mg_coarse.initialize(mg_smoother);

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
      mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
      mg_interface_matrices);

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    PreconditionMG<dim,
                   LinearAlgebra::distributed::Vector<float>,
                   MGTransferMatrixFree<dim, float>>
      preconditioner(dof_handler, mg, mg_transfer);


    SolverControl solver_control(500, 1e-12 * system_rhs.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);

    non_homogeneous_constraints.set_zero(solution);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    non_homogeneous_constraints.distribute(solution);

    const double linfty =
      Utilities::MPI::max(solution.linfty_norm(), MPI_COMM_WORLD);

    pcout << "Solved in " << solver_control.last_step() << " iterations"
          << std::endl
          << "Linfty=" << linfty << std::endl;
  }



  template <int dim>
  void
  LaplaceProblem<dim>::output_results(const unsigned int cycle) const
  {
    if (triangulation.n_global_active_cells() > 1000000)
      return;

    DataOut<dim> data_out;

    solution.update_ghost_values();
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output(
      "solution-" + std::to_string(cycle) + "." +
      std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) +
      ".vtu");
    data_out.write_vtu(output);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
             ++i)
          filenames.emplace_back("solution-" + std::to_string(cycle) + "." +
                                 std::to_string(i) + ".vtu");

        std::string pvtu_filename =
          "solution-" + Utilities::to_string(cycle) + ".pvtu";
        std::ofstream pvtu_output(pvtu_filename);
        data_out.write_pvtu_record(pvtu_output, filenames);
      }

    if (dim == 2 && Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1)
      {
        std::map<types::global_dof_index, Point<dim>> support_points;
        MappingQ<dim> mapping(degree_finite_element);
        DoFTools::map_dofs_to_support_points(mapping,
                                             dof_handler,
                                             support_points);

        const std::string base_filename =
          "grid" + dealii::Utilities::int_to_string(dim) + "_" +
          dealii::Utilities::int_to_string(cycle);
        const std::string filename = base_filename + ".gp";
        std::ofstream     f(filename);

        f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
          << std::endl
          << "set output \"" << base_filename << ".png\"" << std::endl
          << "set size square" << std::endl
          << "set view equal xy" << std::endl
          << "unset xtics" << std::endl
          << "unset ytics" << std::endl
          << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
          << std::endl;
        GridOut().write_gnuplot(triangulation, f);
        f << 'e' << std::endl;

        DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

        f << 'e' << std::endl;
      }
  }



  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 2; ++cycle)
      {
        pcout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            Point<dim> center;
            GridGenerator::hyper_ball(triangulation, center, 12);

            const types::manifold_id            sphere_id = 0;
            static const SphericalManifold<dim> boundary_ball(center);
            triangulation.set_all_manifold_ids_on_boundary(sphere_id);
            triangulation.set_manifold(sphere_id, boundary_ball);

            triangulation.refine_global(dim == 2 ? 1 : 3);
          }
        else
          {
            for (auto &cell : triangulation.active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  if ((cell->center()[0] < -8. && cell->center()[1] < 0) ||
                      (cell->center()[0] > 9. && cell->center()[1] > 5))
                    cell->set_refine_flag();
                }
            triangulation.execute_coarsening_and_refinement();
            // triangulation.refine_global();
          }
        setup_system();
        assemble_rhs();

        if (cycle == 0 && dim == 2)
          {
            pcout << "Constraints:" << std::endl;
            non_homogeneous_constraints.print(std::cout);
          }

        solve();
        output_results(cycle);
        pcout << std::endl;
      };
  }
} // namespace Step37



int
main(int argc, char *argv[])
{
  try
    {
      using namespace Step37;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      LaplaceProblem<dimension> laplace_problem;
      laplace_problem.run();
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
