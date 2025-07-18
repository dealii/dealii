// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// simplified form for step-48 test


#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


namespace Step48
{
  const unsigned int dimension = 2;
  const unsigned int fe_degree = 4;



  template <int dim>
  class SineGordonOperation
  {
  public:
    SineGordonOperation(const MatrixFree<dim, double> &data_in,
                        const double                   time_step);

    void
    apply(LinearAlgebra::distributed::Vector<double>                      &dst,
          const std::vector<LinearAlgebra::distributed::Vector<double> *> &src)
      const;

  private:
    const MatrixFree<dim, double>             &data;
    const VectorizedArray<double>              delta_t_sqr;
    LinearAlgebra::distributed::Vector<double> inv_mass_matrix;

    void
    local_apply(
      const MatrixFree<dim, double>                                   &data,
      LinearAlgebra::distributed::Vector<double>                      &dst,
      const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const;
  };



  template <int dim>
  SineGordonOperation<dim>::SineGordonOperation(
    const MatrixFree<dim, double> &data_in,
    const double                   time_step)
    : data(data_in)
    , delta_t_sqr(make_vectorized_array(time_step * time_step))
  {
    VectorizedArray<double> one = make_vectorized_array(1.);

    data.initialize_dof_vector(inv_mass_matrix);

    FEEvaluation<dim, -1> fe_eval(data);
    const unsigned int    n_q_points = fe_eval.n_q_points;

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        fe_eval.reinit(cell);
        for (unsigned int q = 0; q < n_q_points; ++q)
          fe_eval.submit_value(one, q);
        fe_eval.integrate(EvaluationFlags::values);
        fe_eval.distribute_local_to_global(inv_mass_matrix);
      }

    inv_mass_matrix.compress(VectorOperation::add);
    for (unsigned int k = 0; k < inv_mass_matrix.locally_owned_size(); ++k)
      if (inv_mass_matrix.local_element(k) > 1e-15)
        inv_mass_matrix.local_element(k) =
          1. / inv_mass_matrix.local_element(k);
      else
        inv_mass_matrix.local_element(k) = 0;
  }



  template <int dim>
  void
  SineGordonOperation<dim>::local_apply(
    const MatrixFree<dim>                                           &data,
    LinearAlgebra::distributed::Vector<double>                      &dst,
    const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    AssertDimension(src.size(), 2);
    FEEvaluation<dim, -1> current(data), old(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        current.reinit(cell);
        old.reinit(cell);

        current.read_dof_values(*src[0]);
        old.read_dof_values(*src[1]);

        current.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        old.evaluate(EvaluationFlags::values);

        for (unsigned int q = 0; q < current.n_q_points; ++q)
          {
            const VectorizedArray<double> current_value = current.get_value(q);
            const VectorizedArray<double> old_value     = old.get_value(q);

            current.submit_value(2. * current_value - old_value, q);
            current.submit_gradient(-delta_t_sqr * current.get_gradient(q), q);
          }

        current.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        current.distribute_local_to_global(dst);
      }
  }



  template <int dim>
  void
  SineGordonOperation<dim>::apply(
    LinearAlgebra::distributed::Vector<double>                      &dst,
    const std::vector<LinearAlgebra::distributed::Vector<double> *> &src) const
  {
    dst = 0;
    data.cell_loop(&SineGordonOperation<dim>::local_apply, this, dst, src);
    dst.scale(inv_mass_matrix);
  }



  template <int dim>
  class InitialSolution : public Function<dim>
  {
  public:
    InitialSolution(const unsigned int n_components = 1, const double time = 0.)
      : Function<dim>(n_components, time)
    {}
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
  };

  template <int dim>
  double
  InitialSolution<dim>::value(const Point<dim> &p,
                              const unsigned int /* component */) const
  {
    return 4. * std::exp(-p.square() * 10);
  }



  template <int dim>
  class SineGordonProblem
  {
  public:
    SineGordonProblem();
    void
    run();

  private:
    void
    make_grid_and_dofs();
    void
    output_results(const unsigned int timestep_number);

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif
    FE_Q<dim>                 fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;
    IndexSet                  locally_relevant_dofs;

    MatrixFree<dim, double> matrix_free_data;

    LinearAlgebra::distributed::Vector<double> solution, old_solution,
      old_old_solution;

    const unsigned int n_global_refinements;
    double             time, time_step;
    const double       final_time;
    const double       cfl_number;
    const unsigned int output_timestep_skip;
  };



  template <int dim>
  SineGordonProblem<dim>::SineGordonProblem()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(MPI_COMM_WORLD)
    ,
#endif
    fe(QGaussLobatto<1>(fe_degree + 1))
    , dof_handler(triangulation)
    , n_global_refinements(9 - 2 * dim)
    , time(-10)
    , final_time(-9)
    , cfl_number(.1 / fe_degree)
    , output_timestep_skip(200)
  {}


  template <int dim>
  void
  SineGordonProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -15, 15);
    triangulation.refine_global(n_global_refinements);

    deallog << "   Number of global active cells: "
#ifdef DEAL_II_WITH_P4EST
            << triangulation.n_global_active_cells()
#else
            << triangulation.n_active_cells()
#endif
            << std::endl;

    dof_handler.distribute_dofs(fe);

    deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;


    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    constraints.clear();
    constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    QGaussLobatto<1>                         quadrature(fe_degree + 1);
    typename MatrixFree<dim>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim>::AdditionalData::partition_partition;

    matrix_free_data.reinit(
      MappingQ1<dim>{}, dof_handler, constraints, quadrature, additional_data);

    matrix_free_data.initialize_dof_vector(solution);
    old_solution.reinit(solution);
    old_old_solution.reinit(solution);
  }



  template <int dim>
  void
  SineGordonProblem<dim>::output_results(const unsigned int timestep_number)
  {
    constraints.distribute(solution);
    solution.update_ghost_values();

    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Functions::ZeroFunction<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe_degree + 1),
                                      VectorTools::L2_norm);
    solution.zero_out_ghost_values();
    const double solution_norm =
      std::sqrt(Utilities::MPI::sum(norm_per_cell.norm_sqr(), MPI_COMM_WORLD));

    deallog << "   Time:" << std::setw(8) << std::setprecision(3) << time
            << ", solution norm: " << std::setprecision(5) << std::setw(7)
            << solution_norm << std::endl;
  }



  template <int dim>
  void
  SineGordonProblem<dim>::run()
  {
    make_grid_and_dofs();

    const double local_min_cell_diameter =
      triangulation.last()->diameter() / std::sqrt(dim);
    const double global_min_cell_diameter =
      -Utilities::MPI::max(-local_min_cell_diameter, MPI_COMM_WORLD);
    time_step = cfl_number * global_min_cell_diameter;
    time_step = (final_time - time) / (int((final_time - time) / time_step));
    deallog << "   Time step size: " << time_step
            << ", finest cell: " << global_min_cell_diameter << std::endl
            << std::endl;


    VectorTools::interpolate(dof_handler,
                             InitialSolution<dim>(1, time),
                             solution);
    VectorTools::interpolate(dof_handler,
                             InitialSolution<dim>(1, time - time_step),
                             old_solution);
    output_results(0);

    std::vector<LinearAlgebra::distributed::Vector<double> *>
      previous_solutions;
    previous_solutions.push_back(&old_solution);
    previous_solutions.push_back(&old_old_solution);

    SineGordonOperation<dim> sine_gordon_op(matrix_free_data, time_step);

    unsigned int timestep_number = 1;

    for (time += time_step; time <= final_time;
         time += time_step, ++timestep_number)
      {
        old_old_solution.swap(old_solution);
        old_solution.swap(solution);
        sine_gordon_op.apply(solution, previous_solutions);

        if (timestep_number % output_timestep_skip == 0)
          output_results(timestep_number / output_timestep_skip);
      }
    output_results(timestep_number / output_timestep_skip + 1);

    deallog << std::endl
            << "   Performed " << timestep_number << " time steps." << std::endl
            << std::endl;
  }
} // namespace Step48



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      {
        deallog.push("2d");
        Step48::SineGordonProblem<2> sg_problem;
        sg_problem.run();
        deallog.pop();
      }
      {
        deallog.push("3d");
        Step48::SineGordonProblem<3> sg_problem;
        sg_problem.run();
        deallog.pop();
      }
    }
  else
    {
      {
        Step48::SineGordonProblem<2> sg_problem;
        sg_problem.run();
      }
      {
        Step48::SineGordonProblem<3> sg_problem;
        sg_problem.run();
      }
    }
}
