/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2025 by the deal.II authors
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
 * Description:
 *
 * This test compares the MatrixFree and Portable::MatrixFree
 * infrastructure on the CPU. Considered are the initialization
 * costs and the costs for an operator evaluation.
 * Portable::MatrixFree was written with CUDA and now uses
 * Kokkos as backend and, as consequence, favors GPU hardware. This
 * performance test is meant to track the improvement of
 * the performance of Portable::MatrixFree on the CPU.
 *
 * Status: experimental
 */

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#define ENABLE_MPI

#include "performance_test_driver.h"

using namespace dealii;



template <int dim, int fe_degree, typename Number, typename MemorySpace>
class LaplaceOperator;

template <int dim, int fe_degree, typename Number>
class LaplaceOperator<dim, fe_degree, Number, MemorySpace::Host>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  LaplaceOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature)
  {
    typename MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.cell_loop(&LaplaceOperator::local_apply, this, dst, src);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        phi.read_dof_values_plain(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  MatrixFree<dim, Number> matrix_free;
};



template <int dim, int fe_degree, typename Number>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> *fe_eval,
    const int q_point) const
  {
    fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
  }
};

template <int dim, int fe_degree, typename Number>
class LaplaceOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *gpu_data,
             const Portable::DeviceVector<double>                   &src,
             Portable::DeviceVector<double>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(
      gpu_data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::gradients);
    fe_eval.apply_for_each_quad_point(
      LaplaceOperatorQuad<dim, fe_degree, Number>());
    fe_eval.integrate(EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
  }
  static const unsigned int n_dofs_1d    = fe_degree + 1;
  static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);
};

template <int dim, int fe_degree, typename Number>
class LaplaceOperator<dim, fe_degree, Number, MemorySpace::Default>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  LaplaceOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature)
  {
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_JxW_values | update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    LaplaceOperatorLocal<dim, fe_degree, Number> local_operator;
    matrix_free.cell_loop(local_operator, src, dst);
  }

private:
  Portable::MatrixFree<dim, Number> matrix_free;
};



template <int dim, typename T>
class AnalyticalFunction : public Function<dim, T>
{
public:
  virtual T
  value(const Point<dim, T> &p, const unsigned int component = 0) const override
  {
    (void)component;

    double temp = 0.0;

    for (unsigned int d = 0; d < dim; ++d)
      temp += std::sin(p[d]);

    return temp;
  }
};



template <unsigned int dim, const int degree, typename MemorySpace>
std::vector<double>
run(const unsigned int n_refinements)
{
  ConvergenceTable table;

  const MPI_Comm comm = MPI_COMM_WORLD;

  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  const unsigned n_repetitions_setup = 10;
  const unsigned n_repetitions_vmult = 100;

  parallel::distributed::Triangulation<dim> tria(comm);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  table.add_value("n_levels", tria.n_global_levels());
  table.add_value("degree", degree);

  table.add_value("n_cells", tria.n_global_active_cells());

  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe(degree);
  const QGauss<1>      quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  table.add_value("n_dofs", dof_handler.n_dofs());

  AffineConstraints<Number>                         constraints;
  LaplaceOperator<dim, degree, Number, MemorySpace> laplace_operator;

  std::chrono::time_point<std::chrono::system_clock> now_setup =
    std::chrono::system_clock::now();

  for (unsigned int i = 0; i < n_repetitions_setup; ++i)
    laplace_operator.reinit(mapping, dof_handler, constraints, quadrature);

  double dt_setup = std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::system_clock::now() - now_setup)
                      .count() /
                    1e9;

  dt_setup = Utilities::MPI::sum(dt_setup, MPI_COMM_WORLD) /
             Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  VectorType src, dst;

  laplace_operator.initialize_dof_vector(src);
  laplace_operator.initialize_dof_vector(dst);

  {
    LinearAlgebra::distributed::Vector<Number> src_host(src.get_partitioner());

    VectorTools::interpolate(dof_handler,
                             AnalyticalFunction<dim, Number>(),
                             src_host);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src_host, VectorOperation::insert);
    src.import_elements(rw_vector, VectorOperation::insert);

    dst = 0.0;
  }

  const std::chrono::time_point<std::chrono::system_clock> now_vmult =
    std::chrono::system_clock::now();

  for (unsigned int i = 0; i < n_repetitions_vmult; ++i)
    laplace_operator.vmult(dst, src);

  double dt_vmult = std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::system_clock::now() - now_vmult)
                      .count() /
                    1e9;

  dt_vmult = Utilities::MPI::sum(dt_vmult, MPI_COMM_WORLD) /
             Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);


  table.add_value("time_setup", dt_setup);
  table.set_scientific("time_setup", true);
  table.add_value("time_avg", dt_vmult);
  table.set_scientific("time_avg", true);

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
#if 0
      table.write_text(std::cout);
      std::cout << std::endl;
#endif
    }

  return {dt_setup, dt_vmult};
}


std::tuple<Metric, unsigned int, std::vector<std::string>>
describe_measurements()
{
  return {Metric::timing,
          4,
          {"mf_setup", "mf_vmult", "mf_kokkos_setup", "mf_kokkos_vmult"}};
}

Measurement
perform_single_measurement()
{
  const unsigned int dim           = 3;
  const unsigned int fe_degree     = 4;
  unsigned int       n_refinements = 5;

  switch (get_testing_environment())
    {
      case TestingEnvironment::light:
        break;
      case TestingEnvironment::medium:
        break;
      case TestingEnvironment::heavy:
        n_refinements += 1;
        break;
    }

  const auto result0 = run<dim, fe_degree, MemorySpace::Host>(n_refinements);
  const auto result1 = run<dim, fe_degree, MemorySpace::Default>(n_refinements);

  return {result0[0], result0[1], result1[0], result1[1]};
}
