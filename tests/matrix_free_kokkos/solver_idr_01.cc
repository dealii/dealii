// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test SolverIDR with a portable matrix-free linear elasticity operator.


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_idr.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

#include "../tests.h"

using namespace dealii;

namespace
{
  constexpr double lambda = 1.;
  constexpr double mu     = 1.;

  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution()
      : Function<dim>(dim)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      AssertIndexRange(component, dim);
      double result = 1.;
      for (unsigned int d = 0; d < dim; ++d)
        result *= p[d] * (1. - p[d]);
      return result;
    }
  };

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      AssertIndexRange(component, dim);

      double laplacian = 0.;
      for (unsigned int d = 0; d < dim; ++d)
        {
          double other = 1.;
          for (unsigned int e = 0; e < dim; ++e)
            if (e != d)
              other *= p[e] * (1. - p[e]);
          laplacian += -2. * other;
        }

      double grad_div = 0.;
      for (unsigned int d = 0; d < dim; ++d)
        {
          if (d == component)
            {
              double other = 1.;
              for (unsigned int e = 0; e < dim; ++e)
                if (e != d)
                  other *= p[e] * (1. - p[e]);
              grad_div += -2. * other;
            }
          else
            {
              double other = 1.;
              for (unsigned int e = 0; e < dim; ++e)
                if (e != d && e != component)
                  other *= p[e] * (1. - p[e]);
              grad_div += (1. - 2. * p[d]) * (1. - 2. * p[component]) * other;
            }
        }

      // -div(lambda div(u) I + mu (grad(u) + grad(u)^T)).
      return -mu * laplacian - (lambda + mu) * grad_div;
    }
  };
} // namespace


template <int dim, int fe_degree, typename Number>
class ElasticityCellOperatorQuadrature
{
public:
  DEAL_II_HOST_DEVICE
  ElasticityCellOperatorQuadrature(const double mu, const double lambda)
    : mu(mu)
    , lambda(lambda)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> *fe_eval,
    const int q_point) const
  {
    const Tensor<2, dim, Number> gradient = fe_eval->get_gradient(q_point);
    Tensor<2, dim, Number>       stress = mu * (gradient + transpose(gradient));
    const Number                 volumetric_stress = lambda * trace(gradient);
    for (unsigned int d = 0; d < dim; ++d)
      stress[d][d] += volumetric_stress;
    fe_eval->submit_gradient(stress, q_point);
  }

private:
  double mu;
  double lambda;
};

template <int dim, int degree, typename Number>
class ElasticityCellOperator
{
public:
  static constexpr unsigned int n_q_points = Utilities::pow(degree + 1, dim);

  ElasticityCellOperator(const double mu, const double lambda)
    : mu(mu)
    , lambda(lambda)
  {}

  DEAL_II_HOST_DEVICE
  void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, degree, degree + 1, dim, Number> fe(data, 0);
    fe.read_dof_values(src);
    fe.evaluate(EvaluationFlags::gradients);

    ElasticityCellOperatorQuadrature<dim, degree, Number> quad_operator(mu,
                                                                        lambda);
    data->for_each_quad_point(
      [&](const int q_point) { quad_operator(&fe, q_point); });

    fe.integrate(EvaluationFlags::gradients);
    fe.distribute_local_to_global(dst);
  }

private:
  double mu;
  double lambda;
};

template <int dim, int degree, typename Number = double>
class ElasticityOperator : public EnableObserverPointer
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  ElasticityOperator() = default;

  explicit ElasticityOperator(
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data)
    : data(std::move(data))
  {}

  void
  reinit(std::shared_ptr<Portable::MatrixFree<dim, Number>> data_in)
  {
    data = std::move(data_in);
  }

  void
  initialize_dof_vector(VectorType &v) const
  {
    data->initialize_dof_vector(v, 0);
  }

  types::global_dof_index
  m() const
  {
    return data->get_vector_partitioner(0)->size();
  }

  double
  el(const types::global_dof_index row,
     const types::global_dof_index column) const
  {
    Assert(row == column, ExcNotImplemented());
    Assert(inverse_diagonal.get() != nullptr, ExcNotInitialized());
    return 1. / (*inverse_diagonal)(row, row);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0.;
    data->cell_loop(ElasticityCellOperator<dim, degree, Number>(mu, lambda),
                    src,
                    dst);
    data->copy_constrained_values(src, dst, 0);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_diagonal()
  {
    inverse_diagonal     = std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &diagonal = inverse_diagonal->get_vector();
    data->initialize_dof_vector(diagonal, 0);
    ElasticityCellOperatorQuadrature<dim, degree, Number> quad_operator(mu,
                                                                        lambda);

    MatrixFreeTools::compute_diagonal<dim, degree, degree + 1, dim, Number>(
      *data,
      diagonal,
      quad_operator,
      EvaluationFlags::gradients,
      EvaluationFlags::gradients);

    double *values = diagonal.get_values();
    Kokkos::parallel_for(
      "invert diagonal",
      diagonal.locally_owned_size(),
      KOKKOS_LAMBDA(const int i) {
        values[i] = (values[i] == 0.0) ? 1.0 : 1. / values[i];
      });
  }

  std::shared_ptr<DiagonalMatrix<VectorType>>
  diagonal_inverse() const
  {
    return inverse_diagonal;
  }

private:
  std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
  std::shared_ptr<DiagonalMatrix<VectorType>>        inverse_diagonal;
};

template <int dim, int degree>
void
test(const unsigned int refinement)
{
  using VectorType =
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>;
  using Operator = ElasticityOperator<dim, degree>;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(refinement);
  MappingQ1<dim>  mapping;
  const FE_Q<dim> scalar_fe(degree);
  FESystem<dim>   fe(scalar_fe, dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "refinement: " << refinement << " n_dofs: " << dof_handler.n_dofs()
          << std::endl;

  AffineConstraints<double> constraints;
  constraints.reinit(dof_handler.locally_owned_dofs(),
                     DoFTools::extract_locally_relevant_dofs(dof_handler));
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  auto data = std::make_shared<Portable::MatrixFree<dim, double>>();
  typename Portable::MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  data->reinit(
    mapping, dof_handler, constraints, QGauss<1>(degree + 1), additional_data);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> rhs_host;
  data->initialize_dof_vector(rhs_host, 0);
  FEValues<dim>                        fe_values(mapping,
                          fe,
                          QGauss<dim>(degree + 1),
                          update_values | update_quadrature_points |
                            update_JxW_values);
  Vector<double>                       cell_rhs(fe.n_dofs_per_cell());
  std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
  RightHandSide<dim>                   forcing;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        cell_rhs = 0.;
        for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            cell_rhs[i] +=
              forcing.value(fe_values.quadrature_point(q),
                            fe.system_to_component_index(i).first) *
              fe_values.shape_value(i, q) * fe_values.JxW(q);
        cell->get_dof_indices(indices);
        constraints.distribute_local_to_global(cell_rhs, indices, rhs_host);
      }
  rhs_host.compress(VectorOperation::add);
  constraints.distribute(rhs_host);

  VectorType rhs, solution;
  data->initialize_dof_vector(rhs, 0);
  data->initialize_dof_vector(solution, 0);
  rhs.import_elements(rhs_host, VectorOperation::insert);
  Operator op(data);

  SolverControl         control(100, 1e-10 * rhs.l2_norm());
  SolverIDR<VectorType> solver(control);
  PreconditionIdentity  preconditioner;
  check_solver_within_range(solver.solve(op, solution, rhs, preconditioner),
                            control.last_step(),
                            5,
                            20);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> solution_host;
  data->initialize_dof_vector(solution_host, 0);
  solution_host.import_elements(solution, VectorOperation::insert);
  constraints.distribute(solution_host);
  solution_host.update_ghost_values();
  Vector<double> cellwise_errors(tria.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution_host,
                                    ExactSolution<dim>(),
                                    cellwise_errors,
                                    QGauss<dim>(degree + 1),
                                    VectorTools::L2_norm);
  const double error = VectorTools::compute_global_error(tria,
                                                         cellwise_errors,
                                                         VectorTools::L2_norm);
  deallog << "Error: " << error << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  initlog();
  test<2, 1>(3);
  deallog << "OK" << std::endl;
}
