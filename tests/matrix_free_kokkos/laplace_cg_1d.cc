// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Solve a Laplace Laplace problem in (0,1).

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

constexpr int dim = 1;
using Number      = double;
using VectorType =
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

template <int fe_degree>
class LaplaceOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);
    phi.read_dof_values(src);
    phi.evaluate(EvaluationFlags::gradients);

    data->for_each_quad_point([&](const int &q_point) {
      phi.submit_gradient(phi.get_gradient(q_point), q_point);
    });

    phi.integrate(EvaluationFlags::gradients);
    phi.distribute_local_to_global(dst);
  }

  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);
};



template <int fe_degree>
class LaplaceOperator : public EnableObserverPointer
{
public:
  void
  reinit(const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints)
  {
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_JxW_values | update_gradients;

    const MappingQ1<dim> mapping;
    const QGauss<1>      quadrature_1d(fe_degree + 1);

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature_1d, additional_data);
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0.0;
    LaplaceOperatorLocal<fe_degree> local_operator;
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

private:
  Portable::MatrixFree<dim, Number> matrix_free;
};



// Exact solution of -u'' = 1 on (0,1), u(0) = u(1) = 0.
class ExactSolution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    return 0.5 * p[0] * (1.0 - p[0]);
  }
};


template <int fe_degree>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  triangulation.refine_global(8);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  LaplaceOperator<fe_degree> laplace_operator;
  laplace_operator.reinit(dof_handler, constraints);

  VectorType solution, system_rhs;
  laplace_operator.initialize_dof_vector(solution);
  laplace_operator.initialize_dof_vector(system_rhs);
  solution = 0.0;

  {
    const MappingQ1<dim> mapping;
    const QGauss<dim>    quadrature(fe_degree + 1);

    LinearAlgebra::distributed::Vector<Number> system_rhs_host(
      dof_handler.n_dofs());

    VectorTools::create_right_hand_side<dim, dim>(
      mapping,
      dof_handler,
      quadrature,
      Functions::ConstantFunction<dim>(1.0),
      system_rhs_host,
      constraints);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(dof_handler.n_dofs());
    rw_vector.import_elements(system_rhs_host, VectorOperation::insert);
    system_rhs.import_elements(rw_vector, VectorOperation::insert);
  }

  SolverControl        solver_control(1000, 1e-12 * system_rhs.l2_norm());
  SolverCG<VectorType> solver(solver_control);

  solver.solve(laplace_operator, solution, system_rhs, PreconditionIdentity());

  LinearAlgebra::distributed::Vector<Number> solution_host(
    dof_handler.n_dofs());
  LinearAlgebra::ReadWriteVector<Number> rw_vector(dof_handler.n_dofs());
  rw_vector.import_elements(solution, VectorOperation::insert);
  solution_host.import_elements(rw_vector, VectorOperation::insert);

  LinearAlgebra::distributed::Vector<Number> interpolated_exact(
    dof_handler.n_dofs());
  {
    const MappingQ1<dim> mapping;
    VectorTools::interpolate(mapping,
                             dof_handler,
                             ExactSolution(),
                             interpolated_exact);
  }

  LinearAlgebra::distributed::Vector<Number> diff = solution_host;
  diff -= interpolated_exact;
  const double relative_error = diff.l2_norm() / interpolated_exact.l2_norm();

  AssertThrow(relative_error < 1e-12,
              ExcMessage("Relative error " + std::to_string(relative_error) +
                         " exceeds 1e-12 for degree " +
                         std::to_string(fe_degree)));

  deallog << "degree " << fe_degree << ": relative error < 1e-12: OK"
          << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Kokkos::initialize();


  test<1>();
  test<2>();
  test<3>();

  Kokkos::finalize();

  deallog << "OK" << std::endl;
}
