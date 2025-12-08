// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test Portable::MatrixFree with the Stokes operator using multiple
// DoFHandlers. We create a manufactured solution, apply the Stokes
// operator once, and compute errors.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

#include "../tests.h"


template <int dim>
class VelocityRightHandSide : public Function<dim>
{
public:
  VelocityRightHandSide()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};

template <int dim>
double
VelocityRightHandSide<dim>::value(const Point<dim>  &p,
                                  const unsigned int component) const
{
  const double x   = p[0];
  const double y   = p[1];
  const double pi  = numbers::PI;
  const double pi2 = pi * pi;

  const double sx = std::sin(pi * x);
  const double cx = std::cos(pi * x);
  const double sy = std::sin(pi * y);
  const double cy = std::cos(pi * y);

  if (component == 0)
    return pi * cy * (16.0 * pi2 * sx * sx * sy - 4.0 * pi2 * sy - sx);
  else
    return pi * cx * (-16.0 * pi2 * sx * sy * sy + 4.0 * pi2 * sx - sy);
}

template <int dim>
class VelocitySolution : public Function<dim>
{
public:
  VelocitySolution()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};

template <int dim>
double
VelocitySolution<dim>::value(const Point<dim>  &p,
                             const unsigned int component) const
{
  const double x  = p[0];
  const double y  = p[1];
  const double pi = numbers::PI;

  const double sx  = std::sin(pi * x);
  const double sy  = std::sin(pi * y);
  const double s2x = std::sin(2.0 * pi * x);
  const double s2y = std::sin(2.0 * pi * y);

  if (component == 0)
    return pi * sx * sx * s2y;
  else
    return -pi * s2x * sy * sy;
}

template <int dim>
class PressureSolution : public Function<dim>
{
public:
  PressureSolution()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p,
        const unsigned int /*component*/ = 0) const override
  {
    return std::cos(numbers::PI * p[0]) * std::cos(numbers::PI * p[1]);
  }
};



template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class StokesOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceBlockVector<Number>              &src,
             Portable::DeviceBlockVector<Number>                    &dst) const;
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
DEAL_II_HOST_DEVICE void
StokesOperator<dim, degree_u, degree_p, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *data,
  const Portable::DeviceBlockVector<Number>              &src,
  Portable::DeviceBlockVector<Number>                    &dst) const
{
  Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim> fe_u(data, 0);
  Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1>   fe_p(data, 1);

  fe_u.read_dof_values(src.block(0));
  fe_p.read_dof_values(src.block(1));
  fe_u.evaluate(EvaluationFlags::gradients);
  fe_p.evaluate(EvaluationFlags::values);

  data->for_each_quad_point([&](const int &q_point) {
    const Tensor<2, dim, Number> gradient_u = fe_u.get_gradient(q_point);
    Tensor<2, dim, Number>       vel_term   = gradient_u;
    for (unsigned int d = 0; d < dim; ++d)
      vel_term[d][d] -= fe_p.get_value(q_point);
    fe_u.submit_gradient(vel_term, q_point);

    const Number pressure_term = trace(gradient_u);
    fe_p.submit_value(pressure_term, q_point);
  });

  fe_u.integrate(EvaluationFlags::gradients);
  fe_p.integrate(EvaluationFlags::values);
  fe_u.distribute_local_to_global(dst.block(0));
  fe_p.distribute_local_to_global(dst.block(1));
}

template <
  int dim,
  int degree_u,
  int degree_p,
  typename Number = double,
  typename VectorType =
    LinearAlgebra::distributed::BlockVector<double, MemorySpace::Default>,
  int n_q_points_1d = degree_u + 1>
class PortableMFStokesOperator
{
public:
  PortableMFStokesOperator(const Portable::MatrixFree<dim, double> &data_in)
    : data(data_in)
  {}

  const Portable::MatrixFree<dim, Number> &data;

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = static_cast<Number>(0.);
    StokesOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
      stokes_operator;
    data.cell_loop(stokes_operator, src, dst);

    data.copy_constrained_values(src, dst);
  }
};


template <int dim, int fe_degree>
void
test(unsigned int n_refinements)
{
  using Number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  const unsigned int degree_u = fe_degree + 1;
  const unsigned int degree_p = fe_degree;

  FESystem<dim>   fe_u(FE_Q<dim>(degree_u), dim);
  FE_Q<dim>       fe_p(degree_p);
  DoFHandler<dim> dof_u(tria);
  DoFHandler<dim> dof_p(tria);
  dof_u.distribute_dofs(fe_u);
  dof_p.distribute_dofs(fe_p);

  const IndexSet &owned_set_u = dof_u.locally_owned_dofs();
  const IndexSet  relevant_set_u =
    DoFTools::extract_locally_relevant_dofs(dof_u);
  AffineConstraints<double> constraints_u(owned_set_u, relevant_set_u);
  DoFTools::make_hanging_node_constraints(dof_u, constraints_u);
  VectorTools::interpolate_boundary_values(dof_u,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
                                           constraints_u);
  constraints_u.close();

  const IndexSet &owned_set_p = dof_p.locally_owned_dofs();
  const IndexSet  relevant_set_p =
    DoFTools::extract_locally_relevant_dofs(dof_p);
  AffineConstraints<double> constraints_p(owned_set_p, relevant_set_p);
  DoFTools::make_hanging_node_constraints(dof_p, constraints_p);
  constraints_p.close();

  std::vector<const DoFHandler<dim> *> dof_handlers          = {&dof_u, &dof_p};
  std::vector<const AffineConstraints<double> *> constraints = {&constraints_u,
                                                                &constraints_p};

  MappingQ<dim>                     mapping(fe_degree);
  Portable::MatrixFree<dim, Number> mf_data;
  const QGauss<1>                   quad(fe_degree + 2);
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         update_quadrature_points;
  mf_data.reinit(mapping, dof_handlers, constraints, quad, additional_data);

  PortableMFStokesOperator<dim, degree_u, degree_p> stokes_operator(mf_data);

  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> src;
  mf_data.initialize_dof_vector(src);
  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> dst;
  mf_data.initialize_dof_vector(dst);

  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Host> src_host;
  mf_data.initialize_dof_vector(src_host);

  VectorTools::interpolate(dof_u, VelocitySolution<dim>(), src_host.block(0));
  VectorTools::interpolate(dof_p, PressureSolution<dim>(), src_host.block(1));
  src.block(0).import_elements(src_host.block(0), VectorOperation::insert);
  src.block(1).import_elements(src_host.block(1), VectorOperation::insert);

  stokes_operator.vmult(dst, src);

  src_host.block(0).import_elements(dst.block(0), VectorOperation::insert);
  src_host.block(1).import_elements(dst.block(1), VectorOperation::insert);

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> f_host;
  f_host.reinit(mf_data.get_vector_partitioner(0));
  VectorTools::create_right_hand_side(dof_u,
                                      QGauss<dim>(degree_u + 2),
                                      VelocityRightHandSide<dim>(),
                                      f_host,
                                      constraints_u);

  src_host.block(0) -= f_host;

  deallog << "velocity error: " << src_host.block(0).l2_norm()
          << " pressure error: " << src_host.block(1).l2_norm() << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  MPILogInitAll mpi_inilog;

  test<2, 1>(2);
  test<2, 1>(3);
  test<2, 1>(4);
  test<2, 1>(5);

  deallog << "OK" << std::endl;
}
