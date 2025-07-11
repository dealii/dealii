// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test vector valued FEM for Portable::MatrixFree and Porable::FEEvaluation.

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "Kokkos_Core.hpp"

template <int dim,
          int fe_degree,
          int n_components,
          typename Number,
          typename MemorySpace>
class LaplaceOperator;

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim, fe_degree, n_components, Number, MemorySpace::Host>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  LaplaceOperator(const unsigned int version)
    : version(version)
  {}

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
    matrix_free.cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number> phi(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values_plain(src);

        if (version == 0)
          phi.evaluate(EvaluationFlags::values);
        else if (version == 1)
          phi.evaluate(EvaluationFlags::gradients);
        else if (version == 2)
          phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            if (version == 0 || version == 2)
              phi.submit_value(phi.get_value(q), q);
            if (version == 1 || version == 2)
              phi.submit_gradient(phi.get_gradient(q), q);
          }

        if (version == 0)
          phi.integrate(EvaluationFlags::values);
        else if (version == 1)
          phi.integrate(EvaluationFlags::gradients);
        else if (version == 2)
          phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);

        phi.distribute_local_to_global(dst);
      }
  }

  const unsigned int      version;
  MatrixFree<dim, Number> matrix_free;
};



template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE
  LaplaceOperatorQuad(const unsigned int version)
    : version(version)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
             *phi,
    const int q) const
  {
    if (version == 0 || version == 2)
      phi->submit_value(phi->get_value(q), q);
    if (version == 1 || version == 2)
      phi->submit_gradient(phi->get_gradient(q), q);
  }

private:
  const unsigned int version;
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorLocal
{
public:
  LaplaceOperatorLocal(const unsigned int version)
    : version(version)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
      phi(data);
    phi.read_dof_values(src);

    if (version == 0)
      phi.evaluate(EvaluationFlags::values);
    else if (version == 1)
      phi.evaluate(EvaluationFlags::gradients);
    else if (version == 2)
      phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    phi.apply_for_each_quad_point(
      LaplaceOperatorQuad<dim, fe_degree, n_components, Number>(version));

    if (version == 0)
      phi.integrate(EvaluationFlags::values);
    else if (version == 1)
      phi.integrate(EvaluationFlags::gradients);
    else if (version == 2)
      phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);

    phi.distribute_local_to_global(dst);
  }
  static const unsigned int n_local_dofs =
    Utilities::pow(fe_degree + 1, dim) * n_components;
  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);

private:
  const unsigned int version;
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim,
                      fe_degree,
                      n_components,
                      Number,
                      MemorySpace::Default>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  LaplaceOperator(const unsigned int version)
    : version(version)
  {}

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
    dst = 0.0; // TODO: annoying
    LaplaceOperatorLocal<dim, fe_degree, n_components, Number> local_operator(
      version);
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst); // TODO: annoying
  }

private:
  const unsigned int                version;
  Portable::MatrixFree<dim, Number> matrix_free;
};



template <int dim, typename T>
class AnalyticalFunction : public Function<dim, T>
{
public:
  AnalyticalFunction(const unsigned int n_components)
    : Function<dim, T>(n_components)
  {}

  virtual T
  value(const Point<dim, T> &p, const unsigned int component = 0) const override
  {
    double temp = 0.0;

    for (unsigned int d = 0; d < dim; ++d)
      temp += std::sin(p[d]);

    return temp * (1.0 + component);
  }
};



template <unsigned int dim,
          const int    degree,
          int          n_components,
          typename MemorySpace>
void
run(const unsigned int n_refinements, ConvergenceTable &table)
{
  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe_q(degree);
  const FESystem<dim>  fe(fe_q, n_components);
  const QGauss<dim>    quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  for (unsigned int i = 0; i < 3; ++i)
    {
      LaplaceOperator<dim, degree, n_components, Number, MemorySpace>
        laplace_operator(i);

      laplace_operator.reinit(mapping,
                              dof_handler,
                              constraints,
                              quadrature.get_tensor_basis()[0]);

      VectorType src, dst;

      laplace_operator.initialize_dof_vector(src);
      laplace_operator.initialize_dof_vector(dst);

      {
        LinearAlgebra::distributed::Vector<Number> src_host(
          src.get_partitioner());

        VectorTools::create_right_hand_side<dim, dim>(
          mapping,
          dof_handler,
          quadrature,
          AnalyticalFunction<dim, Number>(n_components),
          src_host,
          constraints);

        LinearAlgebra::ReadWriteVector<Number> rw_vector(
          src.get_partitioner()->locally_owned_range());
        rw_vector.import(src_host, VectorOperation::insert);
        src.import(rw_vector, VectorOperation::insert);

        dst = 0.0;
      }

      laplace_operator.vmult(dst, src);

      table.add_value("fe_degree", degree);
      table.add_value("n_refinements", n_refinements);
      table.add_value("n_components", n_components);
      table.add_value("n_dofs", dof_handler.n_dofs());

      if (std::is_same_v<MemorySpace, dealii::MemorySpace::Host>)
        table.add_value("version", "host");
      else
        table.add_value("version", "default");

      table.add_value("norm", dst.l2_norm());
      table.set_scientific("norm", true);
    }
}

int
main()
{
  initlog();
  Kokkos::initialize();

  const unsigned int dim           = 2;
  const unsigned int fe_degree     = 3;
  unsigned int       n_refinements = 3;

  ConvergenceTable table;

  run<dim, fe_degree, 1, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, 1, MemorySpace::Default>(n_refinements, table);
  run<dim, fe_degree, dim, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim, MemorySpace::Default>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Default>(n_refinements, table);

  table.write_text(deallog.get_file_stream());

  Kokkos::finalize();
}
