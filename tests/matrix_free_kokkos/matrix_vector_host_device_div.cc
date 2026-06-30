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



// Test Portable::FEEvaluation::get_divergence() and submit_divergence() by
// applying the div-div operator (div v, div u) on both the host (FEEvaluation)
// and the device (Portable::FEEvaluation) and comparing the two results. This
// is the portable counterpart of the host test
// matrix_free/matrix_vector_div.cc.

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int fe_degree, typename Number, typename MemorySpace>
class DivOperator;



// Host implementation: uses the regular FEEvaluation together with its
// get_divergence()/submit_divergence() functions.
template <int dim, int fe_degree, typename Number>
class DivOperator<dim, fe_degree, Number, MemorySpace::Host>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature)
  {
    typename MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_JxW_values;

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
    matrix_free.cell_loop(&DivOperator::local_apply, this, dst, src, true);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> phi(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values_plain(src);
        phi.evaluate(EvaluationFlags::gradients);

        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_divergence(phi.get_divergence(q), q);

        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  MatrixFree<dim, Number> matrix_free;
};



// Device version: identical but uses Portable::FEEvaluation
template <int dim, int fe_degree, typename Number>
class DivOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> phi(
      data);
    phi.read_dof_values(src);
    phi.evaluate(EvaluationFlags::gradients);

    data->for_each_quad_point([&](const int &q_point) {
      phi.submit_divergence(phi.get_divergence(q_point), q_point);
    });

    phi.integrate(EvaluationFlags::gradients);
    phi.distribute_local_to_global(dst);
  }

  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);
};



// Device implementation: uses Portable::MatrixFree and the functor above.
template <int dim, int fe_degree, typename Number>
class DivOperator<dim, fe_degree, Number, MemorySpace::Default>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

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
    dst = 0.0;
    DivOperatorLocal<dim, fe_degree, Number> local_operator;
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst);
  }

private:
  Portable::MatrixFree<dim, Number> matrix_free;
};



template <int dim, typename T>
class AnalyticalFunction : public Function<dim, T>
{
public:
  AnalyticalFunction()
    : Function<dim, T>(dim)
  {}

  virtual T
  value(const Point<dim, T> &p, const unsigned int component = 0) const override
  {
    T temp = 0.0;
    for (unsigned int d = 0; d < dim; ++d)
      temp += std::sin(p[d]);

    return temp * (1.0 + component);
  }
};



template <int dim, int fe_degree>
void
test()
{
  using Number = double;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(dim == 2 ? 3 : 2);

  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe_q(fe_degree);
  const FESystem<dim>  fe(fe_q, dim);
  const QGauss<dim>    quadrature(fe_degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  LinearAlgebra::distributed::Vector<Number> src_host;
  {
    DivOperator<dim, fe_degree, Number, MemorySpace::Host> tmp_operator;
    tmp_operator.reinit(mapping,
                        dof_handler,
                        constraints,
                        quadrature.get_tensor_basis()[0]);
    tmp_operator.initialize_dof_vector(src_host);

    VectorTools::create_right_hand_side<dim, dim>(
      mapping,
      dof_handler,
      quadrature,
      AnalyticalFunction<dim, Number>(),
      src_host,
      constraints);
  }

  Number host_norm = 0.;
  {
    DivOperator<dim, fe_degree, Number, MemorySpace::Host> div_operator;
    div_operator.reinit(mapping,
                        dof_handler,
                        constraints,
                        quadrature.get_tensor_basis()[0]);

    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> dst;
    div_operator.initialize_dof_vector(dst);
    div_operator.vmult(dst, src_host);
    host_norm = dst.l2_norm();
  }

  Number device_norm = 0.;
  {
    DivOperator<dim, fe_degree, Number, MemorySpace::Default> div_operator;
    div_operator.reinit(mapping,
                        dof_handler,
                        constraints,
                        quadrature.get_tensor_basis()[0]);

    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src, dst;
    div_operator.initialize_dof_vector(src);
    div_operator.initialize_dof_vector(dst);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src_host, VectorOperation::insert);
    src.import_elements(rw_vector, VectorOperation::insert);

    div_operator.vmult(dst, src);
    device_norm = dst.l2_norm();
  }

  // Report both norms; the host and device results should agree.
  deallog << "dim=" << dim << " fe_degree=" << fe_degree
          << ": host norm = " << host_norm << ", device norm = " << device_norm
          << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(10);

  test<2, 1>();
  test<2, 2>();
  test<3, 1>();
  test<3, 2>();

  return 0;
}
