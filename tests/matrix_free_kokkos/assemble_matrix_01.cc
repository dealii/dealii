// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test Portable::FEEvaluation for assembling the system matrix of elasticity
// problems

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include "../tests.h"



template <int dim, int fe_degree>
class LocalDeviceOperator
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, double>::Data *data,
             const Portable::DeviceVector<double>                   &src,
             Portable::DeviceVector<double>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim, double> fe_eval(
      data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::gradients);

    data->for_each_quad_point([&](const int &q_point) {
      fe_eval.submit_symmetric_gradient(fe_eval.get_symmetric_gradient(q_point),
                                        q_point);
    });

    fe_eval.integrate(EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
  }

  static constexpr unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);
};



template <int dim, int fe_degree>
class DeviceOperator : public EnableObserverPointer
{
public:
  DeviceOperator(std::shared_ptr<const Portable::MatrixFree<dim, double>> mf)
    : mf_data(mf)
  {}

  void
  vmult(LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &dst,
        const LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
          &src) const
  {
    dst = 0.;
    LocalDeviceOperator<dim, fe_degree> test_operator;
    mf_data->cell_loop(test_operator, src, dst);
  }

  void
  initialize_dof_vector(
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &vec) const
  {
    mf_data->initialize_dof_vector(vec);
  }

private:
  std::shared_ptr<const Portable::MatrixFree<dim, double>> mf_data;
};



template <int dim, int fe_degree>
class HostOperator
  : public MatrixFreeOperators::
      Base<dim, LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>
{
public:
  HostOperator()
    : MatrixFreeOperators::Base<
        dim,
        LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>()
  {}

  virtual void
  compute_diagonal() override
  {}

private:
  void
  local_apply(
    const MatrixFree<dim, double>                                       &data,
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>       &dst,
    const LinearAlgebra::distributed::Vector<double, MemorySpace::Host> &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, double> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (const unsigned int q_point : phi.quadrature_point_indices())
          phi.submit_symmetric_gradient(phi.get_symmetric_gradient(q_point),
                                        q_point);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  virtual void
  apply_add(
    LinearAlgebra::distributed::Vector<double>       &dst,
    const LinearAlgebra::distributed::Vector<double> &src) const override
  {
    this->data->cell_loop(&HostOperator::local_apply, this, dst, src);
  }
};



template <int dim, int fe_degree>
void
do_test(const DoFHandler<dim> &dof_handler)
{
  deallog << "Testing " << dof_handler.get_fe().get_name() << std::endl;

  MappingQ1<dim>            mapping;
  AffineConstraints<double> constraints;
  QGauss<1>                 quad(fe_degree + 1);

  // Initialize the host matrix
  typename MatrixFree<dim, double>::AdditionalData additional_data_host;
  additional_data_host.mapping_update_flags = update_gradients;
  MatrixFree<dim, double> mf_host;
  mf_host.reinit(mapping, dof_handler, constraints, quad, additional_data_host);

  HostOperator<dim, fe_degree> matrix_host;
  matrix_host.initialize(
    std::make_shared<const MatrixFree<dim, double>>(mf_host));

  // Initialize the device matrix
  typename Portable::MatrixFree<dim, double>::AdditionalData
    additional_data_dev;
  additional_data_dev.mapping_update_flags = update_gradients;
  Portable::MatrixFree<dim, double> mf_dev;
  mf_dev.reinit(mapping, dof_handler, constraints, quad, additional_data_dev);

  DeviceOperator<dim, fe_degree> matrix_dev(
    std::make_shared<const Portable::MatrixFree<dim, double>>(mf_dev));

  // Set up the src and dst vectors
  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> src_host,
    dst_host;
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> src_dev,
    dst_dev;

  matrix_host.initialize_dof_vector(src_host);
  dst_host.reinit(src_host);

  matrix_dev.initialize_dof_vector(src_dev);
  dst_dev.reinit(src_dev);

  for (unsigned int i = 0; i < src_host.size(); ++i)
    src_host[i] = static_cast<double>(Testing::rand()) / RAND_MAX;
  src_dev.import_elements(src_host, VectorOperation::insert);

  // Do the matrix-vector-multiplication
  matrix_host.vmult(dst_host, src_host);
  matrix_dev.vmult(dst_dev, src_dev);

  // Compute the difference between the host and device results
  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> diff;
  diff.reinit(dst_host);
  diff.import_elements(dst_dev, VectorOperation::insert);
  diff.add(-1.0, dst_host);
  for (unsigned int i = 0; i < diff.size(); ++i)
    deallog << diff[i] << " ";
  deallog << std::endl;
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  FE_Q<dim>       feq(fe_degree);
  FESystem<dim>   fe(feq, dim);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  do_test<dim, fe_degree>(dof);
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Kokkos::initialize();

  deallog << std::fixed << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 2>();
    deallog.pop();
    deallog.push("3d");
    test<3, 2>();
    deallog.pop();
  }

  Kokkos::finalize();

  return 0;
}
