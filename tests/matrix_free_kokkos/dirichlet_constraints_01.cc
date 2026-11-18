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



// Check that Portable::FEEvaluation::read_dof_values() treats DoFs with
// Dirichlet constraints as zero, like FEEvaluation::read_dof_values(): apply
// a Laplace operator to a source vector that is nonzero on the constrained
// DoFs, on a mesh with hanging nodes, with Portable::MatrixFree and with
// MatrixFree, and compare the results.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
      phi(data);
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



template <int dim, int fe_degree, int n_components>
void
test()
{
  using Number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned() && cell->center()[0] < 0.5 &&
        cell->center()[1] < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const FESystem<dim> fe(FE_Q<dim>{fe_degree}, n_components);
  DoFHandler<dim>     dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const MappingQ1<dim> mapping;
  const QGauss<1>      quadrature(fe_degree + 1);

  AffineConstraints<Number> constraints(dof_handler.locally_owned_dofs(),
                                        DoFTools::extract_locally_relevant_dofs(
                                          dof_handler));
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(mapping,
                                           dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim, Number>(
                                             n_components),
                                           constraints);
  constraints.close();

  // MatrixFree on the host
  MatrixFree<dim, Number>                          mf_host;
  typename MatrixFree<dim, Number>::AdditionalData data_host;
  data_host.mapping_update_flags = update_gradients;
  mf_host.reinit(mapping, dof_handler, constraints, quadrature, data_host);

  using HostVector =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;
  HostVector src_host, dst_host;
  mf_host.initialize_dof_vector(src_host);
  mf_host.initialize_dof_vector(dst_host);

  // nonzero on all DoFs, including the constrained ones
  for (unsigned int i = 0; i < src_host.locally_owned_size(); ++i)
    src_host.local_element(i) = random_value<Number>(-1., 1.);

  const std::function<void(const MatrixFree<dim, Number> &,
                           HostVector &,
                           const HostVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    local_apply = [](const MatrixFree<dim, Number>               &data,
                     HostVector                                  &dst,
                     const HostVector                            &src,
                     const std::pair<unsigned int, unsigned int> &range) {
      FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number> phi(
        data);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(EvaluationFlags::gradients);
          for (const unsigned int q : phi.quadrature_point_indices())
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(EvaluationFlags::gradients);
          phi.distribute_local_to_global(dst);
        }
    };
  mf_host.cell_loop(local_apply, dst_host, src_host, true);
  for (const auto i : dof_handler.locally_owned_dofs())
    if (constraints.is_constrained(i))
      dst_host(i) = src_host(i);

  // Portable::MatrixFree
  Portable::MatrixFree<dim, Number>                          mf_device;
  typename Portable::MatrixFree<dim, Number>::AdditionalData data_device;
  data_device.mapping_update_flags = update_JxW_values | update_gradients;
  mf_device.reinit(mapping, dof_handler, constraints, quadrature, data_device);

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src_device,
    dst_device;
  mf_device.initialize_dof_vector(src_device);
  mf_device.initialize_dof_vector(dst_device);

  LinearAlgebra::ReadWriteVector<Number> rw(dof_handler.locally_owned_dofs());
  rw.import_elements(src_host, VectorOperation::insert);
  src_device.import_elements(rw, VectorOperation::insert);

  mf_device.cell_loop(
    LaplaceOperatorLocal<dim, fe_degree, n_components, Number>(),
    src_device,
    dst_device);
  mf_device.copy_constrained_values(src_device, dst_device);

  rw.import_elements(dst_device, VectorOperation::insert);
  for (const auto i : dof_handler.locally_owned_dofs())
    dst_host(i) -= rw(i);

  const double error = dst_host.linfty_norm();
  deallog << "dim=" << dim << " degree=" << fe_degree
          << " n_components=" << n_components << ": "
          << (error < 1e-12 ? std::string("OK") :
                              "error " + std::to_string(error))
          << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  MPILogInitAll                    log;

  test<2, 1, 1>();
  test<2, 2, 1>();
  test<2, 3, 2>();
  test<3, 1, 1>();
  test<3, 2, 2>();
}
