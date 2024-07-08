/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Check that ScratchData returns the correct solution values, gradients, etc.
// - Tensor valued FE
// - hp variant

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim,
          int spacedim        = dim,
          typename NumberType = double,
          typename ExtractorType>
void
run(const ExtractorType &extractor)
{
  LogStream::Prefix prefix("Dim " + Utilities::to_string(dim));
  std::cout << "Dim: " << dim << std::endl;

  hp::FECollection<dim, spacedim> fe(
    FESystem<dim, spacedim>(FE_Q<dim, spacedim>(1),
                            Tensor<2, dim>::n_independent_components),
    FESystem<dim, spacedim>(FE_Q<dim, spacedim>(3),
                            Tensor<2, dim>::n_independent_components));
  hp::QCollection<dim>     qf_cell(QGauss<dim>(fe[0].degree + 1),
                               QGauss<dim>(fe[1].degree + 1));
  hp::QCollection<dim - 1> qf_face(QGauss<dim - 1>(fe[0].degree + 1),
                                   QGauss<dim - 1>(fe[1].degree + 1));

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.begin_active()->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
                           Functions::CosineFunction<spacedim>(
                             fe.n_components()),
                           solution);

  const UpdateFlags update_flags = update_values | update_gradients;
  MeshWorker::ScratchData<dim, spacedim> scratch_data(fe,
                                                      qf_cell,
                                                      update_flags);

  const auto cell = dof_handler.begin_active();
  scratch_data.reinit(cell);
  scratch_data.extract_local_dof_values("solution", solution);

  deallog << "Value: " << scratch_data.get_values("solution", extractor)[0]
          << std::endl;
  deallog << "Gradient: "
          << scratch_data.get_gradients("solution", extractor)[0] << std::endl;
  deallog << "Divergence: "
          << scratch_data.get_divergences("solution", extractor)[0]
          << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  const FEValuesExtractors::Tensor<2> extractor(0);

  run<2>(extractor);
  run<3>(extractor);

  deallog << "OK" << std::endl;
}
