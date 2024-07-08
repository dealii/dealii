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
// - Vector valued FE

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

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

  const FESystem<dim, spacedim> fe(FE_Q<dim, spacedim>(3), dim);
  const QGauss<spacedim>        qf_cell(fe.degree + 1);
  const QGauss<spacedim - 1>    qf_face(fe.degree + 1);

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
                           Functions::CosineFunction<spacedim>(
                             fe.n_components()),
                           solution);

  const UpdateFlags update_flags =
    update_values | update_gradients | update_hessians | update_3rd_derivatives;
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
  deallog << "Symmetric gradient: "
          << scratch_data.get_symmetric_gradients("solution", extractor)[0]
          << std::endl;
  deallog << "Divergence: "
          << scratch_data.get_divergences("solution", extractor)[0]
          << std::endl;
  deallog << "Curl: " << scratch_data.get_curls("solution", extractor)[0]
          << std::endl;
  deallog << "Hessian: " << scratch_data.get_hessians("solution", extractor)[0]
          << std::endl;
  deallog << "Laplacian: "
          << scratch_data.get_laplacians("solution", extractor)[0] << std::endl;
  deallog << "Third derivative: "
          << scratch_data.get_third_derivatives("solution", extractor)[0]
          << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  const FEValuesExtractors::Vector extractor(0);

  // Just run in 3d, so we can extract the curl
  run<3>(extractor);

  deallog << "OK" << std::endl;
}
