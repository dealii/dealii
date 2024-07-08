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

// Check that ScratchData returns the correct jumps in, and averages of,
// solution values, gradients, etc. at interfaces.
// - Vector valued FE
// - hp variant

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
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

template <int dim>
class DiscontinuousFunction : public Function<dim>
{
public:
  DiscontinuousFunction(const unsigned int n_components = 1)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    const double value = (p[0] * p[1] > 0 ? +3 : -1); // Non-trivial average
    return value * (1.0 + std::sin(p[0]) * std::cos(p[1]));
  }
};


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
    FESystem<dim, spacedim>(FE_DGQ<dim, spacedim>(3), dim),
    FESystem<dim, spacedim>(FE_DGQ<dim, spacedim>(4), dim));
  hp::QCollection<dim>     qf_cell(QGauss<dim>(fe.max_degree()),
                               QGauss<dim>(fe.max_degree() + 1));
  hp::QCollection<dim - 1> qf_face(QGauss<dim - 1>(fe.max_degree()),
                                   QGauss<dim - 1>(fe.max_degree() + 1));

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(1);

  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  // Don't interpolate, as points precisely on the interfaces get evaluated and
  // we don't end up with a jump in the values across it.
  AffineConstraints<double> constraints;
  constraints.close();
  VectorTools::project(
    dof_handler,
    constraints,
    hp::QCollection<dim>(QGauss<spacedim>(fe.max_degree() + 2),
                         QGauss<spacedim>(fe.max_degree() + 3)),
    DiscontinuousFunction<spacedim>(fe.n_components()),
    solution);

  const UpdateFlags update_flags =
    update_values | update_gradients | update_hessians | update_3rd_derivatives;
  MeshWorker::ScratchData<dim, spacedim> scratch_data(
    fe, qf_cell, update_flags, qf_face, update_flags);

  const auto   cell = dof_handler.begin_active();
  unsigned int face = 0;
  while (cell->face(face)->at_boundary())
    ++face;

  scratch_data.reinit(cell,
                      face,
                      cell->neighbor(face),
                      cell->neighbor_face_no(face));
  scratch_data.extract_local_dof_values("solution", solution);

  deallog << "Jumps in values: "
          << scratch_data.get_jumps_in_values("solution", extractor)[0]
          << std::endl;
  deallog << "Jumps in gradients: "
          << scratch_data.get_jumps_in_gradients("solution", extractor)[0]
          << std::endl;
  deallog << "Jumps in Hessians: "
          << scratch_data.get_jumps_in_hessians("solution", extractor)[0]
          << std::endl;
  deallog << "Jumps in third derivatives: "
          << scratch_data.get_jumps_in_third_derivatives("solution",
                                                         extractor)[0]
          << std::endl;

  deallog << "Averages of values: "
          << scratch_data.get_averages_of_values("solution", extractor)[0]
          << std::endl;
  deallog << "Averages of gradients: "
          << scratch_data.get_averages_of_gradients("solution", extractor)[0]
          << std::endl;
  deallog << "Averages of Hessians: "
          << scratch_data.get_averages_of_hessians("solution", extractor)[0]
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

  run<2>(extractor);
  run<3>(extractor);

  deallog << "OK" << std::endl;
}
