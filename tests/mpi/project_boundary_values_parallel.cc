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

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools_boundary.h>

#include <sstream>

#include "../tests.h"

template <int dim>
void
check()
{
  const FE_Q<dim>                                           fe(1);
  const Functions::ConstantFunction<dim>                    function(3.);
  const std::map<types::boundary_id, const Function<dim> *> functions = {
    {0, &function}};

  Triangulation<dim> serial;
  GridGenerator::hyper_cube(serial);
  DoFHandler<dim> serial_dofs(serial);
  serial_dofs.distribute_dofs(fe);
  std::map<types::global_dof_index, double> values;
  VectorTools::project_boundary_values(serial_dofs,
                                       functions,
                                       QGauss<dim - 1>(2),
                                       values);
  AssertThrow(!values.empty(), ExcInternalError());
  for (const auto &entry : values)
    AssertThrow(std::abs(entry.second - 3.) < 1e-11, ExcInternalError());

  parallel::distributed::Triangulation<dim> distributed(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(distributed);
  distributed.refine_global(3);
  DoFHandler<dim> distributed_dofs(distributed);
  distributed_dofs.distribute_dofs(fe);
  bool rejected = false;
  try
    {
      VectorTools::project_boundary_values(distributed_dofs,
                                           functions,
                                           QGauss<dim - 1>(2),
                                           values);
    }
  catch (const ExceptionBase &exception)
    {
      std::ostringstream diagnostic;
      exception.print_info(diagnostic);
      rejected = diagnostic.str().find(
                   "not implemented for distributed triangulations") !=
                 std::string::npos;
    }
  AssertThrow(rejected,
              ExcMessage("Distributed projection must be rejected explicitly"));
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << dim
            << "D: serial projection exact; distributed projection rejected"
            << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  deal_II_exceptions::disable_abort_on_exception();
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    initlog();
  check<2>();
  check<3>();
}
