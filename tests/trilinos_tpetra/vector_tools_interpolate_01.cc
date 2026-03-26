// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test distributed VectorTools::interpolate with Tpetra vectors

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(){};

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double f = p[0] * 2.0 + 1.0;
    if (dim > 1)
      f *= p[1] * 3.3 - 1.0;
    if (dim > 2)
      f *= p[2] * 5.0;
    return f;
  };
};

template <int dim, typename VectorType>
void
test()
{
  MyFunction<dim>                           func;
  MappingQ<dim>                             mapping(1);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim>        dof_handler(tr);
  static const FE_Q<dim> fe(1);
  dof_handler.distribute_dofs(fe);

  const auto locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  VectorType solution_distributed(dof_handler.locally_owned_dofs(),
                                  MPI_COMM_WORLD);
  VectorTools::interpolate(mapping, dof_handler, func, solution_distributed);

  solution_distributed.compress(VectorOperation::insert);

  VectorType solution_ghosted(dof_handler.locally_owned_dofs(),
                              locally_relevant_set,
                              MPI_COMM_WORLD);
  solution_ghosted = solution_distributed;

  deallog << "norm: " << solution_distributed.l2_norm() << std::endl;
  dealii::Vector<double> difference(tr.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution_ghosted,
                                    func,
                                    difference,
                                    QGauss<dim>(2),
                                    VectorTools::L2_norm);

  deallog << "error: " << difference.l2_norm() << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

#ifdef DEAL_II_TRILINOS_WITH_EPETRA
  deallog.push("2d:TrilinosWrappers");
  test<2, TrilinosWrappers::MPI::Vector>();
  deallog.pop();

  deallog.push("3d:TrilinosWrappers");
  test<3, TrilinosWrappers::MPI::Vector>();
  deallog.pop();
#endif

  deallog.push("2d:TpetraWrappers::Host");
  test<2, LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host>>();
  deallog.pop();

  deallog.push("3d:TpetraWrappers::Host");
  test<3, LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host>>();
  deallog.pop();

  deallog.push("2d:TpetraWrappers::Default");
  test<2,
       LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>();
  deallog.pop();

  deallog.push("3d:TpetraWrappers::Default");
  test<3,
       LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>();
  deallog.pop();
}
