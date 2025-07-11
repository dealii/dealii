// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test VectorTools::compute_mean_value for parallel computations for
// hp-adaptive applications.
//
// we interpolate a linear function onto the grid with a symmetric mesh.
// the mean value of the interpolation must be the mean of the linear function

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class LinearFunction : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p[0] + 2;
  }
};


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  {
    GridGenerator::hyper_ball(tr);
    tr.refine_global(2);
  }

  hp::MappingCollection<dim> mapping_collection((MappingQ1<dim>()));

  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  q_collection;
  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      fe_collection.push_back(FE_Q<dim>(degree));
      q_collection.push_back(QGauss<dim>(degree + 1));
    }

  DoFHandler<dim> dofh(tr);
  {
    // set fe indices arbitrarily
    unsigned int i = 0;
    for (const auto &cell :
         dofh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      cell->set_active_fe_index(i++ % fe_collection.size());
    dofh.distribute_dofs(fe_collection);
  }

  const IndexSet relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set,
                                      dofh.get_mpi_communicator());
  {
    TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
                                               dofh.get_mpi_communicator());
    VectorTools::interpolate(dofh, LinearFunction<dim>(), interpolated);
    x_rel = interpolated;
  }

  const double mean = VectorTools::compute_mean_value(
    mapping_collection, dofh, q_collection, x_rel, 0);

  deallog << "mean=" << mean << std::endl;

  Assert(std::fabs(mean - 2) < 1e-3, ExcInternalError());
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
