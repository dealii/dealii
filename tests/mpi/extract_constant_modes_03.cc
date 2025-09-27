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



// test DoFTools::extract_constant_modes with FE_Nothing

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test(unsigned int fe_nothing_index)
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  std::vector<unsigned int> sub(2);
  sub[0] = 2 * Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  sub[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(tr,
                                            sub,
                                            Point<2>(0, 0),
                                            Point<2>(1, 1));

  DoFHandler<dim> dofh(tr);

  {
    // set cells to use FE_Q and FE_Nothing alternately
    unsigned int last_index = 1;
    for (const auto &cell :
         dofh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      cell->set_active_fe_index(last_index = 1 - last_index);
  }

  if (fe_nothing_index == 1)
    dofh.distribute_dofs(
      hp::FECollection<dim>(FE_Q<dim>(1), FE_Nothing<dim>(1)));
  else
    dofh.distribute_dofs(
      hp::FECollection<dim>(FE_Nothing<dim>(1), FE_Q<dim>(1)));

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Total dofs=" << dofh.n_dofs() << std::endl;

  // extract constant modes and print

  if (myid == 0)
    {
      ComponentMask mask(1, true);

      const auto constant_modes = DoFTools::extract_constant_modes(dofh, mask);

      for (unsigned int i = 0; i < constant_modes.size(); ++i)
        {
          for (unsigned int j = 0; j < constant_modes[i].size(); ++j)
            deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
          deallog << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));


  if (myid == 0)
    {
      initlog();

      deallog.push("FE_Nothing second");
      test<2>(1);
      deallog.pop();

      deallog.push("FE_Nothing first");
      test<2>(0);
      deallog.pop();
    }
  else
    {
      deallog.push("FE_Nothing second");
      test<2>(1);
      deallog.pop();

      deallog.push("FE_Nothing first");
      test<2>(0);
      deallog.pop();
    }
}
