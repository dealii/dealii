// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double> for a distributed mesh,
// also compare with/without sparse line_cache via IndexSet.
// Refine one corner.

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  for (int i = 0; i < 12; ++i)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "step " << i << std::endl;


      tr.begin_active()->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      DoFHandler<dim> dofh(tr);

      static const FE_Q<dim> fe(1);
      dofh.distribute_dofs(fe);

      const IndexSet dof_set = DoFTools::extract_locally_relevant_dofs(dofh);

      AffineConstraints<double> cm;
      DoFTools::make_hanging_node_constraints(dofh, cm);
      AffineConstraints<double> cm2(dofh.locally_owned_dofs(), dof_set);
      DoFTools::make_hanging_node_constraints(dofh, cm2);

      if (myid == 0)
        {
          std::stringstream s;
          cm.print(s);
          deallog << s.str();

          deallog << "****" << std::endl;

          std::stringstream s2;
          cm2.print(s2);
          deallog << s2.str();

          if (s.str() == s2.str())
            deallog << "ok" << std::endl;
          else
            deallog << "not ok" << std::endl;
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

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
