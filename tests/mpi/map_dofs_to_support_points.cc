// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Test DoFTools::map_dofs_to_support_points for parallel DoFHandlers

#include "../tests.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>



template <int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tr);
  tr.reset_manifold(0);
  tr.refine_global (1);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  std::map<types::global_dof_index, Point<dim> > points;
  DoFTools::map_dofs_to_support_points (MappingQGeneric<dim>(1),
                                        dofh,
                                        points);
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (typename std::map<types::global_dof_index, Point<dim> >::const_iterator
           p = points.begin();
           p != points.end();
           ++p)
        deallog << p->first << " -> " << p->second
                << std::endl;
    }

  // the result of the call above is
  // supposed to be a map that
  // contains exactly the locally
  // relevant dofs, so test this
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  for (unsigned int i=0; i<dofh.n_dofs(); ++i)
    {
      if (relevant_set.is_element(i))
        {
          AssertThrow (points.find(i) != points.end(),
                       ExcInternalError());
        }
      else
        {
          AssertThrow (points.find(i) == points.end(),
                       ExcInternalError());
        }
    }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }

}
