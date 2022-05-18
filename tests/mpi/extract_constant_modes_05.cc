// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// test DoFTools::extract_constant_modes for multi-physical problems like in
// step-46 in parallel

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim>
void
test(const unsigned int fe_degree)
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  Point<dim>                p1, p2;
  std::vector<unsigned int> sub(dim);
  sub[0] = 2 * Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  for (unsigned int i = 0; i < dim; ++i)
    {
      if (i > 0)
        {
          sub[i] = 1;
        }
      p2[i] = 1.0;
    }

  GridGenerator::subdivided_hyper_rectangle(
    static_cast<Triangulation<dim> &>(triangulation), sub, p1, p2);

  DoFHandler<dim> dof_handler(triangulation);

  {
    // set cells to switch between FE_QxFE_Nothing and FE_NothingxFE_Q
    // alternately
    unsigned int last_index = 1;
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            cell->set_active_fe_index(last_index = 1 - last_index);
          }
      }
  }

  FESystem<dim>         fe_1(FE_Q<dim>(fe_degree), 1, FE_Nothing<dim>(), 1);
  FESystem<dim>         fe_2(FE_Nothing<dim>(), 1, FE_Q<dim>(fe_degree), 1);
  hp::FECollection<dim> fe_collection(fe_1, fe_2);

  dof_handler.distribute_dofs(fe_collection);
  DoFRenumbering::component_wise(dof_handler);

  if (myid == 0)
    {
      deallog << "FE degree=" << fe_degree << std::endl;
      deallog << "Total dofs=" << dof_handler.n_dofs() << std::endl;
    }



  // extract constant modes and print
  if (myid == 0)
    {
      for (unsigned int component = 0; component < fe_collection.n_components();
           ++component)
        {
          const FEValuesExtractors::Scalar component_extractor(component);
          const ComponentMask              component_mask(
            fe_collection.component_mask(component_extractor));

          std::vector<std::vector<bool>> constant_modes(
            fe_collection.n_components(),
            std::vector<bool>(dof_handler.n_locally_owned_dofs(), false));

          DoFTools::extract_constant_modes(dof_handler,
                                           component_mask,
                                           constant_modes);



          for (unsigned int i = 0; i < constant_modes.size(); ++i)
            {
              for (unsigned int j = 0; j < constant_modes[i].size(); ++j)
                {
                  deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
                }
              deallog << std::endl;
            }
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

      test<2>(1);
      test<2>(2);
      test<2>(3);
      test<2>(4);
    }
  else
    {
      test<2>(1);
      test<2>(2);
      test<2>(3);
      test<2>(4);
    }
}
