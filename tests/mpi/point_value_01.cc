// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

//
// Checks whether VectorTools::point_value shows the expected behavior on
// distributed meshes. If the point is locally available everything should be
// fine, else the exception ExcPointNotAvailableHere should be thrown.
//

#include "../tests.h"
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/parallel_vector.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  
  FE_Q<dim> fe(1);
  
  DoFHandler<dim> dof_handler (triangulation);
  
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (1);
  
  dof_handler.distribute_dofs (fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  parallel::distributed::Vector<double> locally_owned_solution
    (locally_owned_dofs, MPI_COMM_WORLD);

  locally_owned_solution=1;

  std::vector<Point<dim> > points;

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      if (dim==3)
        for (int k=0; k<2; k++)
          points.push_back(Point<dim>(.25+.5*i,.25+.5*j,.25+.5*k));
      else
        points.push_back(Point<dim>(.25+.5*i,.25+.5*j));

  typename std::vector<Point<dim> >::iterator point_iterator, points_end;
  point_iterator=points.begin();
  points_end=points.end();

  Vector<double> value(1);

  for(; point_iterator!=points_end; point_iterator++)
  {
    try
    {
      VectorTools::point_value (dof_handler, locally_owned_solution,
                                *point_iterator, value);
      if (std::abs(value[0]-1.)>1e-8)
        ExcInternalError();
    }
    catch (const VectorTools::ExcPointNotAvailableHere &)
    {}

    MPI_Barrier(MPI_COMM_WORLD);
  
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << *point_iterator << " OK" << std::endl;
  }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
  {
    std::ofstream logfile("output");
    deallog.attach(logfile, false);
    deallog.depth_console(0);
    deallog.threshold_double(1.e-10);
    
    deallog.push("2d");
    test<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
  else
  {
    test<2>();
    test<3>();
  }

  return 0;
}

