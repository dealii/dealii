// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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


// The same as dof_tools_04 but for a parallel::distributed::parallel::distributed::Triangulation
// instead of a parallel::distributed::Triangulation:
// check
//   DoFTools::extract_hanging_node_constraints
// uses a slightly different refinement and less different FiniteElements

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <string>
#include <algorithm>



template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  types::global_dof_index n_dofs = dof_handler.n_dofs();

  std::vector<bool> is_hanging_node_constrained (n_dofs);
  DoFTools::extract_hanging_node_dofs (dof_handler,
                                       is_hanging_node_constrained);

  MappingQGeneric<dim> mapping(1);
  std::map< types::global_dof_index, Point< dim > > support_point_of_dof;
  DoFTools::map_dofs_to_support_points (mapping, dof_handler, support_point_of_dof);

  std::vector<Point<dim>> constrained_points;

  for (const auto &pair : support_point_of_dof)
    if (is_hanging_node_constrained[pair.first])
      constrained_points.push_back(pair.second);

  const int root = 0;
  const int mylen = constrained_points.size()*dim;
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::vector<int> recvcounts (n_processes);

  MPI_Gather(&mylen, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  int totlen = 0;
  std::vector<int> displs;
  std::vector<Point<dim>> all_points;

  if (my_id==0)
    {
      displs.resize(n_processes);

      displs[0] = 0;
      totlen += recvcounts[0];

      for (unsigned int i=1; i<n_processes; i++)
        {
          totlen += recvcounts[i];
          displs[i] = displs[i-1] + recvcounts[i-1];
        }

      all_points.resize(totlen/dim);
    }

  MPI_Gatherv(constrained_points.data(), mylen, MPI_DOUBLE,
              all_points.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  std::sort(all_points.begin(), all_points.end(),
            [](const Point<dim> &a, const Point<dim> &b)
  {
    for (unsigned int i=0; i<dim; ++i)
      {
        if (a(i) < b(i))
          return true;
        if (a(i) > b(i))
          return false;
      }
    return false;
  });
  all_points.erase(std::unique(all_points.begin(), all_points.end()), all_points.end());

  deallog << all_points.size() << std::endl;
  for (const auto &point: all_points)
    deallog << point << std::endl;
  deallog << std::endl;
}



template <int dim>
void
check (const FiniteElement<dim> &fe,
       const std::string        &name)
{
  deallog << "Checking " << name
          << " in " << dim << "d:"
          << std::endl;

  // create tria and dofhandler
  // objects. set different boundary
  // and sub-domain ids
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global (1);
  for (unsigned int ref=0; ref<2; ++ref)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned() &&
            cell->center()(0)<.5 && cell->center()(1)<.5)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  check_this (dof_handler);
}

#define CHECK(EL,deg,dim)\
  { FE_ ## EL<dim> EL(deg);   \
    check(EL, #EL #deg); }

#define CHECK_ALL(EL,deg)\
  { CHECK(EL,deg,2); \
    CHECK(EL,deg,3); \
  }

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
  CHECK_ALL(Q,1);
  CHECK_ALL(Q,2);

  return 0;
}
