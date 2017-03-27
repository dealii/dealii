// ---------------------------------------------------------------------
// $Id: dof_handler_number_cache.cc 31761 2013-11-22 14:42:37Z heister $
//
// Copyright (C) 2008 - 2015 by the deal.II authors
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



// check number cache for split_tria

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/split_tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <cstdlib>
#include <numeric>


template<int dim>
void test()
{
  const unsigned int nproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  Triangulation<dim> basetria;
  {
    GridGenerator::hyper_L(basetria);
    basetria.refine_global (1);
    typename Triangulation<dim>::active_cell_iterator
    cell = basetria.begin_active();

    for (unsigned int i=0; i < basetria.n_active_cells(); ++i)
      {
        cell->set_subdomain_id((i*nproc/basetria.n_active_cells()) % nproc);
        ++cell;
      }
  }


  /*
  const unsigned int n_refinements[] = { 0, 4, 3, 2 };
  for (unsigned int i=0; i<n_refinements[dim]; ++i)
  {
    // refine one-fifth of cells randomly
    std::vector<bool> flags (triangulation.n_active_cells(), false);
    for (unsigned int k=0; k<flags.size()/5 + 1; ++k)
      flags[Testing::rand() % flags.size()] = true;
    // make sure there's at least one that
    // will be refined
    flags[0] = true;

    // refine triangulation
    unsigned int index=0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        if (flags[index])
          cell->set_refine_flag();
        ++index;
      }

    Assert (index <= triangulation.n_active_cells(), ExcInternalError());

    // flag all other cells for coarsening
    // (this should ensure that at least
    // some of them will actually be
    // coarsened)
    index=0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        if (!flags[index])
          cell->set_coarsen_flag();
        ++index;
      }

    triangulation.execute_coarsening_and_refinement ();*/

  parallel::split::Triangulation<dim> tr(MPI_COMM_WORLD);
  tr.copy_triangulation(basetria);

  GridOut go;
  go.write_mesh_per_processor_as_vtu (tr, "mesh", false, true);


  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tr);

  dof_handler.distribute_dofs (fe);

  int n_coarse = tr.n_active_cells()>0 ? tr.n_cells(0) : 0;
  deallog
      << " locally_owned_subdomain(): " << tr.locally_owned_subdomain() << "\n"
      << " n_active_cells: " << tr.n_active_cells() << "\n"
      << " n_levels: " << tr.n_levels() << "\n"
      << " n_global_levels: " << tr.n_global_levels()  << "\n"
      << " n_locally_owned_active_cells: " << tr.n_locally_owned_active_cells() << "\n"
      << " n_global_active_cells: " << tr.n_global_active_cells() << "\n"
      << " n_coarse_cells: " << n_coarse << "\n"
      << " n_global_coarse_cells: " << tr.n_global_coarse_cells()
      << std::endl;

  deallog << "available_coarse_cells: ";
  tr.available_coarse_cells().print(deallog);

  deallog
      << "n_dofs: " << dof_handler.n_dofs() << std::endl
      << "n_locally_owned_dofs: " << dof_handler.n_locally_owned_dofs() << std::endl;

  deallog << "n_locally_owned_dofs_per_processor: ";
  std::vector<types::global_dof_index> v = dof_handler.n_locally_owned_dofs_per_processor();
  unsigned int sum = 0;
  for (unsigned int i=0; i<v.size(); ++i)
    {
      deallog << v[i] << " ";
      sum += v[i];
    }
  deallog << " sum: " << sum << std::endl;

  //  Assert(dof_handler.n_locally_owned_dofs() == dof_handler.n_locally_owned_dofs_per_processor()[tr.locally_owned_subdomain()], ExcInternalError());
  //  Assert( dof_handler.n_locally_owned_dofs() == dof_handler.locally_owned_dofs().n_elements(), ExcInternalError());

  const unsigned int N = dof_handler.n_dofs();

  Assert (dof_handler.n_locally_owned_dofs() <= N,
          ExcInternalError());
  /*Assert (std::accumulate (dof_handler.n_locally_owned_dofs_per_processor().begin(),
                           dof_handler.n_locally_owned_dofs_per_processor().end(),
                           0U) == N,
          ExcInternalError());
  */
  IndexSet all (N);
  for (unsigned int i=0;
       i<dof_handler.locally_owned_dofs_per_processor().size(); ++i)
    {
      IndexSet intersect = all & dof_handler.locally_owned_dofs_per_processor()[i];
      //  Assert(intersect.n_elements()==0, ExcInternalError());
      all.add_indices(dof_handler.locally_owned_dofs_per_processor()[i]);
    }

  //Assert(all == complete_index_set(N), ExcInternalError());

  std::ofstream out(Utilities::int_to_string(tr.locally_owned_subdomain())
                    + ".gpl");

  out << "plot '-' using 1:2 with lines, '-' with labels point pt 2 offset 1,1" << std::endl;
  GridOut().write_gnuplot (tr, out);

  out << "e" << std::endl;

  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
                                        dof_handler,
                                        support_points);
  DoFTools::write_gnuplot_dof_support_point_info(out,
                                                 support_points);

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  //test<3>();
  deallog.pop();
}
