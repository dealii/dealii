// ---------------------------------------------------------------------
// $Id: 3d_refinement_01.cc 31349 2013-10-20 19:07:06Z maier $
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


// create a split tria mesh

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/split_tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

template <int dim, int spacedim>
void write_mesh (const Triangulation<dim,spacedim> &tria,
                 const char                                *filename_)
{
  DataOut<dim> data_out;
  data_out.attach_triangulation (tria);
  Vector<float> subdomain (tria.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = tria.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");

  data_out.build_patches ();
  const std::string filename = (filename_ +
                                Utilities::int_to_string
                                (tria.locally_owned_subdomain(), 4));
  {
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
  }
}



template<int dim>
void test()
{
  Triangulation<dim> basetria;
  GridGenerator::hyper_L(basetria);
  basetria.refine_global();

  typename Triangulation<dim>::active_cell_iterator it = basetria.begin_active();
  it->set_refine_flag();
  //  ++it;
  it->set_refine_flag();
  basetria.execute_coarsening_and_refinement ();

  typename Triangulation<dim>::active_cell_iterator
  cell = basetria.begin_active();

  for (unsigned int i=0; i < basetria.n_active_cells(); ++i)
    {
      cell->set_subdomain_id((i*3/basetria.n_active_cells()) % 3);
      ++cell;
    }



  parallel::split::Triangulation<dim> tr(MPI_COMM_WORLD);
  tr.copy_triangulation(basetria);

  deallog
      << " locally_owned_subdomain(): " << tr.locally_owned_subdomain() << "\n"
      << " n_active_cells: " << tr.n_active_cells() << "\n"
      << " n_levels: " << tr.n_levels() << "\n"
      << " n_global_levels: " << tr.n_global_levels()  << "\n"
      << " n_locally_owned_active_cells: " << tr.n_locally_owned_active_cells() << "\n"
      << " n_global_active_cells: " << tr.n_global_active_cells() << "\n"
      << std::endl;

  /*deallog << "n_locally_owned_active_cells_per_processor: ";
  std::vector<unsigned int> v = tr.n_locally_owned_active_cells_per_processor();
  for (unsigned int i=0;i<v.size();++i)
    deallog << v[i] << " ";
    deallog << std::endl;*/

  {
    deallog << "subdomains: ";
    typename  Triangulation<dim>::active_cell_iterator it=tr.begin_active();
    for (unsigned int index=0; it!=tr.end(); ++it,++index)
      {
        deallog << it->subdomain_id() << " ";
      }
    deallog << std::endl;
  }

  GridOut go;
  go.write_mesh_per_processor_as_vtu (tr, "mesh", false, true);

  //write_mesh(tr, "mesh");
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  //  test<3>();
  deallog.pop();
}
