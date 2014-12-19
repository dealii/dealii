// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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



// Test Tria<2,3> and DataOutput.

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <deal.II/fe/fe_q.h>

template <int dim, int spacedim>
void write_vtk (const parallel::distributed::Triangulation<dim,spacedim> &tria,
                const char                                *filename)
{
  deallog << "Checksum: "
          << tria.get_checksum ()
          << std::endl;

  tria.write_mesh_vtk (filename);

  // copy the .pvtu and .vtu files
  // into the logstream
  int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid==0)
    {
      cat_file((std::string(filename) + ".pvtu").c_str());
      cat_file((std::string(filename) + "_0000.vtu").c_str());
    }
  else if (myid==1)
      cat_file((std::string(filename) + "_0001.vtu").c_str());
  else if (myid==2)
      cat_file((std::string(filename) + "_0002.vtu").c_str());
  else
    AssertThrow(false, ExcNotImplemented());
  
}

template<int dim, int spacedim>
void test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim,spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::torus(tr, 1, 0.2);
  tr.refine_global();

  tr.begin_active()->set_refine_flag();
  (++(tr.begin_active()))->set_refine_flag();

  tr.execute_coarsening_and_refinement ();

  write_vtk(tr, "file");

  FE_Q<dim,spacedim> fe(1);
  DoFHandler<dim,spacedim> dh(tr);
  dh.distribute_dofs(fe);
  deallog << "dofs " << dh.n_dofs() << std::endl;

  DataOut<dim,DoFHandler<dim,spacedim> > data_out;
    data_out.attach_triangulation (tr);
    Vector<float> subdomain (tr.n_active_cells());
   for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = tr.locally_owned_subdomain();
   data_out.add_data_vector (subdomain, "subdomain");

   std::string name = "f0.vtu";
   name[1] += tr.locally_owned_subdomain();

   {
     std::ofstream file(name.c_str());
     data_out.build_patches ();
     data_out.write_vtu (file);
   }
   cat_file(name.c_str());
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

  deallog.push("2-3");
  test<2,3>(deallog.get_file_stream());
  deallog.pop();
}
