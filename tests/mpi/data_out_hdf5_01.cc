// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

// tests DataOut with HDF5

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>


std::string output_file_name = "output";


template <int dim>
void
test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global (1);
  //  tria.begin_active()->set_refine_flag();
  //tria.execute_coarsening_and_refinement ();

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs());
  for (unsigned int i=0; i<v1.size(); ++i) v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector (dof1, v1, "linear");
  data_out.build_patches (2);

  data_out.write_hdf5_parallel ("out.h5", MPI_COMM_WORLD);
  std::vector<XDMFEntry> xdmf_entries;
  xdmf_entries.push_back(
    data_out.create_xdmf_entry("out.h5", 0, MPI_COMM_WORLD));

  data_out.write_xdmf_file(xdmf_entries, "out.xdmf",
			   MPI_COMM_WORLD);

  deallog <<"ok" << std::endl;
  
  
  // Sadly hdf5 is binary and we can not use hd5dump because it might
  // not be in the path. At least we can look at the xdmf:
  if (0==Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
    cat_file("out.xdmf");  
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  try
    {
      test<2>();
//      test<3>();

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

