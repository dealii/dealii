//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test interaction with p4est with a few simple coarse grids in 2d

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/utilities.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
      GridIn<dim> gi;
      gi.attach_triangulation (tr);
      std::ifstream in ("../deal.II/grid_in_02/2d.xda");
      try
	{
	  gi.read_xda (in);
	}
      catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
	{
					   // ignore distorted cells
	  deallog << distorted_cells.distorted_cells.size()
		  << " distorted cells after creating mesh."
		  << std::endl;
	}

      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "subdomainid = "
		<< tr.begin_active()->subdomain_id()
		<< std::endl;

//      std::vector< unsigned int > cell_subd;
//      cell_subd.resize(tr.n_active_cells());

//      GridTools::get_subdomain_association(tr, cell_subd);
//       for (unsigned int i=0;i<tr.n_active_cells();++i)
// 	deallog << cell_subd[i] << " ";
//       deallog << std::endl;

      if (myid == 0)
	{
	  deallog << "#cells = " << tr.n_global_active_cells() << std::endl;

	  Assert(tr.n_global_active_cells() == tr.n_active_cells(),
		 ExcInternalError() );
	}

      const unsigned int checksum = tr.get_checksum ();
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "Checksum: "
		<< checksum
		<< std::endl;
    }


  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_coarse_01").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
