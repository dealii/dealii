//---------------------------------------------------------------------------
//    $Id: 2d_coarse_grid_01.cc 17444 2008-10-31 19:35:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check existence of ghost layer in 2d with global refinement

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>
#include <base/utilities.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      for (int i=0;i<5;++i)
	{
	  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	    deallog << "refine loop:" << i << std::endl;

	  tr.refine_global(1);

	  if (myid==0)
	    {

	      std::vector< unsigned int > cell_subd;
	      cell_subd.resize(tr.n_active_cells());

	      GridTools::get_subdomain_association(tr, cell_subd);
	      for (unsigned int i=0;i<tr.n_active_cells();++i)
		deallog << cell_subd[i] << " ";
	      deallog << std::endl;
	    }

					   //check that all local
					   //neighbors have the
					   //correct level
	  typename Triangulation<dim,dim>::active_cell_iterator cell;

	  for (cell = tr.begin_active();
	       cell != tr.end();
	       ++cell)
	    {
	      if (cell->subdomain_id() != (unsigned int)myid)
		{
		  Assert (cell->is_ghost() || cell->is_artificial(),
			  ExcInternalError());
		  continue;
		}

	      for (unsigned int n=0;n<GeometryInfo<dim>::faces_per_cell;++n)
		{
		  if (cell->at_boundary(n))
		    continue;
		  Assert (cell->neighbor(n).state() == IteratorState::valid,
			  ExcInternalError());

		  Assert( cell->neighbor(n)->level() == cell->level(),
			  ExcInternalError());

		  Assert(!cell->neighbor(n)->has_children(), ExcInternalError() );
		}
	    }

	  write_vtk (tr, "p4est_2d_ghost_02", "1");

	  Assert (tr.n_global_active_cells() ==
		  static_cast<unsigned int>(std::pow (1.*GeometryInfo<dim>::max_children_per_cell, i+1)),
		  ExcInternalError());
	}
    }

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
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

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  std::cout << myid << ":" << getpid() << std::endl;
				   //system("sleep 20");


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_ghost_02").c_str());
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
