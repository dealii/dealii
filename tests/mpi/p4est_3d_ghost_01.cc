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


// Reproduce a bug in the ghostlayer construction for a simple
// 3d mesh.

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
#include <deal.II/base/utilities.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      unsigned int rep=2;
      unsigned int ref=2;
      std::vector<unsigned int> repetitions;
      repetitions.push_back(rep);
      repetitions.push_back(rep);
      repetitions.push_back(rep);

      Point<3> p(0,0,0), q(1,1,1);
      GridGenerator::subdivided_hyper_rectangle(tr, repetitions, p, q, false);

      tr.refine_global(ref);


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
	    }
	  else
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

      const unsigned int checksum = tr.get_checksum ();
      if (myid==0)
	{
	  deallog << "Checksum: "
		  << checksum
		  << std::endl;

	  std::ofstream file("p4est_3d_ghost_01/1.gpl");
	  GridOut().write_gnuplot (tr, file);
	}

      deallog << "#global cells " << tr.n_global_active_cells() << std::endl;
    }

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

  std::cout << myid << ":" << getpid() << std::endl;
				   //system("sleep 20");


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_3d_ghost_01").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
