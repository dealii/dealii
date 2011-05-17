//---------------------------------------------------------------------------
//    $Id: p4est_3d_refine_02.cc 23327 2011-02-11 03:19:07Z bangerth $
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


// refine a 3d shell once (currently a bug):
//[t500u:16597] [17] /scratch/p4est-0.3.1.55//DEBUG/lib/libsc.so.0(sc_abort_verbosef+0) [0x2afcc77b943b]
//[t500u:16597] [18] /scratch/p4est-0.3.1.55//DEBUG/lib/libp4est.so.0(p8est_quadrant_parent+0x5d) [0x2afcc7548d4b]
//[t500u:16597] [19] /scratch/p4est-0.3.1.55//DEBUG/lib/libp4est.so.0(+0x5046c) [0x2afcc754446c]
//[t500u:16597] [20] /scratch/p4est-0.3.1.55//DEBUG/lib/libp4est.so.0(p8est_partition_ext+0x1295) [0x2afcc7543daf]
//[t500u:16597] [21] /scratch/deal-trunk/deal.II/lib/libdeal_II.g.so.6.4.pre(_ZN6dealii8parallel11distributed13TriangulationILi3ELi3EE33execute_coarsening_and_refinementEv+0x483) [0x2afcbfba1b8f]
//[t500u:16597] [22] ./p4est_3d_refine_02/exe(_Z4testILi3EEvv+0x109) [0x410f1c]

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
#include <ostream>

template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
      GridGenerator::hyper_shell (tr,
			      Point<dim>(),
				  0.5, 1.0,
			      12,
			      true);

      int ind = 0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = tr.begin_active(); cell != tr.end(); ++cell, ++ind)
	if (!cell->is_artificial())
	  {
	    if (myid==0 && (ind==4 || ind==5 || ind==6 || ind== 8))
	      cell->set_refine_flag();
	    if (myid==1 && (ind==0 || ind==2 || ind==10))
	    cell->set_refine_flag();
	  }

      tr.execute_coarsening_and_refinement ();

      unsigned int checksum = tr.get_checksum ();
      if (myid == 0)
	{
	  deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
	  deallog << "Checksum: "
		  << checksum
		  << std::endl;
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


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_3d_refine_02").c_str());
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
