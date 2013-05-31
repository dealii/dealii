//---------------------------------------------------------------------------
//    $Id: p4est_2d_coarse_01.cc 29074 2013-03-27 14:53:45Z heister $
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


// test n_levels and iterators on those levels
// begin(i) should work even if we have no cells on that level

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

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  if (myid == 0)
    deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
  
  deallog << "proc " << myid << ", n levels " << tr.n_levels() << std::endl;
  deallog << "proc " << myid << ", n gobal levels " << tr.n_global_levels() << std::endl;
  
  deallog << "begin().cell_index " << tr.begin()->index() << std::endl;
  deallog << "begin(0).cell_index " << tr.begin(0)->index() << std::endl;
  
  deallog << "begin(1)==end(1)? " << (tr.begin(1)==tr.end(1)) << std::endl;

  
  deallog << "subdomainid = "
	  << tr.begin_active()->subdomain_id()
	  << std::endl;

  //if (myid!=0)
    //   Assert(tr.begin(1)==tr.end(1), ExcInternalError());  
  
  const unsigned int checksum = tr.get_checksum ();
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "Checksum: "
	    << checksum
	    << std::endl;

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log(__FILE__);
  
  deallog.push("2d");
  test<2>();
  deallog.pop();
}
