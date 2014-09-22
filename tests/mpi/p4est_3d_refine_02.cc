// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

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

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();

}
