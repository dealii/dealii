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



// refine  a 3d cell that is not marked once (currently a bug in p4est):
/*
#17 0x00007fffebfde43b in sc_abort_verbose (
    filename=0x7fffec2ef7c8 "/scratch/p4estbuild/p4est-0.3.1.55-67fe1/src/p4est_bits.c", lineno=756, msg=0x7fffec2efc07 "Assertion 'q->level > 0'")
    at /scratch/p4estbuild/p4est-0.3.1.55-67fe1/sc/src/sc.c:603
#18 0x00007fffec2a2d4b in p8est_quadrant_parent (q=0x834670, r=0x833a50)
    at /scratch/p4estbuild/p4est-0.3.1.55-67fe1/src/p4est_bits.c:756
#19 0x00007fffec29e46c in p4est_correct_partition (p4est=0x835b40,
    num_quadrants_in_proc=0x834750)
    at /scratch/p4estbuild/p4est-0.3.1.55-67fe1/src/p4est.c:2430
#20 0x00007fffec29ddaf in p8est_partition_ext (p4est=0x835b40,
    partition_for_coarsening=1, weight_fn=0)
    at /scratch/p4estbuild/p4est-0.3.1.55-67fe1/src/p4est.c:2293
#21 0x00007ffff68a73c7 in dealii::parallel::distributed::Triangulation<3, 3>::execute_coarsening_and_refinement (this=0x7fffffffb7d0)
    at /scratch/deal-trunk/deal.II/source/distributed/tria.cc:2612
#22 0x0000000000410e6b in test<3> () at p4est_3d_refine_03.cc:64
#23 0x000000000040c6e5 in main (argc=1, argv=0x7fffffffdd98)
    at p4est_3d_refine_02.cc:106
    */

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
  GridGenerator::hyper_cube(tr);

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
