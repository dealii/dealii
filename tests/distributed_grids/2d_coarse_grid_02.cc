//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test interaction with p4est with a complicated 2d grid read from file. the
// grid describes a cross-section of an airfoil with flaps at the front and
// back. it has some 30,000 cells

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>



template<int dim>
void test(std::ostream& /*out*/)
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

  write_vtk (tr, "2d_coarse_grid_02", "1");
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  std::ofstream logfile("2d_coarse_grid_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
