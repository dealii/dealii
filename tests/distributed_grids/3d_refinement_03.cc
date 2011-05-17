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


// Like refinement_02, but with a complex grid

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
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test(std::ostream& /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  Triangulation<dim> tr2 (Triangulation<dim>::limit_level_difference_at_vertices);

  {
    GridIn<dim> gi;
    gi.attach_triangulation (tr);
    std::ifstream in ("../deal.II/grid_in_3d/4.in");
    gi.read_xda (in);
  }

  {
    GridIn<dim> gi;
    gi.attach_triangulation (tr2);
    std::ifstream in ("../deal.II/grid_in_3d/4.in");
    gi.read_xda (in);
  }

  Assert (tr.n_active_cells() == tr2.n_active_cells(),
	  ExcInternalError());


  for (unsigned int i=0; i<1; ++i)
    {
      std::vector<bool> flags (tr.n_active_cells());
      for (unsigned int j=0; j<flags.size(); ++j)
	flags[j] = (rand() < 0.2*RAND_MAX);

      InterGridMap<Triangulation<dim> > intergrid_map;
      intergrid_map.make_mapping (tr, tr2);

				       // refine tr and tr2
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = tr.begin_active();
	   cell != tr.end(); ++cell, ++index)
	if (flags[index])
	  {
	    cell->set_refine_flag();
	    intergrid_map[cell]->set_refine_flag();
	  }
      Assert (index == tr.n_active_cells(), ExcInternalError());
      tr.execute_coarsening_and_refinement ();
      tr2.execute_coarsening_and_refinement ();

      write_vtk (tr, "3d_refinement_03", "1");
      deallog << std::endl;

      deallog << i << " Number of cells: "
	      << tr.n_active_cells() << ' '
	      << tr2.n_active_cells()
	      << std::endl;

      assert_tria_equal("3d_refinement_03", tr, tr2);

    }
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  std::ofstream logfile("3d_refinement_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
