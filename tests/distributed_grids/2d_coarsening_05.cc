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


// Pretty much exactly like coarsening_02, except that we go to around 50,000
// cells. this is a similar case to coarsening_03 (where we start with a
// coarse grid of 30,000 cells, however) and that takes a ton of time at the
// time of writing this

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

  GridGenerator::hyper_cube(tr);
  tr.refine_global (1);

  GridGenerator::hyper_cube(tr2);
  tr2.refine_global (1);

  Assert (tr.n_active_cells() == tr2.n_active_cells(),
	  ExcInternalError());


  while (tr.n_active_cells() < 50000)
    {
      InterGridMap<Triangulation<dim> > intergrid_map;
      intergrid_map.make_mapping (tr, tr2);

      std::vector<bool> flags (tr.n_active_cells(), false);

				       // refine one fifth of all cells each
				       // time (but at least one)
      for (unsigned int i=0; i<tr.n_active_cells() / 5 + 1; ++i)
	flags[rand() % flags.size()] = true;

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

				       // flag all other cells for coarsening
				       // (this should ensure that at least
				       // some of them will actually be
				       // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = tr.begin_active();
	   cell != tr.end(); ++cell, ++index)
	if (!flags[index])
	  {
	    cell->set_coarsen_flag();
	    intergrid_map[cell]->set_coarsen_flag();
	  }


      tr.execute_coarsening_and_refinement ();
      tr2.execute_coarsening_and_refinement ();

      deallog << " Number of cells: "
	      << tr.n_active_cells() << ' '
	      << tr2.n_active_cells()
	      << std::endl;

      assert_tria_equal("2d_coarsening_05", tr, tr2);
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

  std::ofstream logfile("2d_coarsening_05/output");
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
