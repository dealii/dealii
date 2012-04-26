//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test DoFTools::map_dofs_to_support_points for parallel DoFHandlers

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tr);
  tr.refine_global (1);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  std::map<unsigned int, Point<dim> > points;
  DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
					dofh,
					points);
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (typename std::map<unsigned int, Point<dim> >::const_iterator
	     p = points.begin();
	   p != points.end();
	   ++p)
	deallog << p->first << " -> " << p->second
		<< std::endl;
    }

				   // the result of the call above is
				   // supposed to be a map that
				   // contains exactly the locally
				   // relevant dofs, so test this
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  for (unsigned int i=0; i<dofh.n_dofs(); ++i)
    if (relevant_set.is_element(i))
      Assert (points.find(i) != points.end(),
	      ExcInternalError())
      else
      Assert (points.find(i) == points.end(),
	      ExcInternalError());
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
      std::ofstream logfile(output_file_for_mpi("map_dofs_to_support_points").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
