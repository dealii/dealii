//---------------------------------------------------------------------------
//    $Id: integrate_difference.cc 24418 2011-09-26 12:52:19Z bangerth $
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


// In Trilinos 10.12.1 interpolate() hangs with a small number of cells (smaller than the number of CPUs). This works fine in 10.4.2 and 10.8.5
// As it turns out, this is because vector creation after interpolate()
// reacts differently in 10.12.1 and causes a hang. The bug is
// that interpolate() does not call compress()!

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>


template <int dim>
class LinearFunction : public Function<dim> 
{
  public:
    double value (const Point<dim> &p,
		  const unsigned int) const
      {
	return 1;
      }
};


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

//tr.refine_global (2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);
  
  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
					     MPI_COMM_WORLD);
  VectorTools::interpolate (dofh,
			    LinearFunction<dim>(),
			    interpolated);
  interpolated.compress();

  double norm = interpolated.l2_norm();
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "norm = " << norm
	    << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("interpolate_03").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
				       //    test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
//      test<3>();
      deallog.pop();
    }

#ifdef DEAL_II_WITH_MPI
  MPI_Finalize();
#endif
}
