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


// Test VectorTools::compute_mean_value for parallel computations. we
// interpolate a linear function onto the grid with a symmetric mesh. the mean
// value of the interpolation must be the mean of the linear function

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
#include <deal.II/numerics/vectors.h>

#include <fstream>


template <int dim>
class LinearFunction : public Function<dim> 
{
  public:
    double value (const Point<dim> &p,
		  const unsigned int) const
      {
	return p[0] + 2;
      }
};


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tr);
  tr.refine_global (2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
					     MPI_COMM_WORLD);
  VectorTools::interpolate (dofh,
			    LinearFunction<dim>(),
			    interpolated);
  
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = interpolated;

  const double mean
    = VectorTools::compute_mean_value (dofh, QGauss<dim>(2), x_rel, 0);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "mean=" << mean
	    << std::endl;
  
  Assert (std::fabs(mean - 2) < 1e-3, ExcInternalError());
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
      std::ofstream logfile(output_file_for_mpi("compute_mean_value").c_str());
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
