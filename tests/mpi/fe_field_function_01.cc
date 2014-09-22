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



// Test VectorTools::fe_field_function_01 for parallel computations. we
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
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

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

  // create a mesh and refine it sufficiently often that only one
  // processor will have the cell we care about
  GridGenerator::hyper_cube(tr);
  tr.refine_global (5);

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

  Functions::FEFieldFunction<dim,DoFHandler<dim>,TrilinosWrappers::MPI::Vector>
  field_function (dofh, x_rel);

  for (unsigned int test=0; test<4; ++test)
    {
      double value;
      Point<dim> p = (dim == 2 ?
                      Point<dim>(test/2+1,test%2+1)/3 :
                      Point<dim>(test/2+1,test/2+1,test%2+1)/3);

      // see if we can find the point on the current processor
      bool point_found = false;
      try
        {
          value = field_function.value (p);
          point_found = true;

          Assert (std::fabs(value - (p[0]+2)) < 1e-8* std::fabs(value + (p[0]+2)),
                  ExcInternalError());
        }
      catch (typename Functions::FEFieldFunction<dim,DoFHandler<dim>,TrilinosWrappers::MPI::Vector>::ExcPointNotAvailableHere &)
        {
          point_found = false;
        }

      Assert (Utilities::MPI::sum(point_found ? 1 : 0, MPI_COMM_WORLD) == 1,
              ExcInternalError());
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

}
