// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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



// test distributed solution_transfer with one processor (several solutions)

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (1);
  DoFHandler<dim> dofh(tr);
  static const FE_Q<dim> fe(2);
  dofh.distribute_dofs (fe);

  parallel::distributed::SolutionTransfer<dim, Vector<double> > soltrans(dofh);
  parallel::distributed::SolutionTransfer<dim, Vector<double> > soltrans2(dofh);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    {
      cell->set_refine_flag();
    }

  tr.prepare_coarsening_and_refinement();

  Vector<double> solution(dofh.n_dofs());

  Vector<double> solution1(dofh.n_dofs());
  solution1=1.0;

  Vector<double> solution2(dofh.n_dofs());
  solution2=2.0;

  std::vector< const Vector<double>* > sols;
  sols.push_back(&solution1);
  sols.push_back(&solution2);

  soltrans.prepare_for_coarsening_and_refinement(solution);
  soltrans2.prepare_for_coarsening_and_refinement(sols);

  tr.execute_coarsening_and_refinement ();

  dofh.distribute_dofs (fe);

  Vector<double> interpolated_solution(dofh.n_dofs());
  Vector<double> interpolated_solution1(dofh.n_dofs());
  Vector<double> interpolated_solution2(dofh.n_dofs());

  std::vector< Vector<double>* > sols_i;
  sols_i.push_back(&interpolated_solution1);
  sols_i.push_back(&interpolated_solution2);

  soltrans.interpolate(interpolated_solution);
  soltrans2.interpolate(sols_i);

  deallog << "#dofs:" << dofh.n_dofs() << std::endl;

  deallog << "norm: " << interpolated_solution1.l2_norm() << " " << interpolated_solution2.l2_norm() << " " << interpolated_solution.l2_norm() << std::endl;

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();

}
