// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check DataOut for complex vectors, using the vector-valued path
//
// this test uses the method that generates the names of vector
// components from one base name by appending _1, _2; a second step
// then adds _re and _im
//
// like the test without _mixed, but output both a real-valued and a
// complex-valued vector

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <string>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FESystem<dim>   fe(FE_Q<dim>(1), 2);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(1. * i, -1. * i);

  Vector<double> w(dof_handler.n_dofs());
  for (unsigned int i = 0; i < w.size(); ++i)
    w(i) = 2 * i;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(v, "node_data");
  data_out.add_data_vector(w, "real_node_data");
  data_out.build_patches();

  data_out.write_gnuplot(deallog.get_file_stream());
}



int
main()
{
  initlog();

  try
    {
      check<1>();
      check<2>();
      check<3>();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
}
