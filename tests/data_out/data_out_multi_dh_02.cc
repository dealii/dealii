// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests DataOut with multiple DoFHandler objects (two times the same element /
// dh)

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs()), v2(dof1.n_dofs());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = i;
  for (unsigned int i = 0; i < v2.size(); ++i)
    v2(i) = -v1(i);

  DataOut<dim> data_out;
  data_out.add_data_vector(dof1, v1, "data1");
  data_out.add_data_vector(dof1, v2, "data2");
  data_out.build_patches();

  data_out.write_vtk(deallog.get_file_stream());
}


int
main()
{
  try
    {
      initlog();
      deallog << std::setprecision(2);

      test<1>();
      test<2>();
      test<3>();

      return 0;
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
    };
}
