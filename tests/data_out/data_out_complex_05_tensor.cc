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


// Like the _05 test, just for tensor-valued quantities.

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

  // Create a vector-valued finite element
  FESystem<dim>   fe(FE_Q<dim>(1), dim * dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Pick a special solution function: constant in space, but one
  // where the individual tensor components have the form c-c*i
  // (counting components from 1, for uniqueness of components). This
  // allows us to check that the ordering of components in the output
  // file is correct.
  //
  // For example, in 2d, the ordering of output components should then
  // be
  //   [1 2 3 4] [-1 -2 -3 -4]
  // where [...] indicates which components jointly form a tensor.
  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(1, -1) * double(i % fe.n_components() + 1);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(
    v,
    std::vector<std::string>(dim * dim, "tensor_field"),
    DataOut<dim>::type_dof_data,
    std::vector<DataComponentInterpretation::DataComponentInterpretation>(
      dim * dim, DataComponentInterpretation::component_is_part_of_tensor));
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
