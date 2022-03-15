// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check DataOut for complex vectors. The finite element used here is
// vector-valued, and we want to output both real and imaginary parts
// of the vector as separate real-valued vector fields. For this,
// DataOut needed to learn that it needs to duplicate not only the
// vector, but keep the dim real parts together, the dim complex parts
// together, and describe everything necessary to DataOutBase.
//
// This test checks the output in .gnuplot format.

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
  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Pick a special solution function: constant in space, but one
  // where the individual vector components have the form c-c*i
  // (counting components from 1, for uniqueness of components). This
  // allows us to check that the ordering of components in the output
  // file is correct.
  //
  // For example, in 2d, the ordering of output components should then
  // be
  //   [1 2] [-1 -2]
  // where [...] indicates which components jointly form a vector.
  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(1, -1) * double(i % fe.n_components() + 1);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(
    v,
    std::vector<std::string>(dim, "vector_field"),
    DataOut<dim>::type_dof_data,
    std::vector<DataComponentInterpretation::DataComponentInterpretation>(
      dim, DataComponentInterpretation::component_is_part_of_vector));
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
  catch (std::exception &exc)
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
