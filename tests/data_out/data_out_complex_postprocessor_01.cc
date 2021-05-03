// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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


// check DataOut for complex vectors, using a postprocessor (indeed,
// a postprocessor similar to the one from step-29)

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

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
class ComputeMagnitude : public DataPostprocessorScalar<dim>
{
public:
  ComputeMagnitude();

  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override;
};

template <int dim>
ComputeMagnitude<dim>::ComputeMagnitude()
  : DataPostprocessorScalar<dim>("Magnitude", update_values)
{}


template <int dim>
void
ComputeMagnitude<dim>::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim> &inputs,
  std::vector<Vector<double>> &               computed_quantities) const
{
  Assert(computed_quantities.size() == inputs.solution_values.size(),
         ExcDimensionMismatch(computed_quantities.size(),
                              inputs.solution_values.size()));

  for (unsigned int i = 0; i < computed_quantities.size(); i++)
    {
      Assert(computed_quantities[i].size() == 1,
             ExcDimensionMismatch(computed_quantities[i].size(), 1));
      Assert(inputs.solution_values[i].size() == 2,
             ExcDimensionMismatch(inputs.solution_values[i].size(), 2));

      const std::complex<double> u(inputs.solution_values[i](0),
                                   inputs.solution_values[i](1));

      computed_quantities[i](0) = std::abs(u);
    }
}



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(1. * i, -1. * i);

  ComputeMagnitude<dim> postprocessor;
  DataOut<dim>          data_out;
  data_out.attach_dof_handler(dof_handler);

  // Output first the solution, then the postprocessed one. This
  // allows us to compare the latter in a visualization program
  data_out.add_data_vector(v, "solution");
  data_out.add_data_vector(v, postprocessor);
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
