// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Check DataOut for complex vectors, using a postprocessor (indeed,
// a postprocessor similar to the one from step-29)
//
// This is a test similar to the _01 test, but for a vector-valued
// element.

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
class ComputeMagnitudes : public DataPostprocessor<dim>
{
public:
  virtual std::vector<std::string>
  get_names() const override;

  virtual UpdateFlags
  get_needed_update_flags() const override;

  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override;
};


template <int dim>
std::vector<std::string>
ComputeMagnitudes<dim>::get_names() const
{
  return {"abs_u", "abs_v"};
}



template <int dim>
UpdateFlags
ComputeMagnitudes<dim>::get_needed_update_flags() const
{
  return update_values;
}


template <int dim>
void
ComputeMagnitudes<dim>::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim> &inputs,
  std::vector<Vector<double>> &               computed_quantities) const
{
  Assert(computed_quantities.size() == inputs.solution_values.size(),
         ExcDimensionMismatch(computed_quantities.size(),
                              inputs.solution_values.size()));

  for (unsigned int i = 0; i < computed_quantities.size(); i++)
    {
      Assert(computed_quantities[i].size() == 2,
             ExcDimensionMismatch(computed_quantities[i].size(), 2));
      Assert(inputs.solution_values[i].size() == 4,
             ExcDimensionMismatch(inputs.solution_values[i].size(), 4));

      const std::complex<double> u(inputs.solution_values[i](0),
                                   inputs.solution_values[i](1));
      const std::complex<double> v(inputs.solution_values[i](2),
                                   inputs.solution_values[i](3));

      computed_quantities[i](0) = std::abs(u);
      computed_quantities[i](1) = std::abs(v);
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

  // Create a system of 2 elements, each of which will have
  // complex-valued solution components
  FESystem<dim>   fe(FE_Q<dim>(1), 2);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Initialize them with more or less random values
  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(1. * i, -1. * i);

  ComputeMagnitudes<dim> postprocessor;
  DataOut<dim>           data_out;
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
