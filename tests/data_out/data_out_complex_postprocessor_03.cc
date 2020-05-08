// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Check DataOut for complex vectors, using two postprocessors. The
// testcase is extracted from the step-58 tutorial. There was an
// indexing error when using more than one postprocessor.

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <string>

#include "../tests.h"


namespace DataPostprocessors
{
  template <int dim>
  class ComplexMagnitude : public DataPostprocessorScalar<dim>
  {
  public:
    ComplexMagnitude();

    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &               computed_quantities) const;
  };

  template <int dim>
  ComplexMagnitude<dim>::ComplexMagnitude()
    : DataPostprocessorScalar<dim>("Magnitude", update_values)
  {}


  template <int dim>
  void
  ComplexMagnitude<dim>::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &               computed_quantities) const
  {
    Assert(computed_quantities.size() == inputs.solution_values.size(),
           ExcDimensionMismatch(computed_quantities.size(),
                                inputs.solution_values.size()));

    for (unsigned int q = 0; q < computed_quantities.size(); ++q)
      {
        Assert(computed_quantities[q].size() == 1,
               ExcDimensionMismatch(computed_quantities[q].size(), 1));
        Assert(inputs.solution_values[q].size() == 2,
               ExcDimensionMismatch(inputs.solution_values[q].size(), 2));

        computed_quantities[q](0) =
          std::abs(std::complex<double>(inputs.solution_values[q](0),
                                        inputs.solution_values[q](1)));
      }
  }



  template <int dim>
  class ComplexPhase : public DataPostprocessorScalar<dim>
  {
  public:
    ComplexPhase();

    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &               computed_quantities) const;
  };

  template <int dim>
  ComplexPhase<dim>::ComplexPhase()
    : DataPostprocessorScalar<dim>("Phase", update_values)
  {}


  template <int dim>
  void
  ComplexPhase<dim>::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &               computed_quantities) const
  {
    Assert(computed_quantities.size() == inputs.solution_values.size(),
           ExcDimensionMismatch(computed_quantities.size(),
                                inputs.solution_values.size()));

    for (unsigned int q = 0; q < computed_quantities.size(); ++q)
      {
        Assert(computed_quantities[q].size() == 1,
               ExcDimensionMismatch(computed_quantities[q].size(), 1));
        Assert(inputs.solution_values[q].size() == 2,
               ExcDimensionMismatch(inputs.solution_values[q].size(), 2));

        computed_quantities[q](0) =
          std::arg(std::complex<double>(inputs.solution_values[q](0),
                                        inputs.solution_values[q](1)));
      }
  }

} // namespace DataPostprocessors



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Initialize the vector with a constant 0+1i. This will lead to
  // easily recognizable magnitudes and phases in the output file.
  Vector<std::complex<double>> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = std::complex<double>(0, 1);

  const DataPostprocessors::ComplexMagnitude<dim> complex_magnitude;
  const DataPostprocessors::ComplexPhase<dim>     complex_phase;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(v, complex_magnitude);
  data_out.add_data_vector(v, complex_phase);
  data_out.add_data_vector(v, "Psi");

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
