// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// tests std::complex DataPostprocessorVector: create a FE field that has two
// components of the kind cos(something) and sin(something) and then have a
// postprocessor that computes the sum of squares. should always be equal to one
//
// this test uses the shortcut class DataPostprocessorVector to make
// writing postprocessors simpler


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



std::ofstream logfile("output");


template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();
  void
  run();

private:
  void
  make_grid_and_dofs();
  void
  solve();
  void
  output_results() const;

  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;

  Vector<std::complex<double>> solution;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
  : fe(FE_Q<dim>(1), 2)
  , dof_handler(triangulation)
{}



template <int dim>
void
LaplaceProblem<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());
}



template <int dim>
class SinesAndCosines : public Function<dim, std::complex<double>>
{
public:
  SinesAndCosines()
    : Function<dim, std::complex<double>>(2)
  {}

  std::complex<double>
  value(const Point<dim> &p, const unsigned int component) const
  {
    switch (component)
      {
        case 0:
          return std::sin(p.norm());
        case 1:
          return std::cos(p.norm());
        default:
          DEAL_II_NOT_IMPLEMENTED();
          return 0;
      }
  }
};



template <int dim>
void
LaplaceProblem<dim>::solve()
{
  // dummy solve. just insert some
  // values as mentioned at the top
  // of the file
  VectorTools::interpolate(dof_handler, SinesAndCosines<dim>(), solution);
}


template <int dim>
class MyPostprocessor : public DataPostprocessorVector<dim>
{
public:
  MyPostprocessor()
    : DataPostprocessorVector<dim>("magnitude_times_d", update_values)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    for (unsigned int q = 0; q < input_data.solution_values.size(); ++q)
      {
        Assert(computed_quantities[q].size() == dim, ExcInternalError());

        for (unsigned int d = 0; d < dim; ++d)
          computed_quantities[q](d) =
            input_data.solution_values[q](0) *
              input_data.solution_values[q](0) +
            input_data.solution_values[q](1) * input_data.solution_values[q](1);
        AssertThrow(std::fabs(computed_quantities[q](0) - 1) < 1e-12,
                    ExcInternalError());
      }
  }
};



template <int dim>
void
LaplaceProblem<dim>::output_results() const
{
  {
    // Save the data in the log file
    MyPostprocessor<dim> p;
    DataOut<dim>         data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, p);
    data_out.build_patches();
    data_out.write_gnuplot(logfile);
  }
}



template <int dim>
void
LaplaceProblem<dim>::run()
{
  make_grid_and_dofs();
  solve();
  output_results();
}



int
main()
{
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run();

  LaplaceProblem<3> laplace_problem_3d;
  laplace_problem_3d.run();

  return 0;
}
