// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Test FEValues::get_function_values()

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <iostream>

#include "../tests.h"

const double eps = 1e-10;

// argument for build_patches()
const unsigned int patches = 10;


// uncomment when debugging
// #define DATA_OUT_FE_ENRICHED

template <int dim>
class EnrichmentFunction : public Function<dim>
{
public:
  EnrichmentFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    return std::exp(-point.norm());
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> res = point;
    Assert(point.norm() > 0,
           dealii::ExcMessage("gradient is not defined at zero"));
    res *= -value(point) / point.norm();
    return res;
  }
};


template <int dim>
void
test3()
{
  deallog << "FEValues.get_function_values()" << std::endl;
  deallog << "for same underlying FEs: f(qp) * N_{fe}(qp) == N_{pou}(qp)"
          << std::endl;
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  EnrichmentFunction<dim> function;
  FE_Enriched<dim>        fe(FE_Q<dim>(1), FE_Q<dim>(1), &function);

  GridGenerator::hyper_cube(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim>   quadrature(2);
  FEValues<dim> fe_values(fe,
                          quadrature,
                          update_values | update_gradients | update_JxW_values);

  Vector<double> solution_fe(dof_handler.n_dofs()),
    solution_pou(dof_handler.n_dofs()), solution(dof_handler.n_dofs());
  solution_fe  = 0;
  solution_pou = 0;

  // groups are one after another
  // thus set to 1.0 the same shape function in the underlying FE
  solution_fe[0]  = 2.0;
  solution_pou[1] = 2.0;
  solution        = solution_fe;
  solution += solution_pou;

  const unsigned int n_q_points = quadrature.size();

  std::vector<double> solution_values_fe(n_q_points),
    solution_values_pou(n_q_points), solution_values(n_q_points);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      const unsigned int                     dofs_per_cell = fe.dofs_per_cell;
      const std::vector<dealii::Point<dim>> &q_points =
        fe_values.get_quadrature_points();
      fe_values.get_function_values(solution_fe, solution_values_fe);
      fe_values.get_function_values(solution_pou, solution_values_pou);
      fe_values.get_function_values(solution, solution_values);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        deallog << " qp=" << q_points[q_point]
                << " f(qp)=" << function.value(q_points[q_point])
                << " U_fe(qp)=" << solution_values_fe[q_point]
                << " U_pou(qp)=" << solution_values_pou[q_point]
                << " U_[fe+pou](qp)=" << solution_values[q_point] << std::endl;
    }
  dof_handler.clear();
}

template <int dim>
void
plot_shape_function()
{
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  EnrichmentFunction<dim> function;
  FE_Enriched<dim>        fe(FE_Q<dim>(1), FE_Q<dim>(1), &function);

  GridGenerator::hyper_cube(triangulation);
  dof_handler.distribute_dofs(fe);

  std::vector<Vector<double>> shape_functions(dof_handler.n_dofs());
  std::vector<std::string>    names;
  for (unsigned int s = 0; s < shape_functions.size(); s++)
    {
      names.push_back(std::string("N_") + dealii::Utilities::int_to_string(s));

      shape_functions[s].reinit(dof_handler.n_dofs());

      // this is ok for 1 cell only!
      const unsigned int group = fe.system_to_base_index(s).first.first;

      // zero everywhere but
      // the DoF corresponding to this shape function
      if (group == 0)
        shape_functions[s][s] = 1.0;
      else
        shape_functions[s][s] = 1.0; // can potentially put another value here
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  for (unsigned int i = 0; i < shape_functions.size(); i++)
    {
      data_out.add_data_vector(shape_functions[i], names[i]);
    }

  data_out.build_patches(patches);

  std::string filename =
    "shape_functions_" + dealii::Utilities::int_to_string(dim) + "D.vtu";
  std::ofstream output(filename.c_str());
  data_out.write_vtu(output);

  dof_handler.clear();
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4) << std::fixed;
  deallog.depth_console(0);

  try
    {
      test3<3>();
#ifdef DATA_OUT_FE_ENRICHED
      plot_shape_function<1>();
      plot_shape_function<2>();
      plot_shape_function<3>();
#endif
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
