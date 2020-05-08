// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// no refinement, just 4 cells.
// Make sure fe_values and shape functions are continuous

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
test5()
{
  deallog << "4 cells:" << std::endl;

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  EnrichmentFunction<dim> function;
  FE_Enriched<dim>        fe(FE_Q<dim>(1), FE_Q<dim>(1), &function);

  GridGenerator::hyper_cube(triangulation);

  triangulation.refine_global();

  dof_handler.distribute_dofs(fe);

  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (unsigned int s = 0; s < dof_handler.n_dofs(); s++)
    {
      names.push_back(std::string("N_") + dealii::Utilities::int_to_string(s));

      Vector<double> shape_function;
      shape_function.reinit(dof_handler.n_dofs());
      shape_function[s] = 1.0;

      shape_functions.push_back(shape_function);
    }

  // output 11th:
  {
    const unsigned int    global_dof = 11;
    const Vector<double> &solution   = shape_functions[global_dof];
    QTrapez<dim>          quadrature;
    FEValues<dim>         fe_values(fe, quadrature, update_values);

    const unsigned int  n_q_points = quadrature.size();
    std::vector<double> solution_values(n_q_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<dealii::types::global_dof_index> local_dof_indices;
    local_dof_indices.resize(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);

        // find out which
        unsigned int local_dof = 0;
        cell->get_dof_indices(local_dof_indices);
        for (; local_dof < dofs_per_cell; local_dof++)
          if (local_dof_indices[local_dof] == global_dof)
            break;

        const std::vector<dealii::Point<dim>> &q_points =
          fe_values.get_quadrature_points();
        fe_values.get_function_values(solution, solution_values);

        deallog << " cell=" << cell->center() << std::endl;
        for (const auto q_point : fe_values.quadrature_point_indices())
          {
            // find non-zero shape_value
            deallog << " qp=" << q_points[q_point]
                    << " f(qp)=" << function.value(q_points[q_point]);

            // if the cell contains our global dof:
            if (local_dof < dofs_per_cell)
              deallog << " N(" << local_dof
                      << ",qp)=" << fe_values.shape_value(local_dof, q_point);

            deallog << " U(qp)=" << solution_values[q_point] << std::endl;
          }
      }
  }

#ifdef DATA_OUT_FE_ENRICHED
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  for (unsigned int i = 0; i < shape_functions.size(); i++)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string filename =
    "4cell_functions_" + dealii::Utilities::int_to_string(dim) + "D.vtu";
  std::ofstream output(filename.c_str());
  data_out.write_vtu(output);
#endif

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
      test5<2>();
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
