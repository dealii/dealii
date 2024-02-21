// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// hp-refinement. Output constraints.

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

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
test6(const bool         do_href,
      const unsigned int p_feq  = 1,
      const unsigned int p_feen = 2)
{
  deallog << "hp: " << do_href << ' ' << p_feq << ' ' << p_feen << std::endl;

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  EnrichmentFunction<dim> function;

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Enriched<dim>(FE_Q<dim>(p_feq)));
  fe_collection.push_back(
    FE_Enriched<dim>(FE_Q<dim>(p_feen), FE_Q<dim>(1), &function));
  // push back to be able to resolve hp-constrains:
  fe_collection.push_back(FE_Enriched<dim>(FE_Q<dim>(p_feen)));

  GridGenerator::hyper_cube(triangulation);

  // one global refinement and 1 refinement of 1 active cell
  // a sketch with material ids:
  // ---------------
  // |      |      |
  // |  1   |  0   |
  // |------|-------
  // |1  |1 |      |
  // |------|  0   |
  // |1  |1 |      |
  // |------|------|
  triangulation.refine_global();
  {
    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active();
    cell->set_active_fe_index(1); // POU
    cell++;
    cell->set_active_fe_index(1); // POU
  }
  dof_handler.distribute_dofs(fe_collection);

  if (do_href)
    {
      triangulation.begin_active()->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
      dof_handler.distribute_dofs(fe_collection);
    }

  deallog << "n_cells: " << triangulation.n_active_cells() << std::endl;

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());

#ifdef DATA_OUT_FE_ENRICHED
  // output to check if all is good:
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (unsigned int s = 0; s < dof_handler.n_dofs(); ++s)
    {
      Vector<double> shape_function;
      shape_function.reinit(dof_handler.n_dofs());
      shape_function[s] = 1.0;

      // if the dof is constrained, first output unconstrained vector
      if (constraints.is_constrained(s))
        {
          names.push_back(std::string("UN_") +
                          dealii::Utilities::int_to_string(s, 2));
          shape_functions.push_back(shape_function);
        }

      names.push_back(std::string("N_") +
                      dealii::Utilities::int_to_string(s, 2));

      // make continuous/constrain:
      constraints.distribute(shape_function);
      shape_functions.push_back(shape_function);
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get material ids:
  Vector<float> fe_index(triangulation.n_active_cells());
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int index = 0; cell != endc; ++cell, ++index)
    {
      fe_index[index] = cell->active_fe_index();
    }
  data_out.add_data_vector(fe_index, "fe_index");

  for (unsigned int i = 0; i < shape_functions.size(); ++i)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string filename = "hp-shape_functions_ref=" + std::to_string(do_href) +
                         "_p_feq=" + std::to_string(p_feq) +
                         "_p_feenr=" + std::to_string(p_feen) + "_" +
                         dealii::Utilities::int_to_string(dim) + "D.vtu";
  std::ofstream output(filename);
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
      test6<2>(false, 1, 2); // 1 vs 2+1
    }
  catch (const std::exception &exc)
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
