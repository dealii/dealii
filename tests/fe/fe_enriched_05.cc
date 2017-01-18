// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test h-refinement with DoFHandler and print constraints.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_collection.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>

const double eps = 1e-10;

// argument for build_patches()
const unsigned int patches = 10;

using namespace dealii;

// uncomment when debugging
// #define DATA_OUT_FE_ENRICHED

template<int dim>
class EnrichmentFunction : public Function<dim>
{
public:
  EnrichmentFunction()
    : Function<dim>(1)
  {}

  virtual double value(const Point<dim> &point,
                       const unsigned int component = 0 ) const
  {
    return std::exp(-point.norm());
  }

  virtual Tensor<1,dim> gradient(const Point<dim> &point,
                                 const unsigned int component = 0) const
  {
    Tensor<1,dim> res = point;
    Assert (point.norm() > 0,
            dealii::ExcMessage("gradient is not defined at zero"));
    res *= -value(point)/point.norm();
    return res;
  }
};


template <int dim>
void test4()
{
  deallog << "h-refinement:"<<std::endl;

  Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler(triangulation);

  EnrichmentFunction<dim> function;
  FE_Enriched<dim> fe(FE_Q<dim>(2),
                      FE_Q<dim>(1),
                      &function);

  GridGenerator::hyper_cube(triangulation);

  // one global refinement and 1 refinement of 1 active cell (local)
  triangulation.refine_global();
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  dof_handler.distribute_dofs(fe);

  ConstraintMatrix constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  constraints.print(deallog.get_file_stream());

#ifdef DATA_OUT_FE_ENRICHED
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string> names;
  for (unsigned int s=0; s < dof_handler.n_dofs(); s++)
    if (!constraints.is_constrained(s))
      {
        names.push_back(std::string("N_") +
                        dealii::Utilities::int_to_string(s));

        Vector<double> shape_function;
        shape_function.reinit(dof_handler.n_dofs());
        shape_function[s] = 1.0;

        // make continuous:
        constraints.distribute(shape_function);
        shape_functions.push_back(shape_function);
      }

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  for (unsigned int i = 0; i < shape_functions.size(); i++)
    data_out.add_data_vector (shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string filename = "h-shape_functions_"+dealii::Utilities::int_to_string(dim)+"D.vtu";
  std::ofstream output (filename.c_str ());
  data_out.write_vtu (output);
#endif

  dof_handler.clear();
}


int main (int argc,char **argv)
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(4);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  try
    {
      test4<2>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
