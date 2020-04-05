// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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



// test hp_constraints when neither_element_dominates.
// 2D and 3D p-refinement only for 2 elements.
// Check continuity across the face by evaluating the
// field from both sides.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

//#define DEBUG_OUTPUT_VTK
unsigned int counter = 0;

const double eps = 1e-10;


template <int dim>
void
test2cells(const unsigned int p1 = 2, const unsigned int p2 = 1)
{
  Assert(dim > 1, ExcInternalError());
  Triangulation<dim> triangulation;
  {
    Point<dim> p1, p2;
    for (unsigned int d = 0; d < dim; d++)
      p1[d] = -1;
    p2[0] = 1.0;
    std::vector<unsigned int> repetitoins(dim, 1);
    repetitoins[0] = 2;
    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitoins,
                                              p1,
                                              p2);
  }

  hp::DoFHandler<dim> dof_handler(triangulation);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(
    FESystem<dim>(FE_Q<dim>(p1), 1, FE_Nothing<dim>(1, true), 1));
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(p2), 1, FE_Q<dim>(1), 1));
  // push back to be able to resolve hp constrains no matter what:
  fe_collection.push_back(
    FESystem<dim>(FE_Q<dim>(p2), 1, FE_Nothing<dim>(1, true), 1));

  hp::QCollection<dim - 1> q_face_collection;
  q_face_collection.push_back(QIterated<dim - 1>(QTrapez<1>(), 2));
  q_face_collection.push_back(QIterated<dim - 1>(QTrapez<1>(), 2));
  q_face_collection.push_back(QIterated<dim - 1>(QTrapez<1>(), 2));

  deallog << "2cells: " << fe_collection[0].get_name() << " vs "
          << fe_collection[1].get_name() << std::endl;

  dof_handler.begin_active()->set_active_fe_index(1); // Q(p2)xQ(1)

  dof_handler.distribute_dofs(fe_collection);

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());

#ifdef DEBUG_OUTPUT_VTK
  // output to check if all is good:
  counter++;
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (unsigned int s = 0; s < dof_handler.n_dofs(); s++)
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

  DataOut<dim, hp::DoFHandler<dim>> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get material ids:
  Vector<float> fe_index(triangulation.n_active_cells());
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler
                                                              .begin_active(),
                                                     endc = dof_handler.end();
  for (unsigned int index = 0; cell != endc; ++cell, ++index)
    {
      fe_index[index] = cell->active_fe_index();
    }
  data_out.add_data_vector(fe_index, "fe_index");

  for (unsigned int i = 0; i < shape_functions.size(); i++)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(0);
  std::string filename = "shape_functions_" +
                         dealii::Utilities::int_to_string(counter, 1) + "_" +
                         dealii::Utilities::int_to_string(dim) + "D.vtu";
  std::ofstream output(filename.c_str());
  data_out.write_vtu(output);
#endif

  // fill some vector
  Vector<double> solution(dof_handler.n_dofs());
  for (unsigned int dof = 0; dof < dof_handler.n_dofs(); dof++)
    solution[dof] = 21.0 * (dof % 2) + 0.5 + dof % 3;

  constraints.distribute(solution);


  // evaluate field at the interface:
  hp::FEFaceValues<dim> fe_face_values_hp(
    fe_collection, q_face_collection, update_values | update_quadrature_points);

  std::vector<unsigned int>   local_face_dof_indices;
  std::vector<Vector<double>> values;
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      const unsigned int fe_index = cell->active_fe_index();
      local_face_dof_indices.resize(fe_collection[fe_index].dofs_per_face);
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (!cell->at_boundary(f)) // that's enough for our simple mesh
          {
            deallog << "cell=" << cell << " face=" << f << std::endl;
            fe_face_values_hp.reinit(cell, f);
            const FEFaceValues<dim> &fe_face_values =
              fe_face_values_hp.get_present_fe_values();

            const unsigned int n_q_points = fe_face_values.n_quadrature_points;
            values.resize(n_q_points, Vector<double>(2));

            fe_face_values.get_function_values(solution, values);

            const std::vector<dealii::Point<dim>> &q_points =
              fe_face_values.get_quadrature_points();

            for (unsigned int q = 0; q < n_q_points; q++)
              deallog << "u[" << q_points[q] << "]={" << values[q][0] << ","
                      << values[q][1] << "}" << std::endl;
          }
    }

  dof_handler.clear();
}
int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);
  deallog << std::fixed;

  try
    {
      test2cells<2>(
        1, 2); // Q(1) x FE_Nothing vs Q(2) x Q(1) => common Q(1) x FE_Nothing
      test2cells<2>(
        2, 1); // Q(2) x FE_Nothing vs Q(1) x Q(1) => common Q(1) x FE_Nothing
      test2cells<3>(
        1, 2); // Q(1) x FE_Nothing vs Q(2) x Q(1) => common Q(1) x FE_Nothing
      test2cells<3>(
        2, 1); // Q(2) x FE_Nothing vs Q(1) x Q(1) => common Q(1) x FE_Nothing
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
