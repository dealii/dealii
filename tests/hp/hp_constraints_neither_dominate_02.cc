// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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
// 2D and 3D hp-refinement.
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
struct less_than_key
{
  inline bool
  operator()(const std::pair<Point<dim>, Vector<double>> &pair1,
             const std::pair<Point<dim>, Vector<double>> &pair2)
  {
    const double      precision = 1e-3;
    const Point<dim> &p1        = pair1.first;
    const Point<dim> &p2        = pair2.first;

    for (unsigned int d = 0; d < dim; d++)
      {
        const bool is_equal = (std::abs(p1[d] - p2[d]) < precision);
        if (!is_equal)
          return p1[d] < p2[d];
      }

    return false; //< They are equal
  }
};

//  ------------
//  |0 |0 |    |
//  |-----| 2  |
//  |1 |1 |    |
//  ------------
template <int dim>
void
test2cells(const FiniteElement<dim> &fe_0,
           const FiniteElement<dim> &fe_1,
           const FiniteElement<dim> &fe_2,
           const FiniteElement<dim> &fe_common)
{
  const unsigned int n_comp = fe_0.n_components();
  Triangulation<dim> triangulation;
  {
    Point<dim> p1, p2;
    for (unsigned int d = 0; d < dim; d++)
      p2[d] = 1.0;
    p1[0] = -1.0;
    std::vector<unsigned int> repetitoins(dim, 1);
    repetitoins[0] = 2;
    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitoins,
                                              p1,
                                              p2);
  }

  hp::DoFHandler<dim> dof_handler(triangulation);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe_0);
  fe_collection.push_back(fe_1);
  fe_collection.push_back(fe_2);
  fe_collection.push_back(fe_common);

  const QIterated<dim - 1> quad_formula(QTrapez<1>(), 2);

  hp::QCollection<dim - 1> q_face_collection;
  q_face_collection.push_back(quad_formula);
  q_face_collection.push_back(quad_formula);
  q_face_collection.push_back(quad_formula);
  q_face_collection.push_back(quad_formula);

  deallog << "hp-ref: " << fe_0.get_name() << " vs " << fe_1.get_name()
          << " vs " << fe_2.get_name() << " with " << fe_common.get_name()
          << std::endl;

  // refine once:
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  for (typename hp::DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      const Point<dim> &center = cell->center();
      if (center[0] > 0.0) // right cell
        cell->set_active_fe_index(2);
      else if (center[1] > 0.5)
        cell->set_active_fe_index(0);
      else
        cell->set_active_fe_index(1);
    }
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
    solution[dof] = 1.5 * (dof % 7) + 0.5 * dim + 2.0 * (dof % 3);

  constraints.distribute(solution);


  // evaluate field at the interface:
  hp::FEFaceValues<dim> fe_face_values_hp(
    fe_collection, q_face_collection, update_values | update_quadrature_points);

  std::vector<unsigned int>   local_face_dof_indices;
  std::vector<Vector<double>> values;

  std::vector<std::pair<Point<dim>, Vector<double>>> pairs_point_value;

  for (typename hp::DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      const unsigned int fe_index = cell->active_fe_index();
      local_face_dof_indices.resize(fe_collection[fe_index].dofs_per_face);
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (std::abs(cell->face(f)->center()[0]) < 0.1)
          {
            // deallog << "cell="<<cell<<" face="<<f<<std::endl;
            fe_face_values_hp.reinit(cell, f);
            const FEFaceValues<dim> &fe_face_values =
              fe_face_values_hp.get_present_fe_values();

            const unsigned int n_q_points = fe_face_values.n_quadrature_points;
            values.resize(n_q_points, Vector<double>(n_comp));

            fe_face_values.get_function_values(solution, values);

            const std::vector<dealii::Point<dim>> &q_points =
              fe_face_values.get_quadrature_points();

            for (unsigned int q = 0; q < n_q_points; q++)
              {
                // since our face is [0,1]^{dim-1}, the quadrature rule
                // will coincide with quadrature points on mother face.
                // Use that to limit output of sub-faces at the same quadrature
                // points only.
                Point<dim - 1> qpt;
                for (unsigned int d = 0; d < dim - 1; d++)
                  qpt[d] = q_points[q][d + 1];

                unsigned int q_found = 0;
                for (; q_found < quad_formula.size(); q_found++)
                  if (quad_formula.point(q_found).distance(qpt) < 1e-5)
                    break;

                if (q_found < quad_formula.size())
                  {
                    pairs_point_value.push_back(
                      std::make_pair(q_points[q], values[q]));
                    //                      deallog << "@"<<q_points[q]<<" u =
                    //                      {"<<values[q][0]; for (unsigned int
                    //                      c = 1; c < n_comp; c++)
                    //                        deallog << ","<<values[q][c];
                    //                      deallog << "}"<<std::endl;
                  }
              }
          }
    }

  // sort
  std::sort(pairs_point_value.begin(),
            pairs_point_value.end(),
            less_than_key<dim>());
  for (unsigned int p = 0; p < pairs_point_value.size(); p++)
    {
      const Point<dim> &    pt  = pairs_point_value[p].first;
      const Vector<double> &val = pairs_point_value[p].second;

      Assert(val.size() == n_comp, ExcInternalError());
      deallog << "@" << pt << " u = {" << val[0];
      for (unsigned int c = 1; c < n_comp; c++)
        deallog << "," << val[c];
      deallog << "}" << std::endl;
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
      {
        const unsigned int dim = 2;
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
      }

      {
        const unsigned int dim = 3;
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
        test2cells<dim>(
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(2), 1),
          FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
      }
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
