// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// Check that we can generate an ABF(0) element in 3d. There used to
// be a number of exceptions that happened when one tried, so this
// test verifies that that no longer happens.
//
// The test projects a function into the space and checks that the
// error is zero (given that the function itself is inside the space
// -- in fact, it is a constant)

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



/*
 * Check the value of the derivative field.
 */

void EvaluateDerivative(DoFHandler<3> &dof_handler, Vector<double> &solution)
{
  // This quadrature rule determines the points, where the
  // derivative will be evaluated.
  QGauss<3>   quad(3);
  FEValues<3> fe_values(dof_handler.get_fe(),
                        quad,
                        UpdateFlags(update_values | update_quadrature_points |
                                    update_gradients | update_JxW_values));

  const unsigned int n_q_points    = quad.size();
  const unsigned int n_components  = dof_handler.get_fe().n_components();
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  // Cell iterators
  DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();

  double err_l2 = 0, err_hdiv = 0;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (; cell != endc; ++cell)
    {
      cell->get_dof_indices(local_dof_indices);

      fe_values.reinit(cell);

      // Get function values
      std::vector<Vector<double>> this_value(n_q_points,
                                             Vector<double>(n_components));
      fe_values.get_function_values(solution, this_value);

      // Get values from solution vector (For Trap.Rule)
      std::vector<std::vector<Tensor<1, 3>>> grads_here(
        n_q_points, std::vector<Tensor<1, 3>>(n_components));
      fe_values.get_function_gradients(solution, grads_here);

      for (const auto q_point : fe_values.quadrature_point_indices())
        {
          double u0 = 0;
          double v0 = 0;
          double w0 = 0;

          // evaluate the values of the solution at the quadrature point
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              u0 += (solution(local_dof_indices[i]) *
                     fe_values.shape_value_component(i, q_point, 0));
              v0 += (solution(local_dof_indices[i]) *
                     fe_values.shape_value_component(i, q_point, 1));
              w0 += (solution(local_dof_indices[i]) *
                     fe_values.shape_value_component(i, q_point, 2));
            }

          // subtract from it the exact solution
          u0 -= 1.0;
          v0 -= 1.0;
          w0 -= 1.0;

          err_l2 += (u0 * u0 + v0 * v0 + w0 * w0) * fe_values.JxW(q_point);

          // now also evaluate the gradient -- which should be zero
          // given that the exact solution is constant
          double dudx = grads_here[q_point][0][0];
          double dvdy = grads_here[q_point][1][1];
          double dwdz = grads_here[q_point][2][2];

          err_hdiv += (dudx + dvdy + dwdz) * (dudx + dvdy + dwdz) *
                      fe_values.JxW(q_point);
        }
    }

  AssertThrow(err_l2 < 1.e-20, ExcInternalError());
  AssertThrow(err_hdiv < 1.e-20, ExcInternalError());

  deallog << "OK" << std::endl;
}


template <int dim>
void
create_mass_matrix(const Mapping<dim> &       mapping,
                   const DoFHandler<dim> &    dof,
                   const Quadrature<dim> &    q,
                   SparseMatrix<double> &     matrix,
                   const Function<dim> &      rhs_function,
                   Vector<double> &           rhs_vector,
                   const Function<dim> *const coefficient = nullptr)
{
  UpdateFlags update_flags =
    UpdateFlags(update_values | update_JxW_values | update_quadrature_points);
  if (coefficient != nullptr)
    update_flags = UpdateFlags(update_flags | update_quadrature_points);

  FEValues<dim> fe_values(mapping, dof.get_fe(), q, update_flags);

  const unsigned int dofs_per_cell       = fe_values.dofs_per_cell,
                     n_q_points          = fe_values.n_quadrature_points;
  const FiniteElement<dim> &fe           = fe_values.get_fe();
  const unsigned int        n_components = fe.n_components();

  Assert(coefficient == nullptr || coefficient->n_components == 1 ||
           coefficient->n_components == n_components,
         ExcInternalError());

  FullMatrix<double>          cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              cell_vector(dofs_per_cell);
  std::vector<double>         coefficient_values(n_q_points);
  std::vector<Vector<double>> coefficient_vector_values(
    n_q_points, Vector<double>(n_components));

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  std::vector<Vector<double>> rhs_values(n_q_points,
                                         Vector<double>(n_components));

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                 endc = dof.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell->get_dof_indices(dof_indices);

      const std::vector<double> &weights = fe_values.get_JxW_values();
      rhs_function.vector_value_list(fe_values.get_quadrature_points(),
                                     rhs_values);
      cell_vector = 0;

      if (coefficient != nullptr)
        {
          if (coefficient->n_components == 1)
            {
              coefficient->value_list(fe_values.get_quadrature_points(),
                                      coefficient_values);
              for (const auto point : fe_values.quadrature_point_indices())
                {
                  const double weight = fe_values.JxW(point);
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const double v = fe_values.shape_value(i, point);
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const double u = fe_values.shape_value(j, point);

                          if ((n_components == 1) ||
                              (fe.system_to_component_index(i).first ==
                               fe.system_to_component_index(j).first))
                            cell_matrix(i, j) +=
                              (u * v * weight * coefficient_values[point]);
                        }
                    }
                }
            }
          else
            {
              coefficient->vector_value_list(fe_values.get_quadrature_points(),
                                             coefficient_vector_values);
              for (const auto point : fe_values.quadrature_point_indices())
                {
                  const double weight = fe_values.JxW(point);
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const double       v = fe_values.shape_value(i, point);
                      const unsigned int component_i =
                        fe.system_to_component_index(i).first;
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const double u = fe_values.shape_value(j, point);
                          if ((n_components == 1) ||
                              (fe.system_to_component_index(j).first ==
                               component_i))
                            cell_matrix(i, j) +=
                              (u * v * weight *
                               coefficient_vector_values[point](component_i));
                        }
                    }
                }
            }
        }
      else
        {
          // Compute eventual sign changes depending on the neighborhood
          // between two faces.
          std::vector<double>       sign_change(dofs_per_cell, 1.0);
          const unsigned int        dofs_per_face = fe.dofs_per_face;
          std::vector<unsigned int> face_dof_indices(dofs_per_face);

          for (const auto point : fe_values.quadrature_point_indices())
            {
              const double weight = fe_values.JxW(point);
              //      const double weight = q.weight(point);

              std::vector<Vector<double>> val_vector(
                dofs_per_cell, Vector<double>(n_components));

              // Precompute the component values
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components();
                     ++comp_i)
                  {
                    val_vector[i](comp_i) =
                      sign_change[i] *
                      fe_values.shape_value_component(i, point, comp_i);
                  }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components();
                     ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] == true)
                    {
                      const double v = val_vector[i](comp_i);
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        for (unsigned int comp_j = 0;
                             comp_j < fe.n_components();
                             ++comp_j)
                          if (fe.get_nonzero_components(j)[comp_j] == true)
                            {
                              const double u = val_vector[j](comp_j);
                              if ((n_components == 1) || (comp_i == comp_j))
                                cell_matrix(i, j) += (u * v * weight);
                            }
                    }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components();
                     ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] == true)
                    {
                      cell_vector(i) += rhs_values[point](comp_i) *
                                        val_vector[i](comp_i) * weights[point];
                    }
            }
        }
      // transfer everything into the global object
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          if ((n_components == 1) || (cell_matrix(i, j) != 0.0))
            matrix.add(dof_indices[i], dof_indices[j], cell_matrix(i, j));

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        rhs_vector(dof_indices[i]) += cell_vector(i);
    }
}


template <int dim>
void
create_right_hand_side(const Mapping<dim> &   mapping,
                       const DoFHandler<dim> &dof_handler,
                       const Quadrature<dim> &quadrature,
                       const Function<dim> &  rhs_function,
                       Vector<double> &       rhs_vector)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  Assert(fe.n_components() == rhs_function.n_components, ExcInternalError());
  Assert(rhs_vector.size() == dof_handler.n_dofs(),
         ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  rhs_vector = 0;

  UpdateFlags update_flags =
    UpdateFlags(update_values | update_quadrature_points | update_JxW_values);
  FEValues<dim> fe_values(mapping, fe, quadrature, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                     n_q_points    = fe_values.n_quadrature_points,
                     n_components  = fe.n_components();

  std::vector<unsigned int> dofs(dofs_per_cell);
  Vector<double>            cell_vector(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  if (n_components == 1)
    {
      std::vector<double> rhs_values(n_q_points);

      for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);

          const std::vector<double> &weights = fe_values.get_JxW_values();
          rhs_function.value_list(fe_values.get_quadrature_points(),
                                  rhs_values);

          cell_vector = 0;
          for (const auto point : fe_values.quadrature_point_indices())
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_vector(i) += rhs_values[point] *
                                fe_values.shape_value(i, point) *
                                weights[point];

          cell->get_dof_indices(dofs);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs_vector(dofs[i]) += cell_vector(i);
        }
    }
  else
    {
      std::vector<Vector<double>> rhs_values(n_q_points,
                                             Vector<double>(n_components));

      for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);

          const std::vector<double> &weights = fe_values.get_JxW_values();
          rhs_function.vector_value_list(fe_values.get_quadrature_points(),
                                         rhs_values);

          cell_vector = 0;
          for (const auto point : fe_values.quadrature_point_indices())
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int comp_i = 0; comp_i < fe.n_components();
                   ++comp_i)
                //        if (fe.get_nonzero_components(i)[comp_i] == true)
                {
                  double det = weights[point] / quadrature.weight(point);

                  cell_vector(i) +=
                    rhs_values[point](comp_i) *
                    fe_values.shape_value_component(i, point, comp_i) *
                    weights[point];
                }

          cell->get_dof_indices(dofs);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs_vector(dofs[i]) += cell_vector(i);
        }
    }
}



//
// This function replaces the deal.II implementation of the projection.
// The purpose is to have more freedom in assembling the matrix.
//

template <int dim>
void
project(const Mapping<dim> &             mapping,
        const DoFHandler<dim> &          dof,
        const AffineConstraints<double> &constraints,
        const Quadrature<dim> &          quadrature,
        const Function<dim> &            function,
        Vector<double> &                 vec,
        const bool                       enforce_zero_boundary = false,
        const Quadrature<dim - 1> &          = QGauss<dim - 1>(2),
        const bool project_to_boundary_first = false)
{
  Assert(dof.get_fe().n_components() == function.n_components,
         ExcInternalError());

  const FiniteElement<dim> &fe = dof.get_fe();

  // make up boundary values
  std::map<types::global_dof_index, double> boundary_values;

  if (enforce_zero_boundary == true)
    // no need to project boundary
    // values, but enforce
    // homogeneous boundary values
    // anyway
    {
      // loop over all boundary faces
      // to get all dof indices of
      // dofs on the boundary. note
      // that in 3d there are cases
      // where a face is not at the
      // boundary, yet one of its
      // lines is, and we should
      // consider the degrees of
      // freedom on it as boundary
      // nodes. likewise, in 2d and
      // 3d there are cases where a
      // cell is only at the boundary
      // by one vertex. nevertheless,
      // since we do not support
      // boundaries with dimension
      // less or equal to dim-2, each
      // such boundary dof is also
      // found from some other face
      // that is actually wholly on
      // the boundary, not only by
      // one line or one vertex
      typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                     endc = dof.end();
      std::vector<types::global_dof_index> face_dof_indices(fe.dofs_per_face);
      for (; cell != endc; ++cell)
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (cell->face(f)->at_boundary())
            {
              cell->face(f)->get_dof_indices(face_dof_indices);
              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                // enter zero boundary values
                // for all boundary nodes
                //
                // we need not care about
                // vector valued elements here,
                // since we set all components
                boundary_values[face_dof_indices[i]] = 0.;
            };
    }


  // set up mass matrix and right hand side
  vec.reinit(dof.n_dofs());
  SparsityPattern sparsity(dof.n_dofs(),
                           dof.n_dofs(),
                           dof.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof, sparsity);
  constraints.condense(sparsity);

  SparseMatrix<double> mass_matrix(sparsity);
  Vector<double>       tmp(mass_matrix.n());

  create_mass_matrix(mapping, dof, quadrature, mass_matrix, function, tmp);

  constraints.condense(mass_matrix);
  constraints.condense(tmp);
  if (boundary_values.size() != 0)
    MatrixTools::apply_boundary_values(
      boundary_values, mass_matrix, vec, tmp, true);

  SolverControl           control(1000, 1e-16, false, false);
  PrimitiveVectorMemory<> memory;
  SolverCG<>              cg(control, memory);

  PreconditionSSOR<> prec;
  prec.initialize(mass_matrix, 1.2);
  // solve
  cg.solve(mass_matrix, vec, tmp, prec);

  // distribute solution
  constraints.distribute(vec);
}



int
main()
{
  initlog();

  Triangulation<3> tria_test;
  GridGenerator::hyper_cube(tria_test);


  for (Triangulation<3>::active_cell_iterator cell = tria_test.begin_active();
       cell != tria_test.end();
       ++cell)
    {
      deallog << "Cell " << cell << std::endl;
      for (unsigned int v = 0; v < 4; ++v)
        deallog << "    " << cell->vertex(v) << std::endl;
    }


  FE_ABF<3> fe(0);
  deallog << "Dofs/cell " << fe.dofs_per_cell << ", Dofs/face "
          << fe.dofs_per_face << std::endl;

  DoFHandler<3> dof_handler(tria_test);
  dof_handler.distribute_dofs(fe);

  deallog << "Dofs total " << dof_handler.n_dofs() << std::endl;

  Vector<double> solution(dof_handler.n_dofs());
  solution = 1;

  // Project solution onto FE field
  AffineConstraints<double> hn_constraints;
  hn_constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, hn_constraints);
  hn_constraints.close();
  MappingQGeneric<3> map_default(1);
  project(map_default,
          dof_handler,
          hn_constraints,
          QGauss<3>(6),
          Functions::ConstantFunction<3>(1., 3),
          solution);

  EvaluateDerivative(dof_handler, solution);
}
