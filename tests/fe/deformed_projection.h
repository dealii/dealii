// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

/*
 * Project the function [1,1] onto a deformed grid and see whether the element
 * elements can represent it exactly. This is the common code for the various
 * *_projection_* tests.
 */



#include "../tests.h"
#include <deal.II/base/logstream.h>

#define PRECISION 2

#include <fstream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_dgq.h>

void test();


template <int dim>
class TestMap1 : public Function<dim>
{
public:
  TestMap1 (const unsigned int n_components) :
    Function<dim> (n_components)
  {}

  virtual ~TestMap1 () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  void vector_value (const Point<dim> &p,
                     Vector<double>   &return_value) const;
};



template <int dim>
double
TestMap1<dim>::value (const Point<dim> &,
                      const unsigned int  ) const
{
  return (1.0);
}


template <int dim>
void TestMap1<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  // Parabolic inflow profile
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}


/*
 * Check the value of the derivative field.
 */

void EvaluateDerivative (DoFHandler<2> *dof_handler,
                         Vector<double> &solution)
{
  // This quadrature rule determines the points, where the
  // derivative will be evaluated.
  QGauss<2> quad (3);
  FEValues<2> fe_values (dof_handler->get_fe (), quad,
                         UpdateFlags(update_values    |
                                     update_q_points  |
                                     update_gradients |
                                     update_JxW_values));

  const unsigned int   n_q_points    = quad.size();
  const unsigned int   n_components   = dof_handler->get_fe().n_components();
  const unsigned int dofs_per_cell = dof_handler->get_fe().dofs_per_cell;

  // Cell iterators
  DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active(),
                                      endc = dof_handler->end();

  double err_l2 = 0,
         err_hdiv = 0;

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  for (; cell!=endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);

      fe_values.reinit (cell);

      // Get function values
      std::vector<Vector<double> > this_value(n_q_points,
                                              Vector<double>(n_components));
      fe_values.get_function_values (solution, this_value);

      // Get values from solution vector (For Trap.Rule)
      std::vector<std::vector<Tensor<1,2> > >
      grads_here (n_q_points, std::vector<Tensor<1,2> > (n_components));
      fe_values.get_function_grads (solution, grads_here);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          //   double u0 = this_value[q_point](0);
          //double v0 = this_value[q_point](1);

          double u0 = 0;
          double v0 = 0;
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              u0 += (solution (local_dof_indices[i]) *
                     fe_values.shape_value_component(i, q_point, 0));
              v0 += (solution (local_dof_indices[i]) *
                     fe_values.shape_value_component(i, q_point, 1));
            }
          u0 -= 1.0;
          v0 -= 1.0;

          double dudx = grads_here[q_point][0][0];
          double dvdy = grads_here[q_point][1][1];

          /*
            sprintf (buf, "QP %i %e %e %e %e %e %e Div %e\n", q_point, u0, v0,
            dudx, dudy, dvdx, dvdy, dudx + dvdy);
          */

          err_l2 += (u0 * u0 + v0 * v0) * fe_values.JxW (q_point);
          err_hdiv += (dudx + dvdy) * (dudx + dvdy) * fe_values.JxW (q_point);
        }
    }

  deallog << "L2-Err " << pow (err_l2, 0.5)
          << "  Hdiv-Err " << pow (err_hdiv, 0.5)
          << std::endl;
}


template <int dim>
void create_mass_matrix (const Mapping<dim>       &mapping,
                         const DoFHandler<dim>    &dof,
                         const Quadrature<dim>    &q,
                         SparseMatrix<double>     &matrix,
                         const Function<dim>   &rhs_function,
                         Vector<double>        &rhs_vector,
                         const Function<dim> *const coefficient = 0)
{
  UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_q_points);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);

  FEValues<dim> fe_values (mapping, dof.get_fe(), q, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert(coefficient == 0 ||
         coefficient->n_components==1 ||
         coefficient->n_components==n_components, ExcInternalError());

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  std::vector<Vector<double> > coefficient_vector_values (n_q_points,
                                                          Vector<double> (n_components));

  std::vector<types::global_dof_index> dof_indices (dofs_per_cell);

  std::vector<Vector<double> > rhs_values(n_q_points, Vector<double>(n_components));

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active (),
                                                 endc = dof.end ();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell->get_dof_indices (dof_indices);

      const std::vector<double> &weights   = fe_values.get_JxW_values ();
      rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
      cell_vector = 0;

      if (coefficient != 0)
        {
          if (coefficient->n_components==1)
            {
              coefficient->value_list (fe_values.get_quadrature_points(),
                                       coefficient_values);
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const double weight = fe_values.JxW(point);
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                      const double v = fe_values.shape_value(i,point);
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                          const double u = fe_values.shape_value(j,point);

                          if ((n_components==1) ||
                              (fe.system_to_component_index(i).first ==
                               fe.system_to_component_index(j).first))
                            cell_matrix(i,j) +=
                              (u * v * weight * coefficient_values[point]);
                        }
                    }
                }
            }
          else
            {
              coefficient->vector_value_list (fe_values.get_quadrature_points(),
                                              coefficient_vector_values);
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const double weight = fe_values.JxW(point);
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                      const double v = fe_values.shape_value(i,point);
                      const unsigned int component_i=
                        fe.system_to_component_index(i).first;
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                          const double u = fe_values.shape_value(j,point);
                          if ((n_components==1) ||
                              (fe.system_to_component_index(j).first == component_i))
                            cell_matrix(i,j) +=
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
          std::vector<double> sign_change (dofs_per_cell, 1.0);
          const unsigned int dofs_per_face = fe.dofs_per_face;
          std::vector<types::global_dof_index> face_dof_indices (dofs_per_face);

          for (unsigned int f = 0; f < 2; ++f)
            {
              typename DoFHandler<dim>::active_face_iterator face = cell->face (f);
              if (!face->at_boundary ())
                {
                  unsigned int nn = cell->neighbor_of_neighbor (f);
                  deallog << "Face " << f << " NeigNeig " << nn << std::endl;
                  if (nn > 1)
                    {
                      face->get_dof_indices (face_dof_indices);
                      for (unsigned int j = 0; j < dofs_per_face; ++j)
                        {
                          sign_change[face_dof_indices[j]] = -1.0;
                          deallog << "DoF " << face_dof_indices[j] << std::endl;
                        }
                    }
                }
            }

          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const double weight = fe_values.JxW(point);
              //     const double weight = q.weight(point);

              std::vector<Vector<double> > val_vector (dofs_per_cell,
                                                       Vector<double> (n_components));

              // Precompute the component values
              for (unsigned int i=0; i < dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components ();
                     ++comp_i)
                  {
                    val_vector[i](comp_i) = sign_change[i] *
                                            fe_values.shape_value_component(i,point,comp_i);
                  }
              // Now eventually switch the sign of some of the ansatzfunctions.
              // TODO

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components ();
                     ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] == true)
                    {
                      const double v = val_vector[i](comp_i);
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        for (unsigned int comp_j = 0;
                             comp_j < fe.n_components (); ++comp_j)
                          if (fe.get_nonzero_components(j)[comp_j] == true)
                            {
                              const double u = val_vector[j](comp_j);
                              if ((n_components==1) ||
                                  (comp_i == comp_j))
                                cell_matrix(i,j) += (u * v * weight);
                            }
                    }


              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int comp_i = 0; comp_i < fe.n_components ();
                     ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] == true)
                    {
                      cell_vector(i) += rhs_values[point](comp_i) *
                                        val_vector[i](comp_i) * weights[point];
                    }
            }
        }
      // transfer everything into the
      // global object. lock the
      // matrix meanwhile
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          if ((n_components==1) ||
              (cell_matrix (i,j) != 0.0))
            /*
              (fe.system_to_component_index(i).first ==
              fe.system_to_component_index(j).first))
            */
            matrix.add (dof_indices[i], dof_indices[j],
                        cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        rhs_vector(dof_indices[i]) += cell_vector(i);
    }
}


template <int dim>
void create_right_hand_side (const Mapping<dim>    &mapping,
                             const DoFHandler<dim> &dof_handler,
                             const Quadrature<dim> &quadrature,
                             const Function<dim>   &rhs_function,
                             Vector<double>        &rhs_vector)
{
  const FiniteElement<dim> &fe  = dof_handler.get_fe();
  Assert (fe.n_components() == rhs_function.n_components,
          ExcInternalError());
  Assert (rhs_vector.size() == dof_handler.n_dofs(),
          ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  rhs_vector = 0;

  UpdateFlags update_flags = UpdateFlags(update_values   |
                                         update_q_points |
                                         update_JxW_values);
  FEValues<dim> fe_values (mapping, fe, quadrature, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                     n_q_points    = fe_values.n_quadrature_points,
                     n_components  = fe.n_components();

  std::vector<types::global_dof_index> dofs (dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  if (n_components==1)
    {
      std::vector<double> rhs_values(n_q_points);

      for (; cell!=endc; ++cell)
        {
          fe_values.reinit(cell);

          const std::vector<double> &weights   = fe_values.get_JxW_values ();
          rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);

          cell_vector = 0;
          for (unsigned int point=0; point<n_q_points; ++point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              cell_vector(i) += rhs_values[point] *
                                fe_values.shape_value(i,point) *
                                weights[point];

          cell->get_dof_indices (dofs);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            rhs_vector(dofs[i]) += cell_vector(i);
        }
    }
  else
    {
      std::vector<Vector<double> > rhs_values(n_q_points, Vector<double>(n_components));

      for (; cell!=endc; ++cell)
        {
          fe_values.reinit(cell);

          const std::vector<double> &weights   = fe_values.get_JxW_values ();
          rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

          cell_vector = 0;
          for (unsigned int point=0; point<n_q_points; ++point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int comp_i = 0; comp_i < fe.n_components ();
                   ++comp_i)
                //       if (fe.get_nonzero_components(i)[comp_i] == true)
                {
                  cell_vector(i) += rhs_values[point](comp_i) *
                                    fe_values.shape_value_component(i,point,comp_i) *
                                    weights[point];
                }

          cell->get_dof_indices (dofs);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            rhs_vector(dofs[i]) += cell_vector(i);
        }
    }
}



//
// This function replaces the deal.II implementation of the projection.
// The purpose is to have more freedom in assembling the matrix.
//

template <int dim>
void project (const Mapping<dim>       &mapping,
              const DoFHandler<dim>    &dof,
              const ConstraintMatrix   &constraints,
              const Quadrature<dim>    &quadrature,
              const Function<dim>      &function,
              Vector<double>           &vec,
              const bool                enforce_zero_boundary = false,
              const Quadrature<dim-1>  & = QGauss<dim-1>(2),
              const bool                project_to_boundary_first = false)
{
  Assert (dof.get_fe().n_components() == function.n_components,
          ExcInternalError());

  const FiniteElement<dim> &fe = dof.get_fe();

  // make up boundary values
  std::map<types::global_dof_index,double> boundary_values;

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
      std::vector<types::global_dof_index> face_dof_indices (fe.dofs_per_face);
      for (; cell!=endc; ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            {
              cell->face(f)->get_dof_indices (face_dof_indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                // enter zero boundary values
                // for all boundary nodes
                //
                // we need not care about
                // vector valued elements here,
                // since we set all components
                boundary_values[face_dof_indices[i]] = 0.;
            };
    }
  else
    // no homogeneous boundary values
    if (project_to_boundary_first == true)
      // boundary projection required
      {
        /*
                 // set up a list of boundary functions for
                          // the different boundary parts. We want the
                                   // @p{function} to hold on all parts of the
                                            // boundary
                                            typename FunctionMap<dim>::type boundary_functions;
                                            for (unsigned char c=0; c<255; ++c)
                                            boundary_functions[c] = &function;
                                            project_boundary_values (dof, boundary_functions, q_boundary,
                                            boundary_values);
        */
      }


  // set up mass matrix and right hand side
  vec.reinit (dof.n_dofs());
  SparsityPattern sparsity(dof.n_dofs(),
                           dof.n_dofs(),
                           dof.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof, sparsity);
  constraints.condense (sparsity);

  SparseMatrix<double> mass_matrix (sparsity);
  Vector<double> tmp (mass_matrix.n());

  create_mass_matrix (mapping, dof, quadrature, mass_matrix, function, tmp);
  deallog << "MM created" << std::endl;
  create_right_hand_side (mapping, dof, quadrature, function, tmp);
  deallog << "RHS created" << std::endl;

  constraints.condense (mass_matrix);
  constraints.condense (tmp);
  if (boundary_values.size() != 0)
    MatrixTools::apply_boundary_values (boundary_values,
                                        mass_matrix, vec, tmp,
                                        true);

  SolverControl           control(1000,1e-16);
  PrimitiveVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionSSOR<> prec;
  prec.initialize(mass_matrix, 1.2);
  // solve
  cg.solve (mass_matrix, vec, tmp, prec);

  // distribute solution
  constraints.distribute (vec);
}


int create_tria (unsigned int elm, Triangulation<2> &tria)
{
  std::vector<Point<2> > points_glob;
  std::vector<Point<2> > points;

  points_glob.push_back (Point<2> (0.0, 0.0));
  points_glob.push_back (Point<2> (1.0, 0.0));
  points_glob.push_back (Point<2> (1.0, 0.5));
  points_glob.push_back (Point<2> (1.0, 1.0));
  points_glob.push_back (Point<2> (0.6, 0.5));
  points_glob.push_back (Point<2> (0.5, 1.0));
  points_glob.push_back (Point<2> (0.0, 1.0));

  // Prepare cell data
  std::vector<CellData<2> > cells (1);

  switch (elm)
    {
    case 0:
      cells[0].vertices[0] = 0;
      cells[0].vertices[1] = 1;
      cells[0].vertices[2] = 4;
      cells[0].vertices[3] = 2;
      cells[0].material_id = 0;
      break;
    case 1:
      cells[0].vertices[0] = 4;
      cells[0].vertices[1] = 2;
      cells[0].vertices[2] = 5;
      cells[0].vertices[3] = 3;
      cells[0].material_id = 0;
      break;
    case 2:
      cells[0].vertices[0] = 0;
      cells[0].vertices[1] = 4;
      cells[0].vertices[2] = 6;
      cells[0].vertices[3] = 5;
      cells[0].material_id = 0;
      break;
    }

  for (unsigned int i = 0; i < 4; ++i)
    {
      points.push_back (points_glob[cells[0].vertices[i]]);
      cells[0].vertices[i] = i;
    }

  tria.create_triangulation (points, cells, SubCellData());

  return (0);
}


void plot_shapes (DoFHandler<2> &dof_handler)
{
  Vector<double> solution (dof_handler.n_dofs ());
  std::set<unsigned int> face_dofs;

  // Create set of all DoFs, which are on the boundary.
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();
  std::vector<types::global_dof_index> face_dof_indices (dof_handler.get_fe().dofs_per_face);
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      {
        cell->face(f)->get_dof_indices (face_dof_indices);
        for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_face; ++i)
          face_dofs.insert (face_dof_indices[i]);
      }

  // Now set a 1 at the place of the different DoFs and
  // output the solution.
  std::set<unsigned int>::iterator face_dof_iter;
  for (face_dof_iter = face_dofs.begin ();
       face_dof_iter != face_dofs.end (); ++face_dof_iter)
    {
      unsigned int dof = *face_dof_iter;

      deallog << dof << std::endl;

      solution (dof) = 1.0;

      // Test the core functionality2
      DataOut<2> *data_out = new DataOut<2>;
      data_out->attach_dof_handler (dof_handler);
      data_out->add_data_vector (solution, "solution");
      data_out->build_patches (4);

      data_out->write_gnuplot (deallog.get_file_stream());

      delete data_out;

      solution (dof) = 0.0;
    }
}


int main (int /*argc*/, char **/*argv*/)
{
  std::ofstream logfile (logname);
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}



void check (const FiniteElement<2> &fe)
{
  Triangulation<2> tria_test;
  DoFHandler<2> *dof_handler;

  deallog << "Dofs/cell " << fe.dofs_per_cell
          << "Dofs/face " << fe.dofs_per_face << std::endl;

  for (unsigned int elm = 0; elm < 3; ++elm)
    {
      create_tria (elm, tria_test);

      // Create a DoFHandler
      dof_handler = new DoFHandler<2> (tria_test);
      dof_handler->distribute_dofs (fe);

      deallog << "Dofs total " << dof_handler->n_dofs () << std::endl;

      // Alloc some DoFs
      Vector<double> solution(dof_handler->n_dofs ());

      // Project solution onto FE field
      ConstraintMatrix     hn_constraints;
      hn_constraints.clear ();
      DoFTools::make_hanging_node_constraints (*dof_handler,
                                               hn_constraints);
      hn_constraints.close ();

      MappingQ1<2> map_default;

      project (map_default, *dof_handler, hn_constraints,
               QGauss<2> (6), TestMap1<2>(2),
               solution);

      // Test the core functionality
      DataOut<2> *data_out = new DataOut<2>;
      data_out->attach_dof_handler (*dof_handler);
      data_out->add_data_vector (solution, "solution");
      data_out->build_patches (4);

      data_out->write_gnuplot (deallog.get_file_stream());

      // Now write face DoFs ..
      for (unsigned int i = 0; i < 4; ++i)
        {
          deallog << solution(i) << std::endl;
        }

      // Clean up ...
      delete data_out;
      delete (dof_handler);
      tria_test.clear ();
    }
}
