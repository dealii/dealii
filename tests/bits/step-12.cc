// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a un-hp-ified version of hp/step-12


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("output");

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/base/timer.h>

#include <iomanip>
#include <fstream>


template <int dim>
class RHS:  public Function<dim>
{
public:
  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<double> &values,
                           const unsigned int component=0) const;
};


template <int dim>
class BoundaryValues:  public Function<dim>
{
public:
  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<double> &values,
                           const unsigned int component=0) const;
};


template <int dim>
class Beta
{
public:
  Beta () {}
  void value_list (const std::vector<Point<dim> > &points,
                   std::vector<Point<dim> > &values) const;
};


template <int dim>
void RHS<dim>::value_list(const std::vector<Point<dim> > &points,
                          std::vector<double> &values,
                          const unsigned int) const
{
  Assert(values.size()==points.size(),
         ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    values[i]=0;
}


template <int dim>
void Beta<dim>::value_list(const std::vector<Point<dim> > &points,
                           std::vector<Point<dim> > &values) const
{
  Assert(values.size()==points.size(),
         ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      const Point<dim> &p=points[i];
      Point<dim> &beta=values[i];

      beta(0) = -p(1);
      beta(1) = p(0);
      beta /= std::sqrt(beta.square());
    }
}


template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
                                     std::vector<double> &values,
                                     const unsigned int) const
{
  Assert(values.size()==points.size(),
         ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    {
      if (points[i](0)<0.5)
        values[i]=1.;
      else
        values[i]=0.;
    }
}


template <int dim>
class DGTransportEquation
{
public:
  DGTransportEquation();

  void assemble_cell_term(const FEValues<dim> &fe_v,
                          FullMatrix<double> &ui_vi_matrix,
                          Vector<double> &cell_vector) const;

  void assemble_boundary_term(const FEFaceValues<dim> &fe_v,
                              FullMatrix<double> &ui_vi_matrix,
                              Vector<double> &cell_vector) const;

  void assemble_face_term1(const FEFaceValuesBase<dim> &fe_v,
                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                           FullMatrix<double> &ui_vi_matrix,
                           FullMatrix<double> &ue_vi_matrix) const;

  void assemble_face_term2(const FEFaceValuesBase<dim> &fe_v,
                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                           FullMatrix<double> &ui_vi_matrix,
                           FullMatrix<double> &ue_vi_matrix,
                           FullMatrix<double> &ui_ve_matrix,
                           FullMatrix<double> &ue_ve_matrix) const;
private:
  const Beta<dim> beta_function;
  const RHS<dim> rhs_function;
  const BoundaryValues<dim> boundary_function;
};


template <int dim>
DGTransportEquation<dim>::DGTransportEquation ()
  :
  beta_function (),
  rhs_function (),
  boundary_function ()
{}


template <int dim>
void DGTransportEquation<dim>::assemble_cell_term(
  const FEValues<dim> &fe_v,
  FullMatrix<double> &ui_vi_matrix,
  Vector<double> &cell_vector) const
{
  const std::vector<double> &JxW = fe_v.get_JxW_values ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> rhs (fe_v.n_quadrature_points);

  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  rhs_function.value_list (fe_v.get_quadrature_points(), rhs);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      {
        for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
          ui_vi_matrix(i,j) -= beta[point]*fe_v.shape_grad(i,point)*
                               fe_v.shape_value(j,point) *
                               JxW[point];

        cell_vector(i) += rhs[point] * fe_v.shape_value(i,point) * JxW[point];
      }
}


template <int dim>
void DGTransportEquation<dim>::assemble_boundary_term(
  const FEFaceValues<dim> &fe_v,
  FullMatrix<double> &ui_vi_matrix,
  Vector<double> &cell_vector) const
{
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> g(fe_v.n_quadrature_points);

  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  boundary_function.value_list (fe_v.get_quadrature_points(), g);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
      if (beta_n>0)
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
            ui_vi_matrix(i,j) += beta_n *
                                 fe_v.shape_value(j,point) *
                                 fe_v.shape_value(i,point) *
                                 JxW[point];
      else
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          cell_vector(i) -= beta_n *
                            g[point] *
                            fe_v.shape_value(i,point) *
                            JxW[point];
    }
}


template <int dim>
void DGTransportEquation<dim>::assemble_face_term1(
  const FEFaceValuesBase<dim> &fe_v,
  const FEFaceValuesBase<dim> &fe_v_neighbor,
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix) const
{
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
      if (beta_n>0)
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
            ui_vi_matrix(i,j) += beta_n *
                                 fe_v.shape_value(j,point) *
                                 fe_v.shape_value(i,point) *
                                 JxW[point];
      else
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            ue_vi_matrix(i,k) += beta_n *
                                 fe_v_neighbor.shape_value(k,point) *
                                 fe_v.shape_value(i,point) *
                                 JxW[point];
    }
}


template <int dim>
void DGTransportEquation<dim>::assemble_face_term2(
  const FEFaceValuesBase<dim> &fe_v,
  const FEFaceValuesBase<dim> &fe_v_neighbor,
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix,
  FullMatrix<double> &ui_ve_matrix,
  FullMatrix<double> &ue_ve_matrix) const
{
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);

  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
      if (beta_n>0)
        {
          for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
              ui_vi_matrix(i,j) += beta_n *
                                   fe_v.shape_value(j,point) *
                                   fe_v.shape_value(i,point) *
                                   JxW[point];

          for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
              ui_ve_matrix(k,j) -= beta_n *
                                   fe_v.shape_value(j,point) *
                                   fe_v_neighbor.shape_value(k,point) *
                                   JxW[point];
        }
      else
        {
          for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
              ue_vi_matrix(i,l) += beta_n *
                                   fe_v_neighbor.shape_value(l,point) *
                                   fe_v.shape_value(i,point) *
                                   JxW[point];

          for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
              ue_ve_matrix(k,l) -= beta_n *
                                   fe_v_neighbor.shape_value(l,point) *
                                   fe_v_neighbor.shape_value(k,point) *
                                   JxW[point];
        }
    }
}


template <int dim>
class DGMethod
{
public:
  DGMethod ();
  ~DGMethod ();

  void run ();

private:
  void setup_system ();
  void assemble_system1 ();
  void assemble_system2 ();
  void solve (Vector<double> &solution);
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  const MappingQ1<dim> mapping;

  FE_DGQ<dim>          fe;
  DoFHandler<dim>      dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  const QGauss<dim>   quadrature;
  const QGauss<dim-1> face_quadrature;

  Vector<double>       solution1;
  Vector<double>       solution2;
  Vector<double>       right_hand_side;

  const DGTransportEquation<dim> dg;
};


template <int dim>
DGMethod<dim>::DGMethod ()
  :
  mapping (),
  fe (1),
  dof_handler (triangulation),
  quadrature (4),
  face_quadrature (4),
  dg ()
{}


template <int dim>
DGMethod<dim>::~DGMethod ()
{
  dof_handler.clear ();
}


template <int dim>
void DGMethod<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           (GeometryInfo<dim>::faces_per_cell
                            *GeometryInfo<dim>::max_children_per_face+1)*fe.dofs_per_cell);

  DoFTools::make_flux_sparsity_pattern (dof_handler, sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution1.reinit (dof_handler.n_dofs());
  solution2.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());
}


template <int dim>
void DGMethod<dim>::assemble_system1 ()
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dofs (dofs_per_cell);
  std::vector<types::global_dof_index> dofs_neighbor (dofs_per_cell);

  const UpdateFlags update_flags = update_values
                                   | update_gradients
                                   | update_q_points
                                   | update_JxW_values;

  const UpdateFlags face_update_flags = update_values
                                        | update_q_points
                                        | update_JxW_values
                                        | update_normal_vectors;

  const UpdateFlags neighbor_face_update_flags = update_values;

  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);

  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);

  FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_vi_matrix (dofs_per_cell, dofs_per_cell);

  Vector<double>  cell_vector (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      ui_vi_matrix = 0;
      cell_vector = 0;

      fe_v.reinit (cell);

      dg.assemble_cell_term(fe_v,
                            ui_vi_matrix,
                            cell_vector);

      cell->get_dof_indices (dofs);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          typename DoFHandler<dim>::face_iterator face=cell->face(face_no);

          ue_vi_matrix = 0;

          if (face->at_boundary())
            {
              fe_v_face.reinit (cell, face_no);

              dg.assemble_boundary_term(fe_v_face,
                                        ui_vi_matrix,
                                        cell_vector);
            }
          else
            {
              typename DoFHandler<dim>::cell_iterator neighbor=
                cell->neighbor(face_no);;

              if (face->has_children())
                {
                  const unsigned int neighbor2=
                    cell->neighbor_of_neighbor(face_no);


                  for (unsigned int subface_no=0;
                       subface_no<face->n_children(); ++subface_no)
                    {
                      typename DoFHandler<dim>::active_cell_iterator
                      neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);

                      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
                              ExcInternalError());
                      Assert (!neighbor_child->has_children(), ExcInternalError());

                      ue_vi_matrix = 0;

                      fe_v_subface.reinit (cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

                      dg.assemble_face_term1(fe_v_subface,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix);

                      neighbor_child->get_dof_indices (dofs_neighbor);

                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int k=0; k<dofs_per_cell; ++k)
                          system_matrix.add(dofs[i], dofs_neighbor[k],
                                            ue_vi_matrix(i,k));
                    }
                }
              else
                {
                  if (neighbor->level() == cell->level())
                    {
                      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

                      fe_v_face.reinit (cell, face_no);
                      fe_v_face_neighbor.reinit (neighbor, neighbor2);

                      dg.assemble_face_term1(fe_v_face,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix);
                    }
                  else
                    {
                      Assert(neighbor->level() < cell->level(), ExcInternalError());

                      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
                        cell->neighbor_of_coarser_neighbor(face_no);
                      const unsigned int neighbor_face_no=faceno_subfaceno.first,
                                         neighbor_subface_no=faceno_subfaceno.second;

                      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                                   neighbor_subface_no)
                              == cell,
                              ExcInternalError());

                      fe_v_face.reinit (cell, face_no);
                      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
                                                    neighbor_subface_no);

                      dg.assemble_face_term1(fe_v_face,
                                             fe_v_subface_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix);
                    }

                  neighbor->get_dof_indices (dofs_neighbor);

                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int k=0; k<dofs_per_cell; ++k)
                      system_matrix.add(dofs[i], dofs_neighbor[k],
                                        ue_vi_matrix(i,k));
                }
            }
        }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        right_hand_side(dofs[i]) += cell_vector(i);
    }
}



template <int dim>
void DGMethod<dim>::assemble_system2 ()
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dofs (dofs_per_cell);
  std::vector<types::global_dof_index> dofs_neighbor (dofs_per_cell);

  const UpdateFlags update_flags = update_values
                                   | update_gradients
                                   | update_q_points
                                   | update_JxW_values;

  const UpdateFlags face_update_flags = update_values
                                        | update_q_points
                                        | update_JxW_values
                                        | update_normal_vectors;

  const UpdateFlags neighbor_face_update_flags = update_values;

  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);


  FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_vi_matrix (dofs_per_cell, dofs_per_cell);

  FullMatrix<double> ui_ve_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_ve_matrix (dofs_per_cell, dofs_per_cell);

  Vector<double>  cell_vector (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      ui_vi_matrix = 0;
      cell_vector = 0;

      fe_v.reinit (cell);

      dg.assemble_cell_term(fe_v,
                            ui_vi_matrix,
                            cell_vector);

      cell->get_dof_indices (dofs);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          typename DoFHandler<dim>::face_iterator face=
            cell->face(face_no);

          if (face->at_boundary())
            {
              fe_v_face.reinit (cell, face_no);

              dg.assemble_boundary_term(fe_v_face,
                                        ui_vi_matrix,
                                        cell_vector);
            }
          else
            {
              Assert (cell->neighbor(face_no).state() == IteratorState::valid,
                      ExcInternalError());
              typename DoFHandler<dim>::cell_iterator neighbor=
                cell->neighbor(face_no);
              if (face->has_children())
                {
                  const unsigned int neighbor2=
                    cell->neighbor_of_neighbor(face_no);

                  for (unsigned int subface_no=0;
                       subface_no<face->n_children(); ++subface_no)
                    {
                      typename DoFHandler<dim>::cell_iterator neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);
                      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
                              ExcInternalError());
                      Assert (!neighbor_child->has_children(), ExcInternalError());

                      ue_vi_matrix = 0;
                      ui_ve_matrix = 0;
                      ue_ve_matrix = 0;

                      fe_v_subface.reinit (cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

                      dg.assemble_face_term2(fe_v_subface,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix,
                                             ui_ve_matrix,
                                             ue_ve_matrix);

                      neighbor_child->get_dof_indices (dofs_neighbor);

                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            system_matrix.add(dofs[i], dofs_neighbor[j],
                                              ue_vi_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs[j],
                                              ui_ve_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
                                              ue_ve_matrix(i,j));
                          }
                    }
                }
              else
                {
                  if (neighbor->level() == cell->level() &&
                      neighbor->index() > cell->index())
                    {
                      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

                      ue_vi_matrix = 0;
                      ui_ve_matrix = 0;
                      ue_ve_matrix = 0;

                      fe_v_face.reinit (cell, face_no);
                      fe_v_face_neighbor.reinit (neighbor, neighbor2);

                      dg.assemble_face_term2(fe_v_face,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix,
                                             ui_ve_matrix,
                                             ue_ve_matrix);

                      neighbor->get_dof_indices (dofs_neighbor);

                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            system_matrix.add(dofs[i], dofs_neighbor[j],
                                              ue_vi_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs[j],
                                              ui_ve_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
                                              ue_ve_matrix(i,j));
                          }
                    }

                }
            }
        }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        right_hand_side(dofs[i]) += cell_vector(i);
    }
}


template <int dim>
void DGMethod<dim>::solve (Vector<double> &solution)
{
  SolverControl           solver_control (1000, 1e-12, false, false);
  SolverRichardson<>      solver (solver_control);

  PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;

  preconditioner.initialize(system_matrix, fe.dofs_per_cell);

  solver.solve (system_matrix, solution, right_hand_side,
                preconditioner);
}


template <int dim>
void DGMethod<dim>::refine_grid ()
{
  Vector<float> gradient_indicator (triangulation.n_active_cells());

  DerivativeApproximation::approximate_gradient (mapping,
                                                 dof_handler,
                                                 solution2,
                                                 gradient_indicator);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   gradient_indicator,
                                                   0.3, 0.1);

  triangulation.execute_coarsening_and_refinement ();
}


template <int dim>
void DGMethod<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());

  filename += ".eps";
  deallog << "Writing grid to <" << filename << ">..." << std::endl;

  GridOut grid_out;
  grid_out.write_eps (triangulation, deallog.get_file_stream());

  filename = "sol-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());

  filename += ".gnuplot";
  deallog << "Writing solution to <" << filename << ">..."
          << std::endl << std::endl;

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution2, "u");

  data_out.build_patches ();

  data_out.write_gnuplot(deallog.get_file_stream());
}


template <int dim>
void DGMethod<dim>::run ()
{
  for (unsigned int cycle=0; cycle<4; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube (triangulation);

          triangulation.refine_global (2);
        }
      else
        refine_grid ();


      deallog << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

      setup_system ();

      deallog << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

      Timer assemble_timer;
      assemble_system1 ();
      solve (solution1);

      system_matrix = 0;
      right_hand_side = 0;
      assemble_timer.reset();

      assemble_timer.start();
      assemble_system2 ();
      solve (solution2);

      solution1-=solution2;
      const double difference=solution1.linfty_norm();
      if (difference>1e-13)
        deallog << "solution1 and solution2 differ!!" << std::endl;
      else
        deallog << "solution1 and solution2 coincide." << std::endl;

      output_results (cycle);
    }
}

int main ()
{
  try
    {
      deallog << std::setprecision(2);
      logfile << std::setprecision(2);

      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      DGMethod<2> dgmethod;
      dgmethod.run ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
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
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}


