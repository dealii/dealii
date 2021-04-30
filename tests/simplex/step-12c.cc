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



// a copy of tests/hp/step-12 testing different reference-cell types


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <iostream>

#include "../tests.h"

#include "./simplex_grids.h"

template <int dim>
struct ScratchData
{
  hp::MappingCollection<dim>            mapping;
  hp::FECollection<dim>                 fe;
  hp::QCollection<dim>                  quadrature;
  std::vector<hp::QCollection<dim - 1>> face_quadrature;

  std::function<void(Triangulation<dim> &)> mesh_generator;
};


template <int dim>
class RHS : public Function<dim>
{
public:
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double> &          values,
             const unsigned int             component = 0) const;
};


template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double> &          values,
             const unsigned int             component = 0) const;
};


template <int dim>
class Beta
{
public:
  Beta()
  {}
  void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<Point<dim>> &      values) const;
};


template <int dim>
void
RHS<dim>::value_list(const std::vector<Point<dim>> &points,
                     std::vector<double> &          values,
                     const unsigned int) const
{
  // Assert(values.size() == points.size(),
  //       ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < std::min(points.size(), values.size()); ++i)
    values[i] = 0;
}


template <int dim>
void
Beta<dim>::value_list(const std::vector<Point<dim>> &points,
                      std::vector<Point<dim>> &      values) const
{
  // Assert(values.size() == points.size(),
  //       ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < std::min(points.size(), values.size()); ++i)
    {
      const Point<dim> &p    = points[i];
      Point<dim> &      beta = values[i];

      beta(0) = -p(1);
      beta(1) = p(0);

      if (beta.norm() > 1e-10)
        beta /= std::sqrt(beta.square());
    }
}


template <int dim>
void
BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                std::vector<double> &          values,
                                const unsigned int) const
{
  // Assert(values.size() == points.size(),
  //       ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < std::min(values.size(), points.size()); ++i)
    {
      if (points[i](0) < 0.5)
        values[i] = 1.;
      else
        values[i] = 0.;
    }
}


template <int dim>
class DGTransportEquation
{
public:
  DGTransportEquation();

  void
  assemble_cell_term(const hp::FEValues<dim> &fe_v,
                     FullMatrix<double> &     ui_vi_matrix,
                     Vector<double> &         cell_vector) const;

  void
  assemble_boundary_term(const hp::FEFaceValues<dim> &fe_v,
                         FullMatrix<double> &         ui_vi_matrix,
                         Vector<double> &             cell_vector) const;

  template <class X, class Y>
  void
  assemble_face_term1(const X &           fe_v,
                      const Y &           fe_v_neighbor,
                      FullMatrix<double> &ui_vi_matrix,
                      FullMatrix<double> &ue_vi_matrix) const;

  template <class X, class Y>
  void
  assemble_face_term2(const X &           fe_v,
                      const Y &           fe_v_neighbor,
                      FullMatrix<double> &ui_vi_matrix,
                      FullMatrix<double> &ue_vi_matrix,
                      FullMatrix<double> &ui_ve_matrix,
                      FullMatrix<double> &ue_ve_matrix) const;

private:
  const Beta<dim>           beta_function;
  const RHS<dim>            rhs_function;
  const BoundaryValues<dim> boundary_function;
};


template <int dim>
DGTransportEquation<dim>::DGTransportEquation()
  : beta_function()
  , rhs_function()
  , boundary_function()
{}


template <int dim>
void
DGTransportEquation<dim>::assemble_cell_term(const hp::FEValues<dim> &fe_v,
                                             FullMatrix<double> &ui_vi_matrix,
                                             Vector<double> &cell_vector) const
{
  const std::vector<double> &JxW =
    fe_v.get_present_fe_values().get_JxW_values();

  std::vector<Point<dim>> beta(
    fe_v.get_present_fe_values().n_quadrature_points);
  std::vector<double> rhs(fe_v.get_present_fe_values().n_quadrature_points);

  beta_function.value_list(fe_v.get_present_fe_values().get_quadrature_points(),
                           beta);
  rhs_function.value_list(fe_v.get_present_fe_values().get_quadrature_points(),
                          rhs);

  for (unsigned int point = 0;
       point < fe_v.get_present_fe_values().n_quadrature_points;
       ++point)
    for (unsigned int i = 0; i < fe_v.get_present_fe_values().dofs_per_cell;
         ++i)
      {
        for (unsigned int j = 0; j < fe_v.get_present_fe_values().dofs_per_cell;
             ++j)
          ui_vi_matrix(i, j) -=
            beta[point] * fe_v.get_present_fe_values().shape_grad(i, point) *
            fe_v.get_present_fe_values().shape_value(j, point) * JxW[point];

        cell_vector(i) += rhs[point] *
                          fe_v.get_present_fe_values().shape_value(i, point) *
                          JxW[point];
      }
}


template <int dim>
void
DGTransportEquation<dim>::assemble_boundary_term(
  const hp::FEFaceValues<dim> &fe_v,
  FullMatrix<double> &         ui_vi_matrix,
  Vector<double> &             cell_vector) const
{
  const std::vector<double> &JxW =
    fe_v.get_present_fe_values().get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals =
    fe_v.get_present_fe_values().get_normal_vectors();

  std::vector<Point<dim>> beta(
    fe_v.get_present_fe_values().n_quadrature_points);
  std::vector<double> g(fe_v.get_present_fe_values().n_quadrature_points);

  beta_function.value_list(fe_v.get_present_fe_values().get_quadrature_points(),
                           beta);
  boundary_function.value_list(
    fe_v.get_present_fe_values().get_quadrature_points(), g);

  for (unsigned int point = 0;
       point < fe_v.get_present_fe_values().n_quadrature_points;
       ++point)
    {
      const double beta_n = beta[point] * normals[point];
      if (beta_n > 0)
        for (unsigned int i = 0; i < fe_v.get_present_fe_values().dofs_per_cell;
             ++i)
          for (unsigned int j = 0;
               j < fe_v.get_present_fe_values().dofs_per_cell;
               ++j)
            ui_vi_matrix(i, j) +=
              beta_n * fe_v.get_present_fe_values().shape_value(j, point) *
              fe_v.get_present_fe_values().shape_value(i, point) * JxW[point];
      else
        for (unsigned int i = 0; i < fe_v.get_present_fe_values().dofs_per_cell;
             ++i)
          cell_vector(i) -= beta_n * g[point] *
                            fe_v.get_present_fe_values().shape_value(i, point) *
                            JxW[point];
    }
}


template <int dim>
template <class X, class Y>
void
DGTransportEquation<dim>::assemble_face_term1(
  const X &           fe_v,
  const Y &           fe_v_neighbor,
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix) const
{
  const std::vector<double> &JxW =
    fe_v.get_present_fe_values().get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals =
    fe_v.get_present_fe_values().get_normal_vectors();

  std::vector<Point<dim>> beta(
    fe_v.get_present_fe_values().n_quadrature_points);
  beta_function.value_list(fe_v.get_present_fe_values().get_quadrature_points(),
                           beta);

  for (unsigned int point = 0;
       point < fe_v.get_present_fe_values().n_quadrature_points;
       ++point)
    {
      const double beta_n = beta[point] * normals[point];
      if (beta_n > 0)
        for (unsigned int i = 0; i < fe_v.get_present_fe_values().dofs_per_cell;
             ++i)
          for (unsigned int j = 0;
               j < fe_v.get_present_fe_values().dofs_per_cell;
               ++j)
            ui_vi_matrix(i, j) +=
              beta_n * fe_v.get_present_fe_values().shape_value(j, point) *
              fe_v.get_present_fe_values().shape_value(i, point) * JxW[point];
      else
        for (unsigned int i = 0; i < fe_v.get_present_fe_values().dofs_per_cell;
             ++i)
          for (unsigned int k = 0;
               k < fe_v_neighbor.get_present_fe_values().dofs_per_cell;
               ++k)
            ue_vi_matrix(i, k) +=
              beta_n *
              fe_v_neighbor.get_present_fe_values().shape_value(k, point) *
              fe_v.get_present_fe_values().shape_value(i, point) * JxW[point];
    }
}


template <int dim>
template <class X, class Y>
void
DGTransportEquation<dim>::assemble_face_term2(
  const X &           fe_v,
  const Y &           fe_v_neighbor,
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix,
  FullMatrix<double> &ui_ve_matrix,
  FullMatrix<double> &ue_ve_matrix) const
{
  const std::vector<double> &JxW =
    fe_v.get_present_fe_values().get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals =
    fe_v.get_present_fe_values().get_normal_vectors();

  std::vector<Point<dim>> beta(
    fe_v.get_present_fe_values().n_quadrature_points);

  beta_function.value_list(fe_v.get_present_fe_values().get_quadrature_points(),
                           beta);

  for (unsigned int point = 0;
       point < fe_v.get_present_fe_values().n_quadrature_points;
       ++point)
    {
      const double beta_n = beta[point] * normals[point];
      if (beta_n > 0)
        {
          for (unsigned int i = 0;
               i < fe_v.get_present_fe_values().dofs_per_cell;
               ++i)
            for (unsigned int j = 0;
                 j < fe_v.get_present_fe_values().dofs_per_cell;
                 ++j)
              ui_vi_matrix(i, j) +=
                beta_n * fe_v.get_present_fe_values().shape_value(j, point) *
                fe_v.get_present_fe_values().shape_value(i, point) * JxW[point];

          for (unsigned int k = 0;
               k < fe_v_neighbor.get_present_fe_values().dofs_per_cell;
               ++k)
            for (unsigned int j = 0;
                 j < fe_v.get_present_fe_values().dofs_per_cell;
                 ++j)
              ui_ve_matrix(k, j) -=
                beta_n * fe_v.get_present_fe_values().shape_value(j, point) *
                fe_v_neighbor.get_present_fe_values().shape_value(k, point) *
                JxW[point];
        }
      else
        {
          for (unsigned int i = 0;
               i < fe_v.get_present_fe_values().dofs_per_cell;
               ++i)
            for (unsigned int l = 0;
                 l < fe_v_neighbor.get_present_fe_values().dofs_per_cell;
                 ++l)
              ue_vi_matrix(i, l) +=
                beta_n *
                fe_v_neighbor.get_present_fe_values().shape_value(l, point) *
                fe_v.get_present_fe_values().shape_value(i, point) * JxW[point];

          for (unsigned int k = 0;
               k < fe_v_neighbor.get_present_fe_values().dofs_per_cell;
               ++k)
            for (unsigned int l = 0;
                 l < fe_v_neighbor.get_present_fe_values().dofs_per_cell;
                 ++l)
              ue_ve_matrix(k, l) -=
                beta_n *
                fe_v_neighbor.get_present_fe_values().shape_value(l, point) *
                fe_v_neighbor.get_present_fe_values().shape_value(k, point) *
                JxW[point];
        }
    }
}


template <int dim>
class DGMethod
{
public:
  DGMethod(const ScratchData<dim> &scratch_data);
  ~DGMethod();

  void
  run();

private:
  void
  setup_system();
  void
  assemble_system2();
  void
  solve(Vector<double> &solution);
  void
  refine_grid();
  void
  output_results(const unsigned int cycle) const;

  Triangulation<dim>               triangulation;
  const hp::MappingCollection<dim> mapping;

  hp::FECollection<dim> fe;
  DoFHandler<dim>       dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  const hp::QCollection<dim>                  quadrature;
  const std::vector<hp::QCollection<dim - 1>> face_quadrature;
  const hp::QCollection<dim - 1>              face_quadrature_dummy;

  Vector<double> solution1;
  Vector<double> solution2;
  Vector<double> right_hand_side;

  const std::function<void(Triangulation<dim> &)> mesh_generator;

  const DGTransportEquation<dim> dg;
};


template <int dim>
DGMethod<dim>::DGMethod(const ScratchData<dim> &scratch_data)
  : mapping(scratch_data.mapping)
  , fe(scratch_data.fe)
  , dof_handler(triangulation)
  , quadrature(scratch_data.quadrature)
  , face_quadrature(scratch_data.face_quadrature)
  , mesh_generator(scratch_data.mesh_generator)
  , dg()
{}


template <int dim>
DGMethod<dim>::~DGMethod()
{
  dof_handler.clear();
}


template <int dim>
void
DGMethod<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  sparsity_pattern.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          (GeometryInfo<dim>::faces_per_cell *
                             GeometryInfo<dim>::max_children_per_face +
                           1) *
                            fe[0].dofs_per_cell);

  DoFTools::make_flux_sparsity_pattern(dof_handler, sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit(sparsity_pattern);

  solution1.reinit(dof_handler.n_dofs());
  solution2.reinit(dof_handler.n_dofs());
  right_hand_side.reinit(dof_handler.n_dofs());
}



template <int dim>
void
DGMethod<dim>::assemble_system2()
{
  const unsigned int dofs_per_cell = dof_handler.get_fe(0).dofs_per_cell;
  std::vector<types::global_dof_index> dofs(dofs_per_cell);
  std::vector<types::global_dof_index> dofs_neighbor(dofs_per_cell);

  const UpdateFlags update_flags = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;

  const UpdateFlags face_update_flags =
    update_values | update_quadrature_points | update_JxW_values |
    update_normal_vectors;

  const UpdateFlags neighbor_face_update_flags = update_values;

  hp::FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags);
  hp::FEFaceValues<dim>    fe_v_face(mapping,
                                  fe,
                                  face_quadrature,
                                  face_update_flags);
  hp::FESubfaceValues<dim> fe_v_subface(mapping,
                                        fe,
                                        face_quadrature_dummy,
                                        face_update_flags);
  hp::FEFaceValues<dim>    fe_v_face_neighbor(mapping,
                                           fe,
                                           face_quadrature,
                                           neighbor_face_update_flags);


  FullMatrix<double> ui_vi_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_vi_matrix(dofs_per_cell, dofs_per_cell);

  FullMatrix<double> ui_ve_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_ve_matrix(dofs_per_cell, dofs_per_cell);

  Vector<double> cell_vector(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      ui_vi_matrix = 0;
      cell_vector  = 0;

      fe_v.reinit(cell);

      dg.assemble_cell_term(fe_v, ui_vi_matrix, cell_vector);

      cell->get_dof_indices(dofs);

      for (const unsigned int face_no : cell->face_indices())
        {
          typename DoFHandler<dim>::face_iterator face = cell->face(face_no);

          if (face->at_boundary())
            {
              fe_v_face.reinit(cell, face_no);

              dg.assemble_boundary_term(fe_v_face, ui_vi_matrix, cell_vector);
            }
          else
            {
              Assert(cell->neighbor(face_no).state() == IteratorState::valid,
                     ExcInternalError());
              typename DoFHandler<dim>::cell_iterator neighbor =
                cell->neighbor(face_no);
              if (face->has_children())
                {
                  const unsigned int neighbor2 =
                    cell->neighbor_of_neighbor(face_no);

                  for (unsigned int subface_no = 0;
                       subface_no < face->n_children();
                       ++subface_no)
                    {
                      typename DoFHandler<dim>::cell_iterator neighbor_child =
                        cell->neighbor_child_on_subface(face_no, subface_no);
                      Assert(neighbor_child->face(neighbor2) ==
                               face->child(subface_no),
                             ExcInternalError());
                      Assert(!neighbor_child->has_children(),
                             ExcInternalError());

                      ue_vi_matrix = 0;
                      ui_ve_matrix = 0;
                      ue_ve_matrix = 0;

                      fe_v_subface.reinit(cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit(neighbor_child, neighbor2);

                      dg.assemble_face_term2(fe_v_subface,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix,
                                             ui_ve_matrix,
                                             ue_ve_matrix);

                      neighbor_child->get_dof_indices(dofs_neighbor);

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          {
                            system_matrix.add(dofs[i],
                                              dofs_neighbor[j],
                                              ue_vi_matrix(i, j));
                            system_matrix.add(dofs_neighbor[i],
                                              dofs[j],
                                              ui_ve_matrix(i, j));
                            system_matrix.add(dofs_neighbor[i],
                                              dofs_neighbor[j],
                                              ue_ve_matrix(i, j));
                          }
                    }
                }
              else
                {
                  if (neighbor->level() == cell->level() &&
                      neighbor->index() > cell->index())
                    {
                      const unsigned int neighbor2 =
                        cell->neighbor_of_neighbor(face_no);

                      ue_vi_matrix = 0;
                      ui_ve_matrix = 0;
                      ue_ve_matrix = 0;

                      fe_v_face.reinit(cell, face_no);
                      fe_v_face_neighbor.reinit(neighbor, neighbor2);

                      dg.assemble_face_term2(fe_v_face,
                                             fe_v_face_neighbor,
                                             ui_vi_matrix,
                                             ue_vi_matrix,
                                             ui_ve_matrix,
                                             ue_ve_matrix);

                      neighbor->get_dof_indices(dofs_neighbor);

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          {
                            system_matrix.add(dofs[i],
                                              dofs_neighbor[j],
                                              ue_vi_matrix(i, j));
                            system_matrix.add(dofs_neighbor[i],
                                              dofs[j],
                                              ui_ve_matrix(i, j));
                            system_matrix.add(dofs_neighbor[i],
                                              dofs_neighbor[j],
                                              ue_ve_matrix(i, j));
                          }
                    }
                }
            }
        }

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i, j));

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        right_hand_side(dofs[i]) += cell_vector(i);
    }
}



template <int dim>
class Operator
{
public:
  using VectorType = Vector<double>;
  using number     = double;

  Operator(const MatrixFree<dim, double> &matrix_free)
    : matrix_free(matrix_free)
  {}

  void
  initialize_dof_vector(VectorType &vec)
  {
    matrix_free.initialize_dof_vector(vec);
  }

  Tensor<1, dim, VectorizedArray<double>>
  beta(const Point<dim, VectorizedArray<double>> &points) const
  {
    Tensor<1, dim, VectorizedArray<double>> betas;

    for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
      {
        Point<dim, double>     p;
        Tensor<1, dim, double> beta;

        for (int d = 0; d < dim; ++d)
          p[d] = points[d][v];

        beta[0] = -p[1];
        beta[1] = +p[0];

        if (beta.norm() > 1e-10)
          beta /= beta.norm();

        for (int d = 0; d < dim; ++d)
          betas[d][v] = beta[d];
      }

    return betas;
  }

  VectorizedArray<double>
  boundary_values(const Point<dim, VectorizedArray<double>> &points) const
  {
    VectorizedArray<double> betas;

    for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
      betas[v] = (points[0][v] < 0.5) ? 1.0 : 0.0;

    return betas;
  }

  void
  rhs(VectorType &vec) const
  {
    const int dummy = 0;

    matrix_free.template loop<VectorType, int>(
      [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
        FEEvaluation<dim, -1, 0, 1, double> phi(data, cell_range);
        for (unsigned int cell = cell_range.first; cell < cell_range.second;
             ++cell)
          {
            phi.reinit(cell);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(0.0, q);

            phi.integrate_scatter(true, false, dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        (void)data;
        (void)dst;
        (void)src;
        (void)face_range;
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, double> phi(data, face_range);
        for (unsigned int face = face_range.first; face < face_range.second;
             ++face)
          {
            phi.reinit(face);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                const auto beta_n = this->beta(phi.quadrature_point(q)) *
                                    phi.get_normal_vector(q);
                const auto beta_n_m = (-std::abs(beta_n) + beta_n) * 0.5;
                phi.submit_value(-beta_n_m * this->boundary_values(
                                               phi.quadrature_point(q)),
                                 q);
              }

            phi.integrate_scatter(true, false, dst);
          }
      },
      vec,
      dummy,
      true);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.template loop<VectorType, VectorType>(
      [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
        FEEvaluation<dim, -1, 0, 1, double> phi(data, cell_range);
        for (unsigned int cell = cell_range.first; cell < cell_range.second;
             ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, true, false);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(-this->beta(phi.quadrature_point(q)) *
                                    phi.get_value(q),
                                  q);
            phi.integrate_scatter(false, true, dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, double> phi_m(data, face_range, true);
        FEFaceEvaluation<dim, -1, 0, 1, double> phi_p(data, face_range, false);
        for (unsigned int cell = face_range.first; cell < face_range.second;
             ++cell)
          {
            phi_m.reinit(cell);
            phi_m.gather_evaluate(src, true, false);
            phi_p.reinit(cell);
            phi_p.gather_evaluate(src, true, false);
            for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
              {
                const auto beta_n = this->beta(phi_m.quadrature_point(q)) *
                                    phi_m.get_normal_vector(q);

                const auto u_m = phi_m.get_value(q);
                const auto u_p = phi_p.get_value(q);
                const auto temp =
                  0.5 * ((u_m + u_p) * beta_n + (u_m - u_p) * std::abs(beta_n));

                phi_m.submit_value(+temp, q);
                phi_p.submit_value(-temp, q);
              }

            phi_m.integrate_scatter(true, false, dst);
            phi_p.integrate_scatter(true, false, dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, double> phi(data, face_range);
        for (unsigned int cell = face_range.first; cell < face_range.second;
             ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, true, false);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                const auto beta_n = this->beta(phi.quadrature_point(q)) *
                                    phi.get_normal_vector(q);
                const auto beta_n_p = (std::abs(beta_n) + beta_n) * 0.5;
                phi.submit_value(beta_n_p * phi.get_value(q), q);
              }

            phi.integrate_scatter(true, false, dst);
          }
      },
      dst,
      src,
      true);
  }

private:
  const MatrixFree<dim, double> &matrix_free;
};


template <int dim>
void
DGMethod<dim>::solve(Vector<double> &solution)
{
  ReductionControl   solver_control(10000, 1e-4, 1e-3, false, false);
  SolverRichardson<> solver(solver_control);

  PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, fe[0].dofs_per_cell);

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_gradients | update_values | update_quadrature_points;
  additional_data.mapping_update_flags_inner_faces =
    update_gradients | update_values | update_quadrature_points;
  additional_data.mapping_update_flags_boundary_faces =
    update_gradients | update_values | update_quadrature_points;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::none;

  AffineConstraints<double> constraints;
  constraints.close();

  MatrixFree<dim, double> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, quadrature, additional_data);

  Operator<dim> op(matrix_free);

  op.initialize_dof_vector(solution);
  op.initialize_dof_vector(right_hand_side);

  op.rhs(right_hand_side);

  deallog << "  " << right_hand_side.l2_norm() << std::endl;
  solver.solve(op, solution, right_hand_side, preconditioner);
  deallog << "  " << solution.l2_norm() << std::endl;
}


template <int dim>
void
DGMethod<dim>::refine_grid()
{
  Vector<float> gradient_indicator(triangulation.n_active_cells());

  DerivativeApproximation::approximate_gradient(mapping[0],
                                                dof_handler,
                                                solution2,
                                                gradient_indicator);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int cell_no = 0; cell != endc; ++cell, ++cell_no)
    gradient_indicator(cell_no) *=
      std::pow(cell->diameter(), 1 + 1.0 * dim / 2);

  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  gradient_indicator,
                                                  0.3,
                                                  0.1);

  triangulation.execute_coarsening_and_refinement();
}


template <int dim>
void
DGMethod<dim>::output_results(const unsigned int cycle) const
{
  return;

  std::string filename = "grid-";
  filename += ('0' + cycle);
  Assert(cycle < 10, ExcInternalError());

  filename += ".eps";
  deallog << "Writing grid to <" << filename << ">..." << std::endl;

  GridOut grid_out;
  // grid_out.write_eps(triangulation, deallog.get_file_stream());

  filename = "sol-";
  filename += ('0' + cycle);
  Assert(cycle < 10, ExcInternalError());

  filename += ".gnuplot";
  deallog << "Writing solution to <" << filename << ">..." << std::endl
          << std::endl;

  DataOut<dim, DoFHandler<dim>> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution1, "u");

  data_out.build_patches(mapping, 2);

#if false
  data_out.write_vtk(deallog.get_file_stream());
#else
  static unsigned int counter = 0;

  std::ofstream output("test." + std::to_string(dim) + "." +
                       std::to_string(counter++) + ".vtk");
  data_out.write_vtk(output);
#endif
}


template <int dim>
void
DGMethod<dim>::run()
{
  deallog << "Cycle " << 0 << ':' << std::endl;
  mesh_generator(triangulation);
  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl;

  setup_system();
  assemble_system2();

  deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;
  solve(solution1);

  output_results(0);
}

int
main()
{
  try
    {
      initlog();
      deallog.get_file_stream().precision(2);

      // triangle
      for (unsigned int i = 1; i <= 2; ++i)
        {
          ScratchData<2> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<2>(MappingFE<2>(FE_SimplexP<2>(1)));
          scratch_data.fe         = hp::FECollection<2>(FE_SimplexDGP<2>(i));
          scratch_data.quadrature = hp::QCollection<2>(QGaussSimplex<2>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<1>>{
            hp::QCollection<1>(QGaussSimplex<1>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<2> &triangulation) -> void {
            GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                                16);
          };

          DGMethod<2> dgmethod(scratch_data);
          dgmethod.run();
        }

      // quadrilateral
      for (unsigned int i = 1; i <= 2; ++i)
        {
          ScratchData<2> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<2>(MappingQGeneric<2>(1));
          scratch_data.fe              = hp::FECollection<2>(FE_DGQ<2>(i));
          scratch_data.quadrature      = hp::QCollection<2>(QGauss<2>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<1>>{
            hp::QCollection<1>(QGauss<1>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<2> &triangulation) -> void {
            GridGenerator::hyper_cube(triangulation);
            triangulation.refine_global(3);
          };

          DGMethod<2> dgmethod(scratch_data);
          dgmethod.run();
        }

      // tetrahedron
      for (unsigned int i = 1; i <= 2; ++i)
        {
          ScratchData<3> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<3>(MappingFE<3>(FE_SimplexP<3>(1)));
          scratch_data.fe         = hp::FECollection<3>(FE_SimplexDGP<3>(i));
          scratch_data.quadrature = hp::QCollection<3>(QGaussSimplex<3>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<2>>{
            hp::QCollection<2>(QGaussSimplex<2>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<3> &triangulation) -> void {
            GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                                8);
          };

          DGMethod<3> dgmethod(scratch_data);
          dgmethod.run();
        }

      // pyramid
      for (unsigned int i = 1; i <= 1; ++i)
        {
          ScratchData<3> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<3>(MappingFE<3>(FE_PyramidP<3>(1)));
          scratch_data.fe         = hp::FECollection<3>(FE_PyramidDGP<3>(i));
          scratch_data.quadrature = hp::QCollection<3>(QGaussPyramid<3>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<2>>{
            hp::QCollection<2>(QGauss<2>(i + 1),
                               QGaussSimplex<2>(i + 1),
                               QGaussSimplex<2>(i + 1),
                               QGaussSimplex<2>(i + 1),
                               QGaussSimplex<2>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<3> &triangulation) -> void {
            GridGenerator::subdivided_hyper_cube_with_pyramids(triangulation,
                                                               8);
          };

          DGMethod<3> dgmethod(scratch_data);
          dgmethod.run();
        }

      // wedge
      for (unsigned int i = 1; i <= 2; ++i)
        {
          ScratchData<3> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<3>(MappingFE<3>(FE_WedgeP<3>(1)));
          scratch_data.fe         = hp::FECollection<3>(FE_WedgeDGP<3>(i));
          scratch_data.quadrature = hp::QCollection<3>(QGaussWedge<3>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<2>>{
            hp::QCollection<2>(QGaussSimplex<2>(i + 1),
                               QGaussSimplex<2>(i + 1),
                               QGauss<2>(i + 1),
                               QGauss<2>(i + 1),
                               QGauss<2>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<3> &triangulation) -> void {
            GridGenerator::subdivided_hyper_cube_with_wedges(triangulation, 8);
          };

          DGMethod<3> dgmethod(scratch_data);
          dgmethod.run();
        }

      // hexahedron
      for (unsigned int i = 1; i <= 2; ++i)
        {
          ScratchData<3> scratch_data;

          scratch_data.mapping =
            hp::MappingCollection<3>(MappingQGeneric<3>(1));
          scratch_data.fe              = hp::FECollection<3>(FE_DGQ<3>(i));
          scratch_data.quadrature      = hp::QCollection<3>(QGauss<3>(i + 1));
          scratch_data.face_quadrature = std::vector<hp::QCollection<2>>{
            hp::QCollection<2>(QGauss<2>(i + 1))};
          scratch_data.mesh_generator =
            [](Triangulation<3> &triangulation) -> void {
            GridGenerator::subdivided_hyper_cube(triangulation, 8);
          };

          DGMethod<3> dgmethod(scratch_data);
          dgmethod.run();
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

  return 0;
}
