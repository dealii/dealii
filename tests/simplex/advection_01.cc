// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Solve advection problem on an adaptive mixed mesh with DG.
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "simplex_grids.h"


using namespace dealii;

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues() = default;
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &values,
             const unsigned int             component = 0) const override;
};

template <int dim>
void
BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                std::vector<double>           &values,
                                const unsigned int             component) const
{
  (void)component;
  AssertIndexRange(component, 1);
  AssertDimension(values.size(), points.size());

  for (unsigned int i = 0; i < values.size(); ++i)
    {
      if (points[i][0] < 2. / 3)
        values[i] = 1.;
      else
        values[i] = 0.;
    }
}


template <int dim>
Tensor<1, dim>
beta(const Point<dim> &p)
{
  Assert(dim >= 2, ExcNotImplemented());

  Tensor<1, dim> wind_field;
  wind_field[0] = -p[1];
  wind_field[1] = p[0];

  if (wind_field.norm() > 1e-10)
    wind_field /= wind_field.norm();

  return wind_field;
}


template <int dim>
struct ScratchData
{
  ScratchData(const hp::MappingCollection<dim> &mapping,
              const hp::FECollection<dim>      &fe,
              const hp::QCollection<dim>       &quadrature,
              const hp::QCollection<dim - 1>   &quadrature_face,
              const UpdateFlags                 update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags interface_update_flags =
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values | update_normal_vectors)
    : fe_values(mapping, fe, quadrature, update_flags)
    , fe_interface_values(mapping, fe, quadrature_face, interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping_collection(),
                scratch_data.fe_values.get_fe_collection(),
                scratch_data.fe_values.get_quadrature_collection(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(
        scratch_data.fe_interface_values.get_mapping_collection(),
        scratch_data.fe_interface_values.get_fe_collection(),
        scratch_data.fe_interface_values.get_quadrature_collection(),
        scratch_data.fe_interface_values.get_update_flags())
  {}

  hp::FEValues<dim>      fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};



struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};



struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.clear();
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(0);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.clear();
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);

    face_data.clear();
  }
};

template <int dim>
class AdvectionProblem
{
public:
  AdvectionProblem();

  void
  run();

private:
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  refine_grid();
  void
  output_results(const unsigned int cycle);

  Triangulation<dim> triangulation;


  const hp::MappingCollection<dim> mapping;

  const hp::FECollection<dim> fe;

  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;

  const hp::QCollection<dim>     quadrature;
  const hp::QCollection<dim - 1> quadrature_face;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  SparseMatrix<double> system_matrix;
  SparsityPattern      sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;
};


template <int dim>
AdvectionProblem<dim>::AdvectionProblem()
  : triangulation(dealii::Triangulation<
                  dim>::MeshSmoothing::limit_level_difference_at_vertices)
  , mapping(MappingFE<dim>(FE_SimplexDGP<dim>(1)), MappingQ1<dim>())
  , fe(FE_SimplexDGP<dim>(3), FE_DGQ<dim>(3))
  , dof_handler(triangulation)
  , quadrature(QGaussSimplex<dim>(fe[0].degree + 1),
               QGauss<dim>(fe[1].degree + 1))
  , quadrature_face(QGaussSimplex<dim - 1>(fe[0].degree + 1),
                    QGauss<dim - 1>(fe[1].degree + 1))
{}


template <int dim>
void
AdvectionProblem<dim>::setup_system()
{
  dof_handler.clear();
  dof_handler.reinit(triangulation);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->reference_cell().is_hyper_cube())
      cell->set_active_fe_index(1);
    else if (cell->reference_cell().is_simplex())
      cell->set_active_fe_index(0);
    else
      DEAL_II_ASSERT_UNREACHABLE();


  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());

  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
}

template <int dim>
void
AdvectionProblem<dim>::assemble_system()
{
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  const BoundaryValues<dim> boundary_function;

  // This is the function that will be executed for each cell.
  const auto cell_worker = [&](const Iterator   &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData         &copy_data) {
    scratch_data.fe_values.reinit(cell);
    const auto &fe_v = scratch_data.fe_values.get_present_fe_values();

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();

    copy_data.reinit(cell, n_dofs);

    const auto &q_points = fe_v.get_quadrature_points();

    const std::vector<double> &JxW = fe_v.get_JxW_values();

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
      {
        auto beta_q = beta(q_points[point]);
        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              copy_data.cell_matrix(i, j) +=
                -beta_q                      // -\beta
                * fe_v.shape_grad(i, point)  // \nabla \phi_i
                * fe_v.shape_value(j, point) // \phi_j
                * JxW[point];                // dx
            }
      }
  };

  const auto boundary_worker = [&](const Iterator     &cell,
                                   const unsigned int &face_no,
                                   ScratchData<dim>   &scratch_data,
                                   CopyData           &copy_data) {
    const unsigned int fe_index    = cell->active_fe_index();
    const unsigned int other_index = dealii::numbers::invalid_unsigned_int;

    scratch_data.fe_interface_values.reinit(
      cell, face_no, other_index, other_index, fe_index);
    const FEFaceValuesBase<dim> &fe_face =
      scratch_data.fe_interface_values.get_fe_face_values(0);

    const auto &q_points = fe_face.get_quadrature_points();

    const unsigned int n_facet_dofs        = fe_face.get_fe().n_dofs_per_cell();
    const std::vector<double>         &JxW = fe_face.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

    std::vector<double> g(q_points.size());
    boundary_function.value_list(q_points, g);

    for (unsigned int point = 0; point < q_points.size(); ++point)
      {
        const double beta_dot_n = beta(q_points[point]) * normals[point];

        if (beta_dot_n > 0)
          {
            for (unsigned int i = 0; i < n_facet_dofs; ++i)
              for (unsigned int j = 0; j < n_facet_dofs; ++j)
                copy_data.cell_matrix(i, j) +=
                  fe_face.shape_value(i, point)   // \phi_i
                  * fe_face.shape_value(j, point) // \phi_j
                  * beta_dot_n                    // \beta . n
                  * JxW[point];                   // dx
          }
        else
          for (unsigned int i = 0; i < n_facet_dofs; ++i)
            copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i
                                     * g[point]                     // g
                                     * beta_dot_n                   // \beta . n
                                     * JxW[point];                  // dx
      }
  };

  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &f,
                               const unsigned int &sf,
                               const Iterator     &ncell,
                               const unsigned int &nf,
                               const unsigned int &nsf,
                               ScratchData<dim>   &scratch_data,
                               CopyData           &copy_data) {
    const unsigned int fe_index          = cell->active_fe_index();
    const unsigned int fe_neighbor_index = ncell->active_fe_index();


    FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;

    typename FEInterfaceValues<dim>::InterfaceData data(f, sf, nf, nsf);
    data.fe_index               = fe_index;
    data.q_index                = fe_index;
    data.mapping_index          = fe_index;
    data.fe_index_neighbor      = fe_neighbor_index;
    data.q_index_neighbor       = fe_neighbor_index;
    data.mapping_index_neighbor = fe_neighbor_index;

    fe_iv.reinit(cell, ncell, data);

    const auto &q_points = fe_iv.get_quadrature_points();

    copy_data.face_data.emplace_back();
    CopyDataFace &copy_data_face = copy_data.face_data.back();

    const unsigned int n_dofs        = fe_iv.n_current_interface_dofs();
    copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

    copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

    const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      {
        const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint];
        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int j = 0; j < n_dofs; ++j)
            copy_data_face.cell_matrix(i, j) +=
              fe_iv.jump_in_shape_values(i, qpoint) // [\phi_i]
              * fe_iv.shape_value((beta_dot_n > 0),
                                  j,
                                  qpoint) // phi_j^{upwind}
              * beta_dot_n                // (\beta . n)
              * JxW[qpoint];              // dx
      }
  };


  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.cell_rhs,
                                           c.local_dof_indices,
                                           system_matrix,
                                           system_rhs);

    for (const auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_matrix,
                                               cdf.joint_dof_indices,
                                               system_matrix);
      }
  };

  ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData         copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_ghost_faces_once |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}


template <int dim>
void
AdvectionProblem<dim>::solve()
{
  SolverControl solver_control(dof_handler.n_dofs(),
                               1e-4 * system_rhs.l2_norm());

  SolverGMRES solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}


template <int dim>
void
AdvectionProblem<dim>::refine_grid()
{
  const double radius = 2. / 3.;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      bool one_inside  = false;
      bool one_outside = false;

      for (unsigned int v = 0; v < cell->n_vertices(); ++v)
        {
          auto vertex = cell->vertex(v);
          if (vertex.norm() > radius)
            one_outside = true;
          if (vertex.norm() < radius)
            one_inside = true;
        }

      if (one_inside && one_outside)
        cell->set_refine_flag();
    }

  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
AdvectionProblem<dim>::output_results(const unsigned int cycle)
{
  deallog << "Cycle: " << cycle << " "
          << "L2 norm: " << solution.l2_norm() << std::endl;
}

template <int dim>
void
AdvectionProblem<dim>::run()
{
  const unsigned int n_cycles = 3;
  for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
    {
      deallog << "Cycle " << cycle << std::endl;

      if (cycle == 0)
        {
          // GridGenerator::subdivided_hyper_cube_with_simplices_mix(triangulation,
          // 2); triangulation.refine_global(2);
          std::vector<Point<dim>>    vertices;
          std::vector<CellData<dim>> cells;
          vertices.emplace_back(dealii::Point<dim>(0., 0.));
          vertices.emplace_back(dealii::Point<dim>(0.5, 0.));
          vertices.emplace_back(dealii::Point<dim>(1., 0.));
          vertices.emplace_back(dealii::Point<dim>(0., 0.5));
          vertices.emplace_back(dealii::Point<dim>(0.5, 0.5));
          vertices.emplace_back(dealii::Point<dim>(1., 0.5));
          vertices.emplace_back(dealii::Point<dim>(0., 1.));
          vertices.emplace_back(dealii::Point<dim>(0.5, 1.));
          vertices.emplace_back(dealii::Point<dim>(1., 1.));


          CellData<dim> tri1;
          tri1.vertices = {0, 1, 4};
          cells.push_back(tri1);

          CellData<dim> tri2;
          tri2.vertices = {0, 4, 3};
          cells.push_back(tri2);

          CellData<dim> quad1;
          quad1.vertices = {1, 2, 4, 5};
          cells.push_back(quad1);

          CellData<dim> quad2;
          quad2.vertices = {3, 4, 6, 7};
          cells.push_back(quad2);

          CellData<dim> tri3;
          tri3.vertices = {4, 5, 8};
          cells.push_back(tri3);

          CellData<dim> tri4;
          tri4.vertices = {4, 8, 7};
          cells.push_back(tri4);

          triangulation.create_triangulation(vertices, cells, SubCellData());
          triangulation.refine_global(1);
        }
      else
        refine_grid();

      deallog << "  Number of active cells:       "
              << triangulation.n_active_cells() << std::endl;

      setup_system();

      deallog << "  Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

      assemble_system();
      solve();

      output_results(cycle);
    }
  deallog << std::endl;
}

int
main()
{
  using namespace dealii;
  initlog();
  AdvectionProblem<2> advection_problem_2d;
  advection_problem_2d.run();

  return 0;
}
