// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests matrix-free read_cell_data for locally refined mesh


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"


template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class Test
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  void
  run(unsigned int fe_degree)
  {
    Triangulation<dim> tria(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::hyper_cube(tria);

    for (unsigned int i = 0; i < 3; ++i)
      {
        tria.begin_active(i)->set_refine_flag();
        tria.execute_coarsening_and_refinement();
      }

    FE_DGQ<dim>     fe(fe_degree);
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    MappingQ<dim> mapping(1);

    QGauss<1> quad(fe_degree + 1);



    unsigned int level = numbers::invalid_unsigned_int;
    {
      AffineConstraints<Number> constraint;

      typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
        additional_data;
      additional_data.mapping_update_flags                = update_values;
      additional_data.mapping_update_flags_inner_faces    = update_values;
      additional_data.mapping_update_flags_boundary_faces = update_values;
      additional_data.mg_level                            = level;

      MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
      matrix_free.reinit(
        mapping, dof_handler, constraint, quad, additional_data);

      VectorType src, dst;

      matrix_free.initialize_dof_vector(src);
      matrix_free.initialize_dof_vector(dst);

      this->setup_vector(matrix_free);

      matrix_free.loop(&Test::cell_operation,
                       &Test::face_operation,
                       &Test::boundary_operation,
                       this,
                       dst,
                       src);

      deallog << std::endl;
    }
  }

private:
  void
  setup_vector(const MatrixFree<dim, Number, VectorizedArrayType> &data)
  {
    unsigned int n_cells = data.n_cell_batches() + data.n_ghost_cell_batches();
    cell_ids.resize(n_cells);

    FEEvaluation<dim, 1> fe_eval(data);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        fe_eval.reinit(cell);

        std::array<CellId, VectorizedArrayType::size()> cell_ids_local;
        for (auto lane = 0u; lane < data.n_active_entries_per_cell_batch(cell);
             lane++)
          cell_ids_local[lane] = data.get_cell_iterator(cell, lane)->id();

        fe_eval.set_cell_data(cell_ids, cell_ids_local);
      }
  }

  void
  cell_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                 VectorType &,
                 const VectorType &,
                 const std::pair<unsigned int, unsigned int> &pair) const
  {
    FEEvaluation<dim, 1> fe_eval(data);
    for (auto cell = pair.first; cell < pair.second; ++cell)
      {
        fe_eval.reinit(cell);
        const auto cell_data = fe_eval.read_cell_data(cell_ids);
        for (auto lane = 0u; lane < data.n_active_entries_per_cell_batch(cell);
             lane++)
          Assert(cell_data[lane] == data.get_cell_iterator(cell, lane)->id(),
                 ExcInternalError());
      }
  }

  void
  face_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                 VectorType &,
                 const VectorType &,
                 const std::pair<unsigned int, unsigned int> &pair) const
  {
    FEFaceEvaluation<dim, 1> fe_eval_m(data, true);
    FEFaceEvaluation<dim, 1> fe_eval_p(data, false);
    for (auto face = pair.first; face < pair.second; ++face)
      {
        fe_eval_m.reinit(face);
        fe_eval_p.reinit(face);
        const auto cell_data_m = fe_eval_m.read_cell_data(cell_ids);
        const auto cell_data_p = fe_eval_p.read_cell_data(cell_ids);
        for (auto lane = 0u; lane < data.n_active_entries_per_face_batch(face);
             lane++)
          {
            Assert(cell_data_m[lane] ==
                     data.get_face_iterator(face, lane, true).first->id(),
                   ExcInternalError());
            Assert(cell_data_p[lane] ==
                     data.get_face_iterator(face, lane, false).first->id(),
                   ExcInternalError());
          }
      }
  }


  void
  boundary_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                     VectorType &,
                     const VectorType &,
                     const std::pair<unsigned int, unsigned int> &pair) const
  {
    FEFaceEvaluation<dim, 1> fe_eval(data, true);
    for (auto face = pair.first; face < pair.second; ++face)
      {
        fe_eval.reinit(face);
        const auto cell_data = fe_eval.read_cell_data(cell_ids);
        for (auto lane = 0u; lane < data.n_active_entries_per_face_batch(face);
             lane++)
          Assert(cell_data[lane] ==
                   data.get_face_iterator(face, lane, true).first->id(),
                 ExcInternalError());
      }
  }

  AlignedVector<std::array<CellId, VectorizedArrayType::size()>> cell_ids;
};


int
main()
{
  initlog();
  {
    deallog.push("2d");
    Test<2> runner;
    runner.run(1);
    deallog.pop();
  }
  {
    deallog.push("3d");
    Test<3> runner;
    runner.run(1);
    deallog.pop();
  }
}
