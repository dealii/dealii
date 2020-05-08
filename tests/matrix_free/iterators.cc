// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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



// tests matrix-free get_cell_iterator and get_face_iterator for
// a globally refined mesh


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <set>

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

    tria.refine_global(2);

    FE_DGQ<dim>     fe(fe_degree);
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    MappingQ<dim> mapping(1);

    QGauss<1> quad(fe_degree + 1);



    for (unsigned int level = numbers::invalid_unsigned_int;
         level < tria.n_global_levels() ||
         level == numbers::invalid_unsigned_int;
         level++)
      {
        if (level == numbers::invalid_unsigned_int)
          deallog << "Active level (with max_level="
                  << tria.n_global_levels() - 1 << "):" << std::endl;
        else
          deallog << "Multigrid level " << level << ":" << std::endl;

        AffineConstraints<Number> constraint;

        typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
          additional_data;
        additional_data.mapping_update_flags                = update_values;
        additional_data.mapping_update_flags_inner_faces    = update_values;
        additional_data.mapping_update_flags_boundary_faces = update_values;
        additional_data.mg_level                            = level;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none;

        MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
        matrix_free.reinit(
          mapping, dof_handler, constraint, quad, additional_data);

        VectorType src, dst;

        matrix_free.initialize_dof_vector(src);
        matrix_free.initialize_dof_vector(dst);

        matrix_free.loop(&Test::cell_operation,
                         &Test::face_operation,
                         &Test::boundary_operation,
                         this,
                         dst,
                         src);

        for (const auto cell : cells)
          deallog << "cell:      " << cell.to_string() << std::endl;

        for (const auto &face : faces)
          deallog << "face:      " << face.first.to_string() << " "
                  << face.second.to_string() << std::endl;

        for (const auto boundary : boundaries)
          deallog << "boundary:  " << boundary.to_string() << std::endl;

        cells.clear();
        faces.clear();
        boundaries.clear();

        deallog << std::endl;
      }
  }

private:
  void
  cell_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                 VectorType &,
                 const VectorType &,
                 const std::pair<unsigned int, unsigned int> &pair) const
  {
    for (auto cell = pair.first; cell < pair.second; cell++)
      for (auto lane = 0u; lane < data.n_active_entries_per_cell_batch(cell);
           lane++)
        cells.emplace(data.get_cell_iterator(cell, lane)->id());
  }

  void
  face_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                 VectorType &,
                 const VectorType &,
                 const std::pair<unsigned int, unsigned int> &pair) const
  {
    for (auto face = pair.first; face < pair.second; face++)
      for (auto lane = 0u; lane < data.n_active_entries_per_face_batch(face);
           lane++)
        {
          const auto face_0 =
            data.get_face_iterator(face, lane, true).first->id();
          const auto face_1 =
            data.get_face_iterator(face, lane, false).first->id();
          if (face_0 < face_1)
            faces.emplace(face_0, face_1);
          else
            faces.emplace(face_1, face_0);
        }
  }


  void
  boundary_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                     VectorType &,
                     const VectorType &,
                     const std::pair<unsigned int, unsigned int> &pair) const
  {
    for (auto face = pair.first; face < pair.second; face++)
      for (auto lane = 0u; lane < data.n_active_entries_per_face_batch(face);
           lane++)
        boundaries.emplace(
          data.get_face_iterator(face, lane, true).first->id());
  }

  mutable std::set<CellId>                    cells;
  mutable std::set<std::pair<CellId, CellId>> faces;
  mutable std::set<CellId>                    boundaries;
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
