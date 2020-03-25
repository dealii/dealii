// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// tests matrix-free read_cell_data for each cell (incl. inner and outer faces,
// which is needed for element-centric-loops)


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

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
    parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tria);
    tria.refine_global(2);

    FE_DGQ<dim>     fe(fe_degree);
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    MappingQ<dim> mapping(1);

    QGauss<1> quad(fe_degree + 1);

    AffineConstraints<Number> constraint;

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags                = update_values;
    additional_data.mapping_update_flags_inner_faces    = update_values;
    additional_data.mapping_update_flags_boundary_faces = update_values;
    additional_data.mapping_update_flags_faces_by_cells = update_values;
    additional_data.hold_all_faces_to_owned_cells       = true;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

    VectorType src, dst;

    matrix_free.initialize_dof_vector(src);
    matrix_free.initialize_dof_vector(dst);

    this->setup_vector(matrix_free);

    matrix_free.cell_loop(&Test::cell_wise_operation, this, dst, src);

    deallog << std::endl;
  }

private:
  void
  setup_vector(const MatrixFree<dim, Number, VectorizedArrayType> &data)
  {
    unsigned int n_cells = data.n_cell_batches() + data.n_ghost_cell_batches();
    cell_ids.resize(n_cells);

    for (unsigned int cell = 0; cell < n_cells; cell++)
      for (auto lane = 0u; lane < data.n_active_entries_per_cell_batch(cell);
           lane++)
        cell_ids[cell][lane] = data.get_cell_iterator(cell, lane)->id();
  }

  void
  cell_wise_operation(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                      VectorType &,
                      const VectorType &,
                      const std::pair<unsigned int, unsigned int> &pair) const
  {
    FEEvaluation<dim, 1>     fe_eval(data);
    FEFaceEvaluation<dim, 1> fe_eval_m(data, true);
    FEFaceEvaluation<dim, 1> fe_eval_p(data, false);

    for (auto cell = pair.first; cell < pair.second; cell++)
      {
        for (auto lane = 0u; lane < data.n_active_entries_per_cell_batch(cell);
             lane++)
          {
            fe_eval.reinit(cell);
            const auto cell_data = fe_eval.read_cell_data(cell_ids);
            Assert(cell_data[lane] == data.get_cell_iterator(cell, lane)->id(),
                   ExcInternalError());
            // deallog << "c " << cell_data[lane] << std::endl;

            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              {
                fe_eval_m.reinit(cell, face);
                fe_eval_p.reinit(cell, face);

                const auto cell_data_m = fe_eval_m.read_cell_data(cell_ids);
                Assert(cell_data_m[lane] ==
                         data.get_cell_iterator(cell, lane)->id(),
                       ExcInternalError());
                // deallog << "m" << cell_data_m[lane] << std::endl;

                const auto bids =
                  data.get_faces_by_cells_boundary_id(cell, face);

                if (bids[lane] != numbers::invalid_unsigned_int)
                  continue;

                const auto cell_data_p = fe_eval_p.read_cell_data(cell_ids);

                Assert(
                  cell_data_p[lane] ==
                    data.get_cell_iterator(cell, lane)->neighbor(face)->id(),
                  ExcInternalError());
                // deallog << "p " << cell_data_p[lane] << std::endl;
              }
          }
      }
  }

  AlignedVector<std::array<CellId, VectorizedArrayType::size()>> cell_ids;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  MPILogInitAll all;

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
