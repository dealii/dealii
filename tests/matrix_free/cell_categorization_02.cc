// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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



// Tests the categorization of cells for vectorization


#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"

#include "create_mesh.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      // guarantee that the mesh also does not change by more than
      // refinement level across vertices that might connect two cells:
      Triangulation<dim>::limit_level_difference_at_vertices),
    typename parallel::distributed::Triangulation<dim>::Settings(
      // needed for GMG:
      parallel::distributed::Triangulation<
        dim>::construct_multigrid_hierarchy));

  // create mesh
  {
    std::vector<unsigned int> repetitions(dim, 2);
    if (dim == 3)
      repetitions[dim - 1] = 1;

    const Point<dim> bottom_left =
      (dim == 3 ? Point<dim>(0.0, 0.0, -0.5) : Point<dim>(0.0, 0.0));
    const Point<dim> top_right =
      (dim == 3 ? Point<dim>(48.0, 44.0, 0.5) : Point<dim>(48.0, 44.0));

    GridGenerator::subdivided_hyper_rectangle(tria,
                                              repetitions,
                                              bottom_left,
                                              top_right);

    for (auto &cell : tria.active_cell_iterators())
      {
        if (cell->center()[0] < 24.)
          cell->set_material_id(0);
        else
          cell->set_material_id(1);
      }

    tria.refine_global(4);
  }

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();
  AffineConstraints<double> constraints;
  constraints.close();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof);
  mg_constrained_dofs.make_zero_boundary_constraints(dof, {0});

  MatrixFree<dim>                          mf_data;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  data.cell_vectorization_categories_strict = true;
  data.cell_vectorization_category.resize(tria.n_active_cells());
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        AssertIndexRange(cell->active_cell_index(), tria.n_active_cells());
        data.cell_vectorization_category[cell->active_cell_index()] =
          cell->material_id();
      }

  mf_data.reinit(dof, constraints, QGauss<1>(2), data);

  const unsigned int max_level = tria.n_global_levels() - 1;
  std::vector<typename MatrixFree<dim, float>::AdditionalData>
    mg_additional_data(max_level + 1);
  for (unsigned int level = 0; level <= max_level; ++level)
    {
      mg_additional_data[level].tasks_parallel_scheme =
        MatrixFree<dim, float>::AdditionalData::none; // partition_color;

      mg_additional_data[level].mapping_update_flags =
        update_gradients | update_JxW_values;

      mg_additional_data[level].cell_vectorization_categories_strict = true;
      mg_additional_data[level].cell_vectorization_category.resize(
        tria.n_cells(level));
      for (const auto &cell : tria.cell_iterators_on_level(level))
        if (cell->is_locally_owned_on_level())
          {
            AssertIndexRange(cell->index(), tria.n_cells(level));
            mg_additional_data[level]
              .cell_vectorization_category[cell->index()] = cell->material_id();
          }

      mg_additional_data[level].mg_level = level;
    }

  std::vector<std::shared_ptr<MatrixFree<dim, float>>> mg_mf_data(max_level +
                                                                  1);

  for (unsigned int level = 0; level <= max_level; ++level)
    {
      AffineConstraints<double> level_constraints;
      IndexSet                  relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof, level, relevant_dofs);
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      mg_mf_data[level] = std::make_shared<MatrixFree<dim, float>>();

      mg_mf_data[level]->reinit(dof,
                                level_constraints,
                                QGauss<1>(2),
                                mg_additional_data[level]);
    }

  for (unsigned int i = 0; i < mf_data.n_macro_cells(); ++i)
    {
      const unsigned int m_id = mf_data.get_cell_iterator(i, 0)->material_id();
      for (unsigned int c = 0; c < mf_data.n_components_filled(i); ++c)
        {
          const unsigned int c_id =
            mf_data.get_cell_iterator(i, c)->material_id();
          AssertThrow(c_id == m_id, ExcInternalError());
        }
    }

  for (unsigned int level = 0; level <= max_level; ++level)
    {
      const auto &level_data = mg_mf_data[level];
      for (unsigned int i = 0; i < level_data->n_macro_cells(); ++i)
        {
          const unsigned int m_id =
            level_data->get_cell_iterator(i, 0)->material_id();
          for (unsigned int c = 0; c < level_data->n_components_filled(i); ++c)
            {
              const unsigned int c_id =
                level_data->get_cell_iterator(i, c)->material_id();
              AssertThrow(m_id == c_id, ExcInternalError());
            }
        }
    }
  deallog << Utilities::int_to_string(dim) << "d OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            testing_max_num_threads());
  MPILogInitAll                    log;
  test<2>();
  test<3>();

  return 0;
}
