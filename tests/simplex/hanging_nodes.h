// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// miscellaneous helper functions for hanging_nodes_* tests


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <vector>

#include "../tests.h"



// ----- setup -----

/**
 * On a Triangulation @p tria, refine those cells which are descendants of the
 * same coarse cell by a specified amount of times.
 *
 * The index of an entry in @p n_refinements corresponds to the globally unique
 * coarse cell id, while the entry itself describes the number of refinements.
 */
template <int dim>
void
refine(const std::vector<unsigned int> &n_refinements, Triangulation<dim> &tria)
{
  AssertDimension(n_refinements.size(), tria.n_cells(0));

  const unsigned int max_level =
    *std::max_element(n_refinements.begin(), n_refinements.end());

  while (tria.n_levels() <= max_level)
    {
      for (const auto &cell : tria.active_cell_iterators())
        {
          const unsigned int coarse_cell_id = cell->id().get_coarse_cell_id();

          if (static_cast<unsigned int>(cell->level()) <
              n_refinements[coarse_cell_id])
            cell->set_refine_flag();
        }

      tria.execute_coarsening_and_refinement();
    }

  if constexpr (running_in_debug_mode())
    {
      for (const auto &cell : tria.active_cell_iterators())
        {
          const unsigned int coarse_cell_id = cell->id().get_coarse_cell_id();

          AssertDimension(cell->level(), n_refinements[coarse_cell_id]);
        }
    }
}


/**
 * On a DoFHandler @p dofh, assign the same active FE index to all cells which
 * are descendants of the same coarse cell.
 *
 * The index of an entry in @p fe_indices corresponds to the globally unique
 * coarse cell id, while the entry itself is the active FE index.
 */
template <int dim>
void
set_active_fe_indices(const std::vector<unsigned int> &fe_indices,
                      DoFHandler<dim>                 &dofh)
{
  AssertDimension(fe_indices.size(), dofh.get_triangulation().n_cells(0));

  for (const auto &cell : dofh.active_cell_iterators())
    {
      const unsigned int coarse_cell_id = cell->id().get_coarse_cell_id();

      cell->set_active_fe_index(fe_indices[coarse_cell_id]);
    }
}



// ----- diagnostics -----

/**
 * Print the indices of degrees of freedom for each face on each cell in
 * deallog.
 *
 * If a face is divided into subfaces on a neighboring cell which is refined,
 * then also print the dof indices on the subfaces from the perspective of the
 * neighboring cell.
 */
template <int dim>
void
print_dof_indices_on_faces(const DoFHandler<dim> &dofh)
{
  std::vector<types::global_dof_index> dof_indices;

  for (const auto &cell : dofh.active_cell_iterators())
    {
      const auto        &fe       = cell->get_fe();
      const unsigned int fe_index = cell->active_fe_index();

      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        {
          const auto &face = cell->face(f);
          Assert(face->fe_index_is_active(fe_index), ExcInternalError());

          const unsigned int n_dofs =
            internal::DoFAccessorImplementation::Implementation::n_dof_indices(
              *face, fe_index);
          dof_indices.resize(n_dofs);
          face->get_dof_indices(dof_indices, fe_index);

          deallog << "cell:" << cell->active_cell_index() << " face:" << f
                  << " dofs:";
          for (const auto &i : dof_indices)
            deallog << i << " ";
          deallog << std::endl;

          // also print dofs on subfaces if neighboring cell is refined.
          // in this case, only one fe should be active on the subface.
          for (unsigned int sf = 0; sf < face->n_children(); ++sf)
            {
              const auto        &subface = face->child(sf);
              const unsigned int subface_fe_index =
                subface->nth_active_fe_index(0);
              Assert(subface->n_active_fe_indices() == 1, ExcInternalError());

              const unsigned int n_dofs = internal::DoFAccessorImplementation::
                Implementation::n_dof_indices(*subface, subface_fe_index);
              dof_indices.resize(n_dofs);
              subface->get_dof_indices(dof_indices, subface_fe_index);

              deallog << "    subface:" << sf << " dofs:";
              for (const auto &i : dof_indices)
                deallog << i << " ";
              deallog << std::endl;
            }
        }
    }
}


/**
 * Print the index and physical coordinate of each degree of freedom in deallog.
 */
template <int dim>
void
print_dof_points(const DoFHandler<dim> &dofh)
{
  hp::MappingCollection<dim> mapping;
  for (unsigned int i = 0; i < dofh.get_fe_collection().size(); ++i)
    mapping.push_back(MappingFE<dim>(dofh.get_fe(i)));

  std::vector<Point<dim>> points(dofh.n_dofs());
  DoFTools::map_dofs_to_support_points(mapping, dofh, points);

  for (unsigned int i = 0; i < dofh.n_dofs(); ++i)
    deallog << "dof:" << i << " point:" << points[i] << std::endl;
}



// ----- tests -----

/**
 * Verify hanging node constraints on locally refined meshes.
 *
 * Creates a Triangulation based on @p grid_generator, refines the coarse cells
 * @p n_refinements times, and assigns @p fe_indices to all descendants of the
 * coarse cells based on the provided @p fe_collection.
 *
 * Makes hanging node constraints and prints them in deallog.
 */
template <int dim>
void
test(const std::vector<unsigned int>                 &n_refinements,
     const std::vector<unsigned int>                 &fe_indices,
     const hp::FECollection<dim>                     &fe_collection,
     const std::function<void(Triangulation<dim> &)> &grid_generator)
{
  // setup grid
  Triangulation<dim> tria;
  grid_generator(tria);

  AssertDimension(n_refinements.size(), tria.n_cells());
  AssertDimension(fe_indices.size(), tria.n_cells());

  refine(n_refinements, tria);

#if false
  GridOut grid_out;
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif

  DoFHandler<dim> dofh(tria);
  set_active_fe_indices(fe_indices, dofh);
  dofh.distribute_dofs(fe_collection);
  deallog << "ndofs: " << dofh.n_dofs() << std::endl;

#if false
  print_dof_points(dofh);
  print_dof_indices_on_faces(dofh);
#endif

  // hanging node constraints
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.close();
  constraints.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}
