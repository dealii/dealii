// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Agglomerate some cells in a grid, and create a finite element space on the
// bounding box of an agglomeration. To check the correctness, compute the area
// of the agglomerated cells using the weights of a custom quadrature rule over
// the agglomerated element.

#include <deal.II/dofs/agglomeration_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/poly_utils.h>

using namespace dealii;

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  MappingQ<dim> mapping(1);
  tria.refine_global(3);
  GridTools::Cache<dim>     cached_tria(tria, mapping);
  AgglomerationHandler<dim> ah(cached_tria);

  std::vector<typename Triangulation<dim>::active_cell_iterator>
    cells; // each cell = an agglomerate
  for (const auto &cell : tria.active_cell_iterators())
    cells.push_back(cell);

  std::vector<types::global_cell_index> flagged_cells;
  const auto                            store_flagged_cells =
    [&flagged_cells](
      const std::vector<types::global_cell_index> &idxs_to_be_agglomerated) {
      for (const int idx : idxs_to_be_agglomerated)
        flagged_cells.push_back(idx);
    };


  if constexpr (dim == 2)
    {
      std::vector<types::global_cell_index> idxs_to_be_agglomerated = {
        3, 6, 9, 12, 13};
      store_flagged_cells(idxs_to_be_agglomerated);

      std::vector<typename Triangulation<dim>::active_cell_iterator>
        cells_to_be_agglomerated;
      PolyUtils::collect_cells_for_agglomeration(tria,
                                                 idxs_to_be_agglomerated,
                                                 cells_to_be_agglomerated);

      std::vector<types::global_cell_index> idxs_to_be_agglomerated2 = {15,
                                                                        36,
                                                                        37};
      store_flagged_cells(idxs_to_be_agglomerated2);

      std::vector<typename Triangulation<dim>::active_cell_iterator>
        cells_to_be_agglomerated2;
      PolyUtils::collect_cells_for_agglomeration(tria,
                                                 idxs_to_be_agglomerated2,
                                                 cells_to_be_agglomerated2);

      std::vector<types::global_cell_index> idxs_to_be_agglomerated3 = {57,
                                                                        60,
                                                                        54};
      store_flagged_cells(idxs_to_be_agglomerated3);
      std::vector<typename Triangulation<dim>::active_cell_iterator>
        cells_to_be_agglomerated3;
      PolyUtils::collect_cells_for_agglomeration(tria,
                                                 idxs_to_be_agglomerated3,
                                                 cells_to_be_agglomerated3);

      std::vector<types::global_cell_index> idxs_to_be_agglomerated4 = {25,
                                                                        19,
                                                                        22};
      store_flagged_cells(idxs_to_be_agglomerated4);

      std::vector<typename Triangulation<dim>::active_cell_iterator>
        cells_to_be_agglomerated4;
      PolyUtils::collect_cells_for_agglomeration(tria,
                                                 idxs_to_be_agglomerated4,
                                                 cells_to_be_agglomerated4);

      // Agglomerate the cells just stored
      ah.define_agglomerate(cells_to_be_agglomerated);
      ah.define_agglomerate(cells_to_be_agglomerated2);
      ah.define_agglomerate(cells_to_be_agglomerated3);
      ah.define_agglomerate(cells_to_be_agglomerated4);
    }
  else if constexpr (dim == 3)
    {
      std::vector<types::global_cell_index> idxs_to_be_agglomerated = {463,
                                                                       459};
      store_flagged_cells(idxs_to_be_agglomerated);

      std::vector<typename Triangulation<dim>::active_cell_iterator>
        cells_to_be_agglomerated;
      PolyUtils::collect_cells_for_agglomeration(tria,
                                                 idxs_to_be_agglomerated,
                                                 cells_to_be_agglomerated);

      ah.define_agglomerate(cells_to_be_agglomerated);
    }

  for (std::size_t i = 0; i < tria.n_active_cells(); ++i)
    {
      // If not present, agglomerate all the singletons
      if (std::find(flagged_cells.begin(),
                    flagged_cells.end(),
                    cells[i]->active_cell_index()) == std::end(flagged_cells))
        ah.define_agglomerate({cells[i]});
    }

  FE_DGQ<dim> fe_dg(1);
  ah.distribute_agglomerated_dofs(fe_dg);
  ah.initialize_fe_values(QGauss<dim>(1), update_JxW_values);

  // prints
  const auto fu = [&](const types::global_cell_index given_index) {
    for (const auto &polytope : ah.polytope_iterators())
      {
        if (polytope->index() <= given_index)
          {
            const auto &fev = ah.reinit(polytope);
            double      sum = 0.;
            for (const double weight : fev.get_JxW_values())
              sum += weight;
            std::cout << "Sum is: " << sum << std::endl;
          }
      }
  };

  if constexpr (dim == 2)
    {
      fu(3);
    }
  else if constexpr (dim == 3)
    {
      fu(0);
    }
}


int
main()
{
  test<2>();
  test<3>();
  return 0;
}
