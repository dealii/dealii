// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// validate algorithms that will flag cells for p-adaptivity:
// - hp::Refinement::p_adaptivity_from_flags


#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
validate(const DoFHandler<dim> &dh)
{
  deallog << " fe_indices:";
  for (const auto &cell : dh.active_cell_iterators())
    deallog << ' ' << cell->future_fe_index();
  deallog << std::endl;
}



template <int dim>
void
setup(Triangulation<dim> &tria, const DoFHandler<dim> &dh)
{
  // Initialize triangulation.
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // Set all active FE indices to 1.
  // Flag first half of cells for refinement, and the other half for coarsening.
  typename DoFHandler<dim>::cell_iterator cell = dh.begin(1), endc = dh.end(1);
  for (unsigned int counter = 0; cell != endc; ++counter, ++cell)
    {
      Assert(!cell->is_active(), ExcInternalError());
      for (unsigned int child_index = 0; child_index < cell->n_children();
           ++child_index)
        {
          const auto &child = cell->child(child_index);
          Assert(child->is_active(), ExcInternalError());

          child->set_active_fe_index(1);

          if (counter < 0.5 * GeometryInfo<dim>::max_children_per_cell)
            child->set_refine_flag();
          else
            child->set_coarsen_flag();
        }
    }
}



template <int dim>
void
test()
{
  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  Triangulation<dim> tria;
  DoFHandler<dim>    dh(tria);
  setup(tria, dh);
  dh.distribute_dofs(fes);

  deallog << "starting situation" << std::endl;
  validate(dh);


  // We flag the first half of all cells to be refined and the last half of all
  // cells to be coarsened for p-adapativity. Ultimately, the first quarter of
  // all cells will be flagged for p-refinement, and the last quarter for p-
  // coarsening.

  const unsigned int n_active = tria.n_active_cells();
  std::vector<bool>  p_flags(n_active, false);
  std::fill(p_flags.begin(), p_flags.begin() + n_active / 4, true);
  std::fill(p_flags.end() - n_active / 4, p_flags.end(), true);

  hp::Refinement::p_adaptivity_from_flags(dh, p_flags);

  deallog << "p-adaptivity from flags" << std::endl;
  validate(dh);


  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
