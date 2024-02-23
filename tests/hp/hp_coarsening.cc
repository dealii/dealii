// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// validate hp-decision algorithms on grid coarsening
// that depend on the composition of h- and p-adaptivity flags


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"



template <int dim>
void
validate(const Triangulation<dim> &tria, const DoFHandler<dim> &dh)
{
  deallog << "ncells: " << tria.n_global_active_cells() << " fe_indices:";
  for (const auto &cell : dh.active_cell_iterators())
    deallog << ' ' << cell->active_fe_index();
  deallog << std::endl;
}



template <int dim>
void
setup(Triangulation<dim> &tria, const DoFHandler<dim> &dh)
{
  // Initialize triangulation.
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  // Set h- and p-flags on all cells.
  for (const auto &cell : dh.active_cell_iterators())
    {
      cell->set_coarsen_flag();
      cell->set_future_fe_index(1);
    }
}



template <int dim>
void
test()
{
  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 2; ++d)
    fes.push_back(FE_Q<dim>(d));

  deallog << "starting situation: ";
  {
    Triangulation<dim> tria;
    DoFHandler<dim>    dh(tria);
    setup(tria, dh);
    dh.distribute_dofs(fes);

    validate(tria, dh);
  }

  deallog << "full h&p flags" << std::endl;
  {
    deallog << " default behavior: ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " force p over h   : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      hp::Refinement::force_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " choose p over h  : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      hp::Refinement::choose_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }
  }


  deallog << "full p flags" << std::endl;
  {
    deallog << " default behavior: ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_coarsen_flag();
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " force p over h   : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_coarsen_flag();
      hp::Refinement::force_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " choose p over h  : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_coarsen_flag();
      hp::Refinement::choose_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }
  }


  deallog << "full h flags" << std::endl;
  {
    deallog << " default behavior: ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_future_fe_index();
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " force p over h   : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_future_fe_index();
      hp::Refinement::force_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }

    deallog << " choose p over h  : ";
    {
      Triangulation<dim> tria;
      DoFHandler<dim>    dh(tria);
      setup(tria, dh);
      dh.distribute_dofs(fes);

      dh.begin_active()->clear_future_fe_index();
      hp::Refinement::choose_p_over_h(dh);
      tria.execute_coarsening_and_refinement();

      validate(tria, dh);
    }
  }

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
