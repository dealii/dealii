// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test by Korosh Taebi: check that we can gave Q(p) and DGQ(r) in the same
// mesh


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


namespace Step27
{
  template <int dim>
  class MixedFECollection
  {
  public:
    MixedFECollection();
    ~MixedFECollection();

    void
    run();

  private:
    Triangulation<dim>    triangulation;
    DoFHandler<dim>       dof_handler;
    hp::FECollection<dim> fe_collection;
  };

  template <int dim>
  MixedFECollection<dim>::MixedFECollection()
    : dof_handler(triangulation)
  {}

  template <int dim>
  MixedFECollection<dim>::~MixedFECollection()
  {
    dof_handler.clear();
  }

  template <int dim>
  void
  MixedFECollection<dim>::run()
  {
    // add two a CG and a DG finite element object to fe_collection
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_DGQ<dim>(1));
    deallog << " fe_collection size = " << fe_collection.size() << std::endl;

    // produce a simple grid with 4 cells
    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(1);

    // looping over all cells and assigning the FE_DG object to the first cell
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (unsigned int counter = 0; cell != endc; ++cell, counter++)
      if (counter == 0)
        {
          cell->set_active_fe_index(cell->active_fe_index() + 1);
        }

    dof_handler.distribute_dofs(fe_collection);

    deallog << "   Number of active cells:       "
            << triangulation.n_active_cells() << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  }
} // namespace Step27



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  Step27::MixedFECollection<1>().run();
  Step27::MixedFECollection<2>().run();
  Step27::MixedFECollection<3>().run();

  deallog << "OK" << std::endl;
}
