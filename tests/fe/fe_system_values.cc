// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that we get the correct values out of FESystem with a strange set of
// base elements to validate the new indexing code.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"

template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof_handler(tr);
  dof_handler.distribute_dofs(fe);

  deallog << "FE = " << fe.get_name() << std::endl;

  const QGauss<dim> quadrature(2);
  MappingQ<dim>     mapping(2);
  FEValues<dim>     fe_values(mapping,
                          fe,
                          quadrature,
                          update_values | update_quadrature_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      deallog << "Cell at center " << cell->center() << std::endl;

      for (unsigned int qp = 0; qp < quadrature.size(); ++qp)
        for (unsigned int d = 0; d < fe.n_dofs_per_cell(); ++d)
          {
            deallog << "d = " << d;
            for (unsigned int c = 0; c < fe.n_components(); ++c)
              deallog << ", " << fe_values.shape_value_component(d, qp, c);
            deallog << std::endl;
          }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(8);

  {
    Triangulation<2> tria;
    GridGenerator::hyper_ball(tria);

    FESystem<2> fe(FE_RaviartThomas<2>(0),
                   1,
                   FE_Q<2>(2),
                   2,
                   FE_RaviartThomas<2>(2),
                   3,
                   FE_Nothing<2>(),
                   4,
                   FE_DGQ<2>(0),
                   1);
    test(tria, fe);
  }

  {
    Triangulation<2> tria;
    GridGenerator::hyper_ball(tria);

    FESystem<2> stokes_fe(FE_RaviartThomas<2>(1), 1, FE_DGQ<2>(1), 1);
    FESystem<2> base_fe(stokes_fe, 1, FE_RaviartThomas<2>(0), 1);
    FESystem<2> fe(base_fe, 2, FE_DGQ<2>(1), 1);
    test(tria, fe);
  }

  deallog << "done..." << std::endl;
}
