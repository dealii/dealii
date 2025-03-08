// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check gradients in support points for FE_Q_iso_Q1

#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

template <int dim>
void check(const FiniteElement<dim> &fe, const Quadrature<dim>& quad, const char* name)
{
  deallog << name << ":" << std::endl;
  
  Triangulation<dim> tr;
  MappingQ<dim> mapping(1);
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  dof.distribute_dofs(fe);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  
    FEValues<dim>  fev(mapping,
                     fe,
                     quad,
                     UpdateFlags(update_gradients));

    fev.reinit(c);
    const unsigned int n_q_points    = quad.size();

    for (unsigned int q=0;q<n_q_points;++q)
      {
      deallog << quad.point(q);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {
	deallog << " " << fev.shape_grad(i, q).norm();
      }
      deallog << std::endl;
      }
 
}

int
main()
{
  initlog();
  
  FE_Q_iso_Q1<2> f1(2);
  FE_Q_iso_Q1<2> f2(QGaussLobatto<1>(2 + 1).get_points());
  QGauss<2> quadrature(3);

  check<>(f1, quadrature, "iso_q1(2)-subd");
  check<>(f2, quadrature, "iso_q1(2)-GL");
}
