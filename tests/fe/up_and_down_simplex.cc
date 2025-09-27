// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// like fe/up_and_down test but use simplex elements instead
// Restrict a vector to a coarse grid and prolongate it again,
// repeat the procedure twice and check if the reslting vectors are the same
// i.e. if prolongation and restriction are successful.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 3


template <int dim>
Point<dim>
transform(const Point<dim> p)
{
  switch (dim)
    {
      case 1:
        return p;
      case 2:
        return Point<dim>(p[0] * (1 + p[1]), p[1] * (1 + p[0]));
      case 3:
        return Point<dim>(p[0] * (1 + p[1]) * (1 + p[2]),
                          p[1] * (1 + p[0]) * (1 + p[2]),
                          p[2] * (1 + p[0]) * (1 + p[1]));
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return Point<dim>();
    };
}


template <int dim>
void
check_element(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof_handler(tr);
  dof_handler.distribute_dofs(fe);

  // create a mostly arbitrary
  // function plus a trend on this
  // grid
  Vector<double> tmp(dof_handler.n_dofs());
  // for (unsigned int i = 0; i < tmp.size(); ++i)
  //   tmp(i) = i; //(i + 13*i%17);
  for (const auto &cell : tr.cell_iterators())
    for (const auto v : cell->vertex_indices())
      {
        const auto v_global = cell->vertex_index(v);
        tmp(v_global)       = cell->vertex(v)(0) + cell->vertex(v)(1);
        if (dim == 3)
          tmp(v_global) = cell->vertex(v)(0) + cell->vertex(v)(2);
      }


  // restrict this function to the
  // next coarser level and
  // distribute it again to the
  // higher level
  Vector<double> x(tmp.size());
  Vector<double> v(fe.dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    if (cell->has_children() && cell->child(0)->is_active())
      {
        // first make sure that what
        // we do is reasonable. for
        // this, _all_ children have
        // to be active, not only
        // some of them
        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          AssertThrow(cell->child(c)->is_active(), ExcInternalError());

        // then restrict and prolongate
        cell->get_interpolated_dof_values(tmp, v);
        cell->set_dof_values_by_interpolation(v, x);
      };

  // now x is a function on the fine
  // grid that is representable on
  // the coarse grid. so another
  // cycle should not alter it any
  // more:
  Vector<double> x2(x.size());
  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    if (cell->has_children() && cell->child(0)->is_active())
      {
        cell->get_interpolated_dof_values(x, v);
        cell->set_dof_values_by_interpolation(v, x2);
      };

  // then check that this is so:
  x2 -= x;
  const double relative_residual = (x2.l2_norm() / x.l2_norm());

  const double threshold = 1e-6;
  deallog << ", dofs_per_cell=" << fe.dofs_per_cell
          << "; relative residual: " << relative_residual << " "
          << (relative_residual < threshold ? "ok" : "botched up!")
          << std::endl;

  // TODO:[WB] Why this exception with a value different from above. Output of
  // the error should be sufficient!
  //  Assert (relative_residual < threshold*x.l2_norm(), ExcInternalError());
}


template <int dim>
void
test()
{
  // make a coarse triangulation as a
  // hypercube. if in more than 1d,
  // distort it so that it is no more
  // an affine image of the
  // hypercube, to make things more
  // difficult. then refine it twice
  Triangulation<dim> temp, tr;
  GridGenerator::hyper_cube(temp, 0., 1.);
  Point<dim> (*p)(Point<dim>) = &transform<dim>;
  GridTools::transform(p, temp);
  GridGenerator::convert_hypercube_to_simplex_mesh(temp, tr);

  tr.refine_global(2);

  // now for a list of finite
  // elements, for which we want to
  // test. we happily waste tons of
  // memory here, but who cares...
  const FiniteElement<dim> *fe_list[] = {
    // first for some scalar
    // elements:

    // FE_SimplexP
    new FE_SimplexP<dim>(1),
    new FE_SimplexP<dim>(2),
    new FE_SimplexP<dim>(3),
    new FE_SimplexDGP<dim>(1),
    new FE_SimplexDGP<dim>(2),
    new FE_SimplexDGP<dim>(3),

    // some composed elements
    // of increasing
    // complexity, to check the
    // logics by which the
    // matrices of the composed
    // elements are assembled
    // from those of the base
    // elements. note that some
    // of the base elements are
    // additive, some not, so
    // the result will be an
    // element that is mixed in
    // this respect

    new FESystem<dim>(FE_SimplexP<dim>(2), 2),
    new FESystem<dim>(FE_SimplexP<dim>(1), 2, FE_SimplexDGP<dim>(2), 2),
    new FESystem<dim>(FE_SimplexP<dim>(1),
                      2,
                      FE_SimplexDGP<dim>(2),
                      2,
                      FE_SimplexDGP<dim>(0),
                      1),
    new FESystem<dim>(FE_SimplexP<dim>(1),
                      2,
                      FESystem<dim>(FE_SimplexP<dim>(1),
                                    2,
                                    FE_SimplexDGP<dim>(2),
                                    2,
                                    FE_SimplexDGP<dim>(2),
                                    1),
                      2,
                      FE_SimplexDGP<dim>(0),
                      1),
    new FESystem<dim>(FE_SimplexP<dim>(1),
                      2,
                      FESystem<dim>(FE_SimplexP<dim>(1),
                                    2,
                                    FE_SimplexDGP<dim>(2),
                                    2,
                                    FESystem<dim>(FE_SimplexDGP<dim>(0), 3),
                                    1),
                      2,
                      FE_SimplexDGP<dim>(0),
                      1),


  };

  for (unsigned int i = 0; i < sizeof(fe_list) / sizeof(fe_list[0]); ++i)
    if (fe_list[i] != nullptr)
      {
        deallog << dim << "d, uniform grid, fe #" << i;
        check_element(tr, *fe_list[i]);
      }


  for (unsigned int i = 0; i < sizeof(fe_list) / sizeof(fe_list[0]); ++i)
    if (fe_list[i] != nullptr)
      delete fe_list[i];
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  test<2>();
  test<3>();

  return 0;
}
