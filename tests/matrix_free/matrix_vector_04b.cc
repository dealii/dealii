// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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



// Same as matrix_vector_04, but using more quadrature points than the FE
// degree would suggest

#include "../tests.h"

#include "matrix_vector_common.h"


template <int dim, int fe_degree>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 96, true);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);
  if (dim == 2)
    tria.refine_global(2);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  do_test<dim, fe_degree, double, fe_degree + 3>(dof, constraints);
}
