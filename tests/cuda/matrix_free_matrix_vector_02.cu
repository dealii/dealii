// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a hypercube mesh with no hanging nodes, but with zero
// Dirichlet conditions.

#include <deal.II/base/function.h>

#include "../tests.h"

#include "matrix_vector_common.h"

template <int dim, int fe_degree, typename Number>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<Number> constraints;
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  do_test<dim,
          fe_degree,
          Number,
          LinearAlgebra::CUDAWrappers::Vector<Number>,
          fe_degree + 1>(dof, constraints, tria.n_active_cells());
}
