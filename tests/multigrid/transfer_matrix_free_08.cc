// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check MGTransferMatrixFree for high polynomial degrees beyond 10 by
// checking a linear function that gets prolongated between the mesh
// levels. Since the function is linear, it should be exactly represented on
// the finest level. Then also look into the result of restriction

#include <deal.II/base/function_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, typename Number>
void
check(const FiniteElement<dim> &fe)
{
  deallog << "FE: " << fe.get_name() << std::endl;

  // run a few different sizes...

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  Triangulation<dim> trcoarse(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(trcoarse);
  DoFHandler<dim> dofcoarse(trcoarse);
  dofcoarse.distribute_dofs(fe);
  Tensor<1, dim> exponents_monomial;
  for (unsigned int d = 0; d < dim; ++d)
    exponents_monomial[d] = 1;
  LinearAlgebra::distributed::Vector<double> vrefcoarse;
  vrefcoarse.reinit(dofcoarse.n_dofs());
  VectorTools::interpolate(dofcoarse,
                           Functions::Monomial<dim>(exponents_monomial),
                           vrefcoarse);

  deallog << "no. cells: " << tr.n_global_active_cells() << std::endl;

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  // build matrix-free transfer
  MGTransferMatrixFree<dim, Number> transfer;
  transfer.build(mgdof);

  LinearAlgebra::distributed::Vector<double> vrefdouble;
  LinearAlgebra::distributed::Vector<Number> vref;
  AssertDimension(mgdof.n_dofs(tr.n_global_levels() - 1), mgdof.n_dofs());
  vrefdouble.reinit(mgdof.n_dofs());
  VectorTools::interpolate(mgdof,
                           Functions::Monomial<dim>(exponents_monomial),
                           vrefdouble);

  vref.reinit(mgdof.n_dofs());
  vref = vrefdouble;
  std::vector<LinearAlgebra::distributed::Vector<Number>> vec(
    tr.n_global_levels());
  for (unsigned int level = 0; level < tr.n_global_levels(); ++level)
    vec[level].reinit(mgdof.n_dofs(level));
  vec.back() = vref;

  // prolongate monomial from coarse to fine, should be exact for monomial
  vec[0] = vrefcoarse;
  for (unsigned int level = 1; level < tr.n_global_levels(); ++level)
    transfer.prolongate(level, vec[level], vec[level - 1]);
  vec.back() -= vref;
  const Number tolerance = 1000. * std::numeric_limits<Number>::epsilon();
  deallog << "Error after prolongation: "
          << filter_out_small_numbers(vec.back().linfty_norm(), tolerance)
          << std::endl;

  // unfortunately, no completely trivial expression of what should happen
  // during restriction
  vec.back() = vref;
  for (unsigned int level = tr.n_global_levels() - 1; level > 0; --level)
    transfer.restrict_and_add(level, vec[level - 1], vec[level]);
  deallog << "Norm after restriction: " << vec[0].l2_norm() << std::endl;
}


int
main(int argc, char **argv)
{
  // no threading in this test...
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  check<2, double>(FE_DGQ<2>(3));
  check<2, double>(FE_DGQ<2>(8));
  check<2, double>(FE_DGQ<2>(16));
  check<2, double>(FE_Q<2>(16));
  check<3, double>(FE_DGQ<3>(11));
  check<2, float>(FE_DGQ<2>(3));
  check<2, float>(FE_DGQ<2>(8));
  check<2, float>(FE_DGQ<2>(16));
  check<2, float>(FE_Q<2>(16));
  check<3, float>(FE_DGQ<3>(11));
}
