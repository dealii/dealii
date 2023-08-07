// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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


// Check that MGTransferMatrixFree::clear works correctly by comparing a
// cleared transfer system with a freshly constructed one
#include <deal.II/base/function_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, typename Number>
void
check(const FiniteElement<dim> &fe)
{
  deallog << "FE: " << fe.get_name() << std::endl;

  MGConstrainedDoFs         mg_constrained_dofs;
  MGTransferMF<dim, Number> transfer;
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(6 - dim);

  // run a few different sizes...
  for (unsigned int c = 0; c < 4; ++c)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        if ((cell->center().norm() < 0.5 &&
             (cell->level() < 5 || cell->center().norm() > 0.45)) ||
            (dim == 2 && cell->center().norm() > 1.2))
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      DoFHandler<dim> mgdof(tr);
      mgdof.distribute_dofs(fe);
      mgdof.distribute_mg_dofs();

      Tensor<1, dim> exponents_monomial;
      for (unsigned int d = 0; d < dim; ++d)
        exponents_monomial[d] = 1;
      LinearAlgebra::distributed::Vector<double> vref;
      vref.reinit(mgdof.n_dofs());
      VectorTools::interpolate(mgdof,
                               Functions::Monomial<dim>(exponents_monomial),
                               vref);

      deallog << "no. cells: " << tr.n_global_active_cells() << std::endl;

      const std::set<types::boundary_id> dirichlet_boundary = {0};
      mg_constrained_dofs.initialize(mgdof);
      mg_constrained_dofs.make_zero_boundary_constraints(mgdof,
                                                         dirichlet_boundary);

      // build matrix-free transfer
      transfer.clear();
      transfer.initialize_constraints(mg_constrained_dofs);
      transfer.build(mgdof);
      MGTransferMF<dim, Number> transfer_ref(mg_constrained_dofs);
      transfer_ref.build(mgdof);
      MGLevelObject<LinearAlgebra::distributed::Vector<Number>> vectors(
        0, tr.n_global_levels() - 1);
      transfer_ref.copy_to_mg(mgdof, vectors, vref);
      for (unsigned int level = vectors.max_level(); level > 0; --level)
        {
          LinearAlgebra::distributed::Vector<Number> vec2(vectors[level - 1]);
          transfer_ref.restrict_and_add(level,
                                        vectors[level - 1],
                                        vectors[level]);
          transfer.restrict_and_add(level, vec2, vectors[level]);
          vec2 -= vectors[level - 1];
          deallog << "Error in restriction:  " << (double)vec2.linfty_norm()
                  << std::endl;
        }

      for (unsigned int level = 1; level < vectors.max_level(); ++level)
        {
          LinearAlgebra::distributed::Vector<Number> vec2(vectors[level + 1]);
          transfer_ref.prolongate(level + 1,
                                  vectors[level + 1],
                                  vectors[level]);
          transfer.prolongate(level + 1, vec2, vectors[level]);
          vec2 -= vectors[level + 1];
          deallog << "Error in prolongation: " << (double)vec2.linfty_norm()
                  << std::endl;
        }
    }
}


int
main(int argc, char **argv)
{
  // no threading in this test...
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  check<2, double>(FE_DGQ<2>(2));
  check<2, double>(FE_Q<2>(2));
  check<3, double>(FE_Q<3>(1));
  check<2, float>(FE_Q<2>(2));
}
