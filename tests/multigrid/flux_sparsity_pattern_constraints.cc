// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that generating a mg flux sparsity pattern without constraints and
// later calling constraints.condense() on it results in the same pattern as
// when creating the condensed pattern right away

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#include "../tests.h"



template <int dim>
void
check()
{
  constexpr unsigned int min_level{0}, max_level{2};

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, -1., 1., true);

  // create periodicity constraints
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;

  GridTools::collect_periodic_faces(tr,
                                    /*b_id1*/ 0,
                                    /*b_id2*/ 1,
                                    /*direction*/ 0,
                                    periodicity_vector);

  tr.add_periodicity(periodicity_vector);
  tr.refine_global(max_level);

  // create a system element composed
  // of one Q1 and one Q2 element
  FESystem<dim>   element(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);
  dof.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof);

  // test sparsity patterns on each MG level
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      //--------------- Regular sparsity pattern checks -----------------
      // first way: directly
      SparsityPattern sparsity_1(dof.n_dofs(l), dof.n_dofs(l));
      MGTools::make_flux_sparsity_pattern(dof, sparsity_1, l);
      mg_constrained_dofs.get_level_constraints(l).condense(sparsity_1);
      sparsity_1.compress();

      // second way: via direct elimination of constraints
      SparsityPattern        sparsity_2;
      DynamicSparsityPattern dsp_2(dof.locally_owned_mg_dofs(l));
      MGTools::make_flux_sparsity_pattern(
        dof, dsp_2, l, mg_constrained_dofs.get_level_constraints(l));
      sparsity_2.copy_from(dsp_2);

      // tests if sparsity_[12] are equal
      deallog << "Check 1:"
              << " -- " << (sparsity_1 == sparsity_2 ? "ok" : "failed")
              << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
