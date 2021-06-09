// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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



// check that generating a mg flux sparsity pattern without constraints and
// later calling constraints.condense() on it results in the same pattern as
// when creating the condensed pattern right away

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"



template <int dim>
void
check()
{
  constexpr unsigned int min_level{0}, max_level{2};

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(max_level);

  // create a system element composed
  // of one Q1 and one Q2 element
  FESystem<dim>   element(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);
  dof.distribute_mg_dofs();

  MGLevelObject<AffineConstraints<double>> mg_constraints(min_level, max_level);

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof);

  std::set<types::boundary_id> set_dirichlet;
  set_dirichlet.insert(0);
  mg_constrained_dofs.make_zero_boundary_constraints(dof, set_dirichlet);

  // test sparsity patterns on each MG level
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      mg_constraints[l].reinit(dof.locally_owned_mg_dofs(l));
      mg_constraints[l].add_lines(
        mg_constrained_dofs.get_refinement_edge_indices(l));
      mg_constraints[l].add_lines(mg_constrained_dofs.get_boundary_indices(l));
      mg_constraints[l].close();

      //--------------- Regular sparsity pattern checks -----------------
      // first way: directly
      SparsityPattern sparsity_1(dof.n_dofs(l), dof.n_dofs(l));
      MGTools::make_flux_sparsity_pattern(dof, sparsity_1, l);
      mg_constraints[l].condense(sparsity_1);
      sparsity_1.compress();

      // second way: via direct elimination of
      // constraints
      SparsityPattern        sparsity_2;
      DynamicSparsityPattern dsp_2(dof.locally_owned_mg_dofs(l));
      MGTools::make_flux_sparsity_pattern(dof, dsp_2, l, mg_constraints[l]);
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
