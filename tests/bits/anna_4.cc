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



// check for an abort in VectorTools::interpolate_boundary_values.
//
// this program is a modified version of one by Anna Schneebeli,
// University of Basel

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// We need a FESystem
#include <deal.II/fe/fe_system.h>

// we need DG-elements
// and Q1-elements
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>

#include <vector>


template <int dim>
class VectorBoundaryValues : public Function<dim>
{
public:
  VectorBoundaryValues();
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template <int dim>
VectorBoundaryValues<dim>::VectorBoundaryValues()
  : Function<dim>(2)
{}

template <int dim>
inline void
VectorBoundaryValues<dim>::vector_value(const Point<dim> &p,
                                        Vector<double>   &values) const
{
  Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2));

  for (unsigned int i = 0; i < 2; ++i)
    values(i) = p[i] * p[i];
}



template <int dim>
class FindBug
{
public:
  FindBug();
  void
  run();

private:
  void
  make_grid_and_dofs();
  void
  dirichlet_conditions();

  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;
  Vector<double>     solution;
};


// Construct FESystem with
// first component: Q1-Element,
// second component: lowest order DG_Element
template <int dim>
FindBug<dim>::FindBug()
  : fe(FE_Q<dim>(1), 1, FE_DGP<dim>(0), 1)
  , dof_handler(triangulation)
{}


template <int dim>
void
FindBug<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;

  deallog << "Total number of cells: " << triangulation.n_cells() << std::endl;


  dof_handler.distribute_dofs(fe);


  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

  solution.reinit(dof_handler.n_dofs());
}


template <int dim>
void
FindBug<dim>::dirichlet_conditions()
{
  // we want to set the Boundary DoFs
  // of the selected component to a
  // given value, say zero.  To do
  // so, we want to use VectorTools::
  // interpolate_boundary_values
  //
  // This works fine if all the
  // components have support on the
  // faces.  (This, of course, has to
  // be requested when fixing the
  // boundary DoFs.)  However,
  // getting the values for the
  // boundary DoFs of a valid
  // component by the function
  // VectorTools::
  // interpolate_boundary_values and
  // an appropriate component_mask
  // (unselecting the DG-component)
  // does not work yet. Or, better:
  // it did not work and this program
  // checks that it has been
  // correctly implemented by now


  std::map<types::global_dof_index, double> dirichlet_dofs;

  // we declare a vector of bools,
  // which tells the
  // VectorTools::interpolate_boundary_values
  // on which components of the
  // FESystem we want to impose
  // Dirichlet BC.
  ComponentMask component_mask(2, false);
  // Dirichlet-BC for the
  // Q1-Component
  component_mask.set(0, true);
  // no Dirichlet-BC for the
  // DG-component
  component_mask.set(1, false);

  // This is just for the final
  // output-test
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    dirichlet_dofs[i] = 1.;


  // Here comes the crucial call....
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(2),
                                           dirichlet_dofs,
                                           component_mask);


  const std::set<types::boundary_id> boundary_ids = {0};

  // get a list of those boundary DoFs which
  // we want to be fixed:
  const IndexSet fixed_dofs =
    DoFTools::extract_boundary_dofs(dof_handler, component_mask, boundary_ids);

  // (Primitive) Check if the DoFs
  // where adjusted correctly (note
  // that we have preset all values
  // to 1, and interpolate_b_v should
  // have overwritten those for
  // component 0 by 0)
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      if (fixed_dofs.is_element(i))
        {
          AssertThrow(dirichlet_dofs[i] == 0, ExcInternalError());
        }
      else
        {
          AssertThrow(dirichlet_dofs[i] == 1, ExcInternalError());
        };
    };

  // check 1 has obviously succeeded,
  // so also check a more complicated
  // boundary value function
  dirichlet_dofs.clear();
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           VectorBoundaryValues<dim>(),
                                           dirichlet_dofs,
                                           component_mask);
  for (const types::global_dof_index dof : fixed_dofs)
    deallog << dof << ' ' << dirichlet_dofs[dof] << std::endl;
}



template <int dim>
void
FindBug<dim>::run()
{
  make_grid_and_dofs();
  dirichlet_conditions();
}



int
main()
{
  initlog();

  FindBug<2>().run();
  FindBug<3>().run();

  return 0;
}
