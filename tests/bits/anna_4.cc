// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// check for an abort in VectorTools::interpolate_boundary_values.
//
// this program is a modified version of one by Anna Schneebeli,
// University of Basel

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

// We need a FESystem
#include <deal.II/fe/fe_system.h>

// we need DG-elements
// and Q1-elements
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>

#include <fstream>
#include <iomanip>
#include <vector>


template <int dim>
class VectorBoundaryValues :  public Function<dim>
{
public:
  VectorBoundaryValues ();
  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &values) const;
};

template <int dim>
VectorBoundaryValues<dim>::VectorBoundaryValues () :
  Function<dim> (2)
{}

template <int dim>
inline
void VectorBoundaryValues<dim>::vector_value (const Point<dim> &p,
                                              Vector<double>   &values) const
{
  Assert (values.size() == 2,
          ExcDimensionMismatch (values.size(), 2));

  for (unsigned int i=0; i<2; ++i)
    values(i) = p(i)*p(i);
}



template <int dim>
class FindBug
{
public:
  FindBug ();
  void run ();
private:
  void make_grid_and_dofs ();
  void dirichlet_conditions ();

  Triangulation<dim>     triangulation;
  FESystem<dim>              fe;
  DoFHandler<dim>        dof_handler;
  Vector<double>          solution;
};


// Construct FESystem with
// first component: Q1-Element,
// second component: lowest order DG_Element
template <int dim>
FindBug<dim>::FindBug () :
  fe (FE_Q<dim>(1), 1,
      FE_DGP<dim>(0), 1),
  dof_handler (triangulation)
{}


template <int dim>
void FindBug<dim>::make_grid_and_dofs ()
{

  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (1);

  deallog << "Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl;

  deallog << "Total number of cells: "
          << triangulation.n_cells()
          << std::endl;


  dof_handler.distribute_dofs (fe);


  deallog << "Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  solution.reinit(dof_handler.n_dofs());
}


template <int dim>
void FindBug<dim>::dirichlet_conditions ()
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


  std::map<types::global_dof_index,double> dirichlet_dofs;

  // we declare a vector of bools,
  // which tells the
  // VectorTools::interpolate_boundary_values
  // on which components of the
  // FESystem we want to impose
  // Dirichlet BC.
  std::vector<bool> component_mask(2);
  // Dirichlet-BC for the
  // Q1-Component
  component_mask[0] = true;
  // no Dirichlet-BC for the
  // DG-component
  component_mask[1] = false;

  // This is just for the final
  // output-test
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    dirichlet_dofs[i] = 1.;


  // Here comes the crucial call....
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim> (2),
                                            dirichlet_dofs,
                                            component_mask);


  std::vector<bool> fixed_dofs (dof_handler.n_dofs());
  std::set<types::boundary_id> boundary_indicators;
  boundary_indicators.insert (0);

  // get a list of those boundary DoFs which
  // we want to be fixed:
  DoFTools::extract_boundary_dofs (dof_handler,
                                   component_mask,
                                   fixed_dofs,
                                   boundary_indicators);

  // (Primitive) Check if the DoFs
  // where adjusted correctly (note
  // that we have preset all values
  // to 1, and interpolate_b_v should
  // have overwritten those for
  // component 0 by 0)
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if (fixed_dofs[i] == true)
        {
          Assert (dirichlet_dofs[i] == 0, ExcInternalError());
        }
      else
        {
          Assert (dirichlet_dofs[i] == 1, ExcInternalError());
        };
    };

  // check 1 has obviously succeeded,
  // so also check a more complicated
  // boundary value function
  dirichlet_dofs.clear ();
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            VectorBoundaryValues<dim> (),
                                            dirichlet_dofs,
                                            component_mask);
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (fixed_dofs[i] == true)
      deallog << i << " " << dirichlet_dofs[i] << std::endl;
}



template <int dim>
void FindBug<dim>::run ()
{
  make_grid_and_dofs ();
  dirichlet_conditions ();
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FindBug<2>().run ();
  FindBug<3>().run ();

  return 0;
}

