// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include "../tests.h"
#include "fe_tools_common.h"
#include <deal.II/lac/sparsity_pattern.h>

// check
//   FETools::interpolation_difference(6)


std::string output_file_name = "output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  // only check if both elements have
  // support points. otherwise,
  // interpolation doesn't really
  // work
  if ((fe1.get_unit_support_points().size() == 0) ||
      (fe2.get_unit_support_points().size() == 0))
    return;
  //  likewise for non-primitive elements
  if (!fe1.is_primitive() || !fe2.is_primitive())
    return;
  // we need to have dof_constraints
  // for this test
  if (!fe1.constraints_are_implemented() ||
      !fe2.constraints_are_implemented())
    return;

  std::auto_ptr<Triangulation<dim> > tria(make_tria<dim>());
  std::auto_ptr<DoFHandler<dim> >    dof1(make_dof_handler (*tria, fe1));
  std::auto_ptr<DoFHandler<dim> >    dof2(make_dof_handler (*tria, fe2));
  ConstraintMatrix cm1, cm2;
  DoFTools::make_hanging_node_constraints (*dof1, cm1);
  DoFTools::make_hanging_node_constraints (*dof2, cm2);
  cm1.close ();
  cm2.close ();

  Vector<double> in (dof1->n_dofs());
  for (unsigned int i=0; i<in.size(); ++i) in(i) = i;
  Vector<double> out (dof1->n_dofs());

  FETools::interpolation_difference (*dof1, cm1, in, *dof2, cm2, out);
  output_vector (out);
}

