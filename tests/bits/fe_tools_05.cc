//----------------------------  fe_tools_05.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_05.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.cc"
#include <lac/sparsity_pattern.h>

// check
//   FETools::interpolate(5)


std::string output_file_name = "fe_tools_05.output";


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
  if (!fe2.constraints_are_implemented())
    return;

  std::auto_ptr<Triangulation<dim> > tria(make_tria<dim>());
  std::auto_ptr<DoFHandler<dim> >    dof1(make_dof_handler (*tria, fe1));
  std::auto_ptr<DoFHandler<dim> >    dof2(make_dof_handler (*tria, fe2));
  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (*dof2, cm);
  cm.close ();
  
  Vector<double> in (dof1->n_dofs());
  for (unsigned int i=0; i<in.size(); ++i) in(i) = i;
  Vector<double> out (dof2->n_dofs());
  
  FETools::interpolate (*dof1, in, *dof2, cm, out);
  output_vector (out);
}

