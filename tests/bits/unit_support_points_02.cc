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



// FESystem::get_unit_support_points would return an empty array if
// one base element had no support points. but that's not necessary
// that way: it should still return something useful if the element
// has no support points but in fact also has no degrees of freedom on
// faces at all. in that case we would simply not care
//
// check this by comparing
//    FESystem(FE_Q, 1,  FE_DGQ, 1)
// with
//    FESystem(FE_Q, 1,  FE_DGP, 1)
// (the latter not having any support points for the second component).

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <fstream>


template <int dim>
void check2 (const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  const std::vector<Point<dim-1> > unit_f_s_p
    = fe.get_unit_face_support_points ();
  for (unsigned int i=0; i<unit_f_s_p.size(); ++i)
    deallog << i << ' ' << unit_f_s_p[i] << std::endl;
}


template <int dim>
void check ()
{
  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGQ<dim> (2), 1));
  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGP<dim> (2), 1));
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
  return 0;
}
