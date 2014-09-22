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



// check equivalence of get_unit_support_points() and
// unit_support_point(), as well as the case where an element does not
// have support points

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_nedelec.h>
#include <fstream>


template <int dim>
void check_cell1 (const FiniteElement<dim> &fe)
{
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    Assert (fe.get_unit_support_points()[i] ==
            fe.unit_support_point(i),
            ExcInternalError());
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void check_face1 (const FiniteElement<dim> &fe)
{
  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
    Assert (fe.get_unit_face_support_points()[i] ==
            fe.unit_face_support_point(i),
            ExcInternalError());
  deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void check_face1 (const FiniteElement<1> &)
{}



template <int dim>
void check1 (const FiniteElement<dim> &fe)
{
  check_cell1 (fe);
  check_face1 (fe);
}


template <int dim>
void check_cell2 (const FiniteElement<dim> &fe,
                  const unsigned int        comp)
{
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    if (fe.system_to_component_index(i).first == comp)
      deallog << i << " " << fe.unit_support_point(i)
              << std::endl;
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void check_face2 (const FiniteElement<dim> &fe,
                  const unsigned int        comp)
{
  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
    if (fe.face_system_to_component_index(i).first == comp)
      deallog << i << " " << fe.unit_face_support_point(i)
              << std::endl;
  deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void check_face2 (const FiniteElement<1> &,
                  const unsigned int)
{}



template <int dim>
void check2 (const FiniteElement<dim> &fe,
             const unsigned int        comp)
{
  check_cell2 (fe, comp);
  check_face2 (fe, comp);
}



template <int dim>
void face_check()
{
  check1 (FE_FaceQ<dim>(2));
  check2 (FE_FaceQ<dim>(2), 0);
}



template <>
void face_check<1>()
{}



template <int dim>
void check ()
{
  check1 (FE_Q<dim>(2));
  check1 (FE_DGQ<dim>(2));
  check1 (FESystem<dim> (FE_Q<dim> (2), 2,
                         FE_DGQ<dim> (2),1));

  check1 (FESystem<dim> (FE_Q<dim> (2), 2,
                         FE_DGQ<dim> (2),1));

  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGQ<dim> (2), 1),
          1);
  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGP<dim> (2), 1),
          0);
  face_check<dim>();
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
  return 0;
}
