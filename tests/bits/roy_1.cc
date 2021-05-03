// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// check has_support_on_face for some elements
//
// this program is a modified version of one by Roy Stogner,
// University of Texas at Austin

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"


template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    for (unsigned int i = 0; i < fe.dofs_per_cell; i++)
      if (fe.has_support_on_face(i, face))
        deallog << "Basis function " << i << " has support on face " << face
                << std::endl;
}


#define check_el(fe)             \
  {                              \
    deallog << #fe << std::endl; \
    check(fe);                   \
  }

template <int dim>
void
check()
{
  deallog << "************ dim = " << dim << std::endl;
  check_el(FE_Q<dim>(1));
  check_el(FE_Q<dim>(2));
  check_el(FE_Q<dim>(3));

  check_el(FE_DGQ<dim>(0));
  check_el(FE_DGQ<dim>(1));
  check_el(FE_DGQ<dim>(2));
  check_el(FE_DGQ<dim>(3));

  check_el(FE_DGP<dim>(0));
  check_el(FE_DGP<dim>(1));
  check_el(FE_DGP<dim>(2));

  if (dim > 1)
    check_el(FE_Nedelec<dim>(0));

  check_el(FESystem<dim>(FE_Q<dim>(1), 2));
  check_el(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));

  if (dim > 1)
    check_el(FESystem<dim>(FE_Nedelec<dim>(0), 2));
}


int
main()
{
  initlog();

  //  check<1> ();
  check<2>();
  //  check<3> ();

  return 0;
}
