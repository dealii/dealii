// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

// test FiniteElement::get_dof_association() with a couple of elements


#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      switch (fe.get_associated_geometry_primitive(i))
        {
          case GeometryPrimitive::vertex:
            deallog << 'v';
            break;
          case GeometryPrimitive::line:
            deallog << 'l';
            break;
          case GeometryPrimitive::quad:
            deallog << 'q';
            break;
          case GeometryPrimitive::hex:
            deallog << 'h';
            break;
          default:
            Assert(false, ExcInternalError());
        }
    }
  deallog << std::endl;
}

int
main()
{
  initlog();

  test(FE_Q<1>(2));
  test(FE_Q<2>(2));
  test(FE_Q<3>(2));

  test(FE_Nedelec<2>(0));
  test(FE_Nedelec<3>(0));

  test(FE_Nedelec<2>(1));
  test(FE_Nedelec<3>(1));
}
