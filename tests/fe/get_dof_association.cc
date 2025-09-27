// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
