// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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



// Test FE_Nothing::operator==(). The base clase operator should suffice

#include <deal.II/fe/fe_nothing.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;
  deallog << std::boolalpha;
  deallog << (FE_Nothing<dim>(1) == FE_Nothing<dim>(1, false)) << std::endl;
  deallog << (FE_Nothing<dim>(1) == FE_Nothing<dim>(2)) << std::endl;
  deallog << (FE_Nothing<dim>(2, true) == FE_Nothing<dim>(2, false))
          << std::endl;
  deallog << (FE_Nothing<dim>(1, true) == FE_Nothing<dim>(2, true))
          << std::endl;
  if (dim == 2)
    {
      deallog << (FE_Nothing<dim>(ReferenceCells::Quadrilateral, 2, true) ==
                  FE_Nothing<dim>(2, true))
              << std::endl;
      deallog << (FE_Nothing<dim>(ReferenceCells::Triangle, 2, true) ==
                  FE_Nothing<dim>(2, true))
              << std::endl;
    }
  if (dim == 3)
    {
      deallog << (FE_Nothing<dim>(ReferenceCells::Hexahedron, 2, true) ==
                  FE_Nothing<dim>(2, true))
              << std::endl;
      deallog << (FE_Nothing<dim>(ReferenceCells::Tetrahedron, 2, true) ==
                  FE_Nothing<dim>(2, true))
              << std::endl;
      deallog << (FE_Nothing<dim>(ReferenceCells::Wedge, 1, false) ==
                  FE_Nothing<dim>(ReferenceCells::Pyramid, 1, false))
              << std::endl;
      deallog << (FE_Nothing<dim>(ReferenceCells::Wedge, 3) ==
                  FE_Nothing<dim>(3))
              << std::endl;
    }
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
