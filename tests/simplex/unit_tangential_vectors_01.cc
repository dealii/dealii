// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test ReferenceCell::unit_tangential_vectors() and ::unit_normal_vectors()
// for all reference-cell types.


#include <deal.II/grid/reference_cell.h>

#include <deal.II/simplex/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  for (const auto face_no :
       internal::Info::get_cell(reference_cell).face_indices())
    {
      deallog << reference_cell.template unit_normal_vectors<dim>(face_no)
              << std::endl;

      for (unsigned int i = 0; i < dim - 1; ++i)
        deallog << reference_cell.template unit_tangential_vectors<dim>(face_no,
                                                                        i)
                << std::endl;
    }
  deallog << std::endl;
}



int
main()
{
  initlog();

  test<2>(ReferenceCell::Tri);
  test<2>(ReferenceCell::Quad);
  test<3>(ReferenceCell::Tet);
  test<3>(ReferenceCell::Pyramid);
  test<3>(ReferenceCell::Wedge);
  test<3>(ReferenceCell::Hex);
}
