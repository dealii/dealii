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

// This is the test function for FiniteElementData with non-local dofs, where
// the size of dofs_per_object is equal to dim+2. The test function output the
// total number of dofs per cell and the number of non-local dofs per cell.


#include <deal.II/fe/fe_base.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test(const std::vector<std::vector<unsigned int>> &dof_per_object_vecs)
{
  for (auto dof_per_object : dof_per_object_vecs)
    {
      FiniteElementData<dim> fe_data(dof_per_object,
                                     1,
                                     1,
                                     FiniteElementData<dim>::H2);
      deallog << "total number of dofs per cell is "
              << fe_data.n_dofs_per_cell()
              << " and number of non-local dofs per cell is "
              << fe_data.n_non_local_dofs_per_cell() << std::endl;
    }
}

int
main()
{
  initlog();

  deallog << "2d elements: " << std::endl;

  std::vector<unsigned int> dof_per_object1{0, 0, 0, 4};
  std::vector<unsigned int> dof_per_object2{0, 0, 4, 10};
  std::vector<unsigned int> dof_per_object3{0, 0, 0, 4};
  std::vector<unsigned int> dof_per_object4{1, 1, 4, 4};

  test<2>({dof_per_object1, dof_per_object2, dof_per_object3, dof_per_object4});

  deallog << "3d elements: " << std::endl;

  std::vector<unsigned int> dof_per_object5{0, 0, 0, 0, 8};
  std::vector<unsigned int> dof_per_object6{0, 0, 4, 1, 10};
  std::vector<unsigned int> dof_per_object7{0, 0, 0, 2, 4};
  std::vector<unsigned int> dof_per_object8{1, 1, 4, 3, 4};

  test<3>({dof_per_object5, dof_per_object6, dof_per_object7, dof_per_object8});
}