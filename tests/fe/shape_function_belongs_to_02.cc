// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test FiniteElement::shape_function_belongs_to(). This test checks
// the case where the finite element is not primitive, which requires
// the use of a more complex code path.


#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // Create an element that has enough components to allow for one
  // scalar, one vector, one symmetric tensor, and one tensor. The
  // vector component is represented by a FE_Nedelec, and the tensor
  // component by 'dim' Nedelec elements. These parts are
  // non-primitive, and therefore so is the whole element.
  const FESystem<dim> fe(FE_Q<dim>(1) ^ 1,
                         FE_Nedelec<dim>(1) ^ 1, // vector component
                         FE_Q<dim>(1) ^
                           (dim * (dim + 1) / 2), // symmetric tensor
                         FE_Nedelec<dim>(1) ^ dim // tensor
  );

  deallog << "Testing " << fe.get_name() << std::endl;

  // Check which shape functions belong to which part of the element
  // (in the ordering outlined above in the creation of the element)
  const FEValuesExtractors::Scalar             scalar(0);
  const FEValuesExtractors::Vector             vector(0 + 1);
  const FEValuesExtractors::SymmetricTensor<2> symm_tensor(0 + 1 + dim);
  const FEValuesExtractors::Tensor<2>          tensor(0 + 1 + dim +
                                             (dim * (dim + 1) / 2));

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    deallog << "phi_" << i << " is part of the scalar part: " << std::boolalpha
            << fe.shape_function_belongs_to(i, scalar) << std::endl;

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    deallog << "phi_" << i << " is part of the vector part: " << std::boolalpha
            << fe.shape_function_belongs_to(i, vector) << std::endl;

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    deallog << "phi_" << i
            << " is part of the symmetric tensor part: " << std::boolalpha
            << fe.shape_function_belongs_to(i, symm_tensor) << std::endl;

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    deallog << "phi_" << i << " is part of the tensor part: " << std::boolalpha
            << fe.shape_function_belongs_to(i, tensor) << std::endl;
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
