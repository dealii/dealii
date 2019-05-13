// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// check PetscScalar

#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"


int
main()
{
  initlog();

  if (typeid(PetscScalar) == typeid(double))
    deallog << "double" << std::endl;
  else if (typeid(PetscScalar) == typeid(float))
    deallog << "float" << std::endl;
  else
    Assert(false, ExcNotImplemented());
}
