// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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



// test internal typetraits used in FEEvaluation

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"


int
main()
{
  initlog();

  deallog << "has_local_element:" << std::endl
          << "LinearAlgebra::distributed::Vector = "
          << internal::has_local_element<
               LinearAlgebra::distributed::Vector<double>>::value
          << std::endl
          << "TrilinosWrappers::MPI::Vector = "
          << internal::has_local_element<TrilinosWrappers::MPI::Vector>::value
          << std::endl;

  deallog << "OK" << std::endl;
}
