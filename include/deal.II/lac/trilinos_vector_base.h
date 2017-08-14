// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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

#ifndef dealii_trilinos_vector_base_h
#define dealii_trilinos_vector_base_h

#warning This file is deprecated. Use <deal.II/lac/trilinos_vector.h> instead.

#include <deal.II/base/config.h>
#include <deal.II/lac/trilinos_vector.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  // All functionality from the former base class VectorBase was moved to
  // MPI::Vector as this ended up to be the only derived class.
  // Hence, using the class VectorBase is deprecated and we only provide
  // an alias for backward compatibility.
  typedef MPI::Vector VectorBase DEAL_II_DEPRECATED;
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_trilinos_vector_base_h
