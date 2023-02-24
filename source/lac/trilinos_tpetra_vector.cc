// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#include <deal.II/lac/trilinos_tpetra_vector.templates.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#  ifdef HAVE_TPETRA_INST_FLOAT
    template class Vector<float>;
#  endif
#  ifdef HAVE_TPETRA_INST_DOUBLE
    template class Vector<double>;
#  endif
#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class Vector<std::complex<float>>;
#    endif
#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class Vector<std::complex<double>>;
#    endif
#  endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
