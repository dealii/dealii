// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2024 by the deal.II authors
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

#include "deal.II/lac/trilinos_tpetra_solver_direct.templates.h"

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_AMESOS2

DEAL_II_NAMESPACE_OPEN

#    ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#      ifdef HAVE_TPETRA_INST_FLOAT
    template class SolverDirectBase<float>;
    template class SolverDirect<float>;
    template class SolverDirectKLU2<float>;
#      endif

#      ifdef HAVE_TPETRA_INST_DOUBLE
    template class SolverDirectBase<double>;
    template class SolverDirect<double>;
    template class SolverDirectKLU2<double>;
#      endif
#      ifdef DEAL_II_WITH_COMPLEX_VALUES
#        ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class SolverDirectBase<std::complex<float>>;
    template class SolverDirect<std::complex<float>>;
    template class SolverDirectKLU2<std::complex<float>>;
#        endif
#        ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class SolverDirectBase<std::complex<double>>;
    template class SolverDirect<std::complex<double>>;
    template class SolverDirectKLU2<std::complex<double>>;
#        endif
#      endif

#    endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_AMESOS2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
