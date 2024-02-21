// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
