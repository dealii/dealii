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
    template class SolverDirectBase<float, MemorySpace::Host>;
    template class SolverDirect<float, MemorySpace::Host>;
    template class SolverDirectKLU2<float, MemorySpace::Host>;

    template class SolverDirectBase<float, MemorySpace::Default>;
    template class SolverDirect<float, MemorySpace::Default>;
    template class SolverDirectKLU2<float, MemorySpace::Default>;
#      endif

#      ifdef HAVE_TPETRA_INST_DOUBLE
    template class SolverDirectBase<double, MemorySpace::Host>;
    template class SolverDirect<double, MemorySpace::Host>;
    template class SolverDirectKLU2<double, MemorySpace::Host>;

    template class SolverDirectBase<double, MemorySpace::Default>;
    template class SolverDirect<double, MemorySpace::Default>;
    template class SolverDirectKLU2<double, MemorySpace::Default>;
#      endif
#      ifdef DEAL_II_WITH_COMPLEX_VALUES
#        ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class SolverDirectBase<std::complex<float>, MemorySpace::Host>;
    template class SolverDirect<std::complex<float>, MemorySpace::Host>;
    template class SolverDirectKLU2<std::complex<float>, MemorySpace::Host>;

    template class SolverDirectBase<std::complex<float>, MemorySpace::Default>;
    template class SolverDirect<std::complex<float>, MemorySpace::Default>;
    template class SolverDirectKLU2<std::complex<float>, MemorySpace::Default>;
#        endif
#        ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class SolverDirectBase<std::complex<double>, MemorySpace::Host>;
    template class SolverDirect<std::complex<double>, MemorySpace::Host>;
    template class SolverDirectKLU2<std::complex<double>, MemorySpace::Host>;

    template class SolverDirectBase<std::complex<double>, MemorySpace::Default>;
    template class SolverDirect<std::complex<double>, MemorySpace::Default>;
    template class SolverDirectKLU2<std::complex<double>, MemorySpace::Default>;
#        endif
#      endif

#    endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_AMESOS2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
