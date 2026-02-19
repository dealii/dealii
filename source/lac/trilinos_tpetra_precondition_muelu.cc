// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/lac/trilinos_tpetra_precondition_muelu.templates.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU

DEAL_II_NAMESPACE_OPEN

#      ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#        ifdef HAVE_TPETRA_INST_FLOAT
    template class PreconditionAMGMueLu<float, MemorySpace::Host>;
    template class PreconditionAMGMueLu<float, MemorySpace::Default>;
#        endif

#        ifdef HAVE_TPETRA_INST_DOUBLE
    template class PreconditionAMGMueLu<double, MemorySpace::Host>;
    template class PreconditionAMGMueLu<double, MemorySpace::Default>;
#        endif

#        ifdef DEAL_II_WITH_COMPLEX_VALUES
#          ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class PreconditionAMGMueLu<std::complex<float>, MemorySpace::Host>;
    template class PreconditionAMGMueLu<std::complex<float>,
                                        MemorySpace::Default>;
#          endif
#          ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class PreconditionAMGMueLu<std::complex<double>,
                                        MemorySpace::Host>;
    template class PreconditionAMGMueLu<std::complex<double>,
                                        MemorySpace::Default>;
#          endif
#        endif


  } // namespace TpetraWrappers

} // namespace LinearAlgebra
#      endif // def DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#    endif
#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
