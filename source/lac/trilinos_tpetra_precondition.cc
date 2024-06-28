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

#include "deal.II/lac/trilinos_tpetra_precondition.templates.h"

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2

DEAL_II_NAMESPACE_OPEN

#    ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#      ifdef HAVE_TPETRA_INST_FLOAT
    template class PreconditionBase<float>;
    template class PreconditionIfpackBase<float>;
    template class PreconditionJacobi<float>;
    template class PreconditionL1Jacobi<float>;
#      endif

#      ifdef HAVE_TPETRA_INST_DOUBLE
    template class PreconditionBase<double>;
    template class PreconditionIfpackBase<double>;
    template class PreconditionJacobi<double>;
    template class PreconditionL1Jacobi<double>;
#      endif
#      ifdef DEAL_II_WITH_COMPLEX_VALUES
#        ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class PreconditionBase<std::complex<float>>;
    template class PreconditionIfpackBase<std::complex<float>>;
    template class PreconditionJacobi<std::complex<float>>;
    template class PreconditionL1Jacobi<std::complex<float>>;
#        endif
#        ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class PreconditonBase<std::complex<double>>;
    template class PreconditionIfpackBase<std::complex<double>>;
    template class PreconditionJacobi<std::complex<double>>;
    template class PreconditionL1Jacobi<std::complex<double>>;
#        endif
#      endif

#    endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
