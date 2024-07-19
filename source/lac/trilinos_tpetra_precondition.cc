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
    template class PreconditionIdentity<float>;
    template class PreconditionIfpackBase<float>;
    template class PreconditionIfpack<float>;
    template class PreconditionJacobi<float>;
    template class PreconditionL1Jacobi<float>;
    template class PreconditionL1GaussSeidel<float>;
    template class PreconditionSOR<float>;
    template class PreconditionSSOR<float>;
    template class PreconditionBlockJacobi<float>;
    template class PreconditionBlockSOR<float>;
    template class PreconditionBlockSSOR<float>;
    template class PreconditionChebyshev<float>;
    template class PreconditionILU<float>;
    template class PreconditionILUT<float>;
#      endif

#      ifdef HAVE_TPETRA_INST_DOUBLE
    template class PreconditionBase<double>;
    template class PreconditionIdentity<double>;
    template class PreconditionIfpackBase<double>;
    template class PreconditionIfpack<double>;
    template class PreconditionJacobi<double>;
    template class PreconditionL1Jacobi<double>;
    template class PreconditionL1GaussSeidel<double>;
    template class PreconditionSOR<double>;
    template class PreconditionSSOR<double>;
    template class PreconditionBlockJacobi<double>;
    template class PreconditionBlockSOR<double>;
    template class PreconditionBlockSSOR<double>;
    template class PreconditionChebyshev<double>;
    template class PreconditionILU<double>;
    template class PreconditionILUT<double>;
#      endif
#      ifdef DEAL_II_WITH_COMPLEX_VALUES
#        ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class PreconditionBase<std::complex<float>>;
    template class PreconditionIdentity<std::complex<float>>;
    template class PreconditionIfpackBase<std::complex<float>>;
    template class PreconditionIfpack<std::complex<float>>;
    template class PreconditionJacobi<std::complex<float>>;
    template class PreconditionL1Jacobi<std::complex<float>>;
    template class PreconditionL1GaussSeidel<std::complex<float>>;
    template class PreconditionSOR<std::complex<float>>;
    template class PreconditionSSOR<std::complex<float>>;
    template class PreconditionBlockJacobi<std::complex<float>>;
    template class PreconditionBlockSOR<std::complex<float>>;
    template class PreconditionBlockSSOR<std::complex<float>>;
    template class PreconditionChebyshev<std::complex<float>>;
    template class PreconditionILU<std::complex<float>>;
    template class PreconditionILUT<std::complex<float>>;
#        endif
#        ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class PreconditionBase<std::complex<double>>;
    template class PreconditionIdentity<std::complex<double>>;
    template class PreconditionIfpackBase<std::complex<double>>;
    template class PreconditionIfpack<std::complex<double>>;
    template class PreconditionJacobi<std::complex<double>>;
    template class PreconditionL1Jacobi<std::complex<double>>;
    template class PreconditionL1GaussSeidel<std::complex<double>>;
    template class PreconditionSOR<std::complex<double>>;
    template class PreconditionSSOR<std::complex<double>>;
    template class PreconditionBlockJacobi<std::complex<double>>;
    template class PreconditionBlockSOR<std::complex<double>>;
    template class PreconditionBlockSSOR<std::complex<double>>;
    template class PreconditionChebyshev<std::complex<double>>;
    template class PreconditionILU<std::complex<double>>;
    template class PreconditionILUT<std::complex<double>>;
#        endif
#      endif

#    endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
