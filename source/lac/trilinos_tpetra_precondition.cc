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
    template class PreconditionBase<float, MemorySpace::Host>;
    template class PreconditionIdentity<float, MemorySpace::Host>;
    template class PreconditionIfpackBase<float, MemorySpace::Host>;
    template class PreconditionIfpack<float, MemorySpace::Host>;
    template class PreconditionJacobi<float, MemorySpace::Host>;
    template class PreconditionL1Jacobi<float, MemorySpace::Host>;
    template class PreconditionL1GaussSeidel<float, MemorySpace::Host>;
    template class PreconditionSOR<float, MemorySpace::Host>;
    template class PreconditionSSOR<float, MemorySpace::Host>;
    template class PreconditionBlockJacobi<float, MemorySpace::Host>;
    template class PreconditionBlockSOR<float, MemorySpace::Host>;
    template class PreconditionBlockSSOR<float, MemorySpace::Host>;
    template class PreconditionChebyshev<float, MemorySpace::Host>;
    template class PreconditionILU<float, MemorySpace::Host>;
    template class PreconditionILUT<float, MemorySpace::Host>;

    template class PreconditionBase<float, MemorySpace::Default>;
    template class PreconditionIdentity<float, MemorySpace::Default>;
    template class PreconditionIfpackBase<float, MemorySpace::Default>;
    template class PreconditionIfpack<float, MemorySpace::Default>;
    template class PreconditionJacobi<float, MemorySpace::Default>;
    template class PreconditionL1Jacobi<float, MemorySpace::Default>;
    template class PreconditionL1GaussSeidel<float, MemorySpace::Default>;
    template class PreconditionSOR<float, MemorySpace::Default>;
    template class PreconditionSSOR<float, MemorySpace::Default>;
    template class PreconditionBlockJacobi<float, MemorySpace::Default>;
    template class PreconditionBlockSOR<float, MemorySpace::Default>;
    template class PreconditionBlockSSOR<float, MemorySpace::Default>;
    template class PreconditionChebyshev<float, MemorySpace::Default>;
    template class PreconditionILU<float, MemorySpace::Default>;
    template class PreconditionILUT<float, MemorySpace::Default>;
#      endif

#      ifdef HAVE_TPETRA_INST_DOUBLE
    template class PreconditionBase<double, MemorySpace::Host>;
    template class PreconditionIdentity<double, MemorySpace::Host>;
    template class PreconditionIfpackBase<double, MemorySpace::Host>;
    template class PreconditionIfpack<double, MemorySpace::Host>;
    template class PreconditionJacobi<double, MemorySpace::Host>;
    template class PreconditionL1Jacobi<double, MemorySpace::Host>;
    template class PreconditionL1GaussSeidel<double, MemorySpace::Host>;
    template class PreconditionSOR<double, MemorySpace::Host>;
    template class PreconditionSSOR<double, MemorySpace::Host>;
    template class PreconditionBlockJacobi<double, MemorySpace::Host>;
    template class PreconditionBlockSOR<double, MemorySpace::Host>;
    template class PreconditionBlockSSOR<double, MemorySpace::Host>;
    template class PreconditionChebyshev<double, MemorySpace::Host>;
    template class PreconditionILU<double, MemorySpace::Host>;
    template class PreconditionILUT<double, MemorySpace::Host>;

    template class PreconditionBase<double, MemorySpace::Default>;
    template class PreconditionIdentity<double, MemorySpace::Default>;
    template class PreconditionIfpackBase<double, MemorySpace::Default>;
    template class PreconditionIfpack<double, MemorySpace::Default>;
    template class PreconditionJacobi<double, MemorySpace::Default>;
    template class PreconditionL1Jacobi<double, MemorySpace::Default>;
    template class PreconditionL1GaussSeidel<double, MemorySpace::Default>;
    template class PreconditionSOR<double, MemorySpace::Default>;
    template class PreconditionSSOR<double, MemorySpace::Default>;
    template class PreconditionBlockJacobi<double, MemorySpace::Default>;
    template class PreconditionBlockSOR<double, MemorySpace::Default>;
    template class PreconditionBlockSSOR<double, MemorySpace::Default>;
    template class PreconditionChebyshev<double, MemorySpace::Default>;
    template class PreconditionILU<double, MemorySpace::Default>;
    template class PreconditionILUT<double, MemorySpace::Default>;
#      endif
#      ifdef DEAL_II_WITH_COMPLEX_VALUES
#        ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class PreconditionBase<std::complex<float>, MemorySpace::Host>;
    template class PreconditionIdentity<std::complex<float>, MemorySpace::Host>;
    template class PreconditionIfpackBase<std::complex<float>,
                                          MemorySpace::Host>;
    template class PreconditionIfpack<std::complex<float>, MemorySpace::Host>;
    template class PreconditionJacobi<std::complex<float>, MemorySpace::Host>;
    template class PreconditionL1Jacobi<std::complex<float>, MemorySpace::Host>;
    template class PreconditionL1GaussSeidel<std::complex<float>,
                                             MemorySpace::Host>;
    template class PreconditionSOR<std::complex<float>, MemorySpace::Host>;
    template class PreconditionSSOR<std::complex<float>, MemorySpace::Host>;
    template class PreconditionBlockJacobi<std::complex<float>,
                                           MemorySpace::Host>;
    template class PreconditionBlockSOR<std::complex<float>, MemorySpace::Host>;
    template class PreconditionBlockSSOR<std::complex<float>,
                                         MemorySpace::Host>;
    template class PreconditionChebyshev<std::complex<float>,
                                         MemorySpace::Host>;
    template class PreconditionILU<std::complex<float>, MemorySpace::Host>;
    template class PreconditionILUT<std::complex<float>, MemorySpace::Host>;

    template class PreconditionBase<std::complex<float>, MemorySpace::Default>;
    template class PreconditionIdentity<std::complex<float>,
                                        MemorySpace::Default>;
    template class PreconditionIfpackBase<std::complex<float>,
                                          MemorySpace::Default>;
    template class PreconditionIfpack<std::complex<float>,
                                      MemorySpace::Default>;
    template class PreconditionJacobi<std::complex<float>,
                                      MemorySpace::Default>;
    template class PreconditionL1Jacobi<std::complex<float>,
                                        MemorySpace::Default>;
    template class PreconditionL1GaussSeidel<std::complex<float>,
                                             MemorySpace::Default>;
    template class PreconditionSOR<std::complex<float>, MemorySpace::Default>;
    template class PreconditionSSOR<std::complex<float>, MemorySpace::Default>;
    template class PreconditionBlockJacobi<std::complex<float>,
                                           MemorySpace::Default>;
    template class PreconditionBlockSOR<std::complex<float>,
                                        MemorySpace::Default>;
    template class PreconditionBlockSSOR<std::complex<float>,
                                         MemorySpace::Default>;
    template class PreconditionChebyshev<std::complex<float>,
                                         MemorySpace::Default>;
    template class PreconditionILU<std::complex<float>, MemorySpace::Default>;
    template class PreconditionILUT<std::complex<float>, MemorySpace::Default>;
#        endif
#        ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class PreconditionBase<std::complex<double>, MemorySpace::Host>;
    template class PreconditionIdentity<std::complex<double>,
                                        MemorySpace::Host>;
    template class PreconditionIfpackBase<std::complex<double>,
                                          MemorySpace::Host>;
    template class PreconditionIfpack<std::complex<double>, MemorySpace::Host>;
    template class PreconditionJacobi<std::complex<double>, MemorySpace::Host>;
    template class PreconditionL1Jacobi<std::complex<double>,
                                        MemorySpace::Host>;
    template class PreconditionL1GaussSeidel<std::complex<double>,
                                             MemorySpace::Host>;
    template class PreconditionSOR<std::complex<double>, MemorySpace::Host>;
    template class PreconditionSSOR<std::complex<double>, MemorySpace::Host>;
    template class PreconditionBlockJacobi<std::complex<double>,
                                           MemorySpace::Host>;
    template class PreconditionBlockSOR<std::complex<double>,
                                        MemorySpace::Host>;
    template class PreconditionBlockSSOR<std::complex<double>,
                                         MemorySpace::Host>;
    template class PreconditionChebyshev<std::complex<double>,
                                         MemorySpace::Host>;
    template class PreconditionILU<std::complex<double>, MemorySpace::Host>;
    template class PreconditionILUT<std::complex<double>, MemorySpace::Host>;

    template class PreconditionBase<std::complex<double>, MemorySpace::Default>;
    template class PreconditionIdentity<std::complex<double>,
                                        MemorySpace::Default>;
    template class PreconditionIfpackBase<std::complex<double>,
                                          MemorySpace::Default>;
    template class PreconditionIfpack<std::complex<double>,
                                      MemorySpace::Default>;
    template class PreconditionJacobi<std::complex<double>,
                                      MemorySpace::Default>;
    template class PreconditionL1Jacobi<std::complex<double>,
                                        MemorySpace::Default>;
    template class PreconditionL1GaussSeidel<std::complex<double>,
                                             MemorySpace::Default>;
    template class PreconditionSOR<std::complex<double>, MemorySpace::Default>;
    template class PreconditionSSOR<std::complex<double>, MemorySpace::Default>;
    template class PreconditionBlockJacobi<std::complex<double>,
                                           MemorySpace::Default>;
    template class PreconditionBlockSOR<std::complex<double>,
                                        MemorySpace::Default>;
    template class PreconditionBlockSSOR<std::complex<double>,
                                         MemorySpace::Default>;
    template class PreconditionChebyshev<std::complex<double>,
                                         MemorySpace::Default>;
    template class PreconditionILU<std::complex<double>, MemorySpace::Default>;
    template class PreconditionILUT<std::complex<double>, MemorySpace::Default>;
#        endif
#      endif

#    endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA
