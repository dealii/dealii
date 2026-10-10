// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/anasazi_solver.templates.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_TRILINOS_WITH_ANASAZI
namespace LinearAlgebra
{
namespace TpetraWrappers
{
  #ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT
  template class AnasaziSolverBlockDavidson<float, MemorySpace::Host>;
  template class AnasaziSolverBlockKrylovSchur<float, MemorySpace::Host>;
  template class AnasaziSolverGeneralizedDavidson<float, MemorySpace::Host>;
  template class AnasaziSolverLOBPCG<float, MemorySpace::Host>;

  template class AnasaziSolverBlockDavidson<float, MemorySpace::Default>;
  template class AnasaziSolverBlockKrylovSchur<float, MemorySpace::Default>;
  template class AnasaziSolverGeneralizedDavidson<float, MemorySpace::Default>;
  template class AnasaziSolverLOBPCG<float, MemorySpace::Default>;
  #endif

  #ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE
  template class AnasaziSolverBlockDavidson<double, MemorySpace::Host>;
  template class AnasaziSolverBlockKrylovSchur<double, MemorySpace::Host>;
  template class AnasaziSolverGeneralizedDavidson<double, MemorySpace::Host>;
  template class AnasaziSolverLOBPCG<double, MemorySpace::Host>;

  template class AnasaziSolverBlockDavidson<double, MemorySpace::Default>;
  template class AnasaziSolverBlockKrylovSchur<double, MemorySpace::Default>;
  template class AnasaziSolverGeneralizedDavidson<double, MemorySpace::Default>;
  template class AnasaziSolverLOBPCG<double, MemorySpace::Default>;
  #endif

#  ifdef DEAL_II_WITH_COMPLEX_VALUES

  #ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_FLOAT
  template class AnasaziSolverBlockDavidson<std::complex<float>>, MemorySpace::Host>;
  template class AnasaziSolverBlockKrylovSchur<std::complex<float>, MemorySpace::Host>;
  template class AnasaziSolverGeneralizedDavidson<std::complex<float>, MemorySpace::Host>;
  template class AnasaziSolverLOBPCG<std::complex<float>, MemorySpace::Host>;

  template class AnasaziSolverBlockDavidson<std::complex<float>, MemorySpace::Default>;
  template class AnasaziSolverBlockKrylovSchur<std::complex<float>,
                                        MemorySpace::Default>;
  template class AnasaziSolverGeneralizedDavidson<std::complex<float>, MemorySpace::Default>;
  template class AnasaziSolverLOBPCG<std::complex<float>, MemorySpace::Default>;

  #endif 
  #ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_COMPLEX_DOUBLE
  template class AnasaziSolverBlockDavidson<std::complex<double>, MemorySpace::Host>;
  template class AnasaziSolverBlockKrylovSchur<std::complex<double>,
                                        MemorySpace::Host>;
  template class AnasaziSolverGeneralizedDavidson<std::complex<double>, MemorySpace::Host>;
  template class AnasaziSolverLOBPCG<std::complex<double>, MemorySpace::Host>;

  template class AnasaziSolverBlockDavidson<std::complex<double>, MemorySpace::Default>;
  template class AnasaziSolverBlockKrylovSchur<std::complex<double>,
                                        MemorySpace::Default>;
  template class AnasaziSolverGeneralizedDavidson<std::complex<double>, MemorySpace::Default>;
  template class AnasaziSolverLOBPCG<std::complex<double>, MemorySpace::Default>;
  #endif
#  endif

} // namespace TrilinosWrappers
}

#endif

DEAL_II_NAMESPACE_CLOSE
