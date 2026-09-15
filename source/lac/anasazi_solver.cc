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
namespace TrilinosWrappers
{

  // template class AnasaziSolverBase<float, MemorySpace::Host>;
  // template class SolverBlockKrylovSchur<float, MemorySpace::Host>;

  // template class AnasaziSolverBase<float, MemorySpace::Default>;
  // template class SolverBlockKrylovSchur<float, MemorySpace::Default>;

  template class AnasaziSolverBase<double, MemorySpace::Host>;
  template class SolverBlockKrylovSchur<double, MemorySpace::Host>;

  template class AnasaziSolverBase<double, MemorySpace::Default>;
  template class SolverBlockKrylovSchur<double, MemorySpace::Default>;

#  ifdef DEAL_II_WITH_COMPLEX_VALUES

  template class AnasaziSolverBase<std::complex<float>, MemorySpace::Host>;
  template class SolverBlockKrylovSchur<std::complex<float>, MemorySpace::Host>;

  template class AnasaziSolverBase<std::complex<float>, MemorySpace::Default>;
  template class SolverBlockKrylovSchur<std::complex<float>,
                                        MemorySpace::Default>;


  template class AnasaziSolverBase<std::complex<double>, MemorySpace::Host>;
  template class SolverBlockKrylovSchur<std::complex<double>,
                                        MemorySpace::Host>;

  template class AnasaziSolverBase<std::complex<double>, MemorySpace::Default>;
  template class SolverBlockKrylovSchur<std::complex<double>,
                                        MemorySpace::Default>;

#  endif

} // namespace TrilinosWrappers

#endif

DEAL_II_NAMESPACE_CLOSE
