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

#include <deal.II/lac/trilinos_tpetra_precondition_frosch.templates.h>

#ifdef DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#    ifdef HAVE_TPETRA_INST_FLOAT
    template class PreconditionFROSch<float>;
    template void
    PreconditionFROSch<float, dealii::MemorySpace::Host>::initialize<2, 2>(
      const SparseMatrix<float, dealii::MemorySpace::Host> &,
      DoFHandler<2, 2> &,
      const PreconditionFROSch<float, dealii::MemorySpace::Host>::AdditionalData
        &);
    template void
    PreconditionFROSch<float, dealii::MemorySpace::Host>::initialize<3, 3>(
      const SparseMatrix<float, dealii::MemorySpace::Host> &,
      DoFHandler<3, 3> &,
      const PreconditionFROSch<float, dealii::MemorySpace::Host>::AdditionalData
        &);
#    endif

#    ifdef HAVE_TPETRA_INST_DOUBLE
    template class PreconditionFROSch<double>;
    template void
    PreconditionFROSch<double, dealii::MemorySpace::Host>::initialize<2, 2>(
      const SparseMatrix<double, dealii::MemorySpace::Host> &,
      DoFHandler<2, 2> &,
      const PreconditionFROSch<double,
                               dealii::MemorySpace::Host>::AdditionalData &);
    template void
    PreconditionFROSch<double, dealii::MemorySpace::Host>::initialize<3, 3>(
      const SparseMatrix<double, dealii::MemorySpace::Host> &,
      DoFHandler<3, 3> &,
      const PreconditionFROSch<double,
                               dealii::MemorySpace::Host>::AdditionalData &);
#    endif

#    ifdef DEAL_II_WITH_COMPLEX_VALUES
#      ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    template class PreconditionFROSch<std::complex<float>>;
    template void
    PreconditionFROSch<std::complex<float>, dealii::MemorySpace::Host>::
      initialize<2, 2>(
        const SparseMatrix<std::complex<float>, dealii::MemorySpace::Host> &,
        DoFHandler<2, 2> &,
        const PreconditionFROSch<std::complex<float>,
                                 dealii::MemorySpace::Host>::AdditionalData &);
    template void
    PreconditionFROSch<std::complex<float>, dealii::MemorySpace::Host>::
      initialize<3, 3>(
        const SparseMatrix<std::complex<float>, dealii::MemorySpace::Host> &,
        DoFHandler<3, 3> &,
        const PreconditionFROSch<std::complex<float>,
                                 dealii::MemorySpace::Host>::AdditionalData &);

#      endif
#      ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    template class PreconditionFROSch<std::complex<double>>;
    template void
    PreconditionFROSch<std::complex<double>, dealii::MemorySpace::Host>::
      initialize<2, 2>(
        const SparseMatrix<std::complex<double>, dealii::MemorySpace::Host> &,
        DoFHandler<2, 2> &,
        const PreconditionFROSch<std::complex<double>,
                                 dealii::MemorySpace::Host>::AdditionalData &);
    template void template class PreconditionFROSch<std::complex<double>>;
    PreconditionFROSch<std::complex<double>, dealii::MemorySpace::Host>::
      initialize<3, 3>(
        const SparseMatrix<std::complex<double>, dealii::MemorySpace::Host> &,
        DoFHandler<3, 3> &,
        const PreconditionFROSch<std::complex<double>,
                                 dealii::MemorySpace::Host>::AdditionalData &);

#      endif
#    endif
  } // namespace TpetraWrappers

} // namespace LinearAlgebra

#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH
