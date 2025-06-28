// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_NOX

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_epetra_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_tpetra_block_vector.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <deal.II/trilinos/nox.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
// We don't build the nox.inst file if Trilinos isn't configured
// with NOX, but doxygen doesn't know that and tries to find that
// file anyway for parsing -- which then of course it fails on. So
// exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "trilinos/nox.inst"
#  endif
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
