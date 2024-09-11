// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_subscriptor_h
#define dealii_subscriptor_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_ref_counting_by_observer_pointer.h>


DEAL_II_NAMESPACE_OPEN

/**
 * A type alias for the EnableRefCountingByObserverPointer class that makes sure
 * the previous name of the class, Subscriptor, continues to be available.
 *
 * @deprecated Use the new name of the class, ObserverPointer, instead.
 */
using Subscriptor DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
  "Use the new name of the class, EnableRefCountingByObserverPointer.") =
  EnableRefCountingByObserverPointer;

DEAL_II_NAMESPACE_CLOSE

#endif
