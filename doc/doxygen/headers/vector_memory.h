// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



/**
 * @defgroup VMemory Vector memory management
 *
 * This page groups a few classes that are used to avoid allocating and
 * deallocating vectors over and over in iterative procedures. These methods
 * all use an object of the base class VectorMemory to get their auxiliary
 * vectors.
 *
 * Some discussion on this topic can be found in the discussion of the
 * InverseMatrix class in step-20.
 *
 * @ingroup LAC
 */
