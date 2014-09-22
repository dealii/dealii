// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



/**
 * @defgroup VMemory Vector memory management
 *
 * This module groups a few classes that are used to avoid allocating and
 * deallocating vectors over and over in iterative procedures. These methods
 * all use an object of the base class VectorMemory to get their auxiliary
 * vectors.
 *
 * Some discussion on this topic can be found in the discussion of the
 * InverseMatrix class in step-20.
 *
 * @ingroup LAC
 */
