//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------


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
