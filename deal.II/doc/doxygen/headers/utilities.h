//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup utilities Utility functions and classes
 *
 * This module simply collects a number of functions and classes that provide
 * general tools for tasks that do not usually have much to do with finite
 * element programs in particular, but happen to be required there just as
 * well.
 */


/**
 * @defgroup data Data storage primitives
 *
 * Here are a few simple classes that help in storage and viewing data. For
 * example, the Table templates allow to use not only arrays of objects (for
 * which one might want to use the std::vector class), but also
 * two-dimensional (rectangular) tables of arbitrary objects, as well as
 * higher-order analogs up to tables with (presently) seven indices.
 *
 * Similarly, the VectorSlice function is a primitive that takes anything that
 * has an interface that resembles a vector (for example the deal.II Vector or
 * the std::vector classes) and presents a view on it as if it were a vector
 * in itself.
 *
 * @ingroup utilities
 */
