//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup Quadrature Quadrature formulæ
 *
 * This module contains the base class Quadrature as well as the
 * quadrature formulæ provided by deal.II.
 *
 * The class QIterated is used, to construct an iterated quadrature
 * formula out of an existing one, thereby increasing the accuracy of
 * the formula without increasing the order.
 *
 * While the usual quadrature formulæ of higher dimensions
 * generate tensor products which are equal in each direction, the
 * class QAnisotropic generates tensor products of possibly different
 * formulæ in each direction.
 *
 * The class QProjector is not actually a quadrature rule by itself,
 * but it provides functions for computing the quadrature on the
 * surfaces of higher dimensional cells.
 *
 * All other classes in this module actually implement quadrature
 * rules of different order and other characteristics.
 */

/**
 * @defgroup IO Input/Output
 *
 * This module collects the classes used for reading and writing
 * meshes and data.
 *
 * The list of supported formats can be found in the
 * description of the classes
 * <ul>
 * <li> GridIn for reading meshes,
 * <li> GridOut for writing meshes,
 * <li> DataOutBase for writing simulation results.
 * </ul>
 *
 * For writing data, you would normally use objects of the class
 * DataOut (DataOutBase only provides the low level output
 * functions). Still, there are some other options for more
 * specialized needs, like DataOutFaces, DataOutRotation and
 * DataOutStack.
 */
