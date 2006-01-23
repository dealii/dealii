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
 * @mainpage
 *
 * This is the main page for class and function documentation in
 * deal.II. Documentation on other aspects of deal.II can be found
 * elsewhere. In particular, there are tutorial programs on the use of this
 * library.
 *
 * Many of the classes in the deal.II library can be grouped into modules (see
 * the Modules entry in the menu at the top of this page). These modules
 * follow roughly the following collaboration diagram that finite element
 * programs follow:
 *
 * @image html collaboration.png "Collaboration diagram for the most important groups of classes in deal.II"
 * @image latex collaboration.eps "Collaboration diagram for the most important groups of classes in deal.II" width=.9\textwidth
 *
 * This classification of groups can be explained as follows:
 * <ol>
 * 
 *   <li> Unit cell: deal.II supports only hypercubes as unit cells, i.e. the
 *   unit cell [0,1] in 1d, the unit square [0,1]^2 in 2d, and the unit cube
 *   [0,1]^3 in 3d. We do not support triangles, tetrahedra, pyramids, or
 *   prisms.
 *   
 *   <li> Triangulation:
 *
 *   <li> Finite Element: Finite element classes describe the
 *   properties of a finite element space as defined on the unit
 *   cell. This includes, for example, how many degrees of freedom are
 *   located at vertices, on lines, or in the interior of cells. In
 *   addition to this, finite element classes of course have to
 *   provide values and gradients of individual shape functions at
 *   points on the unit cell.
 *
 *   The finite element classes are described in the @ref feall module.
 *
 *   <li> Quadrature
 *   <li> DoFHandler:
 *   <li> Mapping:
 *   <li> FEValues:
 *   <li> Linear System:
 *   <li> Linear Solver:
 *   <li> Output:
 * </ol>  
 */
