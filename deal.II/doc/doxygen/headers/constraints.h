//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup constraints Constraints on degrees of freedom
 * @ingroup dofs
 *
 * This module deals with constraints on degrees of
 * freedom. Constraints typically come from several sources, for
 * example:
 * - If you have Dirichlet-type boundary conditions, one usually enforces
 *   them by requiring that that degrees of freedom on the boundary have
 *   particular values, for example $x_12=42$ if the boundary condition
 *   requires that the finite element solution at the location of degree
 *   of freedom 12 has the value 42.
 *
 * This class implements dealing with linear (possibly inhomogeneous)
 * constraints on degrees of freedom. In particular, it handles constraints of
 * the form $x_{i_1} = \sum_{j=2}^M a_{i_j} x_{i_j} + b_i$. In the context of
 * adaptive finite elements, such constraints appear most frequently as
 * "hanging nodes" and for implementing Dirichlet boundary conditions in
 * strong form. The class is meant to deal with a limited number of
 * constraints relative to the total number of degrees of freedom, for example
 * a few per cent up to maybe 30 per cent; and with a linear combination of
 * $M$ other degrees of freedom where $M$ is also relatively small (no larger
 * than at most around the average number of entries per row of a linear
 * system). It is <em>not</em> meant to describe full rank linear systems.
 */
