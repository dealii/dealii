//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_local_integrators_h
#define __deal2__integrators_local_integrators_h

// This file only provides definition and documentation of the
// namespace LocalIntegrators. There is no necessity to include it
// anywhere in a C++ code. Only doxygen will make use of it.

#include <base/config.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @brief Library of integrals over cells and faces
 *
 * This namespace contains application specific local integrals for
 * bilinear forms, forms and error estimates. It is a collection of
 * functions organized into namespaces devoted to certain
 * applications. For instance, the namespace Laplace contains
 * functions for computing cell matrices and cell residuals for the
 * Laplacian operator, as well as functions for the weak boundary
 * conditions by Nitsche or the interior penalty discontinuous
 * Galerkin method. The namespace Maxwell would do the same for
 * curl-curl type problems.
 *
 * There is a namespace Differential containing general
 * first order differential operators like divergence, gradient and
 * curl, and their traces on faces.
 *
 * The namespace L2 contains functions for mass matrices and
 * <i>L<sup>2</sup></i>-inner products.
 *
 * <h3>Signature of functions</h3>
 *
 * Functions in this namespace follow a generic signature. In the
 * simplest case, you have two related functions
 * \begin{code}
 * \end{code} 
 */
namespace LocalIntegrators
{
}

DEAL_II_NAMESPACE_CLOSE

#endif
