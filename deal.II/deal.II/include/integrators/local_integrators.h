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
 * @code
 *   template <int dim>
 *   void
 *   cell_matrix (
 *     FullMatrix<double>& M,
 *     const FEValuesBase<dim>& fe,
 *     const double factor = 1.);
 *
 *   template <int dim>
 *   void
 *   cell_residual (
 *     BlockVector<double>* v,
 *     const FEValuesBase<dim>& fe,
 *     const std::vector<Tensor<1,dim> >& input,
 *     const double factor = 1.);
 * @endcode
 *
 * There is typically a pair of functions for the same operator, the
 * function <tt>cell_residual</tt> implementing the mapping of the
 * operator from the finite element space into its dual, and the
 * function <tt>cell_matrix</tt> generating the bilinear form
 * corresponding to the Frechet derivative of <tt>cell_residual</tt>.
 *
 * The first argument of these functions is the return type, which is
 * <ul>
 * <li> FullMatrix&lt;double&gt; for matrices
 * <li> BlockVector&ltdouble&gt; for vectors
 * </ul>
 *
 * The next argument is the FEValuesBase object representing the
 * finite element for integration. If the integrated operator maps
 * from one finite element space into the dual of another (for
 * instance an off-diagonal matrix in a block system), then first the
 * FEValuesBase for the trial space and after this the one for the
 * test space are specified.
 *
 * This list is followed by the set of required data in the order
 * <ol>
 * <li> Data vectors from finite element functions
 * <li> Data vectors from other objects
 * <li> Additional data
 * <li> A factor which is multiplied with the whole result
 * </ol>
 *
 * <h3>Usage</h3>
 *
 * The local integrators can be used wherever a local integration loop
 * would have been implemented instead. The following example is from
 * the implementation of a Stokes solver, using
 * MeshWorker::Assembler::LocalBlocksToGlobalBlocks. The matrices are
 * <ul>
 * <li> 0: The vector Laplacian for the velocity (here with a vector valued element)
 * <li> 1: The divergence matrix
 * <li> 2: The pressure mass matrix used in the preconditioner
 * </ul>
 *
 * With these matrices, the function called by MeshWorker::loop()
 * could be written like
 * @code
using namespace ::dealii:: LocalIntegrators;

template <int dim>
void MatrixIntegrator<dim>::cell(
  MeshWorker::DoFInfo<dim>& dinfo,
  typename MeshWorker::IntegrationInfo<dim>& info)
{
  Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0));
  Differential::div_primal_matrix(dinfo.matrix(1,false).matrix, info.fe_values(0), info.fe_values(1));
  L2::cell_matrix(dinfo.matrix(2,false).matrix, info.fe_values(1));
}
 * @endcode
 * See step-42 for a worked out example of this code.
 */
namespace LocalIntegrators
{
}

DEAL_II_NAMESPACE_CLOSE

#endif
