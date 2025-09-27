// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2022 by the deal.II authors
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
 * @defgroup UpdateFlags The interplay of UpdateFlags, Mapping, and FiniteElement in FEValues
 *
 * <h2>Introduction</h2>
 *
 * In order to compute local contributions of an individual cell to the global
 * matrix and right hand side, we usually employ two techniques:
 * - First, the integral is transformed from the actual cell $K$ to the
 *   unit/reference cell $\hat K$. For example, for the Laplace equation, we
 *   transform
 *   @f[
 *     A^K_{ij} = \int_K \nabla \varphi_i(\bf x) \cdot \nabla \varphi_j(\bf x) \; dx
 *   @f]
 *   into
 *   @f[
 *     A^K_{ij} =
 *     \int_{\hat K}
 *     \left[ J^{-1}(\hat{\bf x}) \hat \nabla \varphi_i(\hat{\bf x}) \right]
 *     \cdot
 *     \left[ J^{-1}(\hat{\bf x}) \hat \nabla \varphi_j(\hat{\bf x}) \right]
 *     \;
 *     |\textrm{det}\; J(\hat{\bf x})|
 *     \;\;
 *     d\hat x,
 *   @f]
 *   where a hat indicates reference coordinates, and $J(\hat{\bf
 *   x}_q)$ is the Jacobian
 *   $\frac{\partial \bf F_K(\hat{\bf x})}{\partial\bf \hat x}$ of the mapping
 *   $\bf x = \bf F_K(\hat{\bf x})$.
 * - Second, this integral is then approximated through quadrature. This yields
 *   the formula
 *   @f[
 *     A^K_{ij} = \sum_{q}\left[J^{-1}(\hat{\bf x}_q) \hat \nabla \varphi_i(\hat{\bf x}_q)\right] \cdot
 *     \left[J^{-1}(\hat{\bf x}_q) \hat \nabla \varphi_j(\hat{\bf x}_q)\right]\ |\textrm{det}\ J(\hat{\bf x}_q)|
 *     w_q,
 *   @f]
 *   where $q$ indicates the index of the quadrature point, $\hat{\bf x}_q$ its
 *   location on the reference cell, and $w_q$ its weight.
 *
 * In order to evaluate such an expression in an application code, we
 * have to access three different kinds of objects: a quadrature
 * object that describes locations $\hat{\bf x}_q$ and weights $w_q$ of
 * quadrature points on the reference cell; a finite element object that
 * describes the gradients $\hat\nabla \varphi_i(\hat{\bf x}_q)$ of shape
 * functions on the unit cell; and a mapping object that provides the
 * Jacobian as well as its determinant. Dealing with all these
 * objects would be cumbersome and error prone.
 *
 * On the other hand, these three kinds of objects almost always appear together,
 * and it is in fact very rare for deal.II application codes to do anything with
 * quadrature, finite element, or mapping objects besides using them together.
 * For this reason, deal.II uses the FEValues abstraction
 * combining information on the shape functions, the geometry of the actual mesh
 * cell and a quadrature rule on a reference cell. Upon construction it takes one
 * object of each of the three mentioned categories. Later, it can be
 * "re-initialized" for a concrete grid cell and then provides mapped quadrature
 * points and weights, mapped shape function values and derivatives as well as
 * some properties of the transformation from the reference cell to the actual
 * mesh cell.
 *
 * Since computation of any of these values is potentially expensive (for
 * example when using high order mappings with high order elements), the
 * FEValues class only computes what it is explicitly asked for. To this
 * end, it takes a list of flags of type UpdateFlags at construction time
 * specifying which quantities should be updated each time a cell is
 * visited. In the case above, you want the gradients of the shape functions
 * on the real cell, which is encoded by the flag <code>update_gradients</code>,
 * as well as the product of the determinant of the Jacobian times the
 * quadrature weight, which is mnemonically encoded using the
 * term <code>JxW</code> and encoded in the flag <code>update_JxW_values</code>.
 * Because these flags are represented by single bits in integer numbers,
 * producing a <i>set of flags</i> amounts to setting multiple bits
 * in an integer, which is facilitated using the operation
 * <code>update_gradients | update_JxW_values</code> (in other words, and
 * maybe slightly confusingly so, the operation @"this operation <i>and</i> that
 * operation@" is represented by the expression @"single-bit-in-an-integer-for-this-operation
 * <i>binary-or</i> single-bit-in-an-integer-for-that-operation@"). To
 * make operations cheaper, FEValues and the mapping and finite element objects
 * it depends on really only compute those pieces of information that you
 * have specified in the update flags (plus some information necessary to
 * compute what has been specified, see below), and not everything that
 * could possibly be computed on a cell. This optimization makes it much
 * cheaper to iterate over cells for assembly, but it also means that one
 * should take care to provide the minimal set of flags possible.
 *
 * In addition, once you pass a set of flags for what you want, the functions
 * filling the data fields of FEValues are able to distinguish between
 * values that have to be recomputed on each cell (for example mapped
 * gradients) and quantities that do not change from cell to cell (for
 * example the values of shape functions of the usual $Q_p$
 * finite elements at the same quadrature points on different cells; this
 * property does not hold for the shape functions of Raviart-Thomas
 * elements, however, which must be rotated with the local cell).
 * This allows further optimization of the computations underlying assembly.
 *
 *
 * <h2> Tracking dependencies </h2>
 *
 * Let's say you want to compute the Laplace matrix as shown above. In that
 * case, you need to specify the <code>update_gradients</code> flag
 * (for $\nabla\varphi_i(\bf x_q)$) and the <code>update_JxW_values</code>
 * flag (for computing $|\textrm{det}\; J(\bf x_q)|w_q$). Internally, however,
 * the finite element requires the computation of the inverse of the full
 * Jacobian matrix, $J^{-1}(\bf x_q)$ (and not just the determinant of the matrix),
 * and to compute the inverse of the Jacobian, it is also necessary to compute
 * the Jacobian matrix first.
 *
 * Since these are requirements that are not important to the user, it
 * is not necessary to specify this in user code. Rather, given a set
 * of update flags, the FEValues object first asks the finite element
 * object what information it needs to compute in order to satisfy the
 * user's request provided in the update flags. The finite element
 * object may therefore add other flags to the update flags (e.g., in
 * the example above, an FE_Q object would add
 * <code>update_covariant_transformation</code> to the list, since
 * that is the necessary transformation from
 * $\hat\nabla\hat\varphi_i(\hat{\bf x}_q)$ to $\nabla\varphi_i(\bf
 * x_q)$). With these updated flags, FEValues then asks the mapping
 * whether it also wants to add more flags to the list to satisfy the
 * needs of both the user and the finite element object, by calling
 * Mapping::requires_update_flags(). (This procedure of first asking
 * the finite element and then the mapping does not have to be
 * iterated because mappings never require information computed by the
 * finite element classes, while finite element classes typically need
 * information computed by mappings.) Using this final list, the
 * FEValues object then asks both the finite element object and
 * mapping object to create temporary structures into which they can
 * store some temporary information that can be computed once and for
 * all, and these flags will be used when re-computing data on each
 * cell we will visit later on.
 *
 *
 * <h2>Update once or each</h2>
 *
 * As outlined above, we have now determined the final set of things that are
 * necessary to satisfy a user's desired pieces of information as conveyed by
 * the update flags they provided. This information will then typically be queried
 * on every cell the user code visits in a subsequent integration loop.
 *
 * Given that many of the things mappings or finite element classes compute are
 * potentially expensive, FEValues employs a system whereby mappings and finite
 * element objects are encouraged to pre-compute information that can be computed
 * once without reference to a concrete cell, and make use of this when asked
 * to visit a particular cell of the mesh. An example is that the values of
 * the shape functions of the common FE_Q element are defined on the reference
 * cell, and the values on the actual cell are simply exactly the values on
 * the reference cell -- there is consequently no need to evaluate shape functions
 * on every cell, but it is sufficient to do this once at the beginning, store
 * the values somewhere, and when visiting a concrete cell simply copying these
 * values from their temporary location to the output structure. (Note, however,
 * that this is specific to the FE_Q element: this is not the case if we used
 * a FE_RaviartThomas element, since there,
 * computing the values of the shape functions on a cell involves knowing the
 * Jacobian of the mapping which depends on the geometry of the cell we visit;
 * thus, for this element, simply copying pre-computed information is not
 * sufficient to evaluate the values of shape functions on a particular cell.)
 *
 * To accommodate this structure, both mappings and finite element classes
 * may internally split the update flags into two sets commonly referenced as
 * <code>update_once</code> and <code>update_each</code> (though these names
 * do not appear in any public interfaces). The former contains
 * all those pieces of information that can be pre-computed once at the
 * time the FEValues object starts to interact with a mapping or
 * finite element, whereas the latter contains those flags corresponding to
 * things that need to be computed on every cell. For example, if
 * <code>update_flags=update_values</code>, then the FE_Q class will
 * set <code>update_once=update_values</code> and
 * <code>update_each=0</code>, whereas the Raviart-Thomas element will
 * do it the other way around.
 *
 * These sets of flags are intended to be mutually exclusive. There is,
 * on the other hand, nothing that ever provides this decomposition to
 * anything outside the mapping or finite element classes -- it is a purely
 * internal decomposition.
 *
 *
 * <h2>Generation of the actual data</h2>
 *
 * As outlined above, data is computed at two different times: once at
 * the beginning on the reference cell, and once whenever we move to an
 * actual cell. The functions involved in each of these steps are
 * discussed next:
 *
 *
 * <h3>Initialization</h3>
 *
 * Computing data on the reference cell before we even visit the first
 * real cell is a two-step process. First, the constructor of FEValues,
 * FEFaceValues and FESubfaceValues, respectively, need to allow the
 * Mapping and FiniteElement objects to set up internal data
 * structures. These structures are internal in the following sense: the
 * FEValues object asks the finite element and mapping objects to create
 * an object of type FiniteElement::InternalDataBase and
 * Mapping::InternalDataBase each; the actual finite element and mapping
 * class may in fact create objects of a derived type if they wish to
 * store some data beyond what these base classes already provide. The
 * functions involved in this are
 * <ul>
 * <li>Mapping::get_data()
 * <li>Mapping::get_face_data()
 * <li>Mapping::get_subface_data()
 * <li>FiniteElement::get_data()
 * <li>FiniteElement::get_face_data()
 * <li>FiniteElement::get_subface_data()
 * </ul>
 *
 * The FEValues object then takes over ownership of these objects and will
 * destroy them at the end of the FEValues object's lifetime. After this,
 * the FEValues object asks the FiniteElement and Mapping objects to fill
 * these InternalDataBase objects with the data that pertains to what
 * can and needs to be computed on the reference cell. This is done in these
 * functions:
 * <ul>
 * <li>FEValues::initialize()
 * <li>FEFaceValues::initialize()
 * <li>FESubfaceValues::initialize()
 * </ul>
 *
 *
 * <h3>Reinitialization for a mesh cell</h3>
 *
 * Once initialization is over and we call FEValues::reinit, FEFaceValues::reinit
 * or FESubfaceValues::reinit to move to a concrete cell or face, we need
 * to calculate the "update_each" kinds of data. This is done in the following
 * functions:
 * <ul>
 * <li>FEValues::reinit() calls Mapping::fill_fe_values(), then FiniteElement::fill_fe_values()
 * <li>FEFaceValues::reinit() calls Mapping::fill_fe_face_values(), then FiniteElement::fill_fe_face_values()
 * <li>FESubfaceValues::reinit() calls Mapping::fill_fe_subface_values(),
 * then FiniteElement::fill_fe_subface_values()
 * </ul>
 *
 * This is where the actual data fields for FEValues, stored in
 * internal::FEValues::MappingRelatedData and
 * internal::FEValues::FiniteElementRelatedData objects, are
 * computed. These functions call the function in Mapping first, such
 * that all the mapping data required by the finite element is
 * available. Then, the FiniteElement function is called.
 *
 * @ingroup feall
 */
