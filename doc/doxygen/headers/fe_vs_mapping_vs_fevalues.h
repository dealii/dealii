// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**
 * @defgroup FE_vs_Mapping_vs_FEValues How Mapping, FiniteElement, and FEValues work together
 *
 * <h2>Introduction</h2>
 *
 * Most people create finite element (and, potentially, mapping) objects once
 * but then never actually call any member functions on them -- they simply
 * use them for assembly via the FEValues interface. The only other interaction
 * most will have is by reading the FiniteElementData::dofs_per_cell variable,
 * but that is also just set during construction time. In other words, people
 * never observe FiniteElement or Mapping objects actually <i>do</i>
 * anything -- and that is completely by design.
 *
 * This document is therefore for those who are interested in writing finite
 * element or mapping classes and want to understand how FEValues works and
 * interacts with the FiniteElement and Mapping classes. In the following,
 * we will not make a distinction between FEValues (which acts on cells),
 * FEFaceValues (which acts on faces), and FESubfaceValues (which acts on the
 * children of a face of a cell) as they conceptually all work the same.
 * Consequently, the term "FEValues" will be used generally for all three of
 * these classes in the text below.
 *
 *
 * <h2>Who is responsible for what?</h2>
 *
 * Before going into detail about data and control flow, let us define which
 * class is responsible for providing what kind of information.
 *
 * <h3>%FEValues objects</h3>
 *
 * FEValues is an abstraction that derived from the observation that almost
 * everything one ever does in finite element codes only requires the
 * evaluation of finite element shape functions at quadrature points. This
 * could be, for example, the approximation of an integral of the form
 *   @f[
 *     A^K_{ij} = \int_K \nabla \varphi_i(\bf x) \cdot \nabla \varphi_j(\bf x) \; dx
 *   @f]
 * by quadrature
 *   @f[
 *     A^K_{ij} = \sum_q \nabla \varphi_i(\bf x_q) \cdot \nabla \varphi_j(\bf x_q) \;
 *     |\text{det}\; J(\bf x_q)| w_q,
 *   @f]
 * but it is equally valid when wanting to generate graphical output: there we
 * only need to know the values of a finite element field at the vertices
 * of a mesh, and this too can be written as evaluating everything at
 * quadrature points -- these quadrature points are then simply the vertices
 * of the cells (provided, for example, by QTrapez).
 *
 * FEValues' role is to provide a user the values of shape functions, their
 * gradients, etc, at quadrature points. The same is true with some geometric
 * information, e.g., the normal vectors at the quadrature points. To this end,
 * it provides a large number of member functions in the FEValuesBase base
 * class that allow a user to query basically everything one can ask for in
 * regard to shape functions and geometry information, but only at the
 * quadrature points for which the FEValues object was initialized.
 *
 * FEValues does not actually compute this information itself. It really only
 * provides a place to store it, and then orchestrates the interaction
 * between mapping and finite element classes to have them compute what
 * is requested and store the result in the locations provided by
 * FEValues.
 *
 * As a final note, recall that FEValues can provide an incredible array
 * of information, but that almost all of it is not necessary in any given
 * context. For example, to compute the integral above, it is not necessary
 * to know the second derivatives of the shape functions, or to know the
 * normal vectors at quadrature points. To this end, FEValues uses
 * UpdateFlags in its interactions with the Mapping and FiniteElement
 * class to determine what actually needs to be computed. This is discussed
 * in slightly more detail in @ref UpdateFlags.
 *
 *
 * <h3>Mappings</h3>
 *
 * Mappings (i.e., classes derived from the Mapping base class) are responsible
 * for everything that has to do with the mapping from the reference (unit) cell
 * $[0,1]^\text{dim}$ to each of the actual cells
 * $K\subset{\mathbb R}^\text{spacedim}$. This is facilitated by a mapping function
 * $\mathbf F_K:[0,1]^\text{dim} \mapsto K$. The mapping classes therefore
 * implement interfaces that allow evaluating $\mathbf F_K$ to map forward
 * points $\hat{\mathbf x}$ from the reference cell to $K$, and to map backward
 * from the real cell to the reference cell using $\mathbf F_K^{-1}$.
 * Other common operations that mappings provide is to map vectors (which you
 * can think of as vectors attached to a point $\hat{\mathbf x}$ on the
 * reference cell and pointing in certain directions) to their equivalent
 * vectors on the real cell. This is, for example, what one needs to do
 * for the gradients of shape functions: these are vectors defined on the
 * reference cell, and we need to map these gradients to the real cell $K$.
 * Similar operations can also be defined for matrices (tensors of rank 2,
 * as opposed to vectors which are tensors of rank 1) and higher order tensors.
 *
 * Many of these mappings do not only need the map $\mathbf F_K$ itself,
 * but also the gradients of this mapping, typically referred to as the
 * Jacobian $J_K=\hat\nabla \mathbf F_K$, as well as higher derivatives.
 *
 * Since FEValues only ever needs to evaluate these things at quadrature
 * points, mappings do not in general need to provide the ability to
 * evaluate at <i>arbitrary</i> points. Rather, as we will see below, they will
 * be initialized to use a set of quadrature points defined on the
 * reference cell, will then be "re-initialized" for a particular cell,
 * and all further operations will then only require
 * the evaluation of $\mathbf F_K$ at these quadrature points on the
 * real cell.
 *
 * The mapping classes then have the dual role to (i) compute geometric
 * information (e.g., the normal vectors, determinants of the Jacobians, etc)
 * and putting them into the data structures from which FEValues can
 * provide them to the user, and (ii) to provide the support finite
 * elements need to map shape functions and their derivatives from
 * the reference cell to the real cell.
 *
 *
 * <h3>Finite elements</h3>
 *
 * Finite element classes (i.e., classes derived from FiniteElement) are
 * responsible for defining their shape functions, derivatives, and many
 * other aspects on the reference cell, but also for computing the mapped
 * values and derivatives on actual cells (obviously with the help of a
 * mapping object). For the current discussion, only the latter role is
 * important.
 *
 * As with mappings, all that is important for us here is that the finite
 * element classes can provide this information at given quadrature points,
 * and that they can put the computed information into structures provided
 * by FEValues and from which FEValues member functions can then pass
 * it on to the user through the member functions in FEValuesBase.
 *
 *
 * <h2>What to compute?</h2>
 *
 * Let's say a user wants to compute the gradients of shape functions,
 * for example to compute the integral above. Then they would initialize
 * an FEValues object by giving the update_gradients flag (as is done
 * in basically every tutorial program, starting with step-3). What
 * this indicates is that the user expects the FEValues object to be
 * able to provide the gradients of shape functions on the real cell,
 * but expressed no expectation of any other information.
 *
 * FEValues will then first have to find out what the mapping and
 * finite element objects actually require of each other to make this happen.
 * This already happens at the time the FEValues constructor is run.
 * Because the mapping does not depend on the finite element (though the
 * latter does depend on the former), FEValues first asks the finite
 * element via FiniteElement::requires_update_flags() which <i>other</i>
 * pieces of information it also requires to make the user request
 * happen. As an example, if the finite element were of type
 * FE_Q, then it would determine that in order to compute the
 * gradients of the shape functions on the real cell $K$, it will
 * need to compute the gradients of the shape functions on the
 * reference cell (something it can do on its own, without any
 * external help) but that these reference gradients will then have
 * to be multiplied by the inverse of the Jacobian of the mapping,
 * $J^{-1}_K$, at each of the quadrature points. This multiplication
 * is typically referred to as a <i>covariant transformation</i>,
 * and so FE_Q's implementation of FiniteElement::requires_update_flags()
 * function (provided in the intermediate class FE_Poly) will return
 * both the original update_gradients flag as well as
 * update_covariant_transformation.
 *
 * In a second step, the FEValues object will then call the corresponding
 * function in the mapping, Mapping::requires_update_flags() to determine
 * what is required to provide both update_gradients and
 * update_covariant_transformation. The former is not within the realm
 * of the mapping, so is ignored. The latter will typically require
 * the computation of the Jacobian matrix $J_K$ first, which a typical
 * mapping class will indicate by adding update_contravariant_transformation
 * to the list.
 *
 *
 * <h2>Pre-computing things</h2>
 *
 * At this point, the FEValues object has found out the complete
 * set of flags indicating what everyone has to compute to satisfy
 * the user request. The next step, still during the construction
 * of the FEValues object, stems from the realization that
 * many things could be pre-computed once and then re-used every time
 * we move to a real cell. An example would be the fact that to
 * compute the gradients of the shape functions on the real cell,
 * we need to know the gradients of the shape functions on the
 * reference cell (at the quadrature points on the reference cell)
 * and that these will always be the same: every time we visit
 * a new cell, these values will remain the same, so it would be
 * inefficient to re-compute them every time. Similar arguments
 * can be made for some of the information computed by some of
 * the mapping classes.
 *
 * The FEValues object therefore initializes both the mapping and
 * the finite element object it points to, using both the
 * quadrature object and the final set of update flags computed
 * as described in the previous section. This initialization
 * involves pre-computing as much as these classes can already
 * pre-compute given the set of update flags, and then storing
 * this information for later use.
 *
 * The question then arises: where to store this information. In
 * practice, we do not want to store this information in the mapping
 * or finite element object itself, because this would mean that
 * (i) only one FEValues object could use any given mapping or finite
 * element object at a time, and (ii) that these objects could not
 * be used in a multithreaded context.
 *
 * Rather, the approach works like this:
 * - FEValues calls Mapping::get_data() (and FEFaceValues calls
 *   Mapping::get_face_data(), and FESubfaceValues calls
 *   Mapping::get_subface_data()) with the quadrature object and
 *   the final set of update flags. The implementation of these
 *   functions in the classes derived from Mapping will then
 *   allocate an object of a type derived from
 *   Mapping::InternalDataBase where they can store essentially whatever
 *   it is they find useful for later re-use. Mapping::InternalDataBase
 *   itself does not actually provide any member variables of significance,
 *   but it is really left to derived classes what they think they can
 *   usefully pre-compute and store already at this time. If a mapping
 *   has nothing to pre-compute (or the author of the mapping class is
 *   lazy and does not want to think about what could possibly be
 *   pre-computed), then such a class would simply derive its
 *   own InternalData object from Mapping::InternalDataBase without
 *   actually adding any member variables.
 *
 *   The object so produced is then returned to the calling site
 *   in FEValues and stored by the FEValues object. It will be handed
 *   back every time later on the FEValues object wants any information
 *   from the mapping, thereby providing the mapping object the
 *   ability to read the data it had previously stored.
 *
 * - Secondly, FEValues also calls FiniteElement::get_data() (and FEFaceValues
 *   calls Mapping::get_face_data(), and FESubfaceValues calls
 *   Mapping::get_subface_data()), again with the quadrature object and
 *   the final set of update flags. These functions do essentially the
 *   same as their counterparts in the mappings, and again the object
 *   so initialized, this time of a type derived from
 *   FiniteElement::InternalDataBase, will always be given back to the finite
 *   element whenever the FEValues object wants something from the finite
 *   element object at a later time.
 *
 * This approach allows us to use finite element and mapping objects from
 * multiple FEValues objects at the same time, and possibly from multiple
 * threads at the same time. The point is simply that every user of a
 * finite element or mapping object would hold their own, unique, object
 * returned from the <code>get_data()</code> functions, and that everything
 * that ever happens happens on these objects, rather than on the member
 * variables of the mapping or finite element object itself.
 *
 *
 * <h2>Computing on a given cell</h2>
 *
 * All of the previous steps happened at the time the FEValues object
 * was created. Up to this point, all we did was set up data structures,
 * but nothing useful has been computed so far from the perspective of
 * the user. This only happens when FEValues::reinit() is called on
 * a concrete cell $K$.
 *
 * The things FEValues then does are, in this order:
 *
 * - FEValues figures out whether the cell is a translation
 *   or other similarly simple transformation of the previous cell for which
 *   FEValues::reinit() was called. The result of this, stored in a
 *   CellSimilarly::Similarity object will then be passed to mapping and
 *   finite element to potentially simplify some computations. For example,
 *   if the current cell is simply a translation of the previous one, then
 *   there is no need to re-compute the Jacobian matrix $J_K$ of the
 *   mapping (or its inverse) because it will be the same as for the
 *   previous cell.
 *
 * - Next, FEValues::reinit() calls
 *   Mapping::fill_fe_values() (and, obviously,
 *   FEFaceValues calls Mapping::fill_fe_face_values() and
 *   FESubfaceValues calls Mapping::fill_fe_subface_values()). The arguments
 *   to this function include the cell (or face, or subface) which we are
 *   asked to visit, as well as the cell similarity argument from
 *   above, a reference to the object we had previously obtained from
 *   Mapping::get_data(), and a reference to an object of type
 *   internal::FEValues::MappingRelatedData into which the mapping is
 *   supposed to write its results. In particular, it will need to
 *   compute all mapping related information previously specified by
 *   the update flags, and then write them into the output object.
 *   Examples of fields in the output object that the mapping needs
 *   to fill are the computation of JxW values, the computation of
 *   Jacobian matrices and their inverses, and the normal vectors to
 *   cells (if dim is less than spacedim) and faces.
 *
 * - Finally, FEValues::reinit() calls
 *   FiniteElement::fill_fe_values() (and, obviously,
 *   FEFaceValues calls FiniteElement::fill_fe_face_values() and
 *   FESubfaceValues calls FiniteElement::fill_fe_subface_values()). The arguments
 *   to this function include the cell (or face, or subface) which we are
 *   asked to visit, as well as the cell similarity argument from
 *   above, a reference to the object we had previously obtained from
 *   FiniteElement::get_data(), and a reference to an object of type
 *   internal::FEValues::MappingRelatedData into which the mapping is
 *   supposed to write its results.
 *
 *   In addition to these, the FiniteElement::fill_fe_values() function
 *   also receives references to the mapping object in use, as well as the
 *   Mapping::InternalDataBase object we had previously received from
 *   Mapping::get_data(). The reason is that typically, the finite
 *   element wants to map values or gradients of shape functions from the reference
 *   cell to the actual cell, and these mappings are facilitated by the
 *   various Mapping::transform() functions -- which all require a reference
 *   to the internal object that the FEValues object had previously acquired
 *   from the mapping. This is probably best understood by looking at actual code,
 *   and a simple yet instructive example can be found in
 *   FE_Poly::fill_fe_values(), a function that works on general scalar,
 *   polynomial finite element bases.
 *
 *   As with the mapping, the FiniteElement::fill_fe_values() functions then
 *   use whatever information they had previously computed upon construction
 *   of the FEValues object (i.e., when it called FiniteElement::get_data()),
 *   and use this and the functions in the mapping to compute whatever was
 *   requested as specified by the update flags.
 *
 * This all done, we are finally in a position to offer the owner of the
 * FEValues access to the fields originally requested via the update
 * flags.
 *
 * @ingroup feall
 */
