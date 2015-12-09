// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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
 * @defgroup Concepts Concepts, or expectations on template parameters
 *
 * Sometimes imposing constraints on the type of an object without requiring
 * it to belong to a specific inheritance hierarchy is useful. These are
 * usually referred to as <em>concepts</em> in the C++ community. This module
 * lists the concepts commonly used in deal.II with brief descriptions of
 * their intent. The convention in deal.II for listing constraints on a type
 * is to provide the name of the concept as a <code>typename</code> in a
 * template: for example, the type of a Vector depends on the type of the
 * underlying field, but Vector is not a field element, so it is defined as a
 * template:
 * @code
 * template <typename Number>
 * class Vector;
 * @endcode
 * Where @ref ConceptNumber "Number" is understood to be an appropriate field
 * element type.
 *
 * <dl>
 *
 * <dt class="concepts">@anchor ConceptContainerType <b>ContainerType</b></dt>
 *
 * <dd>

 * There are several algorithms (e.g.,
 * GridTools::find_active_cell_around_point) in deal.II that can operate on
 * either a Triangulation or a DoFHandler, as both classes may be considered
 * to be collections of cells: see the @ref GlossMeshAsAContainer
 * "glossary entry" for a further discussion of this idea. %Functions that may
 * be called with either class indicate this by accepting a template parameter
 * like
 * @code
 * template <template <int, int> class Container>
 * @endcode
 * or
 * @code
 * template <typename Container>
 * @endcode
 * which is usually required to have a <code>typedef</code> named
 * <code>active_cell_iterator</code>.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptDoFHandlerType <b>DoFHandlerType</b></dt>
 *
 * <dd>
 * deal.II includes both DoFHandler and hp::DoFHandler as objects which manage
 * degrees of freedom on a mesh. Though the two do not share any sort of
 * inheritance relationship, they are similar enough that many functions just
 * need something which resembles a DoFHandler to work correctly.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptMatrixType <b>MatrixType</b></dt>
 *
 * <dd>
 * Many functions and classes in deal.II require an object which knows how to
 * calculate matrix-vector products (the member function <code>vmult</code>),
 * transposed matrix-vector products (the member function
 * <code>Tvmult</code>), as well as the `multiply and add' equivalents
 * <code>vmult_add</code> and <code>Tvmult_add</code>. Some functions only
 * require <code>vmult</code> and <code>Tvmult</code>, but an object should
 * implement all four member functions if the template requires a MatrixType
 * argument. Writing classes that satisfy these conditions is a sufficiently
 * common occurrence that the LinearOperator class was written to make things
 * easier; see @ref LAOperators for more information.
 *
 * One way to think of <code>MatrixType</code> is to pretend it is a base
 * class with the following signature (this is nearly the interface provided
 * by SparseMatrix):
 *
 * @code
 * class MatrixType
 * {
 * public:
 *   template <typename VectorType>
 *   virtual void vmult(VectorType &u, const VectorType &v) const =0;
 *
 *   template <typename VectorType>
 *   virtual void Tvmult(VectorType &u, const VectorType &v) const =0;
 *
 *   template <typename VectorType>
 *   virtual void vmult_add(VectorType &u, const VectorType &v) const =0;
 *
 *   template <typename VectorType>
 *   virtual void Tvmult_add(VectorType &u, const VectorType &v) const =0;
 * };
 * @endcode
 *
 * Template functions in C++ cannot be virtual (which is the main reason why
 * this approach is not used in deal.II), so implementing this interface with
 * inheritance will not work, but it is still a good way to think about this
 * template concept. One can use the PointerMatrixAux class to implement
 * <code>vmult_add</code> and <code>Tvmult_add</code> instead of implementing
 * them manually.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptNumber <b>Number</b></dt>
 *
 * <dd>
 * This concept describes scalars which make sense as vector or matrix
 * entries, which is usually some finite precision approximation of a field
 * element. The canonical examples are <code>double</code> and
 * <code>float</code>, but deal.II supports <code>std::complex&lt;T&gt;</code>
 * for floating point type <code>T</code> in many places as well.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptPolynomialType <b>PolynomialType</b></dt>
 *
 * <dd>
 * See the description in @ref Polynomials for more information. In some
 * contexts, anything that satisfies the interface resembling
 * @code
 * template <int dim>
 * class PolynomialType
 * {
 *   virtual void compute (const Point<dim>            &unit_point,
 *                         std::vector<Tensor<1,dim> > &values,
 *                         std::vector<Tensor<2,dim> > &grads,
 *                         std::vector<Tensor<3,dim> > &grad_grads) const =0;
 * }
 * @endcode
 *
 * may be considered as a polynomial for the sake of implementing finite
 * elements.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptPreconditionerType <b>PreconditionerType</b></dt>
 *
 * <dd>
 * This is essentially a synonym for <code>MatrixType</code>, but usually only
 * requires that <code>vmult()</code> and <code>Tvmult()</code> be
 * defined. Most of the time defining <code>Tvmult()</code> is not
 * necessary. One should think of <code>vmult()</code> as applying some
 * approximation of the inverse of a linear operator to a vector, instead of
 * the action of a linear operator to a vector, for the preconditioner
 * classes.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptRelaxationType <b>RelaxationType</b></dt>
 *
 * <dd>
 * This is an object capable of relaxation for multigrid methods. One can
 * think of an object satisfying this constraint as having the following
 * interface as well as the constraints required by
 * @ref ConceptMatrixType "MatrixType":
 * @code
 * class RelaxationType
 * {
 * public:
 *   template <typename VectorType>
 *   virtual void step(VectorType &u, const VectorType &v) const =0;
 *
 *   template <typename VectorType>
 *   virtual void Tstep(VectorType &u, const VectorType &v) const =0;
 * };
 * @endcode
 * where these two member functions perform one step (or the transpose of such
 * a step) of the smoothing scheme. In other words, the operation performed by
 * these functions is
 * $dst = dst - P^{-1} (A dst - rhs)$ and $dst = dst - P^{-T} (A dst - rhs)$.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptSparsityPatternType <b>SparsityPatternType</b></dt>
 *
 * <dd>
 * Almost all functions (with the notable exception of
 * SparsityTools::distribute_sparsity_pattern) which take a sparsity pattern
 * as an argument can take either a regular SparsityPattern or a
 * DynamicSparsityPattern, or even one of the block sparsity patterns. See
 * @ref Sparsity for more information.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptStreamType <b>StreamType</b></dt>
 *
 * <dd>
 * Deriving new stream classes in C++ is well-known to be difficult. To get
 * around this, some functions accept a parameter which defines
 * <code>operator&lt;&lt;</code>, which allows for easy output to any kind of
 * output stream.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptVectorType <b>VectorType</b></dt>
 *
 * <dd>
 * deal.II supports many different vector classes, including bindings to
 * vectors in other libraries. These are similar to standard library vectors
 * (i.e., they define <code>begin()</code>, <code>end()</code>,
 * <code>operator[]</code>, and <code>size()</code>) but also define numerical
 * operations like <code>add()</code>. Some examples of VectorType include
 * Vector, TrilinosWrappers::MPI::Vector, and BlockVector.
 * </dd>
 *
 * </dl>
 */
