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
 * @defgroup Concepts Concepts, or expectations on template parameters
 *
 * Sometimes imposing constraints on the type of an object without requiring
 * it to belong to a specific inheritance hierarchy is useful. These are
 * usually referred to as <em>concepts</em> in the C++ community. This module
 * lists the concepts commonly used in deal.II with brief descriptions of
 * their intent. The convention in deal.II for listing constraints on a type
 * is to provide the name of the concept as a <code>typename</code> in a
 * template: for example, the type of a Vector depends on the type of the
 * underlying field, and so it is defined as a
 * template:
 * @code
 * template <typename Number>
 * class Vector;
 * @endcode
 * The point here is that you are creating a vector that can store
 * elements of type @p Number. But there are some underlying
 * assumptions on this. For example, the deal.II Vector class is not
 * intended to be used just as a collection (unlike
 * <code>std::vector</code>) but defines vector space operations such
 * as addition of vectors, or the norm of vectors. Consequently, the
 * data type users can specify for @p Number must satisfy certain
 * conditions (i.e., it must conform to or "model" a "concept"):
 * Specifically, the type must denote objects that represent the
 * elements of what mathematically call a "field" (which you can think
 * of as, well, "numbers": things we can add, multiply, divide, take
 * the absolute value of, etc). The point of a concept is then to
 * describe <em>what conditions a type must satisfy</em> to be a valid
 * template argument in a given context.
 *
 * This page describes these conditions for a number of concepts used
 * throughout deal.II. Specifically, in the example above, the
 * @ref ConceptNumber "Number concept" discussed below describes the
 * types that could be used as argument for the Vector class.
 *
 * Concepts have been proposed as a language
 * extension to C++ for a long time already. They would allow us to
 * describe that a class or function has certain properties in order
 * to be a qualified template argument. For example, it would allow us
 * to express in C++ code that the first argument to, say,
 * GridTools::find_closest_vertex(), must have a type that represents
 * an actual mesh -- which we can currently only describe in words,
 * see below. Using C++ concepts would allow us to describe this in
 * code and trying to call such a function with an object as first
 * argument that is not, in fact, a mesh would yield a compiler error
 * that makes the mismatch clear.
 *
 * Unfortunately, these proposals to C++ have never made it into any
 * official C++ standard; they are proposed for C++20 however. We may
 * start to use them once the vast majority of our users have
 * compilers that support this standard.
 *
 * More information on the topic can be found at
 * <a href="https://en.wikipedia.org/wiki/Concepts_(C%2B%2B)">this wikipedia page</a>.

 *
 * <dl>
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
 * template concept. One can use the LinearOperator class to implement
 * <code>vmult_add</code> and <code>Tvmult_add</code> instead of implementing
 * them manually.
 * </dd>
 *
 * <dt class="concepts">@anchor ConceptMeshType <b>MeshType</b></dt>
 *
 * <dd>
 * Meshes can be thought of as arrays of vertices and connectivities, but a
 * more fruitful view is to consider them as <i>collections of cells</i>. In
 * C++, collections are often called <i>containers</i> (typical containers are
 * std::vector, std::list, etc.) and they are characterized by the ability to
 * iterate over the elements of the collection. The <tt>MeshType</tt> concept
 * refers to any container which defines appropriate methods (such as
 * DoFHandler::begin_active()) and <tt>typedefs</tt> (such as
 * DoFHandler::active_cell_iterator) for managing collections of cells.
 *
 * Instances of Triangulation, DoFHandler, and hp::DoFHandler may
 * all be considered as containers of cells. In fact, the most important parts
 * of the public interface of these classes consists simply of the ability to
 * get iterators to their elements. Since these parts of the interface are
 * generic, i.e., the functions have the same name in all classes, it is
 * possible to write operations that do not actually care whether they work on
 * a triangulation or a DoF handler object. Examples abound, for example, in
 * the GridTools namespace, underlining the power of the abstraction that
 * meshes and DoF handlers can all be considered simply as collections
 * (containers) of cells.
 *
 * On the other hand, meshes are non-standard containers unlike std::vector or
 * std::list in that they can be sliced several ways. For example, one can
 * iterate over the subset of active cells or over all cells; likewise, cells
 * are organized into levels and one can get iterator ranges for only the
 * cells on one level. Generally, however, all classes that implement the
 * containers-of-cells concept use the same function names to provide the same
 * functionality.
 *
 * %Functions that may be called with either class indicate this by accepting
 * a template parameter like
 * @code
 * template <template <int, int> class MeshType>
 * @endcode
 * or
 * @code
 * template <typename MeshType>
 * @endcode
 * The classes that satisfy this concept are collectively referred to as
 * <em>mesh classes</em>. The exact definition of <tt>MeshType</tt> relies a
 * lot on library internals, but it can be summarized as any class with the
 * following properties:
 * <ol>
 *   <li>A <tt>typedef</tt> named <tt>active_cell_iterator</tt>.
 *   </li>
 *   <li>A method <tt>get_triangulation()</tt> which returns a reference to
 *   the underlying geometrical description (one of the Triangulation classes)
 *   of the collection of cells. If the mesh happens to be a Triangulation,
 *   then the mesh just returns a reference to itself.
 *   </li>
 *   <li>A method <tt>begin_active()</tt> which returns an iterator pointing
 *   to the first active cell.
 *   </li>
 *   <li>A static member value <tt>dimension</tt> containing the dimension in
 *       which the object lives.
 *   </li>
 *   <li>A static member value <tt>space_dimension</tt> containing the dimension
 *       of the object (e.g., a 2D surface in a 3D setting would have
 *       <tt>space_dimension = 2</tt>).
 *   </li>
 * </ol>
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
 * a step) of the smoothing scheme. In other words, the operations performed by
 * these functions are
 * $u = u - P^{-1} (A u - v)$ and $u = u - P^{-T} (A u - v)$.
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
