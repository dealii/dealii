// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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
 * @page DEALGlossary Glossary
 *
 * This glossary explains a few terms that are frequently used in the
 * documentation of classes of deal.II. The glossary often only gives
 * a microscopic view of a particular concept; if you struggle with
 * the bigger picture, it may therefore also be worth to consult the
 * global overview of classes on the @ref index page.
 *
 * <dl>
 *
 * <dt class="glossary">@anchor GlossActive <b>Active cells</b></dt>
 * <dd>A cell, face or edge is defined as <i>active</i> if it is not
 * refined any further, i.e., if it does not have children. Unless
 * working with a multigrid algorithm, active cells are the only
 * ones carrying degrees of freedom.</dd>
 *
 *
 *
 * <dt class="glossary">@anchor GlossArtificialCell <b>Artificial cells</b></dt>
 * <dd>
 * If a mesh is distributed across multiple MPI processes using the
 * parallel::distributed::Triangulation class, each processor stores
 * only the cells it owns, one layer of adjacent cells that are owned
 * by other processors (called @ref GlossGhostCell "ghost cells"), all coarse level
 * cells, and all cells that are necessary to maintain the invariant
 * that adjacent cells must differ by at most one refinement
 * level. The cells stored on each process that are not owned by this
 * process and that are not ghost cells are called "artificial cells",
 * and for these cells the predicate
 * <code>cell-@>is_artificial()</code> returns true. Artificial cells
 * are guaranteed to exist in the globally distributed mesh but they
 * may be further refined on other processors. See the
 * @ref distributed_paper "Distributed Computing paper" for more
 * information.
 *
 * The concept of artificial cells has no meaning for triangulations
 * that store the entire mesh on each processor, i.e. the
 * dealii::Triangulation class.  </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBlockLA <b>Block (linear algebra)</b></dt>

 * <dd>It is often convenient to treat a matrix or vector as a collection of
 * individual blocks. For example, in step-20 (and other tutorial
 * programs), we want to consider the global linear system $Ax=b$ in
 * the form
 * @f{eqnarray*}
  \left(\begin{array}{cc}
    M & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),
 * @f}
 * where $U,P$ are the values of velocity and pressure degrees of freedom,
 * respectively, $M$ is the mass matrix on the velocity space, $B$ corresponds to
 * the negative divergence operator, and $B^T$ is its transpose and corresponds
 * to the negative gradient.
 *
 * Using such a decomposition into blocks, one can then define
 * preconditioners that are based on the individual operators that are
 * present in a system of equations (for example the Schur complement,
 * in the case of step-20), rather than the entire matrix. In essence,
 * blocks are used to reflect the structure of a PDE system in linear
 * algebra, in particular allowing for modular solvers for problems
 * with multiple solution components. On the other hand, the matrix
 * and right hand side vector can also treated as a unit, which is
 * convenient for example during assembly of the linear system when
 * one may not want to make a distinction between the individual
 * components, or for an outer Krylov space solver that doesn't care
 * about the block structure (e.g. if only the preconditioner needs
 * the block structure).
 *
 * Splitting matrices and vectors into blocks is supported by the
 * BlockSparseMatrix, BlockVector, and related classes. See the
 * overview of the various linear algebra classes in the @ref LAC
 * module. The objects present two interfaces: one that makes the
 * object look like a matrix or vector with global indexing
 * operations, and one that makes the object look like a collection of
 * sub-blocks that can be individually addressed. Depending on
 * context, one may wish to use one or the other interface.
 *
 * Typically, one defines the sub-structure of a matrix or vector by
 * grouping the degrees of freedom that make up groups of physical
 * quantities (for example all velocities) into individual blocks of
 * the linear system. This is defined in more detail below in the
 * glossary entry on @ref GlossBlock "Block (finite element)".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBlock <b>Block (finite element)</b></dt>
 * <dd>
 * <i>Intent:</i>
 * Blocks are a generalization of @ref GlossComponent "components" in that
 * they group together one or more components of a vector-valued finite
 * element that one would like to consider jointly. One often wants to do this
 * to define operators that correspond to the structure of a (part of a)
 * differential operator acting on the vector-valued solution, such as the
 * Schur complement solver in step-20, or the block solvers and
 * preconditioners of step-22.
 *
 * For the purpose of a discretization, blocks are the better concept to use
 * since it is not always possible to address individual components of a
 * solution. This is, in particular, the case for non-@ref
 * GlossPrimitive "primitive"
 * elements. Take for instance the solution of the mixed Laplacian
 * system with the FE_RaviartThomas element (see step-20). There, the first
 * <tt>dim</tt> components are the directional velocities. Since the shape
 * functions are linear combinations of those, these <tt>dim</tt> components
 * constitute only a single block. On the other hand, the pressure variable is
 * scalar and would form a the second block, but in the <tt>dim+1</tt>st
 * component.
 *
 * The minimal size of each block is dictated by the underlying finite element
 * (a block consists of a single component for scalar elements, but in the
 * case of the FE_RaviartThomas, for example, a block consists of <tt>dim</tt>
 * components). However, several such minimal blocks can be grouped together
 * into user defined blocks at will, and in accordance with the
 * application. For instance, for the
 * <b>Q</b><sub>2</sub><sup><i>d</i></sup>-<b>Q</b><sub>1</sub> (Taylor-Hood) Stokes
 * element, there are <i>d</i>+1 components each of which could in principle form
 * its own block. But we are typically more interested in having only two
 * blocks, one of which consists of all the velocity vector components
 * (i.e. this block would have <i>d</i> components) and the other having only the
 * single pressure component.
 *
 * <i>Implementation:</i>
 * deal.II has a number of different finite element classes, all of which are
 * derived from the FiniteElement base class
 * (see the @ref feall "module on finite element classes").
 * With one exception, whether they are scalar or
 * vector valued, they all define a single block: all vector components the
 * finite element defines through its FiniteElement::n_components() function
 * form a single block, i.e. FiniteElement::n_blocks() returns one.
 *
 * The exception is the FESystem class that takes multiple simpler elements
 * and connects them into more complicated ones. Consequently, it can have
 * more than one block. A FESystem has as many blocks as it has base elements
 * times their multiplicity (see the constructors of FESystem to understand
 * this statement). In other words, it does not care how many blocks each base
 * element has, and consequently you can produce a Stokes element that has
 * only two blocks by creating the object
 * @code
 *    FESystem<dim> (FESystem<dim> (FE_Q<dim>(2), dim), 1,
 *                   FE_Q<dim>(1), 1);
 * @endcode
 * On the other hand, we could have produced a similar object with dim+1
 * blocks using
 * @code
 *    FESystem<dim> (FE_Q<dim>(2), dim,
 *                   FE_Q<dim>(1), 1);
 * @endcode
 * With the exception of the number of blocks, the two objects are the
 * same for all practical purposes, however.
 *
 * <i>Global degrees of freedom:</i>
 * While we have defined blocks above in terms of the vector components of a
 * vector-valued solution function (or, equivalently, in terms of the
 * vector-valued finite element space), every shape function of a finite
 * element is part of one block or another. Consequently, we can partition all
 * degrees of freedom defined on a DoFHandler into individual blocks. Since by
 * default the DoFHandler class enumerates degrees of freedom in a more or
 * less random way, you will first want to call the
 * DoFRenumbering::component_wise function to make sure that all degrees of
 * freedom that correspond to a single block are enumerated consecutively.
 *
 * If you do this, you naturally partition matrices and vectors into blocks as
 * well (see @ref GlossBlockLA "block (linear algebra)).  In most cases, when
 * you subdivide a matrix or vector into blocks, you do so by creating one
 * block for each block defined by the finite element (i.e. in most practical
 * cases the FESystem object). However, this needs not be so: the
 * DoFRenumbering::component_wise function allows to group several vector
 * components or finite element blocks into the same logical block (see, for
 * example, the @ref step_22 "step-22" or step-31 tutorial programs, as
 * opposed to step-20). As a consequence, using this feature, we can achieve
 * the same result, i.e. subdividing matrices into $2\times 2$ blocks and
 * vectors into 2 blocks, for the second way of creating a Stokes element
 * outlined above using an extra argument as we would have using the first way
 * of creating the Stokes element with two blocks right away.
 *
 * More information on this topic can be found in the documentation of
 * FESystem, the @ref vector_valued module and the tutorial programs
 * referenced therein.
 *
 * <i>Selecting blocks:</i>
 * Many functions allow you to restrict their operation to certain
 * vector components or blocks. For example, this is the case for
 * the functions that interpolate boundary values: one may want
 * to only interpolate the boundary values for the velocity block of
 * a finite element field but not the pressure block. The way to do
 * this is by passing a BlockMask argument to such functions, see the
 * @ref GlossBlockMask "block mask entry of this glossary".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBlockMask <b>Block mask</b></dt>
 *
 * <dd>
 * In much the same way as one can think of elements as being composed
 * of physical vector components (see @ref GlossComponent) or logical
 * blocks (see @ref GlossBlock), there is frequently a need to select
 * a set of such blocks for operations that are not intended to be run
 * on <i>all</i> blocks of a finite element space. Selecting which blocks
 * to work on happens using the BlockMask class.
 *
 * Block masks work in much the same way as component masks, including the
 * fact that the BlockMask class has similar semantics to the ComponentMask
 * class. See @ref GlossComponentMask "the glossary entry on component masks"
 * for more information.
 *
 * @note While components and blocks provide two alternate but equally valid
 * viewpoints on finite elements with multiple vector components, the fact
 * is that throughout the library there are far more places where you can
 * pass a ComponentMask argument rather than a BlockMask argument. Fortunately,
 * one can be converted into the other, using the syntax
 * <code>fe.component_mask(block_mask)</code> where <code>block_mask</code>
 * is a variable of type BlockMask. In other words, if you have a block
 * mask but need to call a function that only accepts a component mask, this
 * syntax can be used to obtain the necessary component mask.
 *
 * <b>Creation of block masks:</b>
 * Block masks are typically created by asking the finite element
 * to generate a block mask from certain selected vector components using
 * code such as this where we create a mask that only denotes the
 * velocity components of a Stokes element (see @ref vector_valued):
 * @code
 *   FESystem<dim> stokes_fe (FESystem<dim>(FE_Q<dim>(2), dim), 1,    // Q2 element for the velocities
 *                            FE_Q<dim>(1),                     1);     // Q1 element for the pressure
 *   FEValuesExtractors::Scalar pressure(dim);
 *   BlockMask pressure_mask = stokes_fe.block_mask (pressure);
 * @endcode
 * The result is a block mask that, in 1d as well as 2d and 3d, would have values
 * <code>[false, true]</code>. Similarly, using
 * @code
 *   FEValuesExtractors::Vector velocities(0);
 *   BlockMask velocity_mask = stokes_fe.block_mask (velocities);
 * @endcode
 * would result in a mask <code>[true, false]</code> in any dimension.
 *
 * Note, however, that if we had defined the finite element in the following
 * way:
 * @code
 *   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,    // Q2 element for the velocities
 *                            FE_Q<dim>(1), 1);     // Q1 element for the pressure
 * @endcode
 * then the code
 * @code
 *   FEValuesExtractors::Scalar pressure(dim);
 *   BlockMask pressure_mask = stokes_fe.block_mask (pressure);
 * @endcode
 * would yield a block mask that in 2d has elements
 * <code>[false, false, true]</code> because the element has
 * <code>dim+1</code> components and equally many blocks. See the
 * discussion on what a block represents exactly in the
 * @ref GlossBlock "block entry of this glossary".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBoundaryForm <b>%Boundary form</b></dt>
 *
 * <dd>For a dim-dimensional triangulation in dim-dimensional space,
 * the boundary form is a vector defined on faces. It is the vector
 * product of the image of coordinate vectors on the surface of the
 * unit cell. It is a vector normal to the surface, pointing outwards
 * and having the length of the surface element.
 *
 * A more general definition would be that (at least up to the length
 * of this vector) it is exactly that vector that is necessary when
 * considering integration by parts, i.e. equalities of the form
 * $\int_\Omega \text{div} \vec \phi = -\int_{\partial\Omega} \vec n
 * \cdot \vec \phi$. Using this definition then also explains what
 * this vector should be in the case of domains (and corresponding
 * triangulations) of dimension <code>dim</code> that are embedded in
 * a space <code>spacedim</code>: in that case, the boundary form is
 * still a vector defined on the faces of the triangulation; it is
 * orthogonal to all tangent directions of the boundary and within the
 * tangent plane of the domain. Note that this is compatible with case
 * <code>dim==spacedim</code> since there the tangent plane is the
 * entire space ${\mathbb R}^\text{dim}$.
 *
 * In either case, the length of the vector equals the determinant of
 * the transformation of reference face to the face of the current
 * cell.  </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBoundaryIndicator <b>%Boundary indicator</b></dt>
 *
 * <dd> In a Triangulation object, every part of the boundary may be
 * associated with a unique number (of type types::boundary_id) that
 * is used to determine what kinds of boundary conditions are to be
 * applied to a particular part of a boundary. The boundary is
 * composed of the faces of the cells and, in 3d, the edges of these
 * faces.
 *
 * By default, all boundary indicators of a mesh are zero, unless you are
 * reading from a mesh file that specifically sets them to something different,
 * or unless you use one of the mesh generation functions in namespace GridGenerator
 * that have a 'colorize' option. A typical piece of code that sets the boundary
 * indicator on part of the boundary to something else would look like
 * this, here setting the boundary indicator to 42 for all faces located at
 * $x=-1$:
 * @code
 *   for (typename Triangulation<dim>::active_cell_iterator
 *          cell = triangulation.begin_active();
 *        cell != triangulation.end();
 *        ++cell)
 *     for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
 *       if (cell->face(f)->at_boundary())
 *         if (cell->face(f)->center()[0] == -1)
 *           cell->face(f)->set_boundary_indicator (42);
 * @endcode
 * This calls functions TriaAccessor::set_boundary_indicator. In 3d, it may
 * also be appropriate to call TriaAccessor::set_all_boundary_indicators instead
 * on each of the selected faces. To query the boundary indicator of a particular
 * face or edge, use TriaAccessor::boundary_indicator.
 *
 * In older versions of the library (prior to 8.2), if you wanted also
 * to change the way the Triangulation class treated the boundary for
 * the purposes of mesh refinement, you could call
 * Triangulation::set_boundary to associate a boundary object with a
 * particular boundary indicator. This method is still supported, and
 * it allows the Triangulation object to use a different method of
 * finding new points on faces and edges to be refined; the default is
 * to use a StraightBoundary object for all faces and edges. The
 * results section of step-49 has a worked example that shows all of
 * this in action.
 *
 * The suggested method from version 8.2 onwards, is to split the
 * geometrical description of the boundary from its physical meaning,
 * by using separately manifold_ids and boundary_ids. The former are
 * used to describe how the geometry changes, and the latter are used
 * to identify the boundary conditions.
 *
 * Many of the functions in namespaces DoFTools and VectorTools take
 * arguments that specify which part of the boundary to work on, and
 * they specifically refer to boundary_ids. Examples are
 * DoFTools::make_periodicity_constraints,
 * DoFTools::extract_boundary_dofs,
 * DoFTools::make_zero_boundary_constraints and
 * VectorTools::interpolate_boundary_values,
 * VectorTools::compute_no_normal_flux_constraints.
 *
 * @note Boundary indicators are inherited from mother faces and edges to
 * their children upon mesh refinement. Some more information about boundary
 * indicators is also presented in a section of the documentation of the
 * Triangulation class.
 *
 * @note For parallel triangulations of type parallel::distributed::Triangulation,
 * it is not enough to set boundary indicators only once at the beginning. See
 * the long discussion on this topic in the class documentation of
 * parallel::distributed::Triangulation .
 * </dd>
 *
 * @see @ref boundary "The module on boundaries"
 *
 * <dt class="glossary">@anchor GlossManifoldIndicator <b>%Manifold indicator</b></dt>
 *
 * <dd> Every object that makes up a Triangulation (cells, faces,
 * edges, etc.), is associated with a unique number (of type
 * types::manifol_id) that is used to identify which manifold object
 * is responsible to generate new points when the mesh is refined. 
 *
 * By default, all manifold indicators of a mesh are set to
 * numbers::invalid_manifold_id. A typical piece of code that sets the
 * manifold indicator on a object to something else would look like
 * this, here setting the manifold indicator to 42 for all cells whose
 * center has an $x$ component less than zero:
 *
 * @code
 * for (typename Triangulation<dim>::active_cell_iterator cell =
 *	triangulation.begin_active();
 *	cell != triangulation.end(); ++cell)
 *   if (cell->center()[0] < 0)
 *     cell->set_manifold_id (42);
 * @endcode
 *
 * Here we call the function TriaAccessor::set_manifold_id(). It may
 * also be appropriate to call TriaAccessor::set_all_manifold_ids
 * instead, to set recursively the manifold id on each face (and edge,
 * if in 3d). To query the manifold indicator of a particular object
 * edge, use TriaAccessor::manifold_id().
 *
 * The code above only sets the manifold indicators of a particular
 * part of the Triangulation, but it does not by itself change the way
 * the Triangulation class treats this object for the purposes of mesh
 * refinement. For this, you need to call Triangulation::set_manifold()
 * to associate a manifold object with a particular manifold
 * indicator. This allows the Triangulation objects to use a different
 * method of finding new points on cells, faces or edges to be
 * refined; the default is to use a FlatManifold object for all faces
 * and edges.
 *
 * @note Manifold indicators are inherited from parents to their
 * children upon mesh refinement. Some more information about manifold
 * indicators is also presented in a section of the documentation of
 * the Triangulation class as well as in the
 * @ref manifold "Manifold documentation module".
 * </dd>
 *
 * @see @ref manifold "The module on Manifolds"
 *
 * <dt class="glossary">@anchor GlossComponent <b>Component</b></dt>
 *
 * <dd> When considering systems of equations in which the solution is not
 * just a single scalar function, we say that we have a <i>vector system</i>
 * with a <i>vector-valued solution</i>. For example, the vector solution in
 * the elasticity equation considered in step-8 is $u=(u_x,u_y,u_z)^T$
 * consisting of the displacements in each of the three coordinate
 * directions. The solution then has three elements. Similarly, the 3d Stokes
 * equation considered in step-22 has four elements: $u=(v_x,v_y,v_z,p)^T$. We
 * call the elements of the vector-valued solution <i>components</i> in
 * deal.II. To be well-posed, for the solution to have $n$ components, there
 * need to be $n$ partial differential equations to describe them. This
 * concept is discussed in great detail in the @ref vector_valued module.
 *
 * In finite element programs, one frequently wants to address individual
 * elements (components) of this vector-valued solution, or sets of
 * components. For example, we do this extensively in step-8, and a lot
 * of documentation is also provided in the module on
 * @ref vector_valued "Handling vector valued problems". If you are thinking
 * only in terms of the partial differential equation (not in terms of
 * its discretization), then the concept of <i>components</i> is the natural
 * one.
 *
 * On the other hand, when talking about finite elements and degrees of
 * freedom, <i>components</i> are not always the correct concept because
 * components are not always individually addressable. In particular, this is
 * the case for @ref GlossPrimitive "non-primitive finite elements". Similarly,
 * one may not always <i>want</i> to address individual components but rather
 * sets of components &mdash; e.g. all velocity components together, and
 * separate from the pressure in the Stokes system, without further splitting
 * the velocities into their individual components. In either case, the
 * correct concept to think in is that of a @ref GlossBlock "block".  Since
 * each component, if individually addressable, is also a block, thinking in
 * terms of blocks is most frequently the better strategy.
 *
 * For a given finite element, the number of components can be queried using
 * the FiniteElementData::n_components() function, and you can find out
 * which vector components are nonzero for a given finite element shape
 * function using FiniteElement::get_nonzero_components(). The values and
 * gradients of individual components of a
 * shape function (if the element is primitive) can be queried using the
 * FiniteElement::shape_value_component() and
 * FiniteElement::shape_grad_component() functions on the reference cell. The
 * FEValues::shape_value_component() and FEValues::shape_grad_component()
 * functions do the same on a real cell. See also the documentation of the
 * FiniteElement and FEValues classes.
  *
 * <i>Selecting components:</i>
 * Many functions allow you to restrict their operation to certain
 * vector components or blocks. For example, this is the case for
 * the functions that interpolate boundary values: one may want
 * to only interpolate the boundary values for the velocity components of
 * a finite element field but not the pressure component. The way to do
 * this is by passing a ComponentMask argument to such functions, see the
 * @ref GlossComponentMask "component mask entry of this glossary".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossComponentMask <b>Component mask</b></dt>
 *
 * <dd>
 * When using vector-valued elements (see @ref vector_valued) to solve systems
 * of equations, one frequently wants to restrict some operations to only certain
 * solution variables. For example, when solving the Stokes equations, one may
 * wish to only interpolate boundary values for the velocity components
 * but not the pressure. In deal.II, this is typically done by passing functions
 * a <i>component mask</i>. Component masks are always specified as a
 * ComponentMask object which one can think of as an array with
 * as many entries as the finite element has components (e.g., in the Stokes case, there are
 * <code>dim+1</code> components) and where each entry is either true or false.
 * In the example where we would like to interpolate boundary values only for
 * the velocity components of the Stokes system, this component mask would then
 * be <code>[true, true, false]</code> in 2d and <code>[true, true, true, false]</code>
 * in 3d to indicate that no boundary values shall be set for the pressure variable
 * (the last of the <code>dim+1</code> vector components of the solution.
 *
 * There are many functions that take such component masks, for example
 * DoFTools::make_zero_boundary_values,
 * VectorTools::interpolate_boundary_values,
 * KellyErrorEstimator::estimate, etc. In some cases, there are multiple
 * functions with these names but only some of them have a component mask
 * argument.
 *
 * <b>Semantics of component masks:</b>
 * Many of the functions that take a component mask object that has been default
 * constructed to indicate <i>all components</i>, i.e., as if the vector had the
 * correct length and was filled with only <code>true</code> values. The reason
 * is that default initialized objects can be constructed in place using the code snippet
 * <code>ComponentMask()</code> and can thus be used as a default
 * argument in function signatures.
 *
 * In other words, ComponentMask objects can be in one of two states: They can have
 * been initialized by a vector of booleans with a nonzero length; in that case,
 * they represent a mask of a particular length where some elements may be true
 * and others may be false. Or, the ComponentMask may have been default initialized
 * (using the default constructor) in which case it represents an array of indefinite
 * length (i.e., a length appropriate to the circumstances) in which <i>every entry</i>
 * is true.
 *
 * <b>Creation of component masks:</b>
 * Component masks are typically created by asking the finite element
 * to generate a component mask from certain selected components using
 * code such as this where we create a mask that only denotes the
 * velocity components of a Stokes element (see @ref vector_valued):
 * @code
 *   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,    // Q2 element for the velocities
 *                            FE_Q<dim>(1), 1);     // Q1 element for the pressure
 *   FEValuesExtractors::Scalar pressure(dim);
 *   ComponentMask pressure_mask = stokes_fe.component_mask (pressure);
 * @endcode
 * The result is a component mask that, in 2d, would have values
 * <code>[false, false, true]</code>. Similarly, using
 * @code
 *   FEValuesExtractors::Vector velocities(0);
 *   ComponentMask velocity_mask = stokes_fe.component_mask (velocities);
 * @endcode
 * would result in a mask <code>[true, true, false]</code> in 2d. Of
 * course, in 3d, the result would be <code>[true, true, true, false]</code>.
 *
 * @note Just as one can think of composed elements as being made up of
 * @ref GlossComponent "components" or @ref GlossBlock "blocks", there are
 * component masks (represented by the ComponentMask class) and
 * @ref GlossBlockMask "block masks" (represented by the BlockMask class).
 * The FiniteElement class has functions that convert between the two kinds of
 * objects.
 *
 * @note Not all component masks actually make sense. For example, if you have
 * a FE_RaviartThomas object in 2d, then it doesn't make any sense to have a
 * component mask of the form <code>[true, false]</code> because you try to
 * select individual vector components of a finite element where each shape
 * function has both $x$ and $y$ velocities. In essence, while you can of
 * course create such a component mask, there is nothing you can do with it.
 * </dd>
 *
 *
 *
 * <dt class="glossary">@anchor GlossCompress <b>Compressing distributed
 *                                              vectors and matrices</b></dt>
 *
 * <dd>
 * For %parallel computations, deal.II uses the vector and matrix
 * classes defined in the PETScWrappers and TrilinosWrappers
 * namespaces. When running programs in %parallel using MPI, these
 * classes only store a certain number of rows or elements on the
 * current processor, whereas the rest of the vector or matrix is
 * stored on the other processors that belong to our MPI
 * universe. This presents a certain problem when you assemble linear
 * systems: we add elements to the matrix and right hand side vectors
 * that may or may not be stored locally. Sometimes, we may also want
 * to just <i>set</i> an element, not add to it.
 *
 * Both PETSc and Trilinos allow adding to or setting elements that
 * are not locally stored. In that case, they write the value that we
 * want to store or add into a cache, and we need to call one of the
 * functions TrilinosWrappers::VectorBase::compress(),
 * TrilinosWrappers::SparseMatrix::compress(),
 * PETScWrappers::VectorBase::compress() or
 * PETScWrappers::MatrixBase::compress() which will then ship the
 * values in the cache to the MPI process that owns the element to
 * which it is supposed to be added or written to. Due to the MPI
 * model that only allows to initiate communication from the sender
 * side (i.e. in particular, it is not a remote procedure call), these
 * functions are collective, i.e. they need to be called by all
 * processors.
 *
 * There is one snag, however: both PETSc and Trilinos need to know whether
 * the operation that these <code>compress()</code> functions invoke applies
 * to adding elements or setting them.  In some cases, not all processors may
 * be adding elements, for example if a processor does not own any cells when
 * using a very coarse (initial) mesh.  For this reason, compress() takes an
 * argument of type VectorOperation, which can be either ::%add, or ::%insert.
 * This argument is required for vectors and matrices starting with the 7.3
 * release.
 *
 * In olde releases we also proposed fake add/set operations. Those were the
 * cause of many bugs and deadlocks, so the usage of VectorOperation is now
 * required.
 *
 * In short, you need to call compress() in the following cases (and only in
 * those cases, though calling compress() in other cases just costs some
 * performance):
 *
 * 1. At the end of your assembly loop on matrices and vectors. This needs to
 * be done if you write entries directly or if you use
 * ConstraintMatrix::distribute_local_to_global. Use VectorOperation::add.
 *
 * 2. When you are done setting individual elements in a matrix/vector before
 * any other operations are done (adding to elements, other operations like
 * scaling, solving, reading, etc.). Use VectorOperation::insert.
 *
 * 3. Like in 2., but for adding values to individual elements. Use
 * VectorOperation::add.
 *
 * All other operations like scaling or adding vectors, assignments, calls
 * into deal.II (VectorTools, ConstraintMatrix, ...) or solvers do not require
 * calls to compress().
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDoF <b>Degree of freedom</b></dt>
 *
 * <dd> The term "degree of freedom" (often abbreviated as "DoF") is commonly
 * used in the finite element community to indicate two slightly different,
 * but related things. The first is that we'd like to represent the finite
 * element solution as a linear combination of shape function, in the form
 * $u_h(\mathbf x) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf x)$. Here, $U_j$
 * is a vector of expension coefficients. Because we don't know their values
 * yet (we will compute them as the solution of a linear or nonlinear system),
 * they are called "unknowns" or "degrees of freedom". The second meaning of
 * the term can be explained as follows: A mathematical description of finite
 * element problem is often to say that we are looking for a finite
 * dimensional function $u_h \in V_h$ that satisfies some set of equations
 * (e.g. $a(u_h,\varphi_h)=(f,\varphi_h)$ for all test functions $\varphi_h\in
 * V_h$). In other words, all we say here that the solution needs to lie in
 * some space $V_h$. However, to actually solve this problem on a computer we
 * need to choose a basis of this space; this is the set of shape functions
 * $\varphi_j(\mathbf x)$ we have used above in the expansion of $u_h(\mathbf
 * x)$ with coefficients $U_j$. There are of course many bases of the space
 * $V_h$, but we will specifically choose the one that is described by the
 * finite element functions that are traditionally defined locally on the
 * cells of the mesh. Describing "degrees of freedom" in this context requires
 * us to simply <i>enumerate</i> the basis functions of the space $V_h$. For
 * $Q_1$ elements this means simply enumerating the vertices of the mesh in
 * some way, but for higher elements one also has to enumerate the shape
 * functions that are associated with edges, faces, or cell interiors of the
 * mesh. The class that provides this enumeration of the basis functions of
 * $V_h$ is called DoFHandler.  The process of enumerating degrees of freedom
 * is referred to as "distributing DoFs" in deal.II.
 * </dd>
 *
 * <dt class="glossary">@anchor GlossDirectionFlag <b>Direction flags</b></dt>
 *
 * <dd>The <i>direction flag</i> is used in triangulations embedded in a
 * higher dimensional space to denote the orientation of cells and make the
 * manifold oriented. It is accessed using CellAccessor::direction_flag()
 * and set by the Triangulation class upon creation of a triangulation. You
 * can change all direction flags of a triangulation using the
 * Triangulation::flip_all_direction_flags() function.
 *
 * The flag is necessary to make cases like this work: assume we have a
 * one-dimensional mesh embedded in a two-dimensional space,
 *
 *   @image html direction_flag.png "One dimensional mesh in two dimensions"
 *
 * In one dimensional meshes in one dimensional space, we can always make sure
 * that the location of the left vertex of a cell has a smaller value than the
 * location of the right vertex. However, if we embed a mesh in a higher
 * dimensional space, we can no longer do this. For example, the cells in the
 * mesh above may be described by the following vertex sets: <code>(0,1),
 * (1,2), (3,2), (4,3), (4,5)</code>. (As a side remark, note that here we
 * have vertices -- e.g. vertex 2 -- that are the right end points of more
 * than one cell.) If we define the normal to each cell as that unit vector
 * that is right perpendicular to the vector that connects the first to the
 * second vertex of the line, then we would end up with the following picture:
 *
 *   @image html direction_flag_normals.png "Normal vectors"
 *
 * In other words, this one-dimensional manifold is not oriented. We could in
 * principle revert the order of vertices when creating such a mesh (though
 * there are good reasons not to do so, for example because this mesh may have
 * resulted from extracting the surface mesh of a two dimensional mesh, and we
 * want to preserve the order of vertices of each line segment because they
 * currently match the order of vertices of the faces of the 2d cells). An
 * alternative strategy, chosen in deal.II, is to simply associate with each
 * cell whether the normal should be the left or right normal to the
 * cell. (The default is right normals.) In the example above, the flags for
 * the five cells would be <code>true, true, false, false,
 * true</code>. Multiplying the right normal with plus or minus one, depending
 * on the value of the flag on each cell, yields a set of normal vectors that
 * orient the manifold.
 *
 * Similar issues happen with two-dimensional meshes in three space
 * dimensions. We note that it would not be possible to find consistent
 * direction flags if the two-dimensional manifold is not orientable; such
 * manifolds are not currently supported by deal.II.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDistorted <b>Distorted cells</b></dt>
 *
 * <dd>A <i>distorted cell</i> is a cell for which the mapping from
 * the reference cell to real cell has a Jacobian whose determinant is
 * non-positive somewhere in the cell. Typically, we only check the sign
 * of this determinant at the vertices of the cell. The function
 * GeometryInfo::alternating_form_at_vertices computes these
 * determinants at the vertices.
 *
 * By way of example, if all of the determinants are of roughly equal value
 * and on the order of $h^\text{dim}$ then the cell is well-shaped. For
 * example, a square cell or face has determinants equal to $h^\text{dim}$
 * whereas a strongly sheared parallelogram has a determinant much
 * smaller. Similarly, a cell with very unequal edge lengths will have widely
 * varying determinants. Conversely, a pinched cell in which the location of
 * two or more vertices is collapsed to a single point has a zero determinant
 * at this location. Finally, an inverted or twisted cell in which the
 * location of two vertices is out of order will have negative determinants.
 *
 * The following two images show a well-formed, a pinched, and a twisted cell
 * for both 2d and 3d:
 *
 * @image html distorted_2d.png "A well-formed, a pinched, and a twisted cell in 2d."
 *
 * @image html distorted_3d.png "A well-formed, a pinched, and a twisted cell in 3d."
 * </dd>
 *
 * Distorted cells can appear in two different ways: The original
 * coarse mesh can already contain such cells, or they can be created
 * as the result of mesh refinement if the boundary description in use
 * is sufficiently irregular.
 *
 * If the appropriate flag is given upon creation of a triangulation,
 * the function Triangulation::create_triangulation, which is called
 * by the various functions in GridGenerator and GridIn (but can also
 * be called from user code, see step-14, will signal
 * the creation of coarse meshes with distorted cells by throwing an
 * exception of type Triangulation::DistortedCellList. There are
 * legitimate cases for creating meshes with distorted cells (in
 * particular collapsed/pinched cells) if you don't intend to assemble
 * anything on these cells. For example, consider a case where one
 * would like to simulate the behavior of an elastic material with a
 * fluid-filled crack such as an oil reservoir. If the pressure
 * becomes too large, the crack is closed -- and the cells that
 * discretize the crack volume are collapsed to zero volume. As long
 * as you don't integrate over these cells to simulate the behavior of
 * the fluid (of which there isn't any if the crack has zero volume),
 * such meshes are perfectly legitimate. As a consequence,
 * Triangulation::create_triangulation does not simply abort the
 * program, but throws an exception that contains a list of cells that
 * are distorted; this exception can be caught and, if you believe
 * that you can ignore this condition, you can react by doing nothing
 * with the caught exception.
 *
 * The second case in which distorted cells can appear is through mesh
 * refinement when we have curved boundaries. Consider, for example, the
 * following case where the dashed line shows the exact boundary that the
 * lower edge of the cell is supposed to approximate (let's assume for
 * simplicity that the left, top and right edges are interior edges and
 * therefore will be considered as straight; in fact, for this particular case
 * in 2d where only one side of a cell is at the boundary we have special code
 * that avoids the situation depicted, but you will get the general idea of
 * the problem that holds in 3d or if more than one side of the cell is at the
 * boundary):
 *
 * @image html distorted_2d_refinement_01.png "One cell with an edge approximating a curved boundary"
 *
 * Now, if this cell is refined, we first split all edges and place
 * new mid-points on them. For the left, top and right edge, this is
 * trivial: because they are considered straight, we just take the
 * point in the middle between the two vertices. For the lower edge,
 * the Triangulation class asks the Boundary object associated with
 * this boundary (and in particular the Boundary::new_point_on_line
 * function) where the new point should lie. The four old vertices and
 * the four new points are shown here:
 *
 * @image html distorted_2d_refinement_02.png "Cell after edge refinement"
 *
 * The last step is to compute the location of the new point in the interior
 * of the cell. By default, it is chosen as the average location (arithmetic
 * mean of the coordinates) of these 8 points (in 3d, the 26 surrounding
 * points have different weights, but the idea is the same):
 *
 * @image html distorted_2d_refinement_03.png "Cell after edge refinement"
 *
 * The problem with that is, of course, that the bottom two child cells are
 * twisted, whereas the top two children are well-shaped. While such
 * meshes can happen with sufficiently irregular boundary descriptions
 * (and if the coarse mesh is entirely inadequate to resolve the
 * complexity of the boundary), the Triangulation class does not know
 * what to do in such situations unless one attaches an appropriate manifold
 * object to the cells in question (see the
 * @ref manifold "documentation module on manifolds"). Consequently, absent
 * such a manifold description or if the manifold description does not
 * provide a sufficient description of the geometry, the
 * Triangulation::execute_coarsening_and_refinement function does
 * create such meshes, but it keeps a list of cells whose children are
 * distorted. If this list is non-empty at the end of a refinement
 * step, it will throw an exception of type
 * Triangulation::DistortedCellList that contains those cells that
 * have distorted children. The caller of
 * Triangulation::execute_coarsening_and_refinement can then decide
 * what to do with this situation.
 *
 * One way to deal with this problem is to use the
 * GridTools::fix_up_distorted_child_cells function that attempts to
 * fix up exactly these cells if possible by moving around the node at
 * the center of the cell.
 *
 * Note that the Triangulation class does not test for the presence of
 * distorted cells by default, since the determination whether a cell
 * is distorted or not is not a cheap operation. If you want a
 * Triangulation object to test for distortion of cells, you need to
 * specify this upon creation of the object by passing the appropriate
 * flag.
 *
 *
 * <dt class="glossary">@anchor distributed_paper
 *                           <b>Distributed computing paper</b></dt>
 *
 * <dd>The "distributed computing paper" is a paper by W. Bangerth,
 * C. Burstedde, T. Heister and M. Kronbichler titled "Algorithms and Data
 * Structures for Massively Parallel Generic Finite Element Codes" that
 * describes the implementation of %parallel distributed computing in deal.II,
 * i.e. computations where not only the linear system is split onto different
 * machines as in, for example, step-17, but also the Triangulation and
 * DoFHandler objects. In essence, it is a guide to the parallel::distributed
 * namespace and the techniques used in step-40.
 *
 * The full reference for the paper is as follows:
 * @code
@Article{BBHK11,
  author =       {Wolfgang Bangerth and Carsten Burstedde and Timo Heister
                  and Martin Kronbichler},
  title =        {Algorithms and data structures for massively parallel generic
  adaptive finite element codes},
  journal =      {ACM Trans. Math. Softw.},
  year =         2011,
  volume =       38,
  pages =        {14/1--28}}
 * @endcode
 * It is also available as
 * <a href="http://iamcs.tamu.edu/file_dl.php?type=preprint&preprint_id=237">IAMCS
 * preprint 2011-187</a>.
 *
 * For massively %parallel
 * computations, deal.II builds on the
 * <a href="http://www.p4est.org/" target="_top">p4est</a>
 * library. If you use this functionality, please also cite the
 * p4est paper listed at their website.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossFaceOrientation <b>Face orientation</b></dt>
 * <dd>In a triangulation, the normal vector to a face
 * can be deduced from the face orientation by
 * applying the right hand side rule (x,y -> normal).  We note, that
 * in the standard orientation of faces in 2d, faces 0 and 2 have
 * normals that point into the cell, and faces 1 and 3 have normals
 * pointing outward. In 3d, faces 0, 2, and 4
 * have normals that point into the cell, while the normals of faces
 * 1, 3, and 5 point outward. This information, again, can be queried from
 * GeometryInfo<dim>::unit_normal_orientation.
 *
 * However, it turns out that a significant number of 3d meshes cannot
 * satisfy this convention. This is due to the fact that the face
 * convention for one cell already implies something for the
 * neighbor, since they share a common face and fixing it for the
 * first cell also fixes the normal vectors of the opposite faces of
 * both cells. It is easy to construct cases of loops of cells for
 * which this leads to cases where we cannot find orientations for
 * all faces that are consistent with this convention.
 *
 * For this reason, above convention is only what we call the
 * <em>standard orientation</em>. deal.II actually allows faces in 3d
 * to have either the standard direction, or its opposite, in which
 * case the lines that make up a cell would have reverted orders, and
 * the normal vector would have the opposite direction. You can ask a
 * cell whether a given face has standard orientation by calling
 * <tt>cell->face_orientation(face_no)</tt>: if the result is @p true,
 * then the face has standard orientation, otherwise its normal vector
 * is pointing the other direction. There are not very many places in
 * application programs where you need this information actually, but
 * a few places in the library make use of this. Note that in 2d, the
 * result is always @p true. However, while every face in 2d is always
 * in standard orientation, you can sometimes specify something to
 * assume that this is not so; an example is the function
 * DoFTools::make_periodicity_constraints().
 *
 * There are two other flags that describe the orientation of a face:
 * face_flip and face_rotation. Some documentation for these
 * exists in the GeometryInfo class. An example of their use in user
 * code is given in the DoFTools::make_periodicity_constraints function.
 *
 *
 * <dt class="glossary">@anchor GlossGeneralizedSupport <b>Generalized support points</b></dt>
 * <dd>While @ref GlossSupport "support points" allow very simple interpolation
 * into the finite element space, their concept is restricted to
 * @ref GlossLagrange "Lagrange elements". For other elements, more general
 * interpolation operators can be defined, often relying on integral values
 * or moments. Since these integral values are again computed using a
 * quadrature rule, we consider them a generalization of support
 * points.
 *
 * Note that there is no simple relation between
 * @ref GlossShape "shape functions" and generalized support points, unlike for
 * regular @ref GlossSupport "support points". Instead, FiniteElement defines
 * a couple of interpolation functions doing the actual interpolation.
 *
 * If a finite element is Lagrangian, generalized support points
 * and support points coincide.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossGhostCell <b>Ghost cells</b></dt>
 * <dd>
 * If a mesh is distributed across multiple MPI processes using the
 * parallel::distributed::Triangulation class, each processor stores
 * only the cells it owns, one layer of adjacent cells that are owned
 * by other processors, all coarse level cells, and all cells that are
 * necessary to maintain the invariant that adjacent cells must differ
 * by at most one refinement level. The cells stored on each process
 * that are not owned by this process but that are adjacent to the
 * ones owned by this process are called "ghost cells", and for these
 * cells the predicate <code>cell-@>is_ghost()</code> returns
 * true. Ghost cells are guaranteed to exist in the globally
 * distributed mesh, i.e. these cells are actually owned by another
 * process and are not further refined there. See the
 * @ref distributed_paper "Distributed Computing paper" for more
 * information.
 *
 * The layer of ghost cells consists of all cells that are face, edge, or
 * vertex neighbors of any locally owned cell and that are not locally
 * owned themselves. In other word, the ghost cells completely enclose the
 * subdomain of locally owned cells (with the exception of the boundary of
 * the domain, of course).
 *
 * The concept of ghost cells has no meaning for triangulations that
 * store the entire mesh on each processor, i.e. the
 * dealii::Triangulation class.  </dd>
 *
 *
 * <dt class="glossary">@anchor hp_paper <b>%hp paper</b></dt>
 * <dd>The "hp paper" is a paper by W. Bangerth and O. Kayser-Herold, titled
 * "Data Structures and Requirements for hp Finite Element Software", that
 * describes many of the algorithms and data structures used in the implementation
 * of the hp framework of deal.II. In particular, it summarizes many of the
 * tricky points that have to be considered for %hp finite elements using continuous
 * elements.
 *
 * The full reference for this paper is as follows:
 * @code
Article{BK07,
  author =       {Wolfgang Bangerth and Oliver Kayser-Herold},
  title =        {Data Structures and Requirements for hp Finite Element
                  Software},
  journal =      {ACM Trans. Math. Softw.},
  year =         2009,
  volume =       36,
  number =       1,
  pages =        {4/1--4/31}
}
 * @endcode
 * It is available as Technical Report ISC-07-04-MATH from the
 * <a href="http://www.isc.tamu.edu/publications-reports/technical_reports">Institute
 * for Scientific Computation, Texas A&amp;M University</a>, and also
 * from http://www.math.tamu.edu/~bangerth/publications.html .
 *
 * The numerical examples shown in that paper are generated with a slightly
 * modified version of step-27. The main difference to that
 * tutorial program is that various operations in the program were timed for
 * the paper to compare different options and show that $hp$ methods are
 * really not all that expensive.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossInterpolation <b>Interpolation with finite elements</b></dt>
 * <dd>The purpose of interpolation with finite elements is computing
 * a vector of coefficients representing a finite element function,
 * such that the @ref GlossNodes "node values" of the original
 * function and the finite element function coincide. Therefore, the
 * interpolation process consists of evaluating all @ref GlossNodes
 * "node functionals" <i>N<sub>i</sub></i> for the given function
 * <i>f</i> and store the result as entry <i>i</i> in the coefficient
 * vector.
 *
 *
 * <dt class="glossary">@anchor GlossLagrange <b>Lagrange elements</b></dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyOwnedCell <b>Locally owned cell</b></dt>
 * <dd>This concept identifies a subset of all cells when using
 * distributed meshes, see the @ref distributed module. In such meshes, each
 * cell is owned by exactly one processor. The locally owned ones are those
 * owned by the current processor.
 *
 * Each processor in a parallel computation has a triangulation covering
 * the entire domain that consists of cells that are locally owned, of
 * @ref GlossGhostCell "ghost cells" and of
 * @ref GlossArtificialCell "artificial cells".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyOwnedDof <b>Locally owned degrees of freedom</b></dt>
 * <dd>This concept identifies a subset of all degrees of freedom when using
 * distributed meshes, see the @ref distributed module.  Locally owned degrees
 * of freedom live on locally owned cells. Since degrees of freedom are owned
 * by only one processor, degrees of freedom on interfaces between cells owned
 * by different processors may be owned by one or the other, so not all
 * degrees of freedom on a locally owned cell are also locally owned degrees
 * of freedom.
 *
 * Locally owned DoFs are a subset of the
 * @ref GlossLocallyActiveDof "locally active DoFs".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyActiveDof <b>Locally active degrees of freedom</b></dt>
 * <dd>This concept identifies a subset of all degrees of freedom when using
 * distributed meshes, see the @ref distributed module.  Locally active degrees
 * of freedom are those that live on locally owned cells. Degrees of freedom
 * on interfaces between cells owned by different processors therefore belong
 * to the set of locally active degrees of freedom for more than one processor.
 *
 * Locally active DoFs are a superset of the
 * @ref GlossLocallyOwnedDof "locally owned DoFs" and a subset of the
 * @ref GlossLocallyRelevantDof "locally relevant DoFs".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyRelevantDof <b>Locally relevant degrees of freedom</b></dt>
 * <dd>This concept identifies a subset of all degrees of freedom when using
 * distributed meshes, see the @ref distributed module.  Locally relevant
 * degrees of freedom are those that live on locally owned or ghost cells.
 * Consequently, they may be owned by different processors.
 *
 * Locally relevant DoFs are a superset of the
 * @ref GlossLocallyActiveDof "locally active DoFs."
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossMaterialId <b>Material id</b></dt>
 * <dd>Each cell of a triangulation has associated with it a property called
 * "material id". It is commonly used in problems with heterogeneous
 * coefficients to identify which part of the domain a cell is in and,
 * consequently, which value the coefficient should have on this particular
 * cell. The material id is inherited from mother to child cell upon mesh
 * refinement.
 *
 * The material id is set and queried using the CellAccessor::material_id,
 * CellAccessor::set_material_id and CellAccessor::recursively_set_material_id
 * functions.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossMeshAsAContainer <b>Meshes as containers</b></dt>
 * <dd>
 * Meshes can be thought of as arrays of vertices and connectivities, but a
 * more fruitful view is to consider them as <i>collections of cells</i>. In C++,
 * collections are often called <i>containers</i> (typical containers are std::vector,
 * std::list, etc.) and they are characterized by the ability iterate over the
 * elements of the collection.
 *
 * Triangulations and objects of type DoFHandler or hp::DoFHandler can all be
 * considered as containers of cells. In fact, the most important parts of the
 * public interface of these classes consists simply of the ability to get
 * iterators to their elements, using functions such as Triangulation::begin_active(),
 * Triangulation::end() and their counterparts in DoFHandler and hp::DoFHandler. Since
 * these parts of the interface are generic, i.e., the functions have the same name
 * in all classes, it is possible to write operations that do not actually care whether
 * they work on a triangulation or a DoF handler object. Examples about, for example,
 * in the GridTools namespace, underlining the power of the abstraction that meshes
 * and DoF handlers can all be considered simply as collections (containers) of cells.
 *
 * On the other hand, meshes are non-standard containers unlike std::vector or std::list
 * in that they can be sliced several ways. For example, one can iterate over the
 * subset of active cells or over all cells; likewise, cells are organized into levels
 * and one can get iterator ranges for only the cells on one level. Generally, however,
 * all classes that implement the containers-of-cells concept use the same function
 * names to provide the same functionality.
 * </dd>
 *
 *
 *
 * <dt class="glossary">@anchor mg_paper <b>%Multigrid paper</b></dt>
 * <dd>The "multigrid paper" is a paper by B. Janssen and G. Kanschat, titled
 * "Adaptive multilevel methods with local smoothing", that
 * describes many of the algorithms and data structures used in the implementation
 * of the multigrid framework of deal.II. It underlies the implementation of
 * the classes that are used in step-16 for multigrid
 * methods.
 *
 * The full reference for this paper is as follows:
 * @code
Article{JK10,
  author =       {B. Janssen and G. Kanschat},
  title =        {Adaptive multilevel methods with local smoothing},
  journal =      {submitted},
  year =         2010
}
 * @endcode
 * It is available as Technical Report IAMCS-2009-131 from the
 * <a href="http://iamcs.tamu.edu/research_sub.php?tab_sub=research&cms_id=8">Institute
 * for Applied Mathematics and Computational Science, Texas A&amp;M University</a>.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossNodes <b>Node values or node functionals</b></dt>
 *
 * <dd>It is customary to define a FiniteElement as a pair consisting
 * of a local function space and a set of node values $N_i$ on the
 * mesh cells (usually defined on the @ref GlossReferenceCell
 * "reference cell"). Then, the basis of the local function space is
 * chosen such that $N_i(v_j) = \delta_{ij}$, the Kronecker delta.
 *
 * This splitting has several advantages, concerning analysis as well
 * as implementation. For the analysis, it means that conformity with
 * certain spaces (FiniteElementData::Conformity), e.g. continuity, is
 * up to the node values. In deal.II, it helps simplifying the
 * implementation of more complex elements like FE_RaviartThomas
 * considerably.
 *
 * Examples for node functionals are values in @ref GlossSupport
 * "support points" and moments with respect to Legendre
 * polynomials. Let us give some examples:
 *
 * <table><tr>
 *   <th>Element</th>
 *   <th>Function space</th>
 *   <th>Node values</th></tr>
 *   <tr><th>FE_Q, FE_DGQ</th>
 *     <td><i>Q<sub>k</sub></i></td>
 *     <td>values in support points</td></tr>
 *   <tr><th>FE_DGP</th>
 *     <td><i>P<sub>k</sub></i></td>
 *     <td>moments with respect to Legendre polynomials</td></tr>
 *   <tr><th>FE_RaviartThomas (2d)</th>
 *     <td><i>Q<sub>k+1,k</sub> x Q<sub>k,k+1</sub></i></td>
 *     <td>moments on edges and in the interior</td></tr>
 *   <tr><th>FE_RaviartThomasNodal</th>
 *     <td><i>Q<sub>k+1,k</sub> x Q<sub>k,k+1</sub></i></td>
 *     <td>Gauss points on edges(faces) and anisotropic Gauss points in the interior</td></tr>
 * </table>
 *
 * <dt class="glossary">@anchor GlossPrimitive <b>Primitive finite
 * elements</b></dt>
 * <dd>A finite element (described by its shape functions) is primitive if
 * there is a unique relation from shape function number to vector @ref
 * GlossComponent "component". What this means is that each shape function of
 * a vector-valued element has exactly one nonzero component if an element is
 * primitive. This includes, in particular, all scalar elements as well as
 * vector-valued elements assembled via the FESystem class from other
 * primitive (for example scalar) elements as shown in step-8,
 * step-29, step-22 and several others. On the other hand,
 * the FE_RaviartThomas class used in step-20 and step-21, or the FE_Nedelec
 * class provide non-primitive finite elements because there, each
 * vector-value shape function may have several non-zero components.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossReferenceCell <b>Reference cell</b></dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined. Many properties of the reference
 * cell are described by the GeometryInfo class.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossSerialization <b>Serialization</b></dt>

 * <dd>The term "serialization" refers to the process of writing the state of
 * an object to a stream and later retrieve it again. A typical use case is to
 * save the state of a program to disk for possible later resurrection, often
 * in the context of checkpoint/restart strategies for long running
 * computations or on computers that aren't very reliable (e.g. on very large
 * clusters where individual nodes occasionally fail and then bring down an
 * entire MPI job). In either case, one wants to occasionally save the state
 * of the program so that, upon failure, one can restart it at that point
 * rather than having to run it again from the beginning.
 *
 * deal.II implements serialization facilities by implementing the necessary
 * interfaces for the <a
 * href="http://www.boost.org/doc/libs/1_46_1/libs/serialization/doc/index.html"
 * target="_top">BOOST serialization</a> library. See there for examples on
 * how to save and restore objects. </dd>
 *
 *
 * <dt class="glossary">@anchor GlossShape <b>Shape functions</b></dt>
 * <dd>The restriction of the finite element basis functions to a single
 * grid cell.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossSubdomainId <b>Subdomain id</b></dt>
 * <dd>Each cell of a triangulation has associated with it a property called
 * the "subdomain id" that can be queried using a call like
 * <code>cell-@>subdomain_id()</code> and that can be set for example by using
 * <code>cell-@>set_subdomain_id(13)</code>. (These calls resolve to
 * CellAccessor::subdomain_id() and CellAccessor::set_subdomain_id(), respectively.)
 * While in principle this property
 * can be used in any way application programs deem useful (it is simply an
 * integer associated with each cell that can indicate whatever you want), at
 * least for programs that run in %parallel it usually denotes the processor a
 * cell is associated with.
 *
 * For programs that are parallelized based on MPI but where each processor
 * stores the entire triangulation (as in, for example, step-18, but not in
 * step-32), subdomain ids are assigned to cells by
 * partitioning a mesh, and each MPI process then only works on those cells it
 * "owns", i.e. that belong to a subdomain that the processor is associated with
 * (traditionally, this is the case for the subdomain id whose numerical value
 * coincides with the rank of the MPI process within the MPI
 * communicator). Partitioning is typically done using the
 * GridTools::partition() function, but any other method can also be used to
 * do this.
 *
 * On the other hand, for programs that are parallelized using MPI but
 * where meshes are held distributed across several processors using
 * the parallel::distributed::Triangulation and
 * parallel::distributed::DoFHandler classes, the subdomain id of
 * cells are tied to the processor that owns the cell. In other words,
 * querying the subdomain id of a cell tells you if the cell is owned
 * by the current processor (i.e. if <code>cell-@>subdomain_id() ==
 * triangulation.parallel::distributed::Triangulation::locally_owned_subdomain()</code>)
 * or by another processor. In the %parallel distributed case,
 * subdomain ids are only assigned to cells that the current processor
 * owns as well as the immediately adjacent @ref GlossGhostCell "ghost cells".
 * Cells further away are held on each processor to ensure
 * that every MPI process has access to the full coarse grid as well
 * as to ensure the invariant that neighboring cells differ by at most
 * one refinement level. These cells are called "artificial" (see
 * @ref GlossArtificialCell "here") and have the special subdomain id value
 * types::artificial_subdomain_id.
 *
 * In addition to regular subdomain ids, there is a second, closely related set
 * of flags that are associated with each cell: "level subdomain ids."
 * These exist not only for active cells but, in fact, for every cell in
 * a mesh hierarchy. Their meaning is entirely analogous to the regular
 * subdomain ids, but they are read and written by the
 * CellAccessor::level_subdomain_id() and CellAccessor::set_level_subdomain_id()
 * functions.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossSupport <b>Support points</b></dt> <dd>Support points are
 * by definition those points $p_i$, such that for the shape functions
 * $v_j$ holds $v_j(p_i) = \delta_{ij}$. Therefore, a finite element
 * interpolation can be defined uniquely by the values in the support
 * points.
 *
 * Lagrangian elements fill the vector accessed by
 * FiniteElementBase::get_unit_support_points(), such that the
 * function FiniteElementBase::has_support_points() returns
 * <tt>true</tt>. Naturally, these support points are on the
 * @ref GlossReferenceCell "reference cell".  Then, FEValues can be used
 * (in conjunction with a Mapping) to access support points on the
 * actual grid cells.
 *
 * @note The concept of @ref GlossSupport "support points" is
 * restricted to the finite element families based on Lagrange
 * interpolation. For a more general concept, see
 * @ref GlossGeneralizedSupport "generalized support points".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossTargetComponent <b>Target component</b></dt> <dd>When
 * vectors and matrices are grouped into blocks by component, it is
 * often desirable to collect several of the original components into
 * a single one. This could be for instance, grouping the velocities
 * of a Stokes system into a single block.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitCell <b>Unit cell</b></dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitSupport <b>Unit support points</b></dt>
 * <dd>These are the @ref GlossSupport "support points" on the reference cell, defined in
 * FiniteElementBase. For example, the usual Q1 element in 1d has support
 * points  at <tt>x=0</tt> and <tt>x=1</tt> (and similarly, in higher
 * dimensions at the vertices of the unit square or cube). On the other
 * hand, higher order Lagrangian elements have unit support points also
 * in the interior of the unit line, square, or cube.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossUserFlags <b>User flags</b></dt>
 * <dd>
 *   A triangulation offers one bit per line, quad, etc for user flags.
 *   This field can be
 *   accessed as all other data using iterators, using the syntax
 *   @code
 *      cell->set_user_flag();                // set the user flag of a cell
 *      if (cell->user_flag_set() == false)   // if cell hasn't been flagged yet
 *        {
 *           cell->face(0)->set_user_flag();  // flag its first face
 *        }
 *   @endcode
 *   Typically, this user flag is
 *   used if an algorithm walks over all cells and needs information whether
 *   another cell, e.g. a neighbor, has already been processed. Similarly,
 *   it can be used to flag faces, quads or lines at the boundary for which
 *   some operation has already been performed. The latter is often useful
 *   since a loop such as
 *   @code
 *      // in 3d
 *      for (cell=dof_handler.begin_active();
 *           cell!=dof_handler.end(); ++cell)
 *        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
 *          if (cell->line(l)->at_boundary())
 *            {
 *               do something with this line
 *            }
 *   @endcode
 *   encounters some boundary lines more than once. Consequently, one would
 *   set the user flag of the line in the body of the loop, and only enter the
 *   body if the user flag had not previously been set. There are a number of
 *   additional functions that can be accessed through the iterator interface;
 *   see the TriaAccessor class for more information. Note that there are no
 *   user flags that can be associated with vertices; however, since vertices
 *   are numbered consecutively, this can easily be emulated in user code
 *   using a vector of bools.
 *
 *   There are two functions, Triangulation::save_user_flags and
 *   Triangulation::load_user_flags which
 *   write and read these flags to and from a stream or a vector of bools. Unlike
 *   Triangulation::save_refine_flags and Triangulation::load_refine_flags,
 *   these two functions store
 *   and read the flags of all used lines, quads, etc, i.e., not only of the
 *   active ones.
 *
 *   If you want to store more specific user flags, you can use the functions
 *   Triangulation::save_user_flags_line and Triangulation::load_user_flags_line
 *   and the similarly for quads, etc.
 *
 *   As for the refinement and coarsening flags, there exist two versions of these
 *   functions, one which reads/writes from a stream and one which does so from
 *   a <tt>vector@<bool@></tt>. The latter is used to store flags temporarily, while the
 *   first is used to store them in a file.
 *
 *   It is good practice to clear the user flags using the
 *   Triangulation::clear_user_flags() function before usage, since it is
 *   often necessary to use the flags in more than one function. If the flags may
 *   be in use at the time a function that needs them is called, then this function
 *   should save and restore the flags as described above.
 *
 *   @note If more information than just a single boolean flag needs to be stored
 *   with a cell, line, or face, then see about @ref GlossUserData "user data".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossUserData <b>User pointers and user indices</b></dt>
 * <dd>
 *   Just like the @ref GlossUserFlags "user flags", the Triangulation class offers a
 *   field for each line, quad and hex in which to store more descriptive data than just
 *   a single boolean flag. This is called "user data" and the data that can be stored
 *   in it is either a single unsigned integer or a void pointer. Both are typically
 *   used to index into a bigger array that contains more detailed data an application
 *   wants to attach to a mesh entity.
 *
 *   User data is stored and retrieved in the following manner:
 *   @code
 *      for (cell=dof_handler.begin_active();
 *           cell!=dof_handler.end(); ++cell)
 *        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
 *          if (cell->line(l)->at_boundary())
 *            {
 *              cell->line(l)->set_user_index(42);
 *            }
 *   @endcode
 *   Similarly, there are functions TriaAccessor::set_user_pointer to set a pointer, and
 *   TriaAccessor::user_index and TriaAccessor::user_pointer to retrieve the index
 *   and pointer. To clear all user indices or pointers, use Triangulation::clear_user_data().
 *   As with flags, there are functions that allow to save and restore user data,
 *   either for all entities of the mesh hierarchy or for lines, quads or hexes
 *   separately. There are a number of additional functions that can be accessed
 *   through the iterator interface; see the TriaAccessor class for more information.
 *
 *   @note User pointers and user indices are stored in the same
 *   place. In order to avoid unwanted conversions, Triangulation
 *   checks which one of them is in use and does not allow access to
 *   the other one, until Triangulation::clear_user_data() has been called.
 *
 *   @note The usual warning about the missing type safety of @p void pointers are
 *   obviously in place here; responsibility for correctness of types etc
 *   lies entirely with the user of the pointer.
 *
 *
 * <dt class="glossary">@anchor workstream_paper <b>%WorkStream paper</b></dt>
 * <dd>The "%WorkStream paper" is a paper by B. Turcksin, M. Kronbichler and W. Bangerth
 *   that discusses the design and implementation of WorkStream. WorkStream is, at its
 *   core, a design pattern, i.e., something that is used over and over in finite element
 *   codes and that can, consequently, be implemented generically. In particular, the
 *   paper lays out the motivation for this pattern and then proposes different ways
 *   of implementing it. It also compares the performance of different implementations.
 *
 *   The paper is currently in preparation.
 * </dd>
 *
 * </dl>
 */
