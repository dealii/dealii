// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
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
 * <dd>A cell, face, or edge is defined as <i>active</i> if it is not
 * refined any further, i.e., if it does not have children. Once a cell,
 * face, or edge becomes a parent it is no longer active. Unless
 * working with a multigrid algorithm, active cells are the only
 * ones carrying degrees of freedom.
 * </dd>
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
 * dealii::Triangulation class.
 * </dd>
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
 * respectively, $M$ is the @ref GlossMassMatrix "mass matrix" on the velocity space, $B$ corresponds to
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
 * overview of the various linear algebra classes in the
 * @ref LAC topic. The objects present two interfaces: one that makes the
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
 * (see the @ref feall "topic on finite element classes").
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
 * FESystem, the @ref vector_valued topic and the tutorial programs
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
 * <dt class="glossary">@anchor GlossBoundaryForm <b>Boundary form</b></dt>
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
 * cell.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBoundaryIndicator <b>Boundary indicator</b></dt>
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
 * that have a
 * @ref GlossColorization "colorize"
 * option. A typical piece of code that sets the boundary
 * indicator on part of the boundary to something else would look like
 * this, here setting the boundary indicator to 42 for all faces located at
 * $x=-1$:
 * @code
 *   for (auto &face : triangulation.active_face_iterators())
 *     if (face->at_boundary())
 *       if (face->center()[0] == -1)
 *         face->set_boundary_id (42);
 * @endcode
 * This calls functions TriaAccessor::set_boundary_id. In 3d, it may
 * also be appropriate to call TriaAccessor::set_all_boundary_ids instead
 * on each of the selected faces. To query the boundary indicator of a particular
 * face or edge, use TriaAccessor::boundary_id.
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
 * It is often useful to visualize the boundary indicators used on different
 * parts of the boundary of a domain. For example, when reading in a complex
 * mesh generated by an external program, one typically wants to verify
 * that what is *believed* to be the assignment of boundary indicators is
 * indeed what is in the mesh file, or maybe equally importantly, what deal.II
 * has read in from the mesh file. This can be done by using the
 * DataPostprocessors::BoundaryIds class and the following small piece of code:
 * @code
 *    #include <deal.II/numerics/data_out_faces.h>
 *    #include <deal.II/numerics/data_postprocessor.h>
 *
 *    ...
 *
 *    DataPostprocessors::BoundaryIds<dim> boundary_ids;
 *    DataOutFaces<dim> data_out_faces;
 *    FE_Q<dim>         dummy_fe(1);
 *
 *    DoFHandler<dim>   dummy_dof_handler(triangulation);
 *    dummy_dof_handler.distribute_dofs(dummy_fe);
 *
 *    Vector<double> dummy_solution (dummy_dof_handler.n_dofs());
 *
 *    data_out_faces.attach_dof_handler(dummy_dof_handler);
 *    data_out_faces.add_data_vector(dummy_solution, boundary_ids);
 *    data_out_faces.build_patches();
 *
 *    std::ofstream out("boundary_ids.vtu");
 *    data_out_faces.write_vtu(out);
 * @endcode
 * The code requires setting up a dummy finite element, DoFHandler, and
 * solution vector because that's what DataOutFaces works on. However,
 * while the content of the dummy solution, evaluated on each face of the
 * domain, is passed to the DataPostprocessors::BoundaryIds member functions,
 * they simply ignore this information and only output information obtained
 * from the mesh (namely, the boundary indicators).
 *
 * The example code above uses FE_Q as the finite element. This is appropriate
 * when the mesh consists of quadrilateral or hexahedral meshes. When using
 * simplex meshes (triangular, tetrahedral), FE_SimplexP is the right choice.
 * Mixed meshes will require a bit more work, but the general idea should be
 * clear. In general, which element you choose (include FE_Nothing) does not
 * matter since the actual data at evaluation points is ignored by the
 * postprocessor; the only thing that matters is that the element you choose
 * is compatible with the shape of the cells used in the mesh.
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
 *
 * <dt class="glossary">@anchor GlossCoarseMesh <b>Coarse mesh</b></dt>
 * <dd>
 *   A "coarse mesh" in deal.II is a triangulation object that consists only
 *   of cells that are not refined, i.e., a mesh in which no cell is a child
 *   of another cell. This is generally how triangulations are first
 *   constructed in deal.II, for example using (most of) the functions in
 *   namespace GridGenerator, the functions in class GridIn, or directly
 *   using the function Triangulation::create_triangulation(). One can of
 *   course do computations on such meshes, but most of the time (see, for
 *   example, almost any of the tutorial programs) one first refines the
 *   coarse mesh globally (using Triangulation::refine_global()),
 *   or adaptively (in that case first computing a refinement
 *   criterion, then one of the functions in namespace GridRefinement,
 *   and finally calling
 *   Triangulation::execute_coarsening_and_refinement()). The mesh is
 *   then no longer a "coarse mesh", but a "refined mesh".
 *
 *   In some contexts, we also use the phrase "the coarse mesh of a
 *   triangulation", and by that mean that set of cells that the triangulation
 *   started out with, i.e., from which all the currently
 *   @ref GlossActive "active cells" of the triangulation have been obtained
 *   by mesh refinement. (Some of the coarse mesh cells may of course also
 *   be active if they have never been refined.)
 *
 *   Triangulation objects store cells in <i>levels</i>: in
 *   particular, all cells of a coarse mesh are on level zero. Their
 *   children (if we executed `Triangulation::refine_global(1)` on a
 *   coarse mesh) would then be at level one, etc. The coarse mesh of a
 *   triangulation (in the sense of the previous paragraph) then
 *   consists of exactly the level-zero cells of a triangulation. (Whether
 *   they are active (i.e., have no children) or have been refined is not
 *   important for this definition.)
 *
 *   Most of the triangulation classes in deal.II store the entire coarse
 *   mesh along with at least some of the refined cells. (Both the
 *   dealii::Triangulation and parallel::shared::Triangulation classes
 *   actually store <i>all</i> cells of the entire mesh, whereas some
 *   other classes such as parallel::distributed::Triangulation only
 *   store <i>some</i> of the @ref GlossActive "active cells" on
 *   each process in a parallel computation.) In those cases,
 *   one can query the triangulation for all coarse mesh
 *   cells. Other triangulation classes (e.g.,
 *   parallel::fullydistributed::Triangulation) only store a part
 *   of the coarse mesh. See also
 *   @ref GlossCoarseCellId "the concept of coarse cell ids"
 *   for that case.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossCoarseCellId <b>Coarse cell ID</b></dt>
 * <dd>
 *   Most of the triangulation classes in deal.II, notably
 *   dealii::Triangulation, parallel::shared::Triangulation, and
 *   parallel::distributed::Triangulation, store the entire
 *   @ref GlossCoarseMesh "coarse mesh"
 *   of a triangulation on each process of a parallel computation. On the
 *   other hand, this is not the case for other classes, notably for
 *   parallel::fullydistributed::Triangulation, which is designed for cases
 *   where even the coarse mesh is too large to be stored on each process
 *   and needs to be partitioned.
 *
 *   In those cases, it is often necessary in algorithms to reference a coarse
 *   mesh cell uniquely. Because the triangulation object on the current
 *   process does not actually store the entire coarse mesh, one needs to have
 *   a globally unique identifier for each coarse mesh cell that is independent
 *   of the index within level zero of the triangulation stored locally. This
 *   globally unique ID is called the "coarse cell ID". It can be accessed via
 *   the function call
 *   @code
 *     triangulation.coarse_cell_index_to_coarse_cell_id (coarse_cell->index());
 *   @endcode
 *   where `triangulation` is the triangulation to which the iterator
 *   `coarse_cell` pointing to a cell at level zero belongs. Here,
 *   `coarse_cell->index()` returns the index of that cell within its
 *   refinement level (see TriaAccessor::index()). This is a number
 *   between zero and the number of coarse mesh cells stored on the
 *   current process in a parallel computation; it uniquely identifies
 *   a cell on that parallel process, but different parallel processes may
 *   use that index for different cells located at different coordinates.
 *
 *   For those classes that store all coarse mesh cells on each process,
 *   the Triangulation::coarse_cell_index_to_coarse_cell_id() simply
 *   returns a permutation of the possible argument values. In the
 *   simplest cases, such as for a sequential or a parallel shared
 *   triangulation, the function will in fact simply return the
 *   value of the argument. For others, such as
 *   parallel::distributed::Triangulation, the ordering of
 *   coarse cell IDs is not the same as the ordering of coarse
 *   cell indices. Finally, for classes such as
 *   parallel::fullydistributed::Triangulation, the function returns
 *   the globally unique ID, which is from a larger set of possible
 *   indices than the indices of the coarse cells actually stored on
 *   the current process.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossCollectiveOperation <b>Collective operation</b></dt>
 * <dd>
 *   When running programs in parallel using MPI, a <em>collective
 *   operation</em> is one in which all processes on an
 *   @ref GlossMPICommunicator "MPI communicator"
 *   have to participate. At its core, the concept of collective operations
 *   rests on the mental model that MPI traditionally uses, namely where
 *   processes communicate by sending messages to each other; in this model,
 *   nothing happens if one process sends a message but the receiving process
 *   does not expect or respond to it, or if one process needs access to a
 *   piece of data stored elsewhere, but the storing process does not send it.
 *
 *   Collective operations are then operations that need to be called on all
 *   processes at the same time to execute. An obvious example is calling the
 *   `MPI_Sum` function in which every process provides a number that is then
 *   summed over all processes. If in a program running with 4 processes only
 *   three processes call `MPI_Sum`, the program will hang until the fourth
 *   process eventually also gets to the place where this function is called.
 *   If the fourth process never calls that function, for example because it
 *   calls another MPI function in the belief that the other processes called
 *   that function as well, a "deadlock" results: Every process is now waiting
 *   for something to happen that cannot happen.
 *
 *   Many functions in deal.II are "collective operations" because internally
 *   they call MPI functions that are collective. For some, this is obvious,
 *   such as when you call
 *   parallel::distributed::Triangulation::execute_coarsening_and_refinement()
 *   in step-40, given that this function refines a mesh that is stored
 *   in parallel on all processes in a parallel universe. In some other
 *   cases, the name of the function is a hint that a function is collective,
 *   such as in
 *   GridTools::distributed_compute_point_locations() or
 *   GridTools::build_global_description_tree(), where the "distributed"
 *   and "global" components of the names are an indication; the latter
 *   function also takes an explicit MPI communicator as argument. For yet
 *   other functions, it is perhaps not as obvious that a function is
 *   a collective operation. GridTools::volume() is an example; it is
 *   collective because each process computes the volume of those cells it
 *   locally owns, and these contributions then have to be added up.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossColorization <b>Colorization</b></dt>
 * <dd><em>Colorization</em> is the process of marking certain parts of a
 * Triangulation with different labels. The use of the word <em>color</em>
 * comes from cartography, where countries on a map are made visually distinct
 * from each other by assigning them different colors. Using the same term
 * <em>coloring</em> is common in mathematics, even though we assign integers
 * and not hues to different regions. deal.II refers to two processes as
 * coloring:
 *
 * <ol>
 *   <li> Most of the functions in the GridGenerator namespace take an optional
 *   argument <code>colorize</code>. This argument controls whether or not the
 *   different parts of the boundary will be assigned different
 *   @ref GlossBoundaryIndicator "boundary indicators".
 *   Some functions also assign different
 *   @ref GlossMaterialId "material indicators"
 *   as well.</li>
 *   <li> The function GraphColoring::make_graph_coloring() computes a
 *   decomposition of a Triangulation (more exactly, a range of iterators). No
 *   two adjacent cells are given the same color.</li>
 * </ol>
 * </dd>
 *
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
 * concept is discussed in great detail in the @ref vector_valued topic.
 *
 * In finite element programs, one frequently wants to address individual
 * elements (components) of this vector-valued solution, or sets of
 * components. For example, we do this extensively in step-8, and a lot
 * of documentation is also provided in the topic on
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
 * using a very @ref GlossCoarseMesh "coarse (initial) mesh".
 * For this reason, compress() takes an
 * argument of type VectorOperation, which can be either ::%add, or ::%insert.
 * This argument is required for vectors and matrices starting with the 7.3
 * release.
 *
 * In short, you need to call compress() in the following cases (and only in
 * those cases, though calling compress() in other cases just costs some
 * performance):
 *
 * 1. At the end of your assembly loop on matrices and vectors. This needs to
 * be done if you write entries directly or if you use
 * AffineConstraints::distribute_local_to_global. Use VectorOperation::add.
 *
 * 2. When you are done setting individual elements in a matrix/vector before
 * any other operations are done (adding to elements, other operations like
 * scaling, solving, reading, etc.). Use VectorOperation::insert.
 *
 * 3. Like in 2., but for adding values to individual elements. Use
 * VectorOperation::add.
 *
 * All other operations like scaling or adding vectors, assignments, calls
 * into deal.II (VectorTools, AffineConstraints, ...) or solvers do not
 * require calls to compress().
 * </dd>
 *
 * @note Compressing is an operation that only applies to vectors whose
 * elements are uniquely owned by one and only one processor in a parallel
 * MPI universe. It does not apply to
 * @ref GlossGhostedVector "vectors with ghost elements".
 *
 *
 * <dt class="glossary">@anchor GlossConcept <b>Concepts in deal.II</b></dt>
 *
 * <dd> There are several places in deal.II where we require that a type in a
 * template match a certain interface or behave in a certain way: such
 * constraints are called <em>concepts</em> in C++. See the discussion in
 * @ref Concepts for more information and a list of concepts in deal.II.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDevice <b>Device</b></dt>
 *
 * <dd> We commonly refer to GPUs as "devices" in deal.II. The context is
 * always related to Kokkos that motivated using this term.
 * Occasionally, we also call data corresponding to MemorySpace::Default "device data"
 * (even though it is allocated in CPU memory if Kokkos was configured without
 * a GPU backend) to distinguish between MemorySpace::Default and MemorySpace::Host.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDimension <b>Dimensions `dim` and `spacedim`</b></dt>
 *
 * <dd> Many classes and functions in deal.II have two template parameters,
 * @p dim and @p spacedim. An example is the basic Triangulation class:
 * @code
 *   template <int dim, int spacedim=dim>
 *   class Triangulation {...};
 * @endcode
 * In all of these contexts where you see `dim` and `spacedim` referenced,
 * these arguments have the following meaning:
 *
 * <ul>
 *   <li> @p dim denotes the dimensionality of the mesh. For example, a mesh
 *   that consists of line segments is one-dimensional and consequently
 *   corresponds to `dim==1`. A mesh consisting of quadrilaterals then has
 *   `dim==2` and a mesh of hexahedra has `dim==3`.</li>
 *
 *   <li> @p spacedim denotes the dimensionality of the space in which such a
 *   mesh lives. Generally, one-dimensional meshes live in a one-dimensional
 *   space, and similarly for two-dimensional and three-dimensional meshes
 *   that subdivide two- and three-dimensional domains. Consequently, the @p
 *   spacedim template argument has a default equal to @p dim. But this need
 *   not be the case: For example, we may want to solve an equation for
 *   sediment transport on the surface of the Earth. In this case, the domain
 *   is the two-dimensional surface of the Earth (`dim==2`) that lives in a
 *   three-dimensional coordinate system (`spacedim==3`).</li>
 * </ul>
 *
 * More generally, deal.II can be used to solve partial differential
 * equations on <a href="https://en.wikipedia.org/wiki/Manifold">manifolds</a>
 * that are embedded in higher dimensional space. In other words,
 * these two template arguments need to satisfy `dim <= spacedim`,
 * though in many applications one simply has `dim == spacedim`.
 *
 * Following the convention in geometry, we say that the "codimension" is
 * defined as `spacedim-dim`. In other words, a triangulation consisting of
 * quadrilaterals whose coordinates are three-dimensional (for which we
 * would then use a `Triangulation<2,3>` object) has "codimension one".
 *
 * Examples of uses where these two arguments are not the same are shown in
 * step-34, step-38, step-54.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDoF <b>Degree of freedom</b></dt>
 *
 * <dd> The term "degree of freedom" (often abbreviated as "DoF") is commonly
 * used in the finite element community to indicate two slightly different,
 * but related things. The first is that we'd like to represent the finite
 * element solution as a linear combination of shape functions, in the form
 * $u_h(\mathbf{x}) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf{x})$. Here, $U_j$
 * is a vector of expansion coefficients. Because we don't know their values
 * yet (we will compute them as the solution of a linear or nonlinear system),
 * they are called "unknowns" or "degrees of freedom". The second meaning of
 * the term can be explained as follows: A mathematical description of finite
 * element problem is often to say that we are looking for a finite
 * dimensional function $u_h \in V_h$ that satisfies some set of equations
 * (e.g. $a(u_h,\varphi_h)=(f,\varphi_h)$ for all test functions $\varphi_h\in
 * V_h$). In other words, all we say here that the solution needs to lie in
 * some space $V_h$. However, to actually solve this problem on a computer we
 * need to choose a basis of this space; this is the set of shape functions
 * $\varphi_j(\mathbf{x})$ we have used above in the expansion of $u_h(\mathbf
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
 * The flag is necessary to make cases like this work: Assume we have a
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
 *
 * Finally, the direction flag cannot be used for triangulations where
 * `spacedim>dim+1`, such as for meshes with one-dimensional cells in 3d.
 * In these cases, the normal vector to a cell does not simply point to
 * one side or the other of a cell, but must lie in a two-dimensional sub-space
 * perpendicular to the cell. As a consequence, we cannot make normal vectors
 * consistent between adjacent cells simply by flipping it from one side to
 * the other.
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
 *
 * Distorted cells can appear in two different ways: The original
 * @ref GlossCoarseMesh "coarse mesh" can already contain such cells,
 * or they can be created as the result of moving or distorting a mesh by a
 * relatively large amount.
 *
 * If the appropriate flag is given upon creation of a triangulation,
 * the function Triangulation::create_triangulation, which is called
 * by the various functions in GridGenerator and GridIn (but can also
 * be called from user code, see step-14 and the example at the end of step-49),
 * will signal the creation of coarse meshes with distorted cells by throwing an
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
 * The function GridTools::fix_up_distorted_child_cells can, in some cases,
 * fix distorted cells on refined meshes by moving around the vertices of a
 * distorted child cell that has an undistorted parent.
 *
 * Note that the Triangulation class does not test for the presence of
 * distorted cells by default, since the determination whether a cell
 * is distorted or not is not a cheap operation. If you want a
 * Triangulation object to test for distortion of cells, you need to
 * specify this upon creation of the object by passing the appropriate
 * flag.
 * </dd>
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
 * @code{.bib}
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
 *
 * For massively %parallel
 * computations, deal.II builds on the
 * <a href="https://www.p4est.org/" target="_top">p4est</a>
 * library. If you use this functionality, please also cite the
 * p4est paper listed at their website.
 * </dd>
 *
 * <dt class="glossary">@anchor GlossCombinedOrientation <b>Combined
 * orientation</b></dt>
 * <dd>
 * A Triangulation contains cells as well as lower dimensional objects such as
 * faces (which are either lines in 2d or quadrilaterals or triangles in 3d). In
 * general, the vertices of each cell are numbered in a way which results in a
 * mapping with a positive Jacobian. A consequence of this choice is that the
 * vertices which define a face may be, from the perspective of an arbitrary
 * cell, in a different order than the order given to that face by the
 * neighboring cell. To resolve this inconsistency deal.II stores, for each face
 * and (in 3d) line, a value which may be used to permute the vertices on both
 * faces into a matching configuration. This encoding contains both the number
 * of times a face should be rotated (relative to its neighbor) as well as
 * whether or not the face should be viewed in an opposite orientation (e.g.,
 * for a Quadrilateral, whether or not vertices $1$ and $2$ should be swapped).
 *
 * This value is called the <em>combined_orientation</em> since it combines both
 * the orientation (as defined above) as well as rotations. In some
 * circumstances, to disambiguate between faces and lines, it may alternatively
 * be called either the `combined_face_orientation` or the
 * `combined_line_orientation`. These orientations are represented by
 * types::geometric_orientation, which encodes how the vertices of the canonical
 * definition of a face should be permuted so that they equal the current cell's
 * definition of that face. The binary encoding (which is usually represented as
 * a decomposition of three booleans called orientation, rotation, and flip) of
 * that permutation is an internal library detail and is documented in
 * @ref reordering "the cell reordering page".
 * The default value (which corresponds to the identity permutation) is
 * numbers::default_geometric_orientation. As lines only have two possible
 * orientations (i.e., the vertices are either in the same order as the
 * canonical line or are swapped), the other orientation is encoded as
 * numbers::reverse_line_orientation.
 *
 * These values are taken into consideration by deal.II classes (such as
 * QProjector) to ensure that quantities computed on two adjacent cells use the
 * same quadrature points and shape function orderings. Unless you are working
 * on deal.II internals or new FiniteElement classes, it is not necessary to
 * consider these values. In practice, essentially no applications dependent on
 * deal.II ever need to handle orientation problems.
 * </dd>
 *
 * <dt class="glossary">@anchor GlossGeneralizedSupport <b>Generalized support points</b></dt>
 * <dd>"Generalized support points" are, as the name suggests, a
 * generalization of @ref GlossSupport "support points". The latter
 * are used to describe that a finite element simply <i>interpolates</i>
 * values at individual points (the "support points"). If we call these
 * points $\hat{\mathbf{x}}_i$ (where the hat indicates that these points
 * are defined on the reference cell $\hat{K}$), then one typically defines
 * shape functions $\varphi_j(\mathbf{x})$ in such a way that the
 * <i>nodal functionals</i> $\Psi_i[\cdot]$ simply evaluate the function
 * at the support point, i.e., that $\Psi_i[\varphi]=\varphi(\hat{\mathbf{x}}_i)$,
 * and the basis is chosen so that $\Psi_i[\varphi_j]=\delta_{ij}$ where
 * $\delta_{ij}$ is the Kronecker delta function. This leads to the common
 * @ref GlossLagrange "Lagrange elements".
 *
 * (In the vector valued case, the only other piece of information
 * besides the support points $\hat{\mathbf{x}}_i$ that one needs to provide
 * is the <i>vector component</i> $c(i)$ the $i$th node functional
 * corresponds, so that $\Psi_i[\varphi]=\varphi(\hat{\mathbf{x}}_i)_{c(i)}$.)
 *
 * On the other hand, there are other kinds of elements that are not
 * defined this way. For example, for the lowest order Raviart-Thomas element
 * (see the FE_RaviartThomas class), the node functional evaluates not
 * individual components of a vector-valued finite element function with @p dim
 * components, but the <i>normal component</i> of this vector:
 * $\Psi_i[\varphi]
 *  =
 *  \varphi(\hat{\mathbf{x}}_i) \cdot \mathbf{n}_i
 * $, where the $\mathbf{n}_i$ are the normal vectors to the face of the cell
 * on which $\hat{\mathbf{x}}_i$ is located. In other words, the node functional
 * is a <i>linear combination</i> of the components of $\varphi$ when
 * evaluated at $\hat{\mathbf{x}}_i$. Similar things happen for the BDM,
 * ABF, and Nedelec elements (see the FE_BDM, FE_ABF, FE_Nedelec classes).
 *
 * In these cases, the element does not have <i>support points</i> because
 * it is not purely interpolatory; however, some kind of interpolation
 * is still involved when defining shape functions as the node functionals
 * still require point evaluations at special points $\hat{\mathbf{x}}_i$.
 * In these cases, we call the points <i>generalized support points</i>.
 *
 * Finally, there are elements that still do not fit into this
 * scheme. For example, some hierarchical basis functions
 * (see, for example the FE_Q_Hierarchical element) are defined so
 * that the node functionals are <i>moments</i> of finite element
 * functions,
 * $\Psi_i[\varphi]
 *  =
 *  \int_{\hat{K}} \varphi(\hat{\mathbf{x}})
 *  {\hat{x}_1}^{p_1(i)}
 *  {\hat{x}_2}^{p_2(i)}
 * $ in 2d, and similarly for 3d, where the $p_d(i)$ are the order
 * of the moment described by shape function $i$. Some other elements
 * use moments over edges or faces. In all of these cases, node functionals
 * are not defined through interpolation at all, and these elements then
 * have neither support points, nor generalized support points.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor geometry_paper <b>geometry paper</b></dt>
 * <dd>The "geometry paper" is a paper by L. Heltai, W. Bangerth, M. Kronbichler,
 * and A. Mola, titled
 * "Using exact geometry information in finite element computations", that
 * describes how deal.II describes the geometry of domains. In particular,
 * it discusses the algorithmic foundations on which the Manifold class
 * is based, and what kind of information it needs to provide for mesh
 * refinement, the computation of normal vectors, and the many other places
 * where geometry enters into finite element computations.
 *
 * The paper is currently available on arXiv at https://arxiv.org/abs/1910.09824 .
 * The full reference for this paper is as follows:
 * @code{.bib}
@misc{heltai2019using,
    title={Using exact geometry information in finite element computations},
    author={Luca Heltai and Wolfgang Bangerth and Martin Kronbichler and Andrea Mola},
    year={2019},
    eprint={1910.09824},
    archivePrefix={arXiv},
    primaryClass={math.NA}
}
 * @endcode
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossGhostCell <b>Ghost cells</b></dt>
 * <dd>
 * If a mesh is distributed across multiple MPI processes using the
 * parallel::distributed::Triangulation class, each processor stores
 * only the cells it owns, one layer of adjacent cells that are owned
 * by other processors, all @ref GlossCoarseMesh "coarse level cells",
 * and all cells that are necessary to maintain the invariant that adjacent
 * cells must differ by at most one refinement level. The cells stored on
 * each process that are not owned by this process but that are adjacent to the
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
 * Triangulation and the parallel::shared::Triangulation classes.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossGhostedVector <b>Ghosted vectors</b></dt>
 * <dd>
 * In parallel computations, vectors come in two general kinds:
 * without and with ghost elements. Vectors without ghost
 * elements uniquely partition the vector elements between
 * processors: each vector entry has exactly one processor that
 * owns it, and this is the only processor that stores the value
 * of this entry. In other words, if processor zero stores elements
 * 0...49 of a vector and processor one stores elements 50...99,
 * then processor one is out of luck accessing element 42 of this
 * vector: it is not stored here and the value can not be assessed.
 * This will result in a failed assertion.
 *
 * On the other hand, there are many situations where one *needs* to
 * know vector elements that aren't locally owned, for example to
 * evaluate the solution on a locally owned cell (see
 * @ref GlossLocallyOwnedCell) for which one of the degrees of freedom
 * is at an interface to a cell that we do not own locally (which,
 * in this case must then be a @ref GlossGhostCell "ghost cell")
 * and for which the neighboring cell may be the owner -- in other
 * words, the degree of freedom is not a
 * @ref GlossLocallyOwnedDof "locally owned" but instead only a
 * @ref GlossLocallyActiveDof "locally active" DoF. The values of such
 * degrees of freedom are typically stored on the machine that owns the
 * degree of freedom and, consequently, would not be accessible on the
 * current machine.
 *
 * Because one often needs these values anyway, there is a second kind of
 * vector, often called "ghosted vector". Ghosted vectors store some elements
 * on each processor for which that processor is not the owner.
 * For such vectors, you can read those elements that the
 * processor you are currently on stores but you cannot write into them
 * because to make this work would require propagating the new value to
 * all other processors that have a copy of this value (the list of
 * such processors may be something which the current processor does not
 * know and has no way of finding out efficiently). Since you cannot
 * write into ghosted vectors, the only way to initialize such a vector
 * is by assignment from a non-ghosted vector. This implies having to
 * import those elements we locally want to store from other processors.
 *
 * The way ghosted vectors are actually stored is different between the
 * various implementations of parallel vectors. For PETSc (and the corresponding
 * PETScWrappers::MPI::Vector class), ghosted vectors store the same
 * elements as non-ghosted ones would, plus some additional elements
 * that are owned by other processors. In other words, for each element
 * there is a clear owner among all of the processors and those elements
 * that the current processor stores but does not own (i.e., the
 * "ghost elements") are simply mirror images of a primary value somewhere
 * else -- thus, the name "ghost". This is also the case for the
 * LinearAlgebra::distributed::Vector class.
 *
 * On the other hand, in Trilinos (and consequently in
 * TrilinosWrappers::MPI::Vector), a ghosted vector is simply a view
 * of the parallel vector where the element distributions overlap. The
 * 'ghosted' Trilinos vector in itself has no idea of which entries
 * are ghosted and which are locally owned. In fact, a ghosted vector
 * may not even store all of the elements a non-ghosted vector would
 * store on the current processor. Consequently, for Trilinos vectors,
 * there is no notion of an 'owner' of vector elements in the way we
 * have it in the non-ghost case view (or in the PETSc case) and
 * the name "ghost element" may be misleading since in this view,
 * every element we have available locally may or may not be stored
 * somewhere else as well, but even if it is, the local element is not
 * a mirror value of a primary location as there is no owner of each
 * element.
 *
 * In the end, there are two key take-away messages from the separation between
 * ghosted and non-ghosted vectors:
 * - Ghosted vectors are read-only. You cannot write into them, or add
 *   one such vector into another.
 * - Even if every process participates in storing a vector with
 *   ghost elements, not all elements of the vector may be stored anywhere;
 *   some elements may be stored on multiple processes but no process may
 *   be a designated "owner" of these elements. As a consequence, reduction
 *   operations such as dot products or norms may not be computed on
 *   vectors with ghost elements because in these storage schemes, it is
 *   not possible to ensure that each entry of the vector is counted exactly
 *   once.
 *
 * @note The @ref distributed documentation topic provides a brief
 * overview of where the different kinds of vectors are typically
 * used.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor hp_paper <b>%hp-paper</b></dt>
 * <dd>The "hp-paper" is a paper by W. Bangerth and O. Kayser-Herold, titled
 * "Data Structures and Requirements for hp Finite Element Software", that
 * describes many of the algorithms and data structures used in the implementation
 * of the hp-framework of deal.II. In particular, it summarizes many of the
 * tricky points that have to be considered for %hp-finite elements using continuous
 * elements.
 *
 * The full reference for this paper is as follows:
 * @code{.bib}
@Article{BK07,
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
 * It is available from <a href="https://www.math.colostate.edu/~bangerth/publications.html">https://www.math.colostate.edu/~bangerth/publications.html</a>, also see <a href="https://www.dealii.org/publications.html#details">deal.II publications</a> for details.
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
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossInvalidValue <b>Invalid value</b></dt>
 * <dd>
 * A common problem in software design is what to do if a function needs to
 * return something like "this value does not exist". An example of this
 * could be the function `IndexSet::index_within_set(i)` that returns the
 * how many'th element of the set `i` is. Clearly, the return value of this
 * function should be an unsigned integer type, the result always being
 * a count (zero or positive). The question is what to do if the index
 * `i` is not actually in the set. One *could* consider this a bug: You
 * can't ask an index set for the position of an index that is not in the
 * set, and so an exception should be thrown: The user should first check
 * with `IndexSet::is_element(i)` whether `i` is an element of the set,
 * and only then should they call `IndexSet::index_within_set(i)`.
 * But sometimes there are situations where one simply wants to
 * return a regular count if `i` is in the set, and some kind of
 * "exceptional" value if it is not.
 *
 * Similar questions appear when one writes code of the following kind:
 * @code
 *   unsigned int value;
 *   if (some condition)
 *     value = 13;
 *   [...] // much code
 *   if (some other condition)
 *     value = 42;
 *
 *   launch_the_rocket(value); // something important and expensive
 * @endcode
 * Here, the programmer may know that either `some condition` or
 * `some other condition` is true, and that consequently `value` is
 * always initialized at the end of the block. But there are good
 * reasons not to trust this. First, programmers make mistakes, and
 * so it is conceivable that there are situations where the variable
 * ends up uninitialized, even though that wasn't intended. Second,
 * code changes over time and while the "either `some condition` or
 * `some other condition` is true" situation may hold at the time
 * of development of the software, when code is moved around
 * or undergoes bug fixes and functionality enhancement, things may
 * change and the variable may go uninitialized. A better way
 * to write this code would be like this:
 * @code
 *   unsigned int value = some_invalid_value;
 *   if (some condition)
 *     value = 13;
 *   [...] // much code
 *   if (some other condition)
 *     value = 42;
 *
 *   Assert (value != some_invalid_value, "some error message");
 *   launch_the_rocket(value); // something important and expensive
 * @endcode
 * As before, the issue is what `some_invalid_value` should be.
 *
 * This is such a common problem that many, mostly ad-hoc, solutions
 * are widely used. In some cases, parts of the values of a type can
 * be used. For example, `sqrt` must necessarily return a
 * non-negative value because, well, all square roots are non-negative.
 * As a consequence, this function could return error codes as negative
 * values. (In truth, though,
 * [`sqrt`](https://en.cppreference.com/w/cpp/numeric/math/sqrt)
 * returns `NaN` if one provides
 * a negative input -- here, `NaN` stands for "not a number" and is,
 * just like negative numbers, a stand-in for a value that can not
 * happen as part of the regular operations of this function and can
 * consequently be used to indicate errors.) Similarly, functions such as
 * [`printf()`](https://en.cppreference.com/w/c/io/fprintf)
 * either return the number of characters printed or, in case
 * of an error in the inputs, a negative value. Finally, the
 * [`fopen`](https://en.cppreference.com/w/cpp/io/c/fopen)
 * function returns a pointer to a file descriptor (similar to a
 * `std::iofstream` object) but if the file cannot be opened,
 * for example because the file does not exist or the directory in which
 * it supposedly is does not exist, then the function returns a `nullptr`.
 *
 * All of these examples use that the returned object is of a type whose
 * set of possible values contains values that cannot be legitimate
 * return values (in mathematical language: they are not part of the
 * "range" of the function) and that can consequently be used to indicate
 * errors. This is awkward because the mapping of error codes into the
 * space of possible return values depends very much on the function and
 * what it can and cannot return. For example, `sqrt()` could return `-1`
 * as an error, but `sin()` can not because `-1` is a valid return value.
 * Also, not all functions allow for this. For example, the
 * [`strtol`](https://en.cppreference.com/w/c/string/byte/strtol)
 * ("string to long (integer)") function takes a string as input and
 * returns a long integer as output. But since clearly every possible
 * value of the type `long int` can legitimately be returned, errors in
 * the input (say, if someone provided the string `"nonsense"` as
 * input) cannot be indicated via the return object, and the function
 * needs to indicate errors through another mechanism. The C language
 * does that by letting functions such as `strtol` set the global variable
 * `errno` to a nonzero value to indicate an error (an approach that
 * comes with its own set of problems, among which are that people
 * tend to forget the value of this variable after calling the function).
 *
 * The examples listed above date back to the time when C was first
 * developed, in the late 1960s and 1970s. C++ solves this conundrum
 * in a more systematic way. First, functions can throw exceptions of
 * any type, and one can think of a thrown exception as simply another
 * possible return value of functions that indicates errors without
 * requiring having a part of the value space of the return type of a
 * function occupied for error codes. In fact, the type of an exception
 * is completely decoupled from the usual return type: You can pass as
 * much information through exceptions you throw, even if the function
 * in question returns just a meager `int` in regular operation. This
 * approach is used in a number of deal.II functions: If inputs don't
 * make sense, the program is either aborted (typically via an
 * `Assert` statement) if the inputs are believed to be hard-coded --
 * say, when adding vectors of different length -- or an exception is
 * thrown via C++'s `throw` statement. The function
 * Mapping::transform_real_to_unit_cell() is an example of the latter.
 * Second, in newer C++ standards, one can use the
 * [`std::expected<T,E>`](https://en.cppreference.com/w/cpp/utility/expected)
 * class as the return value that can be thought as "this function
 * returns objects of type `T`, but if an error was detected, then the
 * function instead returns an object of type `E`". You can then ask
 * the returned object whether it stores one or the other. In cases of
 * errors, one would typically store an explanation of the error in `E`,
 * in much the same way as the function could throw an exception of type
 * `E`. For example, a perhaps better design for the the `fopen` function
 * mentioned above could return `std::expected<FILE,std::string>` where
 * if successful, it returns a `FILE` object that identifies the file
 * for writing and reading; if it fails, the function would store a textual
 * description of what went wrong in the second slot of the `std::expected`
 * object (or perhaps an element of an `enum` that simply provides an
 * enumeration of possible reasons for failure). Relatedly, if it is not
 * necessary to provide a reason for the failure, functions could simply
 * return an object of type
 * [`std::optional<T>`](https://en.cppreference.com/w/cpp/utility/optional)
 * that may or may not hold an object of type `T`, and that one can ask
 * about that. This would be the right approach for the
 * `IndexSet::index_within_set(i)` function mentioned above: If `i` is
 * an element of the set, then it returns an object of type
 * `std::optional<IndexSet::size_type>` that contains the requested value;
 * if `i` was not in the set, then it returns an empty
 * `std::optional<IndexSet::size_type>` object.
 *
 * A third approach, widely used in deal.II, is to *explicitly* declare
 * part of the range space as "exceptional". For example, many functions
 * in deal.II deal with indices of degrees of freedom. These are encoded
 * as unsigned integers, except that we explicitly declare the value
 * 4294967295 as an invalid value that indicates an error. (This specific
 * value happens to be the largest unsigned integer; computations are
 * unlikely to be so large that they use this specific value in a legitimate
 * sense.) Many of the data types used in deal.II, such as
 * types::global_dof_index, types::active_fe_index, types::material_id
 * explicitly consider one possible value representable by these
 * types as "invalid" and use it to report errors of uninitialized
 * variables. These values typically have names such as
 * numbers::invalid_unsigned_int, numbers::invalid_material_id, etc.
 *
 * (As a postscript, the `strtol` function mentioned above uses this sort of
 * approach as well. If the input to that function is invalid, it not only
 * sets the global variable `errno`, but *also* returns either `LONG_MAX`
 * or `LONG_MIN`. These are the largest and smallest long integer values.
 * In other words, the function's definition *explicitly* marks these values
 * as "invalid" or "exceptional", even though one could legitimately expect
 * to provide the function with a string for which the conversion to a long
 * integer would result in these values. This is at its core the same
 * approach we use in deal.II with the invalid values mentioned above,
 * except that deal.II uses variable names that reflect the underlying
 * use case (such as whether a value reflects an invalid value
 * for DoF indices or manifold ids), rather than just the type: When
 * using numbers::invalid_material_id, you don't need to know what
 * type is actually used to represent material ids.)
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLagrange <b>Lagrange elements</b></dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points"
 * are called "Lagrange elements". Their node functionals correspond
 * to evaluation of shape functions at these support points.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyOwnedCell <b>Locally owned cell</b></dt>
 * <dd>This concept identifies a subset of all cells when using
 * distributed meshes, see the @ref distributed topic. In such meshes, each
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
 * distributed meshes, see the @ref distributed topic.  Locally owned degrees
 * of freedom live on locally owned cells. Since degrees of freedom are owned
 * by only one processor, degrees of freedom on interfaces between cells owned
 * by different processors may be owned by one or the other, so not all
 * degrees of freedom on a locally owned cell are also locally owned degrees
 * of freedom.
 *
 * Locally owned DoFs are a subset of the
 * @ref GlossLocallyActiveDof "locally active DoFs".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLocallyActiveDof <b>Locally active degrees of freedom</b></dt>
 * <dd>This concept identifies a subset of all degrees of freedom when using
 * distributed meshes, see the @ref distributed topic.  Locally active degrees
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
 * distributed meshes, see the @ref distributed topic.  Locally relevant
 * degrees of freedom are those that live on locally owned or ghost cells.
 * Consequently, they may be owned by different processors.
 *
 * Locally relevant DoFs are a superset of the
 * @ref GlossLocallyActiveDof "locally active DoFs."
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossLumpedMassMatrix <b>Lumped mass matrix</b></dt>
 * <dd>The @ref GlossMassMatrix "mass matrix" is a matrix of the form
 *   @f{align*}{
 *     M_{ij} = \int_\Omega \varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx,
 *   @f}
 * It frequently appears in the solution of time dependent problems where, if
 * one uses an explicit time stepping method, it then leads to the need
 * to solve problems of the form
 *   @f{align*}{
 *     MU^n = MU^{n-1} + k_n BU^{n-1},
 *   @f}
 * in time step $n$, where $U^n$ is the solution to be computed, $U^{n-1}$ is the
 * known solution from the first time step, and $B$ is a matrix related to the
 * differential operator in the PDE. $k_n$ is the size of the time step. A similar
 * linear system of equations also arises out of the discretization of second-order
 * differential equations.
 *
 * The presence of the matrix $M$ on the left side is a nuisance because, even
 * though we have used an explicit time stepping method, we still have to solve a
 * linear system in each time step. It would be much preferable if the matrix were
 * diagonal. "Lumping" the mass matrix is a strategy to replace $M$ by a matrix
 * $M_\text{diagonal}$ that actually is diagonal, yet does not destroy the accuracy
 * of the resulting solution.
 *
 * Historically, mass lumping was performed by adding the elements of a row
 * together and setting the diagonal entries of $M_\text{diagonal}$ to that
 * sum. This works for $Q_1$ and $P_1$ elements, for example, and can be
 * understood mechanically by replacing the continuous medium we are
 * discretizating by one where the continuous mass distribution is replaced by
 * one where (finite amounts of) mass are located only at the nodes. That is,
 * we are "lumping together" the mass of an element at its vertices, thus
 * giving rise to the name "lumped mass matrix". A more mathematical perspective
 * is to compute the integral above for $M_{ij}$ via special quadrature rules;
 * in particular, we replace the computation of
 *   @f{align*}{
 *     M_{ij} = \int_\Omega \varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx
 *            = \sum_K \int_K \varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx,
 *   @f}
 * by quadrature
 *   @f{align*}{
 *     (M_{\text{diagonal}})_{ij} = \sum_K \sum_q \varphi_i(\mathbf x_q^K) \varphi_j(\mathbf x_q^K)
 *     |K| w_q,
 *   @f}
 * where we choose the quadrature points as the *nodes* at which the
 * shape functions are defined. If we order the quadrature points in the
 * same way as the shape functions, then
 *   @f{align*}{
 *     \varphi_i(\mathbf x_q^K) = \delta_{iq},
 *   @f}
 * and consequently
 *   @f{align*}{
 *     (M_{\text{diagonal}})_{ij} = \delta_{ij} \sum_{K, \text{supp}\varphi_i \cap K \neq \emptyset} |K| w_i,
 *   @f}
 * where the sum extends over those cells on which $\varphi_i$ is nonzero.
 * The so-computed mass matrix is therefore diagonal.
 *
 * Whether or not this particular choice of quadrature formula is
 * sufficient to retain the convergence rate of the discretization is
 * a separate question. For the usual $Q_k$ finite elements
 * (implemented by FE_Q and FE_DGQ), the appropriate quadrature
 * formulas are of QGaussLobatto type. Mass lumping can also be done
 * with FE_SimplexP_Bubbles, for example, if appropriate quadrature rules
 * are chosen.
 *
 * For an example of where lumped mass matrices play a role, see step-69.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossManifoldIndicator <b>%Manifold indicator</b></dt>
 *
 * <dd> Every object that makes up a Triangulation (cells, faces,
 * edges, etc.), is associated with a unique number (of type
 * types::manifold_id) which is used to identify which Manifold object
 * describes the coordinate system on that cell, e.g., a PolarManifold will
 * perform calculations in polar coordinates. For example, Manifold objects
 * are responsible for generating new points when a Triangulation is refined,
 * defining the locations of @ref GlossSupport "support points", and defining
 * the locations of quadrature points.
 *
 * By default, all manifold indicators of a mesh are set to
 * numbers::flat_manifold_id. A typical piece of code that sets the
 * manifold indicator on a object to something else would look like
 * this, here setting the manifold indicator to 42 for all cells whose
 * center has an $x$ component less than zero:
 *
 * @code
 * for (auto &cell : triangulation.active_cell_iterators())
 *   if (cell->center()[0] < 0)
 *     cell->set_manifold_id(42);
 * @endcode
 *
 * Here we call the function TriaAccessor::set_manifold_id(). It may
 * also be appropriate to call TriaAccessor::set_all_manifold_ids()
 * instead, to set recursively the manifold id on each face (and edge,
 * if in 3d). To query the manifold indicator of a particular object
 * edge, use TriaAccessor::manifold_id().
 *
 * Every manifold id set on a Triangulation must have an associated Manifold
 * object. This is assigned via Triangulation::set_manifold().
 *
 * @note Manifold indicators are inherited from parents to their
 * children upon mesh refinement. Some more information about manifold
 * indicators is also presented in a section of the documentation of
 * the Triangulation class as well as in the
 * @ref manifold "Manifold documentation topic". Manifold indicators
 * are used in step-53 and step-54.
 * </dd>
 *
 * @see @ref manifold "The topic on Manifolds"
 *
 *
 * <dt class="glossary">@anchor GlossMassMatrix <b>Mass matrix</b></dt>
 * <dd>The "mass matrix" is a matrix of the form
 *   @f{align*}{
 *     M_{ij} = \int_\Omega \varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx,
 *   @f}
 * possibly with a coefficient inside the integral, and
 * where $\varphi_i(\mathbf x)$ are the shape functions of a finite element.
 * The origin of the term refers to the fact that in structural mechanics
 * (where the finite element method originated), one often starts from the
 * elastodynamics (wave) equation
 *   @f{align*}{
 *     \rho \frac{\partial^2 u}{\partial t^2}
 *     -\nabla \cdot C \nabla u = f.
 *   @f}
 * If one multiplies this equation by a test function $\varphi_i$,
 * integrates over $\Omega$, and then discretizes by the substitution
 * $u(\mathbf x,t) \to u_h(\mathbf x)=\sum_j U_j(t) \varphi_j(\mathbf x)$,
 * then the first term above results in
 *   @f{align*}{
 *     \sum_j \left[\int_\Omega \rho \varphi_i \varphi_j \right]
 *     \frac{\partial^2 U_j(t)}{\partial t^2}
 *   @f}
 * which can be written as
 *   @f{align*}{
 *     M
 *     \frac{\partial^2 U(t)}{\partial t^2}
 *   @f}
 * where
 *   @f{align*}{
 *     M_{ij} = \int_\Omega \rho(\mathbf x)\varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx.
 *   @f}
 * Since the matrix entries are a (weighted) integral over a mass density, they
 * have the units of "mass", giving the "mass matrix" its name.
 *
 * In mathematics, where we often consider non-dimensionalized equations, we
 * end up with the case $\rho=1$, and as a consequence the matrix without
 * the coefficient,
 *   @f{align*}{
 *     M_{ij} = \int_\Omega \varphi_i(\mathbf x) \varphi_j(\mathbf x)\; dx,
 *   @f}
 * also carries the name "mass matrix".
 *
 * The mass matrix is almost always written with the symbol $M$. See, for example,
 * step-23, step-26, and a number of the other time dependent equations solved by
 * tutorial programs.
 *
 * The mass matrix is occasionally approximated by a diagonal matrix,
 * see the glossary entry for @ref GlossLumpedMassMatrix "lumped mass matrix".
 * See also the @ref GlossStiffnessMatrix "stiffness matrix"
 * for a related case.
 * </dt>
 *
 *
 * <dt class="glossary">@anchor GlossMaterialId <b>Material id</b></dt>
 * <dd>Each cell of a triangulation has associated with it a property called
 * "material id". It is commonly used in problems with heterogeneous
 * coefficients to identify which part of the domain a cell is in and,
 * consequently, which value the coefficient should have on this particular
 * cell. In practice, the material id of a cell
 * is typically used to identify which cells belong to a particular part of
 * the domain, e.g., when you have different materials (steel, concrete, wood)
 * that are all part of the same domain. One would then usually query the
 * material id associated with a cell during assembly of the bilinear form,
 * and use it to determine (e.g., by table lookup, or a sequence of if-else
 * statements) what the correct material coefficients would be for that cell.
 *
 * This material_id may be set upon construction of a triangulation (through
 * the CellData data structure), or later through use of cell iterators. For a
 * typical use of this functionality, see the step-28 tutorial program. The
 * functions of the GridGenerator namespace typically set the material ID of
 * all cells to zero. When reading a triangulation through the GridIn class,
 * different input file formats have different conventions, but typically
 * either explicitly specify the material id, or if they don't, then GridIn
 * simply sets them to zero. Because the material of a cell is intended
 * to pertain to a particular region of the domain, material ids are inherited
 * by child cells from their parent upon mesh refinement. However, if material
 * ids are being assigned in a parallel distributed computation, any
 * refinement or coarsening step will disregard the id assignments when
 * a cell is moved to another process. So one must be careful when a mesh
 * has more than one material id and one plans on using multiple processors
 * with refinement. It is however safe to assign material ids consistently to
 * all coarse cells on all MPI ranks before executing any refinement.
 *
 * The material id is set and queried using the CellAccessor::material_id,
 * CellAccessor::set_material_id and CellAccessor::recursively_set_material_id
 * functions.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossMPICommunicator <b>MPI Communicator</b></dt>
 * <dd>
 * In the language of the Message Passing Interface (MPI), a communicator
 * can be thought of as a mail system that allows sending messages to
 * other members of the mail system. Within each communicator, each
 * @ref GlossMPIProcess "process" has a
 * @ref GlossMPIRank "rank" (the equivalent of a house number) that
 * allows to identify senders and receivers of messages. It is not
 * possible to send messages via a communicator to receivers that are
 * not part of this communicator/mail service.
 *
 * When starting a parallel program via a command line call such as
 * @code
 *  mpirun -np 32 ./step-17
 * @endcode
 * (or the equivalent used in the batch submission system used on your
 * cluster) the MPI system starts 32 copies of the step-17 executable.
 * Each of these has access to the <code>MPI_COMM_WORLD</code> communicator
 * that then consists of all 32 processors, each with its own rank. A subset
 * of processes within this MPI universe can later agree to create other
 * communicators that allow communication between only a subset of
 * processes.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossMPIProcess <b>MPI Process</b></dt>
 * <dd>
 * When running parallel jobs on distributed memory machines, one
 * almost always uses MPI. There, a command line call such as
 * @code
 *  mpirun -np 32 ./step-17
 * @endcode
 * (or the equivalent used in the batch submission system used on your
 * cluster) starts 32 copies of the step-17 executable. Some of these may actually
 * run on the same machine, but in general they will be running on different
 * machines that do not have direct access to each other's memory space.
 *
 * In the language of the Message Passing Interface (MPI), each of these
 * copies of the same executable running on (possibly different) machines
 * are called <i>processes</i>. The collection of all processes running in
 * parallel is called the "MPI Universe" and is identified by the
 * @ref GlossMPICommunicator "MPI communicator" <code>MPI_COMM_WORLD</code>.
 *
 * Each process has immediate access only to the objects in its own
 * memory space. A process can not read from or write into the memory
 * of other processes. As a consequence, the only way by which
 * processes can communicate is by sending each other messages. That
 * said (and as explained in the introduction to step-17), one
 * typically calls higher level MPI functions in which all processes
 * that are part of a communicator participate. An example would
 * be computing the sum over a set of integers where each process
 * provides one term of the sum.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossMPIRank <b>MPI Rank</b></dt>
 * <dd>
 * In the language of the Message Passing Interface (MPI), the <i>rank</i>
 * of an @ref GlossMPIProcess "MPI process" is the number this process
 * carries within the set <code>MPI_COMM_WORLD</code> of all processes
 * currently running as one parallel job. More correctly, it is the
 * number within an @ref GlossMPICommunicator "MPI communicator" that
 * groups together a subset of all processes with one parallel job
 * (where <code>MPI_COMM_WORLD</code> simply denotes the <i>complete</i>
 * set of processes).
 *
 * Within each communicator, each process has a unique rank, distinct from the
 * all other processes' ranks, that allows
 * identifying one recipient or sender in MPI communication calls. Each
 * process, running on one processor, can inquire about its own rank
 * within a communicator by calling Utilities::MPI::this_mpi_process().
 * The total number of processes participating in a communicator (i.e.,
 * the <i>size</i> of the communicator) can be obtained by calling
 * Utilities::MPI::n_mpi_processes().
 * </dd>
 *
 *
 * <dt class="glossary">@anchor mg_paper <b>%Multigrid paper</b></dt>
 * <dd>The "multigrid paper" is a paper by B. Janssen and G. Kanschat, titled
 * "Adaptive Multilevel Methods with Local Smoothing for H1- and Hcurl-Conforming High Order Finite Element Methods", that
 * describes many of the algorithms and data structures used in the implementation
 * of the multigrid framework of deal.II. It underlies the implementation of
 * the classes that are used in step-16 for multigrid
 * methods.
 *
 * The full reference for this paper is as follows:
 * @code{.bib}
@article{janssen2011adaptive,
  title=    {Adaptive Multilevel Methods with Local Smoothing for H^1- and H^{curl}-Conforming High Order Finite Element Methods},
  author=   {Janssen, B{\"a}rbel and Kanschat, Guido},
  journal=  {SIAM Journal on Scientific Computing},
  volume=   {33},
  number=   {4},
  pages=    {2095--2114},
  year=     {2011},
  publisher={SIAM}}
 * @endcode
 * See
 * <a href="https://dx.doi.org/10.1137/090778523">DOI:10.1137/090778523</a>
 * for the paper and <a href="https://www.dealii.org/publications.html#details">deal.II publications</a> for more details.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossNodes <b>Node values or node functionals</b></dt>
 *
 * <dd>It is customary to define a finite element as a triple
 * $(K,P,\Psi)$ where
 * - $K$ is the cell, where in deal.II this is always a line segment,
 *   quadrilateral, or hexahedron;
 * - $P$ is a finite-dimensional space, e.g., a polynomial space mapped
 *   from the @ref GlossReferenceCell "reference cell" to $K$;
 * - $\Psi$ is a set of "node functionals", i.e., functionals
 *   $\Psi_i : P \rightarrow {\mathbb R}$.
 * The dimension of $P$ must be equal to the number of node functionals.
 * With this definition, we can define a basis of the local function space,
 * i.e., a set of "shape functions" $\varphi_j\in P$, by requiring that
 * $\Psi_i(\varphi_j) = \delta_{ij}$, where $\delta$ is the Kronecker delta.
 *
 * This definition of what a finite element is has several advantages,
 * concerning analysis as well
 * as implementation. For the analysis, it means that conformity with
 * certain spaces (FiniteElementData::Conformity), e.g. continuity, is
 * up to the node functionals. In deal.II, it helps simplifying the
 * implementation of more complex elements like FE_RaviartThomas
 * considerably.
 *
 * Examples for node functionals are values in
 * @ref GlossSupport "support points" and moments with respect to Legendre
 * polynomials. Examples:
 *
 * <table><tr>
 *   <th>Element</th>
 *   <th>%Function space</th>
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
 * The construction of finite elements as outlined above allows writing
 * code that describes a finite element simply by providing a polynomial
 * space (without having to give it any particular basis -- whatever is convenient
 * is entirely sufficient) and the nodal functionals. This is used, for example
 * in the FiniteElement::convert_generalized_support_point_values_to_dof_values()
 * function.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossParallelScaling <b>Parallel scaling</b></dt>
 * <dd>When we say that a parallel program "scales", what we mean is that the
 * program does not become unduly slow (or takes unduly much memory) if we
 * make the problem it solves larger, and that run time and memory consumption
 * decrease proportionally if we keep the problem size the same but increase
 * the number of processors (or cores) that work on it.
 *
 * More specifically, think of a problem whose size is given by a number $N$
 * (which could be the number of cells, the number of unknowns, or some other
 * indicative quantity such as the number of CPU cycles necessary to solve
 * it) and for which you have $P$ processors available for solution. In an
 * ideal world, the program would then require a run time of ${\cal O}(N/P)$,
 * and this would imply that we could reduce the run time to any desired
 * value by just providing more processors. Likewise, for a program to be
 * scalable, its overall memory consumption needs to be ${\cal O}(N)$ and on
 * each involved process needs to be ${\cal O}(N/P)$, again
 * implying that we can fit any problem into the fixed amount of memory
 * computers attach to each processor, by just providing
 * sufficiently many processors.
 *
 * For practical assessments of scalability, we often distinguish between
 * "strong" and "weak" scalability. These assess asymptotic statements
 * such as ${\cal O}(N/P)$ run time in the limits $N\rightarrow \infty$
 * and/or $P\rightarrow \infty$. Specifically, when we say that a program
 * is "strongly scalable", we mean that if we have a problem of fixed
 * size $N$, then we can reduce the run time and memory consumption (on
 * every processor) inversely proportional to $P$ by just throwing more
 * processors at the problem. In particular, strong scalability implies
 * that if we provide twice as many processors, then run time and memory
 * consumption on every process will be reduced by a factor of two. In
 * other words, we can solve the <i>same problem</i> faster and faster
 * by providing more and more processors.
 *
 * Conversely, "weak scalability" means that if we increase the problem
 * size $N$ by a fixed factor, and increase the number of processors
 * $P$ available to solve the problem by the same factor, then the
 * overall run time (and the memory consumption on every processor)
 * remains the same. In other words, we can solve <i>larger and larger
 * problems</i> within the same amount of wallclock time by providing
 * more and more processors.
 *
 * No program is truly scalable in this theoretical sense. Rather, all programs
 * cease to scale once either $N$ or $P$ grows larger than certain limits.
 * We therefore often say things such as "the program scales up to
 * 4,000 cores", or "the program scales up to 100,000,000 unknowns". There are
 * a number of reasons why programs cannot scale without limit; these can
 * all be illustrated by just looking at the (relatively simple) step-17
 * tutorial program:
 * - Sequential sections: Many programs have sections of code that
 *   either cannot or are not parallelized, i.e., where one processor has to do
 *   a certain, fixed amount of work that does not decrease just because
 *   there are a total of $P$ processors around. In step-17, this is
 *   the case when generating graphical output: one processor creates
 *   the graphical output for the entire problem, i.e., it needs to do
 *   ${\cal O}(N)$ work. That means that this function has a run time
 *   of ${\cal O}(N)$, regardless of $P$, and consequently the overall
 *   program will not be able to achieve ${\cal O}(N/P)$ run time but
 *   have a run time that can be described as $c_1N/P + c_2N$ where
 *   the first term comes from scalable operations such as assembling
 *   the linear system, and the latter from generating graphical
 *   output on process 0. If $c_2$ is sufficiently small, then the
 *   program might look like it scales strongly for small numbers of
 *   processors, but eventually strong scalability will cease. In
 *   addition, the program can not scale weakly either because
 *   increasing the size $N$ of the problem while increasing the
 *   number of processors $P$ at the same rate does not keep the
 *   run time of this one function constant.
 * - Duplicated data structures: In step-17, each processor stores the entire
 *   mesh. That is, each processor has to store a data structure of size
 *   ${\cal O}(N)$, regardless of $P$. Eventually, if we make the problem
 *   size large enough, this will overflow each processor's memory space
 *   even if we increase the number of processors. It is thus clear that such
 *   a replicated data structure prevents a program from scaling weakly.
 *   But it also prevents it from scaling strongly because in order to
 *   create an object of size ${\cal O}(N)$, one has to at the very
 *   least write into ${\cal O}(N)$ memory locations, costing
 *   ${\cal O}(N)$ in CPU time. Consequently, there is a component of the
 *   overall algorithm that does not behave as ${\cal O}(N/P)$ if we
 *   provide more and more processors.
 * - Communication: If, to pick just one example, you want to compute
 *   the $l_2$ norm of a vector of which all MPI processes store a few
 *   entries, then every process needs to compute the sum of squares of
 *   its own entries (which will require ${\cal O}(N/P)$ time, and
 *   consequently scale perfectly), but then every process needs to
 *   send their partial sum to one process that adds them all up and takes
 *   the square root. In the very best case, sending a message that
 *   contains a single number takes a constant amount of time,
 *   regardless of the overall number of processes. Thus, again, every
 *   program that does communication cannot scale strongly because
 *   there are parts of the program whose CPU time requirements do
 *   not decrease with the number of processors $P$ you allocate for
 *   a fixed size $N$. In reality, the situation is actually even
 *   worse: the more processes are participating in a communication
 *   step, the longer it will generally take, for example because
 *   the one process that has to add everyone's contributions has
 *   to add everything up, requiring ${\cal O}(P)$ time. In other words,
 *   CPU time <i>increases</i> with the number of processes, therefore
 *   not only preventing a program from scaling strongly, but also from
 *   scaling weakly. (In reality, MPI libraries do not implement $l_2$
 *   norms by sending every message to one process that then adds everything
 *   up; rather, they do pairwise reductions on a tree that doesn't
 *   grow the run time as ${\cal O}(P)$ but as ${\cal O}(\log_2 P)$,
 *   at the expense of more messages sent around. Be that as it may,
 *   the fundamental point is that as you add more processors, the
 *   run time will grow with $P$ regardless of the way the operation
 *   is actually implemented, and it can therefore not scale.)
 *
 * These, and other reasons that prevent programs from scaling perfectly can
 * be summarized in <a href="https://en.wikipedia.org/wiki/Amdahl%27s_law">
 * <i>Amdahl's law</i></a> that states that if a fraction $\alpha$
 * of a program's overall work $W$ can be parallelized, i.e., it can be
 * run in ${\cal O}(\alpha W/P)$ time, and a fraction $1-\alpha$ of the
 * program's work can not be parallelized (i.e., it consists either of
 * work that only one process can do, such as generating graphical output
 * in step-17; or that every process has to execute in a replicated way,
 * such as sending a message with a local contribution to a dedicated
 * process for accumulation), then the overall run time of the program
 * will be
 * @f{align*}
 *   T = {\cal O}\left(\alpha \frac WP + (1-\alpha)W \right).
 * @f}
 * Consequently, the "speedup" you get, i.e., the factor by which your
 * programs run faster on $P$ processors compared to running the program
 * on a single process (assuming this is possible) would be
 * @f{align*}
 *   S = \frac{W}{\alpha \frac WP + (1-\alpha)W}
 *     = \frac{P}{\alpha + (1-\alpha)P}.
 * @f}
 * If $\alpha<1$, which it is for all practically existing programs,
 * then $S\rightarrow \frac{1}{1-\alpha}$ as $P\rightarrow \infty$, implying
 * that there is a point where it does not pay off in any significant way
 * any more to throw more processors at the problem.
 *
 * In practice, what matters is <i>up to which problem size</i> or
 * <i>up to which number of processes</i> or <i>down to which size
 * of local problems ${\cal}(N/P)$</i> a program scales. For deal.II,
 * experience shows that on most clusters with a reasonable fast
 * network, one can solve problems up to a few billion unknowns,
 * up to a few thousand processors, and down to somewhere between
 * 40,000 and 100,000 unknowns per process. The last number is the
 * most relevant: if you have a problem with, say, $10^8$ unknowns,
 * then it makes sense to solve it on 1000-2500 processors since the
 * number of degrees of freedom each process handles remains at more
 * than 40,000. Consequently, there is enough work every process
 * has to do so that the ${\cal O}(1)$ time for communication does
 * not dominate. But it doesn't make sense to solve such a problem with
 * 10,000 or 100,000 processors, since each of these processor's local
 * problem becomes so small that they spend most of their time waiting
 * for communication, rather than doing work on their part of the work.
 * </dd>
 *
 * <dt class="glossary">@anchor GlossPeriodicConstraints <b>Periodic boundary
 * conditions</b></dt>
 * <dd>Periodic boundary condition are often used when only part of the physical
 * relevant domain is modeled. One assumes that the solution simply continues
 * periodically with respect to the boundaries that are considered periodic.
 * In deal.II, support for this is through DoFTools::make_periodicity_constraints()
 * and GridTools::collect_periodic_faces(). As soon as a
 * parallel::distributed::Triangulation is used also
 * parallel::distributed::Triangulation::add_periodicity() has to be called to make
 * sure that all the processes know about relevant parts of the triangulation on both
 * sides of the periodic boundary. A typical process for distributed triangulations would be:
 * -# Create a mesh
 * -# Gather the periodic faces using GridTools::collect_periodic_faces() (Triangulation)
 * -# Add the periodicity information to the mesh
 * using parallel::distributed::Triangulation::add_periodicity()
 * -# Gather the periodic faces using GridTools::collect_periodic_faces() (DoFHandler)
 * -# Add periodicity constraints using DoFTools::make_periodicity_constraints()
 *
 * An example for this can be found in step-45.
 * </dd>
 *
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
 * vector-value shape function may have several non-zero components.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossReferenceCell <b>Reference cell</b></dt>
 * <dd>The finite element method is typically described by providing shape
 * functions on a "reference cell" whose shape, along with the shape functions,
 * is then mapped to the actual cells of the mesh. These reference cells are
 * typically triangles or quadrilaterals for two-dimensional meshes;
 * tetrahedra, hexahedra, wedges, or pyramids for three-dimensional meshes;
 * and simple line segments in the one-dimensional case.
 *
 * Rather than hard-coding properties of the reference cell in all places
 * where one wants to know about, say, the number of vertices of a cell,
 * deal.II uses a single, central place to describe the properties of
 * reference cells: The ReferenceCell class. In loops over the cells of a
 * mesh, one typically asks for properties of these cells using the
 * call `cell->reference_cell()`.
 * </dd>
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
 * href="https://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/index.html"
 * target="_top">BOOST serialization</a> library. See there for examples on
 * how to save and restore objects.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossShape <b>Shape functions</b></dt>
 * <dd>The restriction of the finite element basis functions to a single
 * grid cell.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossStiffnessMatrix <b>Stiffness matrix</b></dt>
 * <dd>The "stiffness matrix" is a matrix of the form
 *   @f{align*}{
 *     A_{ij} = \int_\Omega \nabla\varphi_i(\mathbf x)
 *       \cdot \nabla\varphi_j(\mathbf x)\; dx,
 *   @f}
 * possibly with a coefficient inside the integral, and
 * where $\varphi_i(\mathbf x)$ are the shape functions of a finite element.
 * The term is also used for variations of the case above, for example
 * replacing the gradient by the symmetric gradient in the case where
 * the solution variable is vector-valued (e.g., in elasticity, or the
 * Stokes equations). The key feature is that in the integral, first
 * derivatives are applied to both the test and trial functions,
 * $\varphi_i,\varphi_j$.
 *
 * The origin of the term refers to the fact that in structural mechanics
 * (where the finite element method originated), one often starts from the
 * elastostatics equation
 *   @f{align*}{
 *     -\nabla \cdot C \nabla u = f.
 *   @f}
 * In this equation, $C$ is the stress-strain tensor that, informally
 * speaking, relates how much force one has to apply to obtain a
 * unit displacement. In other words, it encodes the "stiffness" of
 * the material: A large $C$, i.e., a large stiffness, means a large
 * required force for a desired displacement and the other way around.
 *
 * If one multiplies this equation by a test function $\varphi_i$,
 * integrates over $\Omega$, and then discretizes by the substitution
 * $u(\mathbf x,t) \to u_h(\mathbf x)=\sum_j U_j(t) \varphi_j(\mathbf x)$,
 * then after integration by parts one ends up with
 *   @f{align*}{
 *     \sum_j \left[\int_\Omega \nabla \varphi_i \cdot C \varphi_j \right]
 *     U_j
 *   @f}
 * which can be written as
 *   @f{align*}{
 *     AU
 *   @f}
 * where
 *   @f{align*}{
 *     A_{ij} = \int_\Omega \nabla\varphi_i(\mathbf x) \cdot C \nabla \varphi_j(\mathbf x)\; dx.
 *   @f}
 * Since the matrix entries are (weighted) integrals of the stiffness
 * coefficient, the resulting matrix is called the "stiffness matrix".
 *
 * In mathematics, where we often consider non-dimensionalized equations,
 * we end up with the case $C=1$, and as a consequence the matrix without
 * the coefficient,
 *   @f{align*}{
 *     A_{ij} = \int_\Omega \nabla\varphi_i(\mathbf x) \cdot \nabla\varphi_j(\mathbf x)\; dx
 *   @f}
 * which corresponds to the Laplace or Poisson equation,
 *   @f{align*}{
 *     -\Delta u = f,
 *   @f}
 * also carries the name "stiffness matrix".
 *
 * The stiffness matrix is almost always denotes by the symbol $A$. See, for example,
 * step-4, step-6, as well as a number of the time dependent equations considered in
 * programs such as step-23 or step-26.
 *
 * See also the @ref GlossStiffnessMatrix "stiffness matrix"
 * for a related case.
 * </dt>
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
 * least for programs that run in %parallel it usually denotes the
 * @ref GlossMPIRank "MPI rank" of the processor that "owns" this cell.
 *
 * For programs that are parallelized based on MPI but where each processor
 * stores the entire triangulation (as in, for example, step-17 and step-18,
 * but not in step-40), subdomain ids are assigned to cells by
 * partitioning a mesh, and each MPI process then only works on those cells it
 * "owns", i.e., that belong to a subdomain the processor owns
 * (traditionally, this is the case for the subdomain id whose numerical value
 * coincides with the rank of the MPI process within the MPI
 * communicator). Partitioning is typically done using the
 * GridTools::partition() function, but any other method can also be used to
 * do this. (Alternatively, the parallel::shared::Triangulation class can
 * partition the mesh automatically using a similar approach.)
 *
 * On the other hand, for programs that are parallelized using MPI but
 * where meshes are held distributed across several processors using
 * the parallel::distributed::Triangulation class, the subdomain id of
 * cells is tied to the processor that owns the cell. In other words,
 * querying the subdomain id of a cell tells you if the cell is owned
 * by the current processor (i.e. if <code>cell-@>subdomain_id() ==
 * triangulation.parallel::distributed::Triangulation::locally_owned_subdomain()</code>)
 * or by another processor. In the %parallel distributed case,
 * subdomain ids are only assigned to cells that the current processor
 * owns as well as the immediately adjacent @ref GlossGhostCell "ghost cells".
 * Cells further away are held on each processor to ensure
 * that every MPI process has access to the full
 * @ref GlossCoarseMesh "coarse grid" as well
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
 * <dt class="glossary">@anchor GlossSupport <b>Support points</b></dt>
 * <dd>Support points are by definition those points $p_i$, such that for the
 * shape functions $v_j$ holds $v_j(p_i) = \delta_{ij}$. Therefore, a finite
 * element interpolation can be defined uniquely by the values in the support
 * points.
 *
 * Lagrangian elements fill the vector accessed by
 * FiniteElement::get_unit_support_points(), such that the
 * function FiniteElement::has_support_points() returns
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
 * <dt class="glossary">@anchor GlossTargetComponent <b>Target component</b></dt>
 * <dd>
 * When vectors and matrices are grouped into blocks by component, it is
 * often desirable to collect several of the original components into
 * a single one. This could be for instance, grouping the velocities
 * of a Stokes system into a single block.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitCell <b>Unit cell</b></dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitSupport <b>Unit support points</b></dt>
 * <dd>These are the @ref GlossSupport "support points" on the reference cell, defined in
 * FiniteElement. For example, the usual Q1 element in 1d has support
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
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossUserProvidedCallBack <b>User provided callbacks</b></dt>
 * <dd>
 *   Much functionality in deal.II under the hood uses external libraries that
 *   operate by calling back into user-provided functions. Examples are
 *   ODE solvers that solve differential equations of the form
 *   @f[
 *     \mathbf x'(t) = \mathbf f(t,\mathbf x(t)),
 *   @f]
 *   where users need to provide a function that, for a given time $t$ and vector
 *   $\mathbf x$ returns the value of the right hand side $\mathbf f(t,\mathbf x)$.
 *   Other examples are solvers for nonlinear systems
 *   @f[
 *     \mathbf F(\mathbf U) = 0,
 *   @f]
 *   where users need to provide functions that for a given vector $\mathbf U$ returns
 *   the vector $\mathbf F(\mathbf U)$ and, in many cases, also the Jacobian
 *   $\nabla \mathbf F(\mathbf U)$ as well as possibly other information about
 *   the problem such as the relative scaling of solution variables within the
 *   vector $\mathbf U$ that can be used to make the problem better conditioned.
 *
 *   These functions are often called "callbacks" because the ODE or nonlinear
 *   solver libraries "call back" into user code. In code written in the C programming
 *   language, these callbacks would often be described by pointers to user
 *   functions handed to the solver library. Since deal.II is written in C++, we
 *   typically instead use
 *   [`std::function`](https://en.cppreference.com/w/cpp/utility/functional/function)
 *   objects. Examples of classes that use this approach are SUNDIALS::KINSOL,
 *   TrilinosWrappers::NOXSolver, PETScWrapper::NonlinearSolver, and
 *   PETScWrappers::TimeStepper. step-77 illustrates how this can be used with
 *   SUNDIALS::KINSOL.
 *
 *   Many of these libraries use a convention that comes from their origin in the
 *   C programming language: User callbacks for the purposes of these libraries
 *   return the data they are asked through one of their arguments, and indicate
 *   success or failure by returning an `int` that needs to be zero if the function
 *   returned successfully, and a nonzero value if the function failed for whatever
 *   reason. (Examples of failures could include if the right hand side function
 *   $\mathbf f(t,\mathbf x)$ contains a square root of one of the $x_i$, but that
 *   $x_i$ is negative; or perhaps if a function requires table lookup for values
 *   that are outside the range that is tabulated.)
 *
 *   The approach to return integer values is decidedly not in the spirit
 *   of C++: If a function cannot complete the work it is asked to do, the C++ way
 *   to deal with this is to throw an exception. As a consequence, for all of the
 *   places where deal.II wraps external libraries and where user codes need to
 *   provide callbacks, we adopt the following conventions that user callbacks
 *   should follow:
 *   - If a function successfully completes its operations, then it simply
 *     returns as expected. If it is supposed to provide specific information,
 *     then it should do so either via a regular return value, or via a non-`const`
 *     argument, as appropriate for the specific callback.
 *   - If a function cannot successfully complete its operation, it should
 *     throw an exception as appropriate. If possible the classes wrapping
 *     the external library (such as the KINSOL, NOX, or PETSc libraries mentioned
 *     above) will then capture the exception, propagate an appropriate failure
 *     code to the underlying library (say, a nonzero error code), which will
 *     typically lead to some clean-up operations inside that external library,
 *     and a return to the wrapper code. There, the originally thrown exception
 *     will then be re-thrown and become visible again in the place where the
 *     wrappers were thrown. In other words, for all practical purposes, it
 *     looks like the exception thrown in the callback had simply propagated
 *     all the way back to user code.
 *   - There are some libraries that allow callbacks to indicate "recoverable
 *     errors". For example, KINSOL can solve nonlinear systems
 *     $\mathbf F(\mathbf U)=0$ and deal with situations where a function
 *     evaluation for $\mathbf F$ is not possible -- for example the case above
 *     one tries to take the square root of a negative value -- but where it
 *     could then try again for a modified $\mathbf U$ for which evaluation
 *     of the square root is possible. (This is often possible in
 *     [line search algorithms](https://en.wikipedia.org/wiki/Line_search)
 *     where using a shorter step length might actually succeed.) In such
 *     cases, a user-provided callback function should throw an exception
 *     of type StandardExceptions::RecoverableUserCallbackError,
 *     which will then internally be
 *     translated into an appropriate code understandable by the underlying
 *     library. It is worthwhile pointing out that a user callback throwing
 *     a "recoverable" exception does not actually guarantee that the
 *     underlying library can actually recover: For example, KINSOL will
 *     eventually give up if for several shorter and shorter step lengths
 *     the residual computation throws a "recoverable" exception, and
 *     the nonlinear solver will then return to user space with an error
 *     that the deal.II wrappers translate into an exception.
 *
 *   The purpose of these conventions is to provide a unified approach to user
 *   callbacks that is independent of how the underlying library likes to
 *   have errors reported. (That is, independent of whether the underlying
 *   library uses nonzero return values, exceptions, or any other mechanism.)
 *   As a consequence, all deal.II classes that require user callbacks try
 *   to follow the convention above.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor workstream_paper <b>%WorkStream paper</b></dt>
 * <dd>The "WorkStream paper" is a paper by B. Turcksin, M. Kronbichler and W. Bangerth
 *   that discusses the design and implementation of WorkStream. WorkStream is, at its
 *   core, a design pattern, i.e., something that is used over and over in finite element
 *   codes and that can, consequently, be implemented generically. In particular, the
 *   paper lays out the motivation for this pattern and then proposes different ways
 *   of implementing it. It also compares the performance of different implementations.
 *
 * The full reference for this paper is as follows:
 * @code{.bib}
@Article{TKB16,
  author =       {Bruno Turcksin and Martin Kronbichler and Wolfgang Bangerth},
  title =        {\textit{WorkStream} -- a design pattern for multicore-enabled finite element computations},
  journal =      {accepted for publication in the ACM Trans. Math. Softw.},
  year =         2016
}
 * @endcode
 * It is available from <a href="https://www.math.colostate.edu/~bangerth/publications.html">https://www.math.colostate.edu/~bangerth/publications.html</a>, also see <a href="https://www.dealii.org/publications.html#details">deal.II publications</a> for details.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossZOrder <b>Z order</b></dt>
 * <dd>
 *  The "Z order" of cells describes an order in which cells are traversed.
 *
 *  By default, if you write a loop over all cells in deal.II, the cells
 *  will be traversed in an order where coarser cells (i.e., cells that were
 *  obtained from
 *  @ref GlossCoarseMesh "coarse mesh" cells with fewer refinement steps) come
 *  before cells that are finer (i.e., cells that were obtained with more refinement
 *  steps). Within each refinement level, cells are traversed in an order
 *  that has something to do with the order in which they were created;
 *  in essence, however, this order is best of thought of as "unspecified":
 *  you will visit each cell on a given refinement level exactly once, in
 *  some order, but you should not make any assumptions about this order.
 *
 *  Because the order in which cells are created factors into the order
 *  of cells, it can happen that the order in which you traverse cells is
 *  different for two identical meshes. For example, think of a 1d (coarse)
 *  mesh with two cells: If you first refine the first of these cells and then
 *  the other, then you will traverse the four cells on refinement level 1
 *  in a different order than if you had first refined the second coarse
 *  cell and then the first coarse cell.
 *
 *  This order is entirely practical for almost all applications because
 *  in most cases, it does not actually matter in which order one traverses
 *  cells. Furthermore, it allows using data structures that lead to
 *  particularly low cache miss frequencies and are therefore efficient
 *  for high performance computing applications.
 *
 *  On the other hand, there are cases where one would want to traverse
 *  cells in a particular, specified and reproducible order that only
 *  depends on the mesh itself, not its creation history or any other
 *  seemingly arbitrary design decisions. The "Z order" is one way
 *  to achieve this goal.
 *
 *  To explain the concept of the Z order, consider the following sequence
 *  of meshes (with each cell numbered using the "level.index" notation,
 *  where "level" is the number of refinements necessary to get from a
 *  @ref GlossCoarseMesh "coarse mesh" cell to a particular cell, and "index" the index of this
 *  cell within a particular refinement level):
 *
 *  @image html simple-mesh-0.png "A coarse mesh"
 *  @image html simple-mesh-1.png "The mesh after one refinement cycle"
 *  @image html simple-mesh-2.png "The mesh after two refinement cycles"
 *  @image html simple-mesh-3.png "The mesh after three refinement cycles"
 *
 *  Note how the cells on level 2 are ordered in the order in which they
 *  were created. (Which is not always the case: if cells had been removed
 *  in between, then newly created cells would have filled in the holes
 *  so created.)
 *
 *  The "natural" order in which deal.II traverses cells would then be
 *  0.0 -> 1.0 -> 1.1 -> 1.2 -> 1.3 -> 2.0 -> 2.1 -> 2.2 -> 2.3 -> 2.4 ->
 *  2.5 -> 2.6 -> 2.7. (If you want to traverse only over the
 *  @ref GlossActive "active cells", then omit all cells from this
 *  list that have children.)
 *  This can be thought of as the "lexicographic"
 *  order on the pairs of numbers "level.index", but because the index
 *  within each level is not well defined, this is not a particularly useful
 *  notion. Alternatively, one can also think of it as one possible breadth-first
 *  traversal of the tree that corresponds to this mesh and that represents
 *  the parent-child relationship between cells:
 *
 *  @image html simple-mesh-tree.png "The tree that corresponds to the mesh after three refinement cycles"
 *
 *  On the other hand, the Z order corresponds to a particular
 *  depth-first traversal of the tree. Namely: start with a cell, and if it
 *  has children then iterate over these cell's children; this rule is
 *  recursively applied as long as a child has children.
 *
 *  For the given mesh above, this yields the following order: 0.0 -> 1.0 -> 2.4
 *  -> 2.5 -> 2.6 -> 2.7 -> 1.1 -> 1.2 -> 1.3 -> 1.4 -> 2.0 -> 2.1 -> 2.2 -> 2.3.
 *  (Again, if you only care about active cells, then remove 0.0, 1.0, and 1.3
 *  from this list.) Because the order of children of a cell is well defined
 *  (as opposed to the order of cells within each level), this "hierarchical"
 *  traversal makes sense and is, in particular, independent of the history
 *  of a triangulation.
 *
 *  In practice, it is easily implemented using a recursive function:
 *  @code
 *    template <int dim>
 *    void visit_cells_hierarchically (const typename Triangulation<dim>::cell_iterator &cell)
 *    {
 *      if (cell->has_children())
 *        for (unsigned int c=0; c<cell->n_children(); ++c)
 *          visit_cells_hierarchically (cell->child(c));
 *      else
 *        {
 *          ... do whatever you wanted to do on each cell ...;
 *        }
 *    }
 *  @endcode
 *  This function is then called as follows:
 *  @code
 *    // loop over all coarse mesh cells
 *    for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin(0);
 *         cell != triangulation.end(); ++cell)
 *      visit_cells_hierarchically (cell);
 *  @endcode
 *
 *  Finally, as an explanation of the term "Z" order: if you draw a line through
 *  all cells in the order in which they appear in this hierarchical fashion,
 *  then it will look like a left-right inverted Z on each refined cell. Indeed,
 *  the curve so defined can be thought of a space-filling curve and is also
 *  sometimes called "Morton ordering", see
 *  https://en.wikipedia.org/wiki/Z-order_curve .
 * </dd>
 *
 *
 *
 * </dl>
 */
