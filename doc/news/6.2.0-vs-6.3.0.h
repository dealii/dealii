// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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
 * @page changes_between_6_2_and_6_3 Changes between Version 6.2 and 6.3

<p>
This is the list of changes made between the deal.II releases listed above.
made to the three sub-libraries <a href="#base">base</a>,
<a href="#lac">lac</a>, and <a href="#deal.II">deal.II</a>, as well as
changes to the <a href="#general">general infrastructure,
documentation, etc</a>.
</p>

<p>
All entries are signed with the names of the author. Regular
contributor's names are abbreviated by WB (Wolfgang Bangerth), GK
(Guido Kanschat), RH (Ralf Hartmann).
</p>


<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>

<ol>
  <li>
  <p>
  Removed: The class TrilinosWrappers::SolverBlock and
  TrilinosWrappers::PreconditionBlock have been deleted for compatibility
  reasons with Trilinos 10.2. These solvers can be exchanged by deal.II's
  own iterative solvers with block matrices and vectors without loss in
  performance.
  <br>
  (Martin Kronbichler 2010/06/24)
  </p>

  <li>
  <p>
  Changed: The SparseMatrix::operator() used to always check that an entry
  exists when accessing it. If it didn't, it would throw an exception that
  could be caught in a user program. Since these accesses are very frequent,
  this check now only happens in debug mode when the program aborts if a
  nonexistent element is accessed. If you access a nonexistent element in
  optimized mode, anything might happen (as in many other functions if
  input arguments do not satisfy their constraints).
  <br>
  (WB 2010/06/04)
  </p>

  <li>
  <p>
  Removed: The interface to PETSc has been simplified to better handle
  incremental changes in PETSc versions and accommodate changes in
  functionality between versions. As a part of this process, the
  deal.II configure script no longer handles PETSc versions
  &lt;2.3.0. Attempting to configure deal.II with PETSc versions that are
  not supported will result in the error message, "Unknown PETSc
  version". The usage of the PETScWrappers are otherwise not affected
  by this change.
  <br>
  (Toby D. Young 2010/03/06)
  </p>

  <li>
  <p>
  Removed: The class TrilinosWrappers::PreconditionStokes has been deleted
  because it did not work properly, and it is too specific to be part of
  the general library. A preconditioner that has the same properties is
  explained in the @ref step_31 step-31 tutorial program.
  <br>
  (Martin Kronbichler 2009/11/23)
  </p>

  <li>
  <p>
  Changed: The class MGSmootherRelaxation now instead of a
  preconditioner takes a relaxation method with the functions
  <code>step</code> and <code>Tstep</code>. These perform a complete
  relaxation step, thus saving an auxiliary vector and computational
  effort for Gauss-Seidel type methods.

  While all relaxation preconditioners have been provided with the new
  two functions, you may have used MGSmootherRelaxation with your own
  preconditioner. In that case, you have two options:
  <ol>
    <li> Configure deal.II with <tt>--enable-mgcompatibility</tt>,
    which restores the old behavior of MGSmootherRelaxation.
    <li> Use MGSmootherPrecondition, which does what
    MGSmootherRelaxation did before.
  </ol>
  <br>
  (GK 2009/08/04)
  </p>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
   <li>
   <p>
   New: The new tutorial program step-45,
   contributed by Markus B&uuml;rg, shows how to implement periodic boundary
   conditions.
   <br>
   (Markus B&uuml;rg, 2010/06/23)
   </p>
   </li>

   <li>
   <p>
   Improved: Exception classes declared locally to another class using
   DeclException0 through DeclException5 were previously shown as local
   classes in the class overview of the library. This was annoying since
   they were not really of interest but made up most of the list.
   The documentation now shows them in one central place in the
   @ref Exceptions module.
   <br>
   (WB 2010/04/26)
   </p>
   </li>

   <li>
   <p>
   Improved: We now compile all files with their full path name on the command
   line. This makes it simpler for tools like debuggers or profilers (e.g.
   valgrind) to find the source files that corresponds to an executable.
   <br>
   (WB 2010/04/15)
   </p>
   </li>

   <li>
   <p>
   Improved:
          Over the last few months, the multigrid implementation has seen
          significant rewrites, with much of the work done by B&auml;rbel
          Janssen. The goal &mdash; now achieved &mdash; was to finally fully
          support multigrid also for continuous finite elements on adaptively
          refined meshes (uniformly refined meshes and discontinuous elements
          have worked for a long time). As part of this process,
	  step-16 has
          been rewritten and now solves the same problem
	  step-6 solves, just
          with a multigrid solver.
   <br>
   (B&auml;rbel Janssen, WB 2010/02/13)
   </p>
   </li>


   <li>
   <p>
   New: The version of <a href="http://www.boost.org/">boost</a>
   included in the <code>contrib/</code> directory has been updated
   to 1.41.0.
   <br>
   (WB 2009/12/10)
   </p>
   </li>

   <li>
   <p>
   New: There is now a new tutorial program, step-35,
   contributed by Abner Salgado-Gonzalez, that implements a solver
   for the Navier-Stokes equations using a decoupled projection
   scheme.
   <br>
   (Abner Salgado-Gonzalez 2009/10/07)
   </p>
   </li>

   <li>
   <p>
   New: A report has been added on the
   <a href="http://www.dealii.org/reports/codimension-one/desimone-heltai-manigrasso.pdf"
   target="body">codimension one</a> capabilities
   of the library (by Antonio DeSimone, Luca Heltai
   and Cataldo Manigrasso, SISSA, Trieste, Italy). It
   explains in detail how to use the
   library for the solution of problems defined on codimension
   one manifolds, such as, for example, %Boundary Element Methods.
   <br>
   (Luca Heltai, 2009/09/23)
   </p>
   </li>

   <li>
   <p>
   New: The configure switch <code>--with-cpu=...</code> now allows the value
   <code>native</code>, indicating that we would like the compiler to figure
   out which CPU we are running on and optimize for it. The resulting
   libraries may not work on other (previous generation) processors, however.
   <br>
   (WB 2009/08/20)
   </p>
   </li>

   <li>
   <p>
   Fixed: step-31 had a bug in the computation of the
   global scaling parameter in the function that evaluates the artificial
   viscosity: we computed
   $c(\mathbf{u},T) =
    c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
    \frac{1}{|\mathrm{diam}(\Omega)|^{\alpha-2}}$
   when it should have been
   $c(\mathbf{u},T) =
    c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
    \ |\mathrm{diam}(\Omega)|^{\alpha-2}$. This didn't matter much in this
   program because $\mathrm{diam}(\Omega)=2^{1/\textrm{dim}}$ and so is close
   to one. It would matter, however, if the domain had been different, as
   it is, for example, in the future step 32.
   <br>
   (WB 2009/08/19)
   </p>
   </li>

   <li>
   <p>
   Changed: When using Trilinos wrapper objects in %parallel through MPI, each
   object now uses a separate and distinct MPI communicator object. This
   ensures that different objects (such as different matrices, or different
   vectors) communicate on separate channels, thereby simplifying debugging
   and possibly the parallelization of programs.
   <br>
   (WB 2009/08/10)
   </p>
   </li>

   <li>
   <p>
   New: There is now a new tutorial program, step-36,
   contributed by Toby D. Young and Wolfgang Bangerth, that demonstrates
   solving eigenvalue problems.
   <br>
   (Toby D. Young, WB 2009/07/29)
   </p>
   </li>

   <li>
   <p>
   Changed: When configuring to use METIS for partitioning meshes in %parallel,
   the METIS header files had to be modified by hand. In addition, with some
   MPI implementations one would get into trouble if <code>mpi.h</code>
   included <code>mpicxx.h</code>. These two problems have now been
   worked around.
   <br>
   (WB 2009/07/06)
   </p>
   </li>

   <li>
   <p>
   New: As a primary means of parallelizing programs, deal.II now uses
   a task-based, rather than thread-based approach, in which one
   uses a high-level description of <i>what</i> needs to be done,
   rather than how these jobs have to be mapped onto threads. We then
   use the <a href="http://www.threadingbuildingblocks.org">Threading
   Building Blocks (TBB) library</a> to schedule tasks onto available
   hardware resources. This new scheme of describing parallism and
   various abstractions to make programming in this framework easier
   are described in great detail in the
   @ref threads "Parallel computing with multiple processors" module.
   In addition, most of the parallelism already used within deal.II
   has been converted to use tasks, rather than threads, and so have
   some of the tutorial programs.
   <br>
   (WB 2009/01/09)
   </p>
   </li>

   <li>
   <p>
   Changed: The support for threading has been completely re-written. In
   particular, the Threads::spawn functions have been deprecated, and
   new functions Threads::new_thread have been introduced.
   Threading is now discussed in a lot of detail in the
   @ref threads "Parallel computing with multiple processors" module.
   <br>
   (WB 2009/01/09)
   </p>
   </li>

   <li>
   <p>
   Changed: Previously, one had to give the two flags
   <code>--enable-multithreading --with-multithreading</code> to
   <code>./configure</code> to enable thread usage throughout the library.
   This has now been simplified: only the flag <code>--enable-threads</code>
   is now necessary. Furthermore, since most current machines have multiple
   cores these days, the default is now to use threads. This can be switched
   off using <code>--disable-threads</code>, however.
   <br>
   (WB 2008/09/29)
   </p>
   </li>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li><p>New: The Timer class can now accumulate and average run times of
  pieces of code across multiple MPI processes.
  <br>
  (Timo Heister 2010/06/07)
  </p></li>

  <li><p>New: The Utilities::System::compute_point_to_point_communication_pattern
  function can be used to compute who wants to send messages to the
  current processor in unstructured point-to-point MPI communications.
  <br>
  (WB 2010/06/07)
  </p></li>

  <li><p>New: The DataOutBase class (and all derived classes such as DataOut,
  MatrixOut, etc) can now produce the XML-based version of the VTK file format
  (the so-called VTU format). Furthermore, the
  DataOutInterfaces::write_pvtu_record function can be used to describe a set
  of %parallel VTU files as part of a single visualization set.
  <br>
  (Scott Miller 2010/06/01)
  </p></li>

  <li><p>Changed: The Function::vector_gradient_list function was previously
  implemented by calling Function::gradient on each point and each component.
  It has been changed to now call Function::vector_gradient on each point
  only, and derived classes should implement this function accordingly.
  <br>
  (WB 2010/02/10)
  </p></li>

  <li><p>Fixed: The file <code>data_out_base.cc</code> could not be compiled
  when Tecplot was available. This should now be fixed.
  <br>
  (WB 2010/02/10)
  </p></li>

  <li><p>Fixed: The PolynomialsBDM, PolynomialsABF and Functions::FlowFunction
  classes had a race condition in multithreaded programs. This now fixed.
  <br>
  (WB 2010/01/27)
  </p></li>

  <li><p>New: The new class IndexSet can represent sets and ranges of indices.
  <br>
  (WB 2009/10/09)
  </p></li>

  <li><p>New: There is now an <code>operator @<@<</code> for the
  TableIndices class.
  <br>
  (WB 2009/09/24)
  </p></li>

  <li>
  <p>
  Fixed: If anything had been put into a LogStream object without flushing
  it with std::endl before the destruction of the log stream, it was lost.
  This is now fixed.
  <br>
  (WB 2009/09/23)
  </p>
  </li>

  <li><p>New: SymmetricTensor::component_to_unrolled_index() and
  SymmetricTensor::unrolled_to_component_indices() allow to convert between the
  indices of an element of a symmetric tensor and the index within an unrolled vector.
  <br>
  (WB 2009/09/23)
  </p></li>

  <li><p>New: classes NamedData and NamedSelection provide an interface to store and
  retrieve data objects with name identifiers.
  <br>
  (GK 2009/09/13)
  </p></li>

  <li>
  <p>
  New: The Utilities::System::job_supports_mpi() can be used to query whether
  the current job runs under MPI or not.
  <br>
  (WB 2009/08/14)
  </p>
  </li>

  <li>
  <p>
  New: The Utilities::Trilinos::comm_self function return an MPI
  communicator that consists only of the current processor.
  <br>
  (WB 2009/08/07)
  </p>
  </li>

  <li>
  <p>
  New: The Utilities::Trilinos::duplicate_communicator function allows to duplicate
  an Epetra_Comm object to get a unique %parallel MPI communicator out of an
  existing one. Utilities::Trilinos::duplicate_map creates a map that has
  the same members as the given template but uses a separate communicator.
  <br>
  (WB 2009/08/06)
  </p>
  </li>

  <li>
  <p>
  New: There is now a specialization Tensor<0,dim> of tensors of rank 0. Since rank-0
  tensors are scalars, this class essentially acts like a scalar, but it allows for
  some neat template tricks that involve tensors of arbitrary rank.
  <br>
  (WB 2009/07/15)
  </p>
  </li>

  <li>
  <p>
  New: The GeometryInfo::alternating_form_at_vertices can be used
  to investigate the degree of distortion of cells.
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  New: The GeometryInfo::d_linear_shape_function and
  GeometryInfo::d_linear_shape_function_gradient functions can be used
  to represent the $d$-linear shape functions that are frequently
  used to map the reference cell to real cells (though the
  Mapping class hierarchy also allows to use higher order mappings).
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  New: The determinant() function is now implemented for rank-2 Tensor
  arguments of all sizes. The implementation is not efficient for very large
  matrix sizes, however.
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  Improved: The QGaussLobatto::gamma function now returns a long double
  instead of an unsigned int, otherwise we will get an overflow and thus
  meaningless weights for higher QGaussLobatto quadrature rules.
  <br>
  (Tobias Leicht, RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: The new function Utilities::duplicate_communicator can be used
  to duplicate an MPI communicator to produce a unique duplicate.
  <br>
  (WB 2009/05/13)
  </p>
  </li>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li><p>New: The ConstraintMatrix class can now handle storing only
  a subset of all constraints, for example only for degrees of
  freedom that are relevant for the subdomain that is owned by one
  process in an MPI universe.
  <br>
  (Timo Heister, Martin Kronbichler 2010/06/07)
  </p></li>

  <li><p>New: The PETScWrappers::MPI::Vector and TrilinosWrappers::MPI::Vector
  classes can now handle ghost elements, i.e. elements that are not
  owned by the current processor but are available for reading
  anyway. The simplest form of ghosting would be to simply import
  an entire vector to local memory, but the new function allow to
  select the elements we need to support the case of computations
  where importing all elements of even a single vector would
  exceed available memory.
  <br>
  (Timo Heister 2010/06/07)
  </p></li>

  <li>
    <p>
    New: A class SparseDirectMumps that provides an interface to
    the MUltifrontal Massively Parallel sparse direct %Solver (MUMPS).
    </p>
  <br>
  (Markus Buerg 2010/05/10)
  </li>

  <li>
    <p>
    Fixed: BlockSparsityPattern::copy_from accidentally only copied
    n_block_rows times n_block_rows blocks, instead of n_block_rows
    times n_block_cols. This is now fixed.
    </p>
  <br>
  (WB 2010/01/06)
  </li>

  <li>
    <p>
    Fixed: SparsityPattern::copy_from crashed whenever a compressed sparsity
    pattern was copied that had either zero rows or zero columns. This is now
    fixed.
    </p>
  <br>
  (WB 2010/01/06)
  </li>

  <li>
    <p>
    New: The function Householder::least_squares can handle BlockVectors as
    well now. Note that in one place we still have to copy to a Vector to use
    the function backward from FullMatrix.
    </p>
  <br>
  (B&auml;rbel Janssen 2010/01/05)
  </li>

  <li>
    <p>
    New: There are now two ConstraintMatrix::add_lines functions that can
    add several constraints at once.
    </p>
  <br>
  (WB 2010/01/05)
  </li>

  <li>
    <p>
    Improved: The Vector class has been equipped with an improved way to
    calculate sums in inner products and norms. This reduces the accumulation
    of round-off errors. Especially the solution with float vectors should
    profit from the new implementation.
    </p>
  <br>
  (Martin Kronbichler 2009/11/05)
  </li>

  <li>
    <p>
    Improved: The ConstraintMatrix class now uses a cache for random access to
    the constraint lines. This considerably increases performance of the
    *_local_to_global functions, where such an access pattern is usual. Moreover,
    the ConstraintMatrix class has now a function get_dof_values that can import
    data from a global vector to a cell vector with respecting the constraints.
    </p>
  <br>
  (Martin Kronbichler 2009/09/30)
  </li>

  <li>
    <p>
    Fixed: SparsityTools::reorder_Cuthill_McKee would produce an error if the
    input graph had disconnected components. This is now fixed.
    </p>
  <br>
  (WB 2009/09/25)
  </li>

  <li>
    <p>
    Fixed: When using the TrilinosWrappers::MPI::Vector::reinit() function with a %parallel
    vector, and if the vector initialized and the vector given had a local range on one of
    the processors that exactly matched, the program would freeze if the local ranges on
    the other processors did not also match exactly. This is now fixed.
    </p>
  <br>
  (WB 2009/09/02)
  </li>

  <li><p> Improved: BlockVector and several of the block sparsity
  patterns can now be initialized with BlockIndices
  objects. Therefore, if an application needs such an object, it does
  not have to store a vector of block sizes separately.
  <br>
  (GK 2009/08/26)
  </p>
  </li>

  <li><p> New: The class PreconditionChebyshev implements a
  preconditioner based on Chebyshev polynomials. It is based on matrix-vector
  products together with some vector updates.
  <br>
  (Martin Kronbichler 2009/08/25)
  </p>
  </li>

  <li>
  <p>
  Fixed: Crash or strange behaviour (wrong matrix entries written) in
  PETScWrappers::MPI::BlockSparseMatrix when adding or setting elements
  through any of the set() and add() routines. This happened when different
  CPUs access different blocks at the start of assembly or when switching
  between adding and setting.
  <br>
  (Timo Heister 2009/08/05)
  </p>
  </li>

  <li> <p>New: The relaxation preconditioners PreconditionJacobi, PreconditionSOR and
  PreconditionSSOR, as well as their blocked versions PreconditionBlockJacobi,
  PreconditionBlockSOR and PreconditionBlockSSOR now have functions <code>step</code>
  and <tt>Tstep</tt> performing one complete step of these methods.
  <br>
  (GK 2009/08/04)
  </p>
  </li>

  <li>
  <p>
  New: There are new functions FullMatrix::cholesky and
  FullMatrix::outer_product.  FullMatrix::cholesky finds the Cholesky
  decomposition of a matrix in lower triangular form.
  FullMatrix::outer_product calculates <tt>*this</tt> $= VW^T$ where $V$
  and $W$ are vectors.
  <br>
  (Jean Marie Linhart 2009/07/27)
  </p>
  </li>

  <li>
  <p>
  Fixed: The TrilinosWrappers::MPI::BlockVector class declares an assignment
  operator from the non-Trilinos BlockVector class but it could not be
  compiled due to an oversight. This is now fixed.
  <br>
  (WB 2009/06/29)
  </p>
  </li>

  <li>
  <p>
  New: Based on work by Francisco Alvaro, the existing SLEPcWrappers now
  have a handle on the generalized eigenvalue problem where B=I.
  <br>
  (Toby D. Young 2009/06/25)
  </p>
  </li>

  <li>
  <p>
  New: Based on work with Francisco Alvaro and Jose
  E. Roman, the new SLEPcWrappers give a handle on some of the
  features of SLEPc (Scalable Library for Eigenvalue Problem
  Computations): (1) The SLEPcWrappers::SolverBase class can be used
  for specifying an eigenvalue problem, either in standard or
  generalized form, on serial or %parallel architectures with support
  for a few solver types; and (2) The
  SLEPcWrappers::TransformationBase class encapsulates a variety of
  spectral transformations providing some functionality required for
  acceleration techniques based on the transformation of the spectrum.
  <br>
  (Toby D. Young 2009/06/25)
  </p>
  </li>

  <li>
  <p>
  Fixed: The TrilinosWrappers::BlockVector class declares an assignment
  operator from the non-Trilinos BlockVector class but it wasn't implemented.
  This is now fixed.
  <br>
  (WB 2009/06/24)
  </p>
  </li>

  <li>
  <p>
  New: The SparseMatrix class has now a function SparseMatrix::mmult that
  can multiply two sparse matrices with each other.
  <br>
  (Martin Kronbichler 2009/05/04)
  </p>
  </li>
</ol>


<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>

<li> <p> New: The namespace MeshWorker contains a generic MeshWorker::loop() over all cells and faces as well as auxiliary classes, which allow to program integrals over mesh cells and faces in a very generic and local way. In particular, an application programmer will not have to distinguish between regular faces and faces with hanging nodes anymore. Two tutorial programs (step-12 and step-39) highlight the functionality of this framework.
<br>
(GK 2010/06/24)
</p></li>

  <li>
  <p>New: The FE_Q_Hierarchical class now has functions
  FE_Q_Hierarchical::hp_constraints_are_implemented and
  FE_Q_Hierarchical::hp_vertex_dof_identities.
  <br>
  (Markus B&uuml;rg 2010/06/08)
  </p></li>

  <li>
  <p>New: The FEValuesViews::Vector class now has functions
  FEValuesViews::Vector::curl and FEValuesViews::Vector::get_function_curls.
  <br>
  (Markus B&uuml;rg 2010/05/13)
  </p></li>

  <li>
  <p>New: TriaAccessor::extent_in_direction() returns the length
  of an object in a given direction.
  <br>
  (James Avery 2010/05/10)
  </p></li>

  <li>
  <p>New: There is a new function DoFTools::extract_dofs_with_support_on_boundary().
  <br>
  (WB 2010/05/07)
  </p></li>

  <li>
  <p>Fixed: FE_DGQ::has_support_on_face() returned the wrong value in 1d if the
  polynomial degree of the finite element equals zero (i.e. for piecewise
  constants) where the lone shape function is nonzero on all faces. This is now
  fixed.
  <br>
  (WB 2010/05/07)
  </p></li>

  <li>
  <p>Fixed: VectorTools::interpolate_boundary_values inadvertently produces
  an exception when used with hp::DoFHandler objects in 1d. This is now fixed.
  <br>
  (WB 2010/05/04)
  </p></li>

  <li>
  <p>Fixed: The GridIn::read_msh function got into trouble if the mesh file
  had a section that listed physical names for variables. This is now fixed.
  <br>
  (WB 2010/05/03)
  </p></li>

  <li>
  <p> Improved: DoFHandler iterators now can be assigned from a Triangulation iterator
  after the dof handler was set once.
  <br>
  (GK 2010/03/25)
  </p></li>

  <li>
  <p>
  New: The function DoFRenumbering::downstream has now an additional bool
  argument. If enabled, the downstream comparison is performed on a DoF
  basis, as opposed to the cell-based comparison that is used with a false
  argument.
  <br>
  (Martin Kronbichler 2010/03/19)
  </p>
  </li>

  <li><p> Improved: DoFHandler now has a BlockInfo object, automatically
  updated after DoFHandler::distribute_dofs() and accessible by
  DoFHandler::block_info(). This object can be used to initialize block
  vectors and obliterates the necessity to count dofs per block (or dofs
  per component in legacy code) in application programs.
  <br>
  (GK 2010/03/18)
  </p>
  </li>

  <li>
  <p>
  New: The class MGTransferSelect is prepared for use on adaptively refined meshes.
  <br>
  (B&auml;rbel Janssen 2010/02/05)
  </p>
  </li>

  <li>
  <p>
  New: The function Cuthill_McKee in namespace DoFRenumbering is now also compiled for
  MGDoFHandler as well as the make_sparsity_pattern functions in DoFTools.
  <br>
  (B&auml;rbel Janssen 2010/01/08)
  </p>
  </li>

  <li>
  <p>
  New: Constructor of FESystem now also exists for four base elements.
  <br>
  (Thomas Wick 2010/01/08)
  </p>
  </li>

  <li>
  <p>
  New: The functions in namespace DoFRenumbering::boost are now also compiled for
  hp::DoFHandler arguments.
  <br>
  (WB 2009/12/14)
  </p>
  </li>

  <li>
  <p>
  Fixed: The functions DoFTools::count_dofs_per_component and DoFTools::extract_dofs
  produced wrong results for elements that consist of two or more nested FESystems.
  Moreoever, a bug in DoFTools::extract_constant_modes has been corrected and
  DoFTools::distribute_cell_to_dof_vector now works according to the documentation,
  namely leaving unselected components in the result vector unchanged, instead of
  setting these to zero as was done before.
  <br>
  (Martin Kronbichler 2009/12/14)
  </p>
  </li>

  <li>
  <p>
  Fixed: The function Triangulation::n_levels() accidentally turned out to be
  quite expensive, in particular if the mesh has been coarsened significantly
  in the past. Since Triangulation::n_active_cells() calls Triangulation::n_levels()
  repeatedly, by consequence it also is surprisingly expensive, which is
  particularly annoying since this is a frequently called function. This has
  now been fixed by caching the result of both functions.
  <br>
  (Wolfgang Bangerth 2009/11/22)
  </p>
  </li>

  <li>
  <p>
  Fixed: The function DataOut class got confused when its DataOut::first_cell() and
  DataOut::next_cell() were overloaded and cell data was given. This is now fixed.
  <br>
  (Wolfgang Bangerth 2009/11/05)
  </p>
  </li>

  <li>
  <p>
  New: The function DoFRenumbering::component_wise is now also implemented
  for arguments of type hp::DoFHandler.
  <br>
  (Markus B&uuml;rg 2009/10/01)
  </p>
  </li>

  <li><p>New: The class BlockInfo stores and computes all BlockIndices objects related to
  DoFHandler and MGDoFHandler based on FESystem.
  <br>
  (GK 2009/09/13)
  </p>
  </li>

  <li><p>Improved: MGDoFHandler now also has a typedef for MGDoFHandler::Container
  <br>
  (GK 2009/09/13)
  </p>
  </li>

  <li><p> FETools::compute_block_renumbering() can nor return block sizes instead of
  start indices.
  <br>
  (GK 2009/08/26)
  </p>
  </li>

  <li>
  <p>
  New: The function GridGenerator::truncated_cone() and the class ConeBoundary
  can now be used to describe conical objects.
  <br>
  (Markus B&uuml;rg 2009/08/17)
  </p>
  </li>

  <li>
  <p>
  New: Instead of asking face or edge iterators for their boundary indicator
  using TriaAccessor::boundary_indicator() and then the triangulation for
  the boundary object using Triangulation::get_boundary(), you can now directly
  ask the iterator for the boundary object using TriaAccessor::get_boundary().
  <br>
  (WB 2009/07/31)
  </p>
  </li>

  <li>
  <p>
  Fixed: The CellAccessor::recursively_set_material_id function did not
  set the material id for all children, but only for the first two, which
  is obviously a bug. This should now be fixed.
  <br>
  (WB 2009/07/14)
  </p>
  </li>

  <li>
  <p>
  Fixed: The GridIn class sometimes had problems with input files that had
  whitespace at the end of lines. This should now be fixed.
  <br>
  (WB 2009/07/10)
  </p>
  </li>

  <li>
  <p>
  New: Previously, the Triangulation::create_triangulation
  function silently accepted input meshes with inverted cells
  (i.e. cells with a zero or negative determinant of the Jacobian of
  the mapping from the reference cell to the real cell). This can been
  changed now: By passing the appropriate flag to the constructor of
  the Triangulation class, the Triangulation::create_triangulation
  function checks whether cells are distorted or
  inverted, and may throw an exception containing a list of cells
  for which this is the case. If you know that this is harmless, for
  example if you have cells with collapsed vertices in your mesh but
  you do not intend to integrate on them, then you can catch and
  ignore this message. In all other cases, the output of your
  computations are likely to be wrong anyway.
  <br>
  The same is true for the Triangulation::execute_coarsening_and_refinement
  function: if it creates cells that are distorted, it throws a list of cells
  whose children are distorted.
  <br>
  The whole issue is described in some detail in the entry on
  @ref GlossDistorted "distorted cells" in the glossary.
  <br>
  (WB 2009/06/29)
  </p>
  </li>


  <li>
  <p>
  New: The new hp::DoFHandler::set_active_fe_indices function allows
  to distribute all active FE indices at once based on a given
  vector. This might be useful if this information is stored
  somewhere and has to be reconstructed or else if two DoFHandler
  objects with the same FE index distribution should be created.
  There is now also a corresponding
  hp::DoFHandler::get_active_fe_indices function.
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>

  <li>
  <p>
  Fixed: The projection of quadrature points to subfaces in
  MappingQ in case of 3d anisotropic refinement did not respect
  non-standard face orientation/flip/rotation cases. This
  has now been fixed.
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>

  <li>
  <p>
  New: The new Triangulation::n_raw_faces() function forwards
  to Triangulation::n_raw_lines() in 2d and
  Triangulation::n_raw_quads() in 3d.
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>

  <li>
  <p>
  New: There is now a new DataOutFaces::build_patches function which
  takes a Mapping argument. For higher order mappings this allows to
  represent curved boundaries by using more subdivisions. This function
  is also useful in the context of MappingQ1Eulerian.
  <br>
  (Tobias Leicht, RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: For empty triangulations the new Triangulation::set_mesh_smoothing
  function allows to override the MeshSmoothing given to the constructor.
  <br>
  (RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: The new function TriaAccessor::is_translation_of computes
  whether a cell, face, or edge is a translation of another.
  <br>
  (Martin Kronbichler, WB 2009/05/19)
  </p>
  </li>

  <li>
  <p>
  New: The DoFTools::make_sparsity_pattern functions have acquired a
  new paramater <code>subdomain_id</code>. If a value other than the
  default value is passed for it, the function only assembles the
  sparsity pattern on those cells that have the given subdomain id.
  This is useful, for example, in conjunction with the
  TrilinosWrappers::SparsityPattern class that can hold a sparsity
  pattern distributed across several MPI processes; in that case, it is
  not necessary that each process builds the entire sparsity pattern.
  <br>
  (WB 2009/04/29)
  </p>
  </li>

   <li>
   <p>
   Fixed: The DoFRenumbering::component_wise function for MGDoFHandler objects
   did a few things in %parallel that weren't thread-safe. This is now fixed.
   <br>
   (WB, 2009/01/20)
   </p>

   <li>
   <p>
   Changed: The two DataOut::build_patches, DataOutFaces::build_patches, and
   DataOutRotation::build_patches functions have lost the argument
   that indicated the number of threads with which they should build the
   intermediate representation. This is something that now happens
   transparently in the background and doesn't need caller input any more.
   <br>
   (WB 2008/12/16)
   </p>

   <li>
   <p>
   Changed: The KellyErrorEstimator::estimate functions had a parameter
   that indicates the number of threads to be used in the computation.
   This parameter continues to exist for compatibility, but is now ignored.
   Rather, the number of threads is determined automatically by scheduling
   the requested computations on available compute resources.
   <br>
   (WB, 2008/12/29)
   </p>

   <li>
   <p>
   New: The new function internal::hp::FEValuesBase::get_fe_collection function
   allows to query the finite element collection currently in used in an hp::FEValues,
   hp::FEFaceValues, or hp::FESubfaceValues object.
   <br>
   (WB 2008/09/30)
   </p>

   <li>
   <p>
   New: The new function FEValuesBase::get_update_flags allows to query
   the update flags that are currently set on an FEValues, FEFaceValues, or
   FESubfaceValues object.
   <br>
   (WB 2008/09/29)
   </p>
</ol>


*/
