// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2024 by the deal.II authors
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
 * @page changes_between_8_2_1_and_8_3 Changes between Version 8.2.1 and 8.3

<p>
This is the list of changes made between the release of deal.II version
8.2.1 and that of 8.3.0. All entries are signed with the names of the
author.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="821-830-incompatible"></a>
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

  <li>Changed: GridTools::distort_random is now deterministic (gives the same
  distorted mesh when called with the same input Triangulation).
  <br>
  (Timo Heister, 2015/07/09)
  </li>

  <li>Changed: MatrixCreator::create_mass_matrix() took matrix and vector
  objects where the scalar type of the matrix was a template argument
  but the scalar type of the vector was always <code>double</code>. This
  has been changed so that the two always need to match.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2015/06/17)
  </li>

  <li>Changed: Functions such as FEValuesBase::get_function_values() that
  extracted the values of functions at the quadrature points of a cell
  implicitly always assumed that <code>double</code> is a reasonable
  type to store the result in. This, however, is not true if the solution
  vector's underlying scalar type was, for example,
  <code>std::complex@<double@></code>. All of the functions of
  FEValuesBase as well as of the various FEValuesViews class that extract
  values from solution vectors have been changed so that they now return
  their results in vectors that use the same underlying scalar type
  (float, double, or std::complex) as was used in the solution vector.
  <br>
  Most user codes will be entirely unaffected by this because they simply
  use the default vector types which all store their data as doubles.
  You may have to adjust your code, though, if you use non-standard
  types such as Vector@<float@> for solution vectors, or use
  complex-valued data types. This includes compiling PETSc with
  complex scalars.
  <br>
  (Denis Davydov, 2015/05/08)
  </li>

  <li>Removed: The generic, templated vmult, Tvmult, etc. -interfaces of
  LAPACKFullMatrix - they were never implemented.
  <br>
  (Matthias Maier, 2015/04/10)
  </li>

  <li>Removed: SparseDirectUMFPACK::vmult_add and
  SparseDirectUMFPACK::Tvmult_add - they were never implemented.
  <br>
  (Matthias Maier, 2015/04/10)
  </li>

  <li>Removed: The class NamedData has been removed after it had been
  superseded by AnyData a while ago. This affects the use of classes
  in Algorithms and MeshWorker
  <br>
  (Guido Kanschat, 2015/04/02)
  </li>

  <li> Removed: The CMake configuration does not use the variable
  <code>DEAL_II_CMAKE_MACROS_RELDIR</code> any more. Instead, the fixed
  location <code>\${DEAL_II_SHARE_RELDIR}/macros</code> is used
  unconditionally
  <br>
  (Matthias Maier, 2015/03/26)
  </li>

  <li> Changed: The TrilinosWrappers::SparseMatrix::clear_row() function used
  to call TrilinosWrappers::SparseMatrix::compress() before doing its work,
  but this is neither efficient nor safe. You will now have to do this
  yourself after assembling a matrix and before clearing rows.
  <br>
  The changes to the function above also affect the
  MatrixTools::apply_boundary_values() variants that operate on Trilinos
  matrices.
  <br>
  (Wolfgang Bangerth, 2015/03/09)
  </li>

  <li> Changed: Implicit conversion from Tensor@<1,dim@> to Point@<dim@> was
  previously possible. This has now been prohibited (but you can still
  do the conversion with an explicit cast) as such conversions are
  likely incorrect uses of class Point (which should represent only
  points in space, i.e., vectors anchored at the origin) whereas Tensor
  should be used for vectors anchored elsewhere (such as normal vectors,
  directions, differences between points, etc). The difference in
  usage between Point and Tensor have now been clarified in the documentation
  of class Point.
  <br>
  (Wolfgang Bangerth, 2015/01/12)
  </li>

  <li> Changed: The project configuration no longer exports
  <code>[...]/include/deal.II</code>. Thus it is now mandatory to prefix
  all includes of deal.II headers with <code>deal.II/</code>, i.e.
  <code>\#include &lt;deal.II/[...]&gt;</code>.
  <br>
  (Matthias Maier, 2015/01/19)
  </li>

  <li> Changed: ParameterHandler::leave_subsection() no longer returns a bool
  indicating if there was a subsection to leave. This never worked in the
  first place, because an exception was thrown.
  <br>
  (Timo Heister, 2015/01/19)
  </li>

  <li> Changed: Make.global_options was completely redesigned. It still
  contains Makefile sourcable information, but they now closely mimic the
  declarative style of deal.IIConfig.cmake. Thus, projects that still use
  Makefiles that source Make.global_options have to be ported to the new
  layout.
  <br>
  (Matthias Maier, 2015/01/13)
  </li>

  <li> Removed: The Component <code>compat_files</code> was removed entirely. deal.II
  now always configures and installs with a somewhat FSHS compliant
  directory structure. Further, the ancient make_dependencies binary was
  removed. Either migrate your project to CMake, or port your build system
  to the new (incompatible) Make.global_options found at
  <code>\${DEAL_II_SHARE_RELDIR}</code>.
  <br>
  (Matthias Maier, 2015/01/13)
  </li>

  <li> Changed: The two-argument call to the MPI_InitFinalize constructor
  used to imply that the user wanted only one thread per MPI process. This
  has been changed and now means that every processor core on the system is
  used. If you run as many MPI processes as there are processor cores, then
  this means one thread per MPI process, as before. On the other hand, if you
  start fewer MPI processes than there are cores, then your program will now
  be allowed to use more than one thread. You can get the old behavior by
  setting the third (optional) argument to one.
  <br>
  (Wolfgang Bangerth, 2015/01/14)
  </li>

  <li> Removed: TrilinosWrappers::SparseMatrix copy constructor got removed
  to be in line with PETSc and dealii::SparseMatrix. You can use reinit()
  and copy_from().
  <br>
  (Timo Heister, 2015/01/12)
  </li>

  <li> Removed: The following compatibility definitions were removed from
  <code>include/deal.II/base/config.h.in</code> (replacement in brackets):
  - DEAL_II_CAN_USE_CXX11 (new: DEAL_II_WITH_CXX11)
  - DEAL_II_CAN_USE_CXX1X (new: DEAL_II_WITH_CXX11)
  - DEAL_II_COMPILER_SUPPORTS_MPI (new: DEAL_II_WITH_MPI)
  - DEAL_II_MAJOR (new: DEAL_II_VERSION_MAJOR)
  - DEAL_II_MINOR (new: DEAL_II_VERSION_MINOR)
  - DEAL_II_USE_ARPACK (new: DEAL_II_WITH_ARPACK)
  - DEAL_II_USE_CXX11 (new: DEAL_II_WITH_CXX11)
  - DEAL_II_USE_METIS (new: DEAL_II_WITH_METIS)
  - DEAL_II_USE_MT (new: DEAL_II_WITH_THREADS)
  - DEAL_II_USE_P4EST (new: DEAL_II_WITH_P4EST)
  - DEAL_II_USE_PETSC (new: DEAL_II_WITH_PETSC)
  - DEAL_II_USE_SLEPC (new: DEAL_II_WITH_SLEPC)
  - DEAL_II_USE_TRILINOS (new: DEAL_II_WITH_TRILINOS)
  <br>
  (Matthias Maier, 2015/01/12)
  </li>

  <li> Removed: The direct Mumps interface through
  <code>SparseDirectMUMPS</code> has been removed. The MUMPS solver is
  still available through the Trilinos or PETSc interfaces. Alternatively,
  there is <code>SparseDirectUMFPACK</code>, which has a similar interface.
  <br>
  (Matthias Maier, 2015/01/11)
  </li>

  <li> Removed: This release removes a number of functions that have long
  been deprecated and that were previously already marked as
  deprecated (i.e., they would have yielded warnings by the compiler whenever
  you tried to use them). In almost all cases, there is a function with same
  name but different argument list that should be used instead.
  Specifically, the removed functions and classes are:
  <br>
  <em>With headers in <code>deal.II/base/</code>:</em>
  - ThreadManagement::spawn.
  - Threads::ThreadCondition and Threads::ThreadMutex.
  - DataOutBase::create_xdmf_entry with 3 arguments.
  - DataOutBase::write_hdf5_parallel with 2 arguments.
  - The versions of FunctionParser::initialize that took a
    <code>use_degrees</code> or <code>constants</code> argument.
    The implementation as it is now no longer supports either of
    these two concepts (since we switched from the FunctionParser
    library to the muparser library after the deal.II 8.1 release).
  - GridOutFlags::XFig::level_color.
  - class BlockList.
  - The MPI support functions in namespace Utilities and Utilities::System.
  - Deprecated members of namespace types.
  - Namespace deal_II_numbers.
  - MultithreadInfo::n_default_threads.
  - Table::data.

  <br>
  <em>With headers in <code>deal.II/lac/</code>:</em>
  - The deprecated constructors of SparseMIC,
    SparseILU, and SparseLUDecomposition.
  - SparseMIC::decompose and SparseILU::decompose.
  - SparseMIC::reinit and SparseLUDecomposition::reinit.
  - SparseILU::apply_decomposition.
  - SparseLUDecomposition::decompose and SparseLUDecomposition::is_decomposed.
  - The compress() functions without argument in the various vector
    classes. You should use the versions with a VectorOperation
    argument instead.
  - Vector::scale.
  - TrilinosWrappers::*Vector*::%compress with an Epetra_CombineMode
    argument.
  - SparsityPattern and ChunkSparsityPattern functions that take an
    <code>optimize_diagonal</code> argument.
  - SparsityPattern::partition.
  - SparsityPattern::get_rowstart_indices and
    SparsityPattern::get_column_numbers.
  - SparsityPattern::row_iterator and corresponding row_begin() and row_end()
    functions.
  - CompressedSparsityPattern::row_iterator and corresponding row_begin()
    and row_end() functions.
  - The typedef CompressedBlockSparsityPattern.
  - The deprecated constructors of SparsityPattern iterator classes.
  - The deprecated variants of DoFTools::make_periodicity_constraints.
  - BlockMatrixArray and BlockTrianglePreconditioner functions that
    take an explicit VectorMemory object.
  - The SolverSelector constructor that takes a VectorMemory argument.
  - The version of parallel::distributed::Vector::compress_finish
    function that takes a boolean as argument.
  - The version of BlockVector::scale and
    parallel::distributed::Vector::scale,
    parallel::distributed::BlockVector::scale
    function that takes a scalar as argument.
  - PreconditionBlock::size.
  - Classes PreconditionedMatrix and PreconditionLACSolver.
  - PETScWrappers::VectorBase::update_ghost_values.
  - PETScWrappers::MPI::Vector constructors and reinit variants.
  - SparseMatrixIterators::Accessor and SparseMatrixIterators::Iterator
    constructors.
  - SparseMatrix::raw_entry and SparseMatrix::global_entry.
  - The ConstraintMatrix functions that transform a matrix, vector, or
    linear system into a smaller by not just setting the corresponding
    rows and columns to zero, but actually shrinking the size of the
    linear system.

  <br>
  <em>With headers in <code>deal.II/deal.II/</code>:</em>
  - GridGenerator::laplace_transformation.
  - The version of GridGenerator::parallelogram where the corners are given
    as a rank-2 tensor rather than as an array of points.
  - GridTools::create_union_triangulation.
  - GridTools::extract_boundary_mesh.
  - Triangulation::distort_random.
  - Triangulation::clear_user_pointers.
  - The refinement listener concept of the Triangulation class. This
    approach to getting notified about what happens to triangulations
    has been superseded by the signals defined by the triangulation
    class.

  <br>
  <em>With headers in <code>deal.II/fe/</code>:</em>
  - In FEValues and related classes, the functions that contain the
    term <code>2nd_derivatives</code> were removed in favor of those
    with names containing <code>hessian</code>. Similarly, functions
    with names including <code>function_grads</code> were removed in
    favor of those called <code>function_gradients</code>. Finally,
    the <code>cell_normal_vector</code> functions were replaced by
    <code>normal_vector</code> ones. In all cases, the new functions
    have been around for a while.
  - Mapping::transform_covariant and Mapping::transform_contravariant.

  <br>
  <em>With headers in <code>deal.II/dofs/</code>:</em>
  - DoFRenumbering::downstream_dg.
  - DoFTools::count_dofs_per_component.
  - DoFTools::make_sparsity_pattern with a vector-of-vector mask.

  <br>
  <em>With headers in <code>deal.II/multigrid/</code>:</em>
  - The constructors of classes MGSmoother, MGSmootherRelaxation and
    MGSmootherPrecondition that take a VectorMemory object.
  - MGLevelObject::get_minlevel and MGLevelObject::get_maxlevel.
  - MGConstrainedDoFs::non_refinement_edge_index
  - MGConstrainedDoFs::at_refinement_edge_boundary
  - MGTools::count_dofs_per_component.
  - MGTools::apply_boundary_values.
  - MGTools::extract_inner_interface_dofs.
  - Class MGMatrix.
  - Multigrid::vmult and friends.

  <br>
  <em>With headers in <code>deal.II/matrix_free/</code>:</em>
  - Classes FEEvaluationDGP, FEEvaluationGeneral and FEEvaluationGL.

  <br>
  <em>With headers in <code>deal.II/mesh_worker/</code>:</em>
  - Deprecated variants of MeshWorker::loop and MeshWorker::integration_loop.

  <br>
  <em>With headers in <code>deal.II/algorithm/</code>:</em>
  - Algorithms::ThetaTimestepping::operator().
  - Algorithms::ThetaTimestepping::initialize.
  - Algorithms::Newton::initialize.

  <br>
  <em>With headers in <code>deal.II/numerics/</code>:</em>
  - TimeDependent::end_sweep (with an argument).
  - PointValueHistory::mark_locations.
  - The DataPostprocessor::compute_derived_quantities_scalar and
    DataPostprocessor::compute_derived_quantities_vector functions without
    evaluation points. If you have
    data postprocessor classes implemented in your program that overload these
    functions, you will have to change it in a way that they overload the
    functions of same name but with the evaluation point argument instead.
  <br>
  This release also removes the deprecated class MGDoFHandler. The
  functionality of this class had previously been incorporated into
  the DoFHandler class. Unlike the changes above, if you were still
  using this class, you will need to do the following changes to
  your code:
  - Where you called <code>mg_dof_handler.distribute_dofs()</code>
    you now also need to explicitly call
    <code>mg_dof_handler.distribute_mg_dofs()</code>.
  - If you called <code>mg_dof_handler.begin(level)</code>, you
    will now have to write this as
    <code>mg_dof_handler.begin_mg(level)</code> to make clear that
    you are not just interested in an iterator to a cell on a given
    level, but in fact to a cell that can access the degrees of
    freedom on a particular level of a multigrid hierarchy.
  - The type previously referred to as
    <code>MGDoFHandler::cell_iterator</code> now corresponds to
    <code>MGDoFHandler::level_cell_iterator</code>.
  - Where you previously called DoFRenumbering::component_wise
    for the entire MGDoFHandler object, you now need to call
    this function for the DoFHandler object, and then call the
    same function with the <code>level</code> argument for each
    of the levels of the triangulation individually.
  <br>
  (Wolfgang Bangerth, 2014/12/29-2015/01/22)
  </li>

  <li> Removed: The config.h file no longer exports HAVE_* definitions.
  Those are either entirely removed (for the blas/lapack symbols) or
  renamed to DEAL_II_HAVE_*. This change is done in order to avoid clashes
  with external projects also exporting HAVE_* definitions in their header
  files.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="821-830-general"></a>
<h3>General</h3>


<ol>
  <li> New: IndexSet now can be constructed using Epetra_Map.
  <br>
  (Luca Heltai, 2015/07/25)
  </li>

  <li> New: Added the class Functions::Polynomial for representation of polynomials.
  The new class is derived from the Function class.
  <br>
  (Angel Rodriguez, 2015/07/01)
  </li>

  <li> New: deal.II now supports compilation in C++14 mode, which may be
  enabled with the CMake option <code>DEAL_II_WITH_CXX14</code>.
  <br>
  (David Wells, 2015/06/21)
  </li>

  <li> New: Implement a modified version of the Kelly error estimator, which
  effectively provides the boundary residual term for the hp-FEM error estimators.
  <br>
  (Denis Davydov, 2015/06/17)
  </li>

  <li> New: DerivativeForm() now takes an additional optional
  template argument specifying the type, similarly to Tensor() classes.
  <br>
  (Luca Heltai, 2015/05/16)
  </li>

  <li> New: Utilities::MPI::min() functions.
  <br>
  (Timo Heister, 2015/05/12)
  </li>

  <li> New: A PackagedOperation class that stores computation expressions.
  The primary purpose of this class is to provide syntactic sugar for vector
  operations (vector space addition, scalar multiplication, application of a
  linear operator) while avoiding intermediate storage.
  <br>
  (Matthias Maier, 2015/05/10)
  </li>

  <li> Fixed: Utilities::generate_normal_random_number() will now
  produce a deterministic sequence of numbers on every thread.
  Furthermore, it will also work on systems that do not have the
  <code>rand_r</code> library function, such as Cygwin.
  <br>
  (Wolfgang Bangerth, 2015/04/19)
  </li>

  <li> New: A LinearOperator class that stores the abstract concept of a
  linear operator. The class is fully compatible with the solver and
  preconditioner interfaces. The primary purpose of this class is to
  provide syntactic sugar for complex matrix-vector operations and free the
  user from having to create, set up and handle intermediate storage
  locations by hand.
  <br>
  (Matthias Maier, 2015/04/08)
  </li>

  <li> Fixed: There were a number of places in the library where we unconditionally
  called functions <code>_mm_malloc()/_mm_free()</code> to allocate and free
  memory with a known alignment. This function, however, is only available on
  systems with x86 or x64_64 compatible processors. These places have now been
  replaced by calling <code>posix_memalign()</code> instead, a function that
  should be more widely available.
  <br>
  (Wolfgang Bangerth, 2015/04/15)
  </li>

  <li> Deprecated: The library uses functions such as CellAccessor::subdomain_id(),
  TriaAccessor::manifold_id(), etc, but used the deviant spelling
  TriaAccessor::boundary_indicator(), TriaAccessor::set_boundary_indicator(),
  TriaAccessor::set_all_boundary_indicators(). These last three functions are now
  deprecated and have been replaced by TriaAccessor::boundary_id(),
  TriaAccessor::set_boundary_id(), and TriaAccessor::set_all_boundary_ids() for
  consistency. Similar, Triangulation::get_boundary_indicators() has been
  deprecated in favor of Triangulation::get_boundary_ids().
  <br>
  (Wolfgang Bangerth, 2015/04/11)
  </li>

  <li> Changed: All example programs used to have calls to set_boundary()
  methods to deal with curved boundaries. These have been replaced with
  the corresponding set_manifold() equivalent.
  <br>
  (Luca Heltai, 2015/04/06)
  </li>

  <li> New: A new flag no_automatic_repartitioning in
  parallel::distributed::Triangulation will disable the automatic
  repartitioning when calling execute_coarsening_and_refinement() (or things
  like refine_global(), ...), resulting in all cells staying on the processor
  they were before. The new function repartition() will execute the
  repartitioning as done automatically before.
  <br>
  (Timo Heister, 2015/03/22)
  </li>

  <li> Changed: All (Block)Compressed*SparsityPattern classes got
  replaced by DynamicSparsityPattern and
  BlockDynamicSparsityPattern, respectively and all examples now
  teach the dynamic way of creating dynamic sparsity patterns.
  <br>
  (Timo Heister, 2015/03/22)
  </li>

  <li> Improved: We have traditionally had a large number of exceptions
  that did not output any useful error message other than the name
  of the exception class. This name was suggestive of the error that
  had happened, but did not convey a sufficient amount of information
  to what happened in many of the places where these kinds of exceptions
  were used, nor what may have caused the exception, or how it could
  be avoided. We have gone through many of these places and changed
  the exception to be much more verbose in what they state about the
  problem, its origin, and how it may be solved.
  <br>
  (Wolfgang Bangerth, 2015/02/28-2015/03/31)
  </li>

  <li> Changed: We have traditionally used Point@<dim@> to represent points
  in physical space, i.e., vectors that are anchored at the origin, whereas
  for vectors anchored elsewhere (e.g., differences between points, normal
  vectors, gradients, etc) we have used Tensor@<1,dim@>. This has now be
  made more formal in the documentation but also in the return types of
  <code>operator-()</code> for Point objects: The difference between two
  points, <code>p1-p2</code> now returns a Tensor@<1,dim@>. On the other
  hand, subtracting a Tensor@<1,dim@> object from a Point, <code>p-t</code>,
  results in a Point@<dim@>.
  <br>
  (Wolfgang Bangerth, 2015/02/05)
  </li>

  <li> New: Examples from 1 to 16 now use the Manifold interface
  instead of the old Boundary interface to describe curved boundaries
  and domains.
  <br>
  (Luca Heltai, 2015/01/15)
  </li>

  <li> New: The build system now queries for git branch name and
  revision sha1 (and automatically reconfigures if necessary). This
  information is used to annotate summary.log and detailed.log with the
  current revision sha1. Further, a header file <deal.II/base/revision.h>
  is now available that exports the macros: DEAL_II_GIT_BRANCH,
  DEAL_II_GIT_REVISION, DEAL_II_GIT_REVISION_SHORT.
  <br>
  (Matthias Maier, 2015/01/02)
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="821-830-specific"></a>
<h3>Specific improvements</h3>



<ol>
  <li> New: VectorTools::get_position_vector now works with arbitrary
  FESystems, provided that the geometrical components are primitive,
  and that you provide a component mask to select what components of
  the finite element to use for the geometrical interpolation.
  <br>
  (Luca Heltai, 2015/07/25)
  </li>

  <li> Fixed: parallel::distributed::refine_and_coarsen_fixed_fraction()
  in rare circumstances decided to not refine any cells at all, even
  if the refinement threshold was nonzero. This is now fixed.
  <br>
  (Wolfgang Bangerth, Andrea Bonito, 2015/07/24)
  </li>

  <li> Fixed: Bug in DynamicSparsityPattern::iterator would cause invalid
  stl::vector::iterator comparison.
  <br>
  (Timo Heister, 2015/07/22)
  </li>

  <li>New: parallel::distributed::Triangulation::add_periodicity
  now allows for arbitrary orientations between matching faces.
  <br>
  (Daniel Arndt, 2015/07/12)
  </li>

  <li> New: Utilities::trim() function removes trailing and leading spaces.
  <br>
  (Timo Heister, 2015/07/11)
  </li>

  <li>Changed: The IsBlockMatrix class is now declared in
  <code>constraint_matrix.h</code> instead of its former home in
  <code>block_indices.h</code>.
  <br>
  (Wolfgang Bangerth, 2015/07/10)
  </li>

  <li>New: CellId::to_string() returns a string representation of a CellId object.
  <br>
  (Timo Heister, 2015/07/05)
  </li>

  <li>New: Utilities::replace_in_string().
  <br>
  (Timo Heister, 2015/07/05)
  </li>

  <li>Improved: GridOut::write_vtk() and GridOut::write_vtu() now
  output material id, level and subdomain ids of the cells.
  <br>
  (Guido Kanschat, 2015/07/05)
  </li>

  <li>Improved: The font scaling in GridOut::write_svg() was broken,
  since the units were missing. It has been fixed and an additional
  parameter GridOutFlags::Svg::cell_font_scaling has been introduced
  for tuning.
  <br>
  (Guido Kanschat, 2015/07/04)
  </li>

  <li> New: VectorizedArray now provides two methods
  vectorized_load_and_transpose() and vectorized_transpose_and_store() that
  perform vectorized reads or writes and convert from array-of-struct into
  struct-of-array or the other way around.
  <br>
  (Martin Kronbichler, 2015/07/02)
  </li>

  <li>New: GridGenerator::cheese() for a mesh with many holes;
  GridGenerator::simplex() for simplices in 2 and 3 dimensions;
  GridGenerator::hyper_cross() for crosses in 2 and 3 dimensions.
  <br>
  (Guido Kanschat, 2015/07/02)
  </li>

  <li> Fixed: The specialization of DoFAccessor for zero-dimensional objects,
  i.e., for vertices as created by accessing the faces of one-dimensional
  cells, had a member function DoFAccessor::child() that was declared but not
  implemented. This is now fixed.
  <br>
  (Wolfgang Bangerth,  2015/07/01)
  </li>

  <li> Improved: Functions::Monomial::gradient function now works when both base and exponent
  are equal to zero for one or more components of the monomial.
  Also, an assertion is added to avoid exponentiation of negative base numbers with real exponents.
  <br>
  (Angel Rodriguez,  2015/06/29)
  </li>

  <li> Fixed: The function numbers::is_finite() produced incorrect results when
  called with a NaN number (specifically, it produces an uncatchable floating
  point exception when called with a signaling NaN). This was clearly not
  intended since such values are definitely not finite.
  <br>
  (Wolfgang Bangerth, 2015/06/29)
  </li>

  <li> Improved: The SparseMatrix class can now also use <code>std::complex</code>
  scalars for its elements.
  <br>
  (Wolfgang Bangerth, 2015/06/26)
  </li>

  <li> Improved: FE_DGQArbitraryNodes::get_name() now also detects if
  the quadrature rule was Gauss points.
  <br>
  (Guido Kanschat, 2015/06/22)
  </li>

  <li> Improved: DoFRenumbering::Cuthill_McKee() can now also
  use starting indices for parallel triangulations.
  <br>
  (Wolfgang Bangerth, 2015/06/11)
  </li>

  <li> Improved: VectorTools::interpolate now works with FE_Nothing.
  <br>
  (Angel Rodriguez, 2015/06/03)
  </li>

  <li> Improved: deal.II now uses a variety of strategies to silence compiler
  warnings about unused variables and unused parameters.
  <br>
  (David Wells, 2015/04/13)
  </li>

  <li> New: Add a clear function to the PETSc::Vector
  and PETSc::MPI::Vector classes similar to the Trilinos vector classes.
  <br>
  (Martin Steigemann 2015/05/22)
  </li>

  <li> New: Three new quadrature formulas in quadrature_lib, based on
  Chebyshev quadrature rules. See functions QGaussChebyshev,
  QGaussRadauChebyshev and QGaussLobattoChebyshev.
  <br>
  (Giuseppe Pitton, Luca Heltai 2015/05/11)
  </li>

  <li> Fixed: MatrixOut now also works with Trilinos and PETSc matrices.
  <br>
  (Wolfgang Bangerth, 2015/05/11)
  </li>

  <li> Changed: TrilinosWrappers::Vector, TrilinosWrappers::BlockVector,
  PETScWrappers::Vector, and PETScWrappers::BlockVector are deprecated. Either
  use the MPI or the deal.II version of the Vector/BlockVector.
  <br>
  (Bruno Turcksin, 2015/05/04)
  </li>

  <li> Fixed: GridGenerator::half_hyper_shell can now be colorized.
  <br>
  (Daniel Arndt, 2015/05/05)
  </li>

  <li> New: The function VectorTools::point_gradient has been added to compute
  the gradient of a given FE function.
  <br>
  (Daniel Arndt, 2015/05/03)
  </li>

  <li> New: dealii:Vector, dealii::BlockVector,
  TrilinosWrappers::MPI::Vector, TrilinosWrappers::MPI::BlockVector,
  PETScWrappers::MPI::Vector and PETScWrappers::MPI::BlockVector now have
  move constructors and move assignment operators in C++11 mode.
  <br>
  (Matthias Maier, 2015/05/01)
  </li>

  <li> New: Introduce DoFRenumbering::block_wise for multigrid computation.
  <br>
  (Timo Heister, Florian Sonner, 2015/04/30)
  </li>

  <li> New: There are now MPI sum functions for Tensors and SymmetricTensors
  in the Utilities::MPI namespace.
  <br>
  (Ian Rose, 2015/04/24)
  </li>

  <li> Fixed: project_boundary_values_curl_conforming_l2() produced incorrect
  results for non-uniform grids in 2D. Adjustment to the way 2D tangents to edges are
  computed fixes this.
  <br>
  (Ross Kynch, 2015/04/23)
  </li>

  <li> Fixed: The TimerOutput class reported abnormally large cpu time when run
  with more than one process with MPI. Now this is fixed.
  <br>
  (Lei Qiao, 2015/04/19)
  </li>

  <li> New: The VectorTools::integrate_difference() function can now
  also compute the $H_\text{div}$ seminorm, using the
  VectorTools::Hdiv_seminorm argument.
  <br>
  (Zhen Tao, Arezou Ghesmati, Wolfgang Bangerth, 2015/04/17)
  </li>

  <li> Fixed: The class SymmetricTensor<2,dim> is now usable also for dim>3.
  <br>
  (Martin Kronbichler, 2015/04/14)
  </li>

  <li> New: The DynamicSparsityPattern class (formerly called
  CompressedSparsityPattern) now has an iterator class that allows to
  walk over the nonzero elements of a matrix represented by this class.
  <br>
  (Wolfgang Bangerth, 2015/04/13)
  </li>

  <li> New: The GridGenerator::subdivided_hyper_cube() and
  GridGenerator::subdivided_hyper_rectangle() now work also for codimension
  one and two Triangulation;
  <br>
  (Luca Heltai, 2015/04/12)
  </li>

  <li> New: A new VectorTools::get_position_vector() function has been
  added to the library that allows one to interpolate the Geometry of
  a (possibly curved) triangulation to vector finite element fields
  of at least spacedim components.
  <br>
  (Luca Heltai, 2015/04/11)
  </li>

  <li> New: TrilinosWrappers::BlockSparseMatrix now has member functions
  TrilinosWrappers::BlockSparseMatrix::domain_paritioner() and
  TrilinosWrappers::BlockSparseMatrix::range_partitioner() that return a
  vector of the underlying block Epetra_Map.
  <br>
  (Matthias Maier, 2015/04/08)
  </li>

  <li> New: A new MappingFEField() class has been added to the library
  that generalizes MappingQEulerian to allow arbitrary FiniteElements.
  <br>
  (Marco Tezzele, Luca Heltai, 2015/04/06)
  </li>

  <li> Changed: The cells of coarse meshes in
  parallel::distributed::Triangulation used to be ordered in a Cuthill-McKee
  numbering, which yields very high surface-to-volume ratios of the individual
  processors' partition in case the coarse mesh consists of many cells, in
  particular in 3D. The algorithm now uses SparsityTools::reorder_hierarchical
  in order to get more compact partitions, similarly to the z-ordering applied
  by p4est.
  <br>
  (Martin Kronbichler, 2015/04/10)
  </li>

  <li> New: There is now a new method
  GridTools::get_vertex_connectivity_of_cells.
  <br>
  (Martin Kronbichler, 2015/04/10)
  </li>

  <li> New: There is now a new method SparsityTools::reorder_hierarchical to
  sort nodes of a graph (sparsity pattern) in a z-like way by hierarchical
  grouping of neighboring nodes.
  <br>
  (Martin Kronbichler, 2015/04/10)
  </li>

  <li> Changed: The methods SparsityTools::reorder_Cuthill_McKee and
  GridTools::get_face_connectivity_of_cells used to take a SparsityPattern as
  argument. The data type has been changed to DynamicSparsityPattern in order
  to avoid copying things around. The old interface is still available but
  marked as deprecated.
  <br>
  (Martin Kronbichler, 2015/04/10)
  </li>

  <li> New: One can now ask a cell the how many-th active cell it is,
  using CellAccessor::active_cell_index() (which is called as
  <code>cell-@>active_cell_index()</code>).
  <br>
  (Wolfgang Bangerth, 2015/04/10)
  </li>

  <li> Fixed: Added missing hp-related functions to FE_Q_Hierarchical together with
  a couple of unit tests. Improved code comments.
  <br>
  (Denis Davydov, 2015/04/10)
  </li>

  <li> Fixed: deal.II did not compile on 32-bit systems when using newer
  p4est versions (1.0 and later) due to a type mismatch. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2015/04/08)
  </li>

  <li> Changed: In the spirit of the changes made to the distinction
  between Point and Tensor objects discussed above, the first argument
  to GridTools::shift() has been changed from a Point to a Tensor@<1,dim@>.
  <br>
  (Wolfgang Bangerth, 2015/04/02)
  </li>

  <li> New: There is now a new quadrature formula in quadrature_lib. It is
  now possible to use Telles' quadrature rules through the function QTelles
  to integrate singular integrals
  <br>
  (Nicola Giuliani, 2015/04/01)
  </li>

  <li> New: Added FE_Bernstein: a scalar finite element based on Bernstein basis polynomials.
  <br>
  (Marco Tezzele, Luca Heltai, 2015/03/31)
  </li>

  <li> New: A function to get a map with all vertices at boundaries has
  been added at GridTools::get_all_vertices_at_boundary(). This function
  will return a map which can be used in functions like
  GridTools::laplace_transform().
  <br>
  (Fernando Posada, 2015/03/31)
  </li>

  <li> Fixed: TrilinosWrappers::SparseMatrix::local_range() erroneously
  threw an exception in 64-bit mode. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2015/03/24)
  </li>

  <li> New: The GridOut::write_gnuplot() function produced output
  for 1d meshes embedded in higher dimensional spaces that was
  invalid in that the lines showing individual cells were connected.
  While this is not wrong for singly connected 1d meshes, it leads to wrong
  results if the domain is not singly connected and not every cell is the
  right neighbor of the previous cell.
  <br>
  (Wolfgang Bangerth, 2015/03/23)
  </li>

  <li> New: The various file format writers of class GridOut were not
  instantiated for 1d meshes in 3d space. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2015/03/23)
  </li>

  <li> New: ParameterHandler::declare_alias() allows to define
  alternate names for parameters. This is primarily intended to allow
  for backward compatible ways of changing the names of parameters
  to applications.
  <br>
  (Wolfgang Bangerth, 2015/03/22)
  </li>

  <li> New: GridGenerator::create_triangulation_with_removed_cells() creates
  a new mesh out of an existing one by dropping individual cells.
  <br>
  (Wolfgang Bangerth, 2015/03/13)
  </li>

  <li> New: Add MueLu preconditioner from Trilinos through the class
  TrilinosWrappers::PreconditionAMGMueLu. This is a new algebraic
  multigrid package. The input parameters are almost the same as the ones
  from ML so that the two preconditioners can be easily swapped.
  <br>
  (Bruno Turcksin, 2015/03/11)
  </li>

  <li> Fixed: Iterating over the elements of a TrilinosWrappers::SparseMatrix
  object previously led to errors if the matrix was in fact stored in
  parallel across multiple MPI processes. This is now fixed: rows not
  stored locally on the processor where you run the iteration simply look
  like they're empty.
  <br>
  (Wolfgang Bangerth, 2015/03/08)
  </li>

  <li> New: There is now a new macro DeclExceptionMsg that allows to
  declare an exception that does not take any run-time arguments
  yet still allows to specify an error message.
  <br>
  (Wolfgang Bangerth, 2015/02/27)
  </li>

  <li> New: There is now a class std_cxx11::unique_ptr that provides
  the functionality of std::unique_ptr in C++11 mode, and that
  provides an emulation for older compilers.
  <br>
  (Wolfgang Bangerth, 2015/02/22)
  </li>

  <li> New: IndexSet now has a local typedef IndexSet::value_type.
  <br>
  (Wolfgang Bangerth, 2015/02/22)
  </li>

  <li> New: FE_TraceQ can now also be used in 1D.
  <br>
  (Martin Kronbichler, 2015/02/21)
  </li>

  <li> New: FE_TraceQ now implements get_face_interpolation_matrix and
  get_subface_interpolation_matrix, enabling
  DoFTools::make_hanging_node_constraints on this element.
  <br>
  (Anton Evgrafov, 2015/02/21)
  </li>

  <li> Fixed: MappingQEulerian would previously not move interior points
  in 1D for higher order mappings. This has been fixed by removing a few
  specializations of MappingQ for 1D that are no longer necessary.
  <br>
  (Martin Kronbichler, 2015/02/19)
  </li>

  <li> Fixed: The implementation of the class GrowingVectorMemory has been
  moved from source/lac/vector_memory.cc to the new file
  include/deal.II/lac/vector_memory.templates.h. This allows users to
  create instantiations of GrowingVectorMemory for their own vector classes
  in case they intend to use them for the deal.II solvers.
  <br>
  (Martin Kronbichler, 2015/02/18)
  </li>

  <li> Changed: All members of MultithreadInfo are now static so it is no
  longer necessary to use the global instance multithread_info (now
  deprecated) or create your own instance (which does not work correctly
  anyway).
  <br>
  (Timo Heister, 2015/02/13)
  </li>

  <li> Fixed: There was a bug in the energy source term of step-33 whereby
  the term was erroneously multiplied by the density. This is now fixed.
  <br>
  (Praveen C, Lei Qiao, 2015/02/13)
  </li>

  <li> Changed: If you take the product of a Tensor and a scalar number,
  you previously got a Tensor back that stored its elements in the same
  data type as the original tensor. This leads to problems if you
  multiply a <code>Tensor@<1,dim,double@></code> by a
  <code>std::complex@<double@></code> because the result clearly needs
  to store its elements as complex numbers, rather than as double
  variables. This is now changed: The result of the product of a Tensor
  and a scalar number is now a Tensor that stores its elements in a data
  type appropriate for this product. The same approach is taken for the
  SymmetricTensor class.
  <br>
  (Wolfgang Bangerth, 2015/02/11)
  </li>

  <li> New: There is now a new class ProductType that can be used
  to infer the type of the product of two objects. There is now also
  a class EnableIfScalar that helps restrict some templates to only
  cases where a type is a scalar.
  <br>
  (Wolfgang Bangerth, 2015/02/04)
  </li>

  <li> New: The Tensor classes now have copy constructors and copy
  operators that allow assignment from other tensors with different
  underlying scalar types.
  <br>
  (Denis Davydov, 2015/02/03)
  </li>

  <li> New: Class hp::DoFHandler can now also be serialized.
  <br>
  (Lukas Korous, 2015/01/31)
  </li>

  <li> Bugfix: deal.II now correctly links against librt in case of bundled
  boost being used.
  <br>
  (Matthias Maier, 2015/01/27)
  </li>

  <li> New: A new macro <code>DEAL_II_QUERY_GIT_INFORMATION</code> is
  provided to query user projects for git repository information similarly
  to those exported by deal.II.
  <br>
  (Matthias Maier, 2015/01/21)
  </li>

  <li> Fixed: FEFaceValues and FESubfaceValues did not fill the
  jacobians and inverse jacobians if requested via the update flags.
  This is now fixed.
  <br>
  (Martin Kronbichler, 2015/01/23)
  </li>

  <li> Fixed: ParameterHandler::read_input() now checks that
  'subsection'/'end' are balanced in the input.
  <br>
  (Timo Heister, 2015/01/19)
  </li>

  <li> Fixed: In 3d, when you set the <code>colorize</code> flag of
  GridGenerator::hyper_shell(), the faces of the domain were colored but
  the edges were not. This was an oversight because to refine correctly,
  the edges also have to have the appropriate boundary indicator set.
  <br>
  (Wolfgang Bangerth, 2015/01/16)
  </li>

  <li> New: MultithreadInfo::n_cpus() returns the correct number of CPU
  on FreeBSD.
  <br>
  (Bruno Turcksin, 2015/01/14)
  </li>

  <li> Improved: MPI collective operations such as MPI::sum, MPI::max now
  check for job_supports_mpi() internally, which allows running them also
  without a call to MPI_Init.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> Changed: The method job_supports_mpi() now resides in the namespace
  Utilities::MPI instead of Utilities::System for consistency with other MPI
  methods. The old method has been marked deprecated and will be removed in
  a future version.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> Fixed: The update of ghost values in parallel::distributed::Vector when
  calling the assignment operator is now active when one of the two vector had
  its ghost values updated before or when the layout of the right hand side
  vector is one-to-one, more consistent with parallel PETSc and Trilinos
  vectors.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> New: PETScWrappers::MPI::SparseMatrix::reinit(other) copies
  the layout of another matrix. TrilinosWrappers::SparseMatrix
  operator= and copy constructor are now disabled. This brings
  functionality between PETSc and Trilinos in line.
  <br>
  (Timo Heister, 2015/01/12)
  </li>

  <li> New: Triangulation::set_all_manifold_ids() and
  Triangulation::set_all_manifold_ids_on_boundaries()
  set all manifold ids on every object or on every
  boundary object respectively.
  <br>
  (Luca Heltai, 2015/01/12)
  </li>

  <li> New: GridTools::copy_boundary_to_manifold_id() and
  GridTools::copy_material_to_manifold_id() copy
  boundary_ids and material_ids to manifold_ids for
  faces on the boundary and for cells respectively.
  <br>
  (Luca Heltai, 2015/01/09)
  </li>

  <li> Fixed: Utilities::int_to_string() produced wrong results if
  the number of digits specified was ten or greater.
  <br>
  (David Wells, 2015/01/08)
  </li>

  <li> Fixed: VectorTools::interpolate_to_different_mesh() was accidentally
  only instantiated for dealii::Vector arguments, rather than all vector
  classes. This is now fixed.
  <br>
  (Benjamin Brands, Wolfgang Bangerth, 2014/12/29)
  </li>

  <li> Fixed: Use CASROOT environment variable as additional hint for
  opencasacade.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>

  <li> Fixed: Update the run_testsuite.cmake script to also pick up
  muparser and opencascade configuration.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>

  <li> Fixed: Update several places in the documentation that were not
  updated from functionparser to muparser. Add several forgotten
  DEAL_II_WITH_* variables to certain places in the documentation.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>
</ol>

*/
