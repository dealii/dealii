// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
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
@page changes_between_9_4_0_and_9_5_0 Changes between Version 9.4.0 and 9.5.0

<p>
This is the list of changes made between the release of deal.II version
9.4.0 and that of 9.5.0. All entries are signed with the names of the
author.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="940-950-incompatible"></a>
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
  Changed: The SUNDIALS::KINSOL, SUNDIALS::ARKode, and SUNDIALS::IDA
  interfaces used to use callbacks through which users provide
  information to these solvers. If user code encountered errors, this
  was indicated via integer return codes -- a zero return value
  indicating success.
  This is not in line with how C++ typically operates; in C++ errors are
  typically indicated via exceptions. The interfaces to the libraries
  above have been changed to now use exceptions instead, and this
  convention has been documented in a glossary entry on
  @ref GlossUserProvidedCallBack "user provided callback".
  <br>
  (Wolfgang Bangerth, 2023/05/31)
 </li>

 <li>
  Changed: The integrate() function of FEPointEvaluation now multiplies the submitted
  entities internally with a JxW value to represent a proper integration like in FEEvaluation.
  The old behavior (where the entities are only multiplied with the test functions and summed up)
  is available through the new function test_and_sum().
  <br>
  (Maximilian Bergbauer, 2023/05/29)
 </li>

 <li>
  Deprecated: All vector classes had member functions called `import()`
  to copy elements from a source into the current vector. Unfortunately,
  `import` is a
  [keyword (of sorts) in C++20](https://en.cppreference.com/w/cpp/keyword/import).
  While not strictly necessary because we do not use the name in a
  context where it would be recognized as a keyword, it is useful to
  avoid the name nonetheless, if only to avoid confusing readers and
  IDEs. As a consequence, these functions have been renamed to
  `import_elements()`. The old name remains for now, but is deprecated.
  <br>
  (Wolfgang Bangerth, 2023/05/23)
 </li>

 <li>
  Extensive rework of CUDAWrappers::MatrixFree and CUDAWrappers::FEEvaluation.
  <br>
  (Bruno Turcksin, 2023/05/11)
 </li>

 <li>
  Changed: The get_mpi_communicator() methods no longer return references to MPI_Comm.
  They all now return the MPI_Comm by value.
  <br>
  (Stefano Zampini, 2023/04/12)
 </li>

 <li>
  Removed: The docker image with root user has been removed.
  <br>
  If you would like to use the docker image in the context of github actions,
  you can use the regular image and override the default user like this:
  ```
  container:
    image: dealii/dealii:master-focal
    options: --user root
  ```
  (Marc Fehling, 2023/04/09)
 </li>

 <li>
  Removed: The CUDAWrappers::MatrixFee::AdditionalData::ParallelizationScheme
  parameter has been removed.
  <br>
  (Bruno Turcksin, 2023/03/30)
 </li>

 <li>
  Updated: deal.II no longer exports the <code>-fPIC/-fpic</code> compiler
  flags for compiling relocatable, position independent code. This should
  have minimal impact on user projects because CMake is handling the
  <code>-fPIC/-fpic</code> automatically. But for some less common
  scenarios in user projects (e.g., object targets) it might be necessary
  to set the POSITION_INDEPENDENT_CODE target property in CMake.
  <br>
  (Matthias Maier, 2023/03/29)
 </li>

 <li>
  Changed: The default behavior of the hdf5 output is changed from "not compressed" to "compressed" (default compression_level: best_speed).
  <br>
  (Christoph Schmidt, 2023/03/24)
 </li>

 <li>
  Updated: deal.II now requires CMake version 3.13.4 or later.
  <br>
  (Matthias Maier, 2023/01/26)
 </li>

 <li>
   Updated: The minimum version for Trilinos has been bumped to 12.14.1 if Trilinos bundles Kokkkos.
   <br>
   (Daniel Arndt, 2022/12/30)
 </li>

 <li>
  Fixed: The function ParameterHandler::add_parameter() used to
  call the internal action. Within that step, the action
  converts the default value to a string and back afterwards.
  This can lead to round-off errors so that the default
  values might change in the case of floating-point numbers.
  The action is not called any more during ParameterHandler::add_parameter(),
  fixing the problem.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/12/11)
 </li>

 <li>
  Removed: A number of obscure CMake configuration options have been removed:
  <code>DEAL_II_PREFER_STATIC_LINKAGE</code>,
  <code>DEAL_II_STATIC_EXECUTABLE</code>,
  <code>DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS</code>.
  <br>
  (Matthias Maier, 2022/12/02)
 </li>

 <li>
  Removed: The build system no longer generates pkg-config files.
  <br>
  (Matthias Maier, 2022/11/30)
 </li>

 <li>
  Removed: The <code>deal.IIFeatureConfig.cmake</code> configuration file
  has been removed. This file used to contain a record of all
  <code>FEATURE_[suffix]</code> variables for a given external library.
  Said information is now exported via interface targets in CMake.
  <br>
  (Matthias Maier, 2022/11/30)
 </li>

 <li>
  Changed: We introduced a stronger datatype for FE indices:
  <ul>
    <li> Instead of unsigned int, we now use types::fe_index for FE indices throughout the library.
    <li> For invalid FE indices, we no longer use numbers::invalid_unsigned_int, but numbers::invalid_fe_index.
  </ul>
  This affects interfaces to a lot of functions. Most of them are backwards compatible via implicit conversion.
  Functions that do not fall in that category are listed below and might need your attention.
  Some of them have been deprecated, some are replaced because we imagine them only being used internally.
  <ul>
    <li> DoFHandler::set_active_fe_indices() (deprecated)
    <li> DoFHandler::get_active_fe_indices() (deprecated)
  </ul>
  Further, the serialization of parallel::distributed::Triangulation objects is also affected by this change
  as we write FE indices to disk. Thus, the version in the serialization metafile has been bumped to `5`.
  <br>
  (Marc Fehling, 2022/10/27)
 </li>

 <li>
  Removed: The base class of SparsityPattern, SparsityPatternBase, has been removed.
  <br>
  (David Wells, 2022/09/22)
 </li>

 <li>
  Updated: The version of the vtu output has been increased to 2.2
  if high-order output is requested.
  This is necessary to fix visualization issues in Paraview.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/09/18)
 </li>

 <li>
  Changed: The function DoFHandler::get_active_fe_indices() now returns
  the result vector. The previous version that writes into the argument
  as well as DoFTools::get_active_fe_indices() have been deprecated.
  <br>
  (Marc Fehling, 2022/08/20)
 </li>

 <li>
  Moved: The functions in Utilities::Trilinos have been moved to a separate header,
  base/trilinos_utilities.h.
  <br>
  (David Wells, 2022/08/16)
 </li>

 <li>
  Changed: The oldest supported version of PETSc has been increased from 3.3.0
  to 3.7.0.
  <br>
  (David Wells, 2022/07/27)
 </li>

 <li>
  Removed: The deprecated QTrapez, ParticleHandler::locally_relevant_ids(),
  and Arkode::reinit_vector have been removed.
  <br>
  (Daniel Arndt, 2022/06/29)
 </li>

 <li>
  Removed: The deprecated member function CellId::to_cell() has been removed.
  <br>
  (Daniel Arndt, 2022/06/28)
 </li>

 <li>
  Removed: The deprecated header file lac/parallel_block_vector.h
  has been removed.
  <br>
  (Daniel Arndt, 2022/06/27)
 </li>

 <li>
  Removed: The deprecated classes ConstantFunction and ZeroFunction
  have been removed.
  <br>
  (Daniel Arndt, 2022/06/24)
 </li>

 <li>
  Removed: The deprecated member functions DoFHandler::initialize() and
  DoFHandler::set_fe() have been removed.
  <br>
  (Daniel Arndt, 2022/06/24)
 </li>

 <li>
  Removed: The deprecated member functions of QProjector not using a ReferenceCell
  object have been removed.
  <br>
  (Daniel Arndt, 2022/06/24)
 </li>

 <li>
  Removed: The deprecated functions MatrixFree::n_macro_cells(),
  MatrixFree::get_hp_cell_iterator(), MatrixFree::n_components_filled(),
  MatrixFree::get_dof_handler(), and
  MatrixFree::reinit() with default Mapping argument have been removed.
  <br>
  (Daniel Arndt, 2022/06/22)
 </li>

 <li>
  Removed: The deprecated constructors for SparseVanka and SparseBlockVanka
  have been removed.
  <br>
  (Daniel Arndt, 2022/06/22)
 </li>

 <li>
  Removed: The deprecated LinearAlgebra::CUDAWrappers::atomicAdd_wrapper
  has been removed.
  <br>
  (Daniel Arndt, 2022/06/22)
 </li>

 <li>
  Rotated: The list of incompatible changes has been rotated.
  <br>
  (Matthias Maier, 2022/06/16)
 </li>

 <li>
  Removed: The previously deprecated class hp::DoFHandler has been
  removed, and with it all functions and classes that required the
  template parameter DoFHandlerType.
  <br>
  From now on, use the standard DoFHandler. All hp-functionality has been
  integrated to this one.
  <br>
  (Marc Fehling, Peter Munch, 2022/05/24)
 </li>

 <li>
  Changed: All GridIn functions now remove unused vertices and will attempt to fix
  pyramids and wedges with negative volumes.
  <br>
  (David Wells, 2022/05/25)
 </li>

 <li>
  Removed: Deprecated constructors for PETScWrappers::MPI::SparseMatrix,
  a deprecated overload for PETScWrappers::MatrixBase::add and a deprecated
  constructor for PETScWrappers::MPI::Vector have been removed.
  <br>
  (Daniel Arndt, 2020/04/21)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="940-950-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: Added support for the VTK library (https://vtk.org/).
  <br>
  (Pasquale Africa, 2023/05/24)
 </li>

 <li>
  New: The new function distributed_compute_intersection_locations() computes
  intersections on parallel::distributed::Triangulations using CGAL. It
  returns a new class DistributedComputeIntersectionLocationsInternal which
  stores the found intersections and relevant information about communication.
  The class contains a helper function
  convert_to_distributed_compute_point_locations_internal() to produce
  DistributedComputePointLocationsInternal without an additional
  (possibly expensive) search of points and can therefore be used to
  initialize RemotePointEvaluation.
  <br>
  (Johannes Heinz, 2023/05/18)
 </li>

 <li>
  Changed: Several parts of the library involve interfacing with
  external libraries by way of user-defined callback functions. Specific
  examples are the interfaces to the SUNDIALS solvers (e.g., the
  SUNDIALS::KINSOL class). These interfaces typically required using the
  convention for error reporting defined by the underlying library -- in
  the case of SUNDIALS, for example, callbacks needed to return zero in
  case of success, a negative value for an irrecoverable error, and a
  positive value for a recoverable error.
  This approach does not scale across the many interfaces we have. As a
  consequence, we standardized how callbacks should behave, as
  documented in @ref GlossUserProvidedCallBack "this glossary entry".
  The interfaces in SUNDIALS::KINSOL and SUNDIALS::ARKode have been
  changed correspondingly.
  <br>
  (Wolfgang Bangerth, 2021/05/08)
 </li>

 <li>
  New: Add a two level transfer operator between non-nested multigrid levels for CG and DG.
  <br>
  (Marco Feder, Johannes Heinz, Peter Munch, Martin Kronbichler, 2023/04/26)
 </li>

 <li>
  New: Add wrappers to the non-linear solver SNES (PETScWrappers::NonLinearSolver ) and
  the ode integrator TS (PETScWrappers::TimeStepper) from the PETSc library.
  <br>
  (Stefano Zampini, 2023/04/22)
 </li>

 <li>
  New: Kokkos is now a required dependency.
  <br>
  (Daniel Arndt, 2022/10/26)
 </li>

 <li>
  New: Add a wrapper to the non-linear solver NOX from the
  Trilinos library.
  <br>
  (Peter Munch, Vladimir Ivannikov, 2022/10/02)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="940-950-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: An out-of-bounds write due to an insufficiently sized buffer in
  the Scalapack wrappers has been fixed.
  <br>
  (Matthias Maier, 2023/07/03)
 </li>

 <li>
  Fixed: With some compilers, calling Threads::Task::Task() with a
  function object that ends with an exception lead to a segmentation
  fault. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2023/06/26)
 </li>

 <li>
  Fixed: Originally, the LAPACK wrappers in
  source/lac/lapack_full_matrix.cc directly printed errors to the console.
  This is fixed by using the AssertThrow mechanism to explicitly display
  various error messages based on the error values that are consistent with
  the LAPACK manual.
  <br>
  (Tao Jin, 2023/06/20)
 </li>

 <li>
  Fixed: The numbering applied within FEPointEvaluation did not work correctly
  for multiple components and when used with first_selected_component!=0 for
  certain FiniteElement types, leading to wrong results. This is now fixed.
  <br>
  (Johannes Heinz, Peter Munch, 2023/06/11)
 </li>

 <li>
  Fixed: DoFCellAccessor::distribute_local_to_global() now works with std:: iterators.
  <br>
  (Johannes Heinz, 2022/06/06)
 </li>

 <li>
  New: All parallel BlockVector classes now have a reinit function that takes
  collections of Utilities::MPI::Partitioner objects as an argument. This
  affects TrilinosWrappers::MPI::BlockVector, PETScWrappers::MPI::BlockVector,
  and LinearAlgebra::distributed::BlockVector. The interface is compatible.
  <br>
  (Marc Fehling, 2023/06/02)
 </li>

 <li>
  Improved: Make it possible to fill JxW in NonMatching::MappingInfo from quadratures which contain JxW as weights.
  <br>
  (Johannes Heinz, 2023/06/01)
 </li>

 <li>
  Fixed: Patterns::pattern_factory now parses lists of lists correctly.
  <br>
  (Daniel Arndt, 2023/05/16)
 </li>

 <li>
  Changed: The KINSOL wrapper used to copy vectors to/from deal.II format
  to/from a native SUNDIALS format. These copies are no longer made.
  <br>
  (Sebastian Proell, 2023/05/15)
 </li>

 <li>
  Fixed: step-77 contained a bug where the residual was not
  zeroed out before assembly.
  <br>
  (Sebastian Proell, 2023/05/14)
 </li>

 <li>
  Changed: MPI communicators (of type `MPI_Comm`) are now treated as
  trivially copyable POD types throughout the library. This means that
  we store and pass them by value instead of by reference.
  <br>
  (Timo Heister, 2023/05/08)
 </li>

 <li>
  New: The new class ScopeExit is useful in ensuring that cleanup code
  is always run, regardless of how a function is left.
  <br>
  (Wolfgang Bangerth, 2023/05/05)
 </li>

 <li>
  New: Added new nonlinear solver class NonlinearSolverSelector.
  <br>
  (Sean Ingimarson, 2023/05/01)
 </li>

 <li>
  New: FEInterfaceValues::reinit() sometimes had trouble automatically
  determining which quadrature or mapping object to use in the hp
  context. There are now two more rules that help determine which
  objects to use, as described in detail in the documentation.
  <br>
  (Wolfgang Bangerth, 2023/04/25)
 </li>

 <li>
  Fixed: Consider the following scenario: You initialize a DoFHandler with a
  parallel::distributed::Triangulation, and afterwards you want to create the
  Triangulation by calling Triangulation::create_triangulation() that takes a
  TriangulationDescription::Description as an argument. This caused an issue
  as the active FE tables of the DoFHandler have not been properly updated
  during the creation of the Triangulation. This is fixed now by treating the
  resulting Triangultion as a 'fresh' one by triggering Triangulation::Signals::create
  at the end of Triangulation::create_triangulation(), which in turn
  reinitializes the DoFHandler.
  <br>
  (Marc Fehling, 2023/04/21)
 </li>

 <li>
  Fixed: ParameterAcceptor::clear() had issues with 'use after free'.
  <br>
  (Luca Heltai, 2023/04/20)
 </li>

 <li>
  Fixed: Mesh files in COMSOL's `.mphtxt` format are often created on
  Windows and then read using deal.II on Unix. Windows uses a different
  line ending style than Linux, and this led to error messages that were
  impossible to understand because the different line ending style is
  not visible in most editors. This is now fixed: The
  GridIn::read_comsol_mphtxt() function also accepts Windows line
  endings.
  <br>
  (Wolfgang Bangerth, 2023/04/20)
 </li>

 <li>
  Fixed: In FEInterfaceValues::reinit(), if `q_index` is not explicitly
  given but none of the two adjacent elements dominates the other, one
  used to get a rather obscure error. The error now produced points at
  the actual problem, and a much expanded discussion of what is going
  wrong is provided in the documentation of the function.
  <br>
  (Simon Sticko, Wolfgang Bangerth, 2023/04/13)
 </li>

 <li>
  Deprecated: The DoFTools::map_dofs_to_support_points() function family
  had functions that returned what they computed via a non-const
  reference argument `support_points` of type `std::map`. These
  overloads are now deprecated and replaced by functions that return the
  map via a regular `return` value.
  <br>
  (Wolfgang Bangerth, 2023/04/07)
 </li>

 <li>
  New: Added Functions::SignedDistance::ZalesakDisk()
  to compute the signed distance function of Zalesak's disk.
  <br>
  (Magdalena Schreter, Simon Sticko, Peter Munch, 2023/04/06)
 </li>

 <li>
  Fixed: In grid_out.cc VtkFlags with default settings were passed to the
  write_vtu method instead of the available vtu_flags. Thus, a
  compression_level that has been set was simply ignored. This has been
  fixed in #14983 by passing the available vtu_flags to the method.
  <br>
  (Christoph Schmidt, 2023/04/05)
 </li>

 <li>
  New: Added Functions::SignedDistance::Rectangle() and BoundingBox::signed_distance()
  to compute the signed distance function of a rectangle.
  <br>
  (Magdalena Schreter, Peter Munch, 2023/03/31)
 </li>

 <li>
  Changed: The default `compression_level` for vtu output is changed from `best_compression` to `best_speed`.
  <br>
  (Christoph Schmidt, 2023/03/28)
 </li>

 <li>
  New: Added an ArrayView::empty() function.
  <br>
  (Maximilian Bergbauer, 2023/03/21)
 </li>

 <li>
  New: Added functionality to the matrix-free CellwiseInverseMassMatrix
  operator for coupling (dyadic) coefficients. `fill_inverse_JxW_values`
  and `CellwiseInverseMassMatrixImplFlexible` have been extended.
  Internal interfaces for the coefficients now use ArrayView instead
  of AlignedVector.
  <br>
  (Bugrahan Temur, 2023/03/15)
 </li>

 <li>
  Fixed: There was no instance of GridTools::collect_periodic_faces() with
  MeshType = parallel::shared::Triangulation<1,1>. But the function works with the
  template parameter as well. Thus, it should be instantiated.
  <br>
  (Nils Schween, 2023/03/14)
 </li>

 <li>
  New: LAPACKFullMatrix now provides two new functions,
  LAPACKFullMatrix::get_right_eigenvectors() and
  LAPACKFullMatrix::get_left_eigenvectors(), to return eigenvectors after calls
  to LAPACKFullMatrix::compute_eigenvalues().
  <br>
  (Martin Kronbichler, 2023/03/09)
 </li>

 <li>
  New: Provide overloads for `make_array_view` for `AlignedVector`.
  <br>
  (Sebastian Proell, 2023/03/06)
 </li>

 <li>
  Fixed: It was previously possible to assign scalar values to
  Tensor and VectorizedArray objects that were temporaries -- say in expressions
  such as
  ```
    VectorizedArray<...> my_function();
    ...
    my_function() = 1.234;
  ```
  This does not make any sense: What `my_function()` returns is a
  temporary object, and assigning a value to it has no consequences
  because the temporary object dies at the end of the line. Whatever the
  programmer intended to do here was almost certainly a mistake.
  As a consequence, this is now prohibited and the compiler will produce
  an error when trying to do this.
  <br>
  (Wolfgang Bangerth, 2023/03/05)
 </li>

 <li>
  Fixed: SolverCG used with matrices and preconditioners that support
  interleaving of the vector updates with the matrix-vector product (e.g. with
  certain matrix-free operators) sometimes performed the final update to the
  result vector `x` with a wrong coefficient. This could affect the returned
  solution for coarse tolerances far away from convergence. This is now fixed.
  <br>
  (Martin Kronbichler, 2023/02/26)
 </li>

 <li>
  Fixed: When used with the linear-polynomial MappingQCache(1),
  FEValues::reinit() could previously trigger CellSimilarity in case the
  underlying mesh was affine, but the deformed state described via MappingQCache
  was not. This is now fixed.
  <br>
  (Martin Kronbichler, 2023/02/22)
 </li>

 <li>
  New: Implement `VectorTools::add_constant` that adds a constant
  to a given component of a finite element solution.
  <br>
  (Timo Heister, 2023/02/17)
 </li>

 <li>
  New: The function ReferenceCell::face_vertex_location() allows
  translating between the vertices of faces of cells and the coordinate
  system of the cell itself.
  <br>
  (Wolfgang Bangerth, 2023/02/09)
 </li>

 <li>
  Fixed: Fix GridTools::find_all_active_cells_around_point() to return the correct cells
  if the requested point is located on a cell edge adjacent to neighbor cells at different
  refinement levels.
  <br>
  (Peter Munch, Magdalena Schreter, 2023/02/07)
 </li>

 <li>
  Fixed: The DerivativeForm class had a constructor that took an array
  of elements via the type
  `Tensor<order, spacedim, Tensor<1, dim,Number>>`.
  This should have been
  `Tensor<1, spacedim, Tensor<order, dim, Number>>`, and it is now fixed.
  <br>
  (Robin Goermer, Wolfgang Bangerth, 2023/02/01)
 </li>

 <li>
  New: The function IndexSet::get_index_vector() returns a vector
  with all indices that are part of an index set.
  <br>
  (Wolfgang Bangerth, 2023/01/26)
 </li>

 <li>
  Improved: Mapping can now return vertices on a face of a cell.
  <br>
  (Johannes Heinz, 2023/01/25)
 </li>

 <li>
  Fixed: When outputting hexahedral meshes in AVS Explorer's UCD file format,
  cells were written with an inverted vertex order. While they were shown
  correctly, these cells had negative volume. This has now been fixed.
  <br>
  (Nils Margenberg, 2023/01/22)
 </li>

 <li>
  Improved: Introduce a function that only checks if Bounding Boxes are overlapping.
  Fixed: If a 1D BBox is completely contained in another one, it was not recognized as mergeable neighbor.
  <br>
  (Johannes Heinz, 2022/07/07)
 </li>

 <li>
  New: The function ReferenceCell::contains_points() is now implemented
  for all reference cell kinds (it had been missing for wedges and
  pyramids before).
  <br>
  (Wolfgang Bangerth, 2023/01/11)
 </li>

 <li>
  New: There are now overloads of the make_array_view() function for objects of
  type std::array.
  <br>
  (Wolfgang Bangerth, 2023/01/09)
 </li>

 <li>
  Deprecated: The functions ReferenceCell::compute_orientation() and
  ReferenceCell::permute_according_orientation() have been deprecated. Use
  ReferenceCell::get_orientation_index() and
  ReferenceCell::reorient_based_on_orientation_index() instead.
  <br>
  (Wolfgang Bangerth, 2023/01/09)
 </li>

 <li>
  New: Writing graphical output in VTU format now uses a workflow that
  uses multiple threads where possible, substantially increasing the
  speed with which functions such as DataOut::write_vtu() work.
  <br>
  (Wolfgang Bangerth, 2022/12/30)
 </li>

 <li>
  Improved: The FEInterfaceValues values public interface has been extended
  further to support hp-FEM. Functions have been added to report if the class has
  been initialized with hp support, and to get the underlying
  hp::MappingCollection, hp::FECollection and hp::QCollection. The various
  FEInterfaceValues::reinit() methods now accept indices to indicate which
  quadrature rule, mapping (and, in some cases, finite elements) should be used to
  compute values across the interface.
  <br>
  (Jean-Paul Pelteret, 2022/12/28)
 </li>

 <li>
  New: ScratchData now has support for hp finite elements on interfaces.
  <br>
  (Jean-Paul Pelteret, 2022/12/23)
 </li>

 <li>
  New: Allow the construction of a deal.II PETSc matrices from a PETSc Mat.
  <br>
  (Stefano Zampini, 2022/12/21)
 </li>

 <li>
  New: Expose PETScWrappers::MPI::BlockVector as a PETSc VecNest vector, and make sure the mpi communicator is always
  queried from the internal vec object.
  <br>
  (Stefano Zampini, 2022/12/21)
 </li>

 <li>
  New: PETScWrappers:BlockSparseMatrix now is also a PETSc MATNEST type.
  <br>
  (Stefano Zampini, 2022/12/20)
 </li>

 <li>
  Fixed: Small fixes to communicator handling in PETSc classes. Move from
  GetArray to GetArrayRead (threadsafe version).
  <br>
  (Stefano Zampini, 2022/12/20)
 </li>

 <li>
  New: LinearAlgebra::distributed::BlockVector::reinit now accepts
  a vector of Utilities::MPI::Partitioner objects.
  <br>
  (Marc Fehling, 2022/12/20)
 </li>

 <li>
  Improved: FEInterfaceValues objects can now also be constructed using hp-collections.
  <br>
  (Marco Feder, 2022/12/20)
 </li>

 <li>
  Fixed: Add default MPI getter for PETSc Mat objects, that queries the underlying PETSc type.
  <br>
  (Stefano Zampini, 2022/12/13)
 </li>

 <li>
  Fixed: PETSc has no concept of ownership, but only shared ownership. Vec and Mat objects are
  reference counted, and automatically cleaned when no longer used. We do not need to keep track
  of ownership manually.
  <br>
  (Stefano Zampini, 2022/12/13)
 </li>

 <li>
  Improved: The function GridTools::rotate()
  can now also handle Triangulation objects with
  dim=1 and spacedim=2.
  <br>
  (Peter Munch, 2022/12/08)
 </li>

 <li>
  New: SolverGMRES and SolverFGMRES now also support classical
  Gram-Schmidt orthonormalization
  alongside to the existing modified Gram-Schmidt algorithm. This
  algorithm allows to reduce the cost of vector operations in terms of
  communication latency and memory transfer significantly.
  <br>
  (Peter Munch, Martin Kronbichler, 2022/12/08)
 </li>

 <li>
  Fixed: The residual assembly in step-77 contained a bug that prevented
  the correct residual from being computed. This is now fixed.
  <br>
  (Stefano Zampini, 2022/12/06)
 </li>

 <li>
  Improved: The history of linear residuals norms inside classes SolverControl and
  ReductionControl is now cleared at the first step of the solution.
  <br>
  (Vladimir Ivannikov, 2022/12/05)
 </li>

 <li>
  New: TrilinosWrappers::MPI::Vector and PETScWrappers::MPI::Vector now both have
  reinit functions that take a Utilities::MPI::Partitioner as an argument, so
  their interface is compatible to LinearAlgebra::distributed::Vector.
  <br>
  (Marc Fehling, 2022/12/05)
 </li>

 <li>
  Improved: The quick_tests mechanism has been redesigned. Quick tests are
  now part of the regular deal.II testsuite. This means they can be
  configured via the <code>setup_tests_quick_tests</code> target, and run via
  invoking ctest from the build directory. The <code>test</code> target will
  now ensure that the library is fully compiled and quick tests are
  configured prior to running all quick tests.
  <br>
  (Matthias Maier, 2022/11/29, 2022/12/10)
 </li>

 <li>
  New: There is now a new function Threads::TaskGroup::return_values() that waits
  for all tasks in a task group and returns a vector of their returned objects.
  <br>
  (Wolfgang Bangerth, 2022/11/29)
 </li>

 <li>
  Changed: All SparsityPattern classes now inherit from a new class SparsityPatternBase
  which provides an API for efficiently adding new entries. As a consequence of this change,
  the majority of templates over the sparsity pattern type have been removed.
  <br>
  (David Wells, 2022/11/25)
 </li>

 <li>
  New: Add utility function Utilities::MPI::scatter.
  <br>
  (Magdalena Schreter, Peter Munch, 2022/11/24)
 </li>

 <li>
  Improved: The deal.II CMake project configuration and the
  deal_II_invoke_autopilot() CMake macro that is used in user projects now
  always export the <code>compile_commands.json</code> file used for various
  tooling (such as completion backends like clangd).
  <br>
  (Matthias Maier, 2022/11/24)
 </li>

 <li>
  Improved: Added capability in the particle handler to generate ghost particles
  in ghost cells close to periodic boundary conditions.
  <br>
  (Bruno Blais, Audrey Collard-Daigneault, 2022-11-18)
 </li>

 <li>
  Changed: The constructors of the various solver classes in the PETScWrappers namespace
  used to take an `MPI_Comm` argument. This has now been changed: There
  is a new set of constructors that no longer take this argument, and
  the old constructors have been deprecated. In practice, the solvers
  now use the same MPI communicator as the matrix give to the `solve()`
  function.
  <br>
  (Wolfgang Bangerth, 2022/11/16)
 </li>

 <li>
  New: PreconditionChebyshev now also supports
  James Lottes
  s fourth-kind Chebyshev.
  <br>
  (Martin Kronbichler, Peter Munch, 2022/11/12)
 </li>

 <li>
  Workaround: The testsuite CMake macros now generate a top-level target
  compile_test_executables in user projects that can be used to force the
  compilation of all test executables. This works around an issue with the
  ninja build system that concurrent invocation of tests
  (via <code>ctest -jN</code>) that build executables will fail when calling
  back into cmake. As a workaround it is now possible to first build all
  executables first via <code>ninja compile_test_executables</code> and then
  run tests in parallel with <code>ctest -jN</code>.
  <br>
  (Matthias Maier, 2022/11/11)
 </li>

 <li>
  Improved: FEImmersedSurfaceValues can now be constructed using a MappingFEField.
  <br>
  (Marco Feder, 2022/10/31)
 </li>

 <li>
  New: A new overload for MatrixFree::loop allows
  to specify a pre- and post-operation, similar
  to MatrixFree::cell_loop.
  <br>
  (Sebastian Proell, 2022/10/11)
 </li>

 <li>
  New: We now instruct doxygen to build the documentation with multiple
  threads, substantially reducing the amount of time it takes to
  complete this step.
  <br>
  (Wolfgang Bangerth, 2022/10/04)
 </li>

 <li>
  New: The GridTools::MarchingCubeAlgorithm supports now also 1D.
  <br>
  (Magdalena Schreter, Peter Munch, 2022/09/23)
 </li>

 <li>
  Fixed: On very coarse meshes, with very particular manifolds, e.g., polar manifolds,
  cylindrical manifolds, Manifold<3, 3>::normal_vector may return a zero vector.
  This is fixed now.
  <br>
  (Timo Heister, Luca Heltai, Jiaqi Zhang, 2022/09/20)
 </li>

 <li>
  New: The new function CGALWrappers::compute_intersection_of_cells computes
  the intersection between two (affine) cells starting from the location of the
  vertices. The intersection is described by a vector of arrays, where each array
  identifies a simplex of the sub-tassellation of the intersection.
  <br>
  (Marco Feder, Johannes Heinz, 2022/09/20)
 </li>

 <li>
  Fixed: The function GridTools::internal::distributed_compute_point_locations()
  now projects reference points outside of a cell (but within a tolerance) onto
  the unit cell. This enables the use of FE_Q_iso_Q1 in
  Utilities::MPI::RemotePointEvaluation.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/09/12)
 </li>

 <li>
  Improved: The functions VectorTools::point_values()/VectorTools::point_gradients() and
  the class Utilities::MPI::RemotePointEvaluation now allow to specify the first
  component to be selected. The feature is used in the class DataOutResample
  to output multi-component solutions.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/09/12)
 </li>

 <li>
  Improved: The vector operations in SolverIDR::solve() have been revised for
  performance.  The new implementation should show a speedup in case the
  matrix-vector product and preconditioner are not much more expensive than the
  vector operations. Furthermore, the usage of temporary vectors has been
  reduced by two.
  <br>
  (Martin Kronbichler, 2022/09/12)
 </li>

 <li>
  New: Functions DoFHandler::get_future_fe_indices() and
  DoFHandler::set_future_fe_indices() that store and recover future FE
  indices on a DoFHandler. They can also be used to transfer future FE
  indices between different DoFHandler objects.
  <br>
  (Marc Fehling, 2022/08/20)
 </li>

 <li>
  New: The CellAccessor::as_dof_handler_iterator() function can be used
  to convert from a Triangulation active cell iterator to a DoFHandler active cell
  iterator, or from an active cell iterator of one DoFHandler to that of another
  DoFHandler.
  <br>
  (Jean-Paul Pelteret, 2022/07/15)
 </li>

 <li>
  Improved: Refinement and coarsening flags are now
  communicated in parallel::shared::Triangulation.
  <br>
  (Peter Munch, 2022/07/12)
 </li>

 <li>
  New: The new class MatrixFreeTools::ElementActivationAndDeactivationMatrixFree
  is a wrapper around MatrixFree designed to deal with
  DoFHandler objects involving cells without degrees of freedom, i.e.,
  cells using FE_Nothing as element type. The class helps to implement the
  "element birth and death technique".
  <br>
  (Peter Munch, Sebastian Proell, 2022/07/09)
 </li>

 <li>
  New: Implement `DataOut::write_deal_II_intermediate_in_parallel()`
  that writes a combined file from all MPI ranks using MPI I/O of the
  internal patches and a corresponding
  `DataOutReader::read_whole_parallel_file()` to read them back in.
  <br>
  (Timo Heister, 2022/07/07)
 </li>

 <li>
  Improved: The function GridGenerator::merge_triangulations() can now copy boundary ids to remaining boundary faces.
  <br>
  (Johannes Heinz, 2022/07/07)
 </li>

 <li>
  New: There are now functions Utilities::MPI::isend() and
  Utilities::MPI::irecv() that can send and receive arbitrary objects to
  and from other processes, and that handle all of the setting up and
  tearing down of MPI objects and temporary buffers.
  <br>
  (Wolfgang Bangerth, 2022/07/07)
 </li>

 <li>
  Improved: The function VectorTools::point_values() can now handle
  cell-data vectors.
  <br>
  (Peter Munch, 2022/07/06)
 </li>

 <li>
  Added the missing explicit instantiations of DoFTools::distribute_cell_to_dof_vector
  for problems with non-zero codim (dim != spacedim)
  <br>
  (Ahmad Shahba, 2022/07/06)
 </li>

 <li>
  Fixed: MappingCartesian used to throw exceptions in debug mode when checking if
  small cells in a large domain size are actually aligned to coordinate axes.
  This is fixed now.
  <br>
  (Rene Gassmoeller, 2022/07/04)
 </li>

 <li>
  Fixed: The function find_active_cell_around_point() now skips cells without any marked vertices.
  <br>
  (Johannes Heinz, 2022/06/28)
 </li>

 <li>
  New: The PetscWrappers::PreconditionBDDC preconditioner has been implemented. It requires a special matrix type (IS), which is initialized with the new function PetscWrappers::reinit.
  <br>
  (Nicolas Barnafi, 2022/05/10)
 </li>

</ol>

*/
