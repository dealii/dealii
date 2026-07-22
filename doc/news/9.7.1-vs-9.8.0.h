// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
@page changes_between_9_7_1_and_9_8_0 Changes between Version 9.7.1 and 9.8.0

<p>
This is the list of changes made between the release of deal.II version
9.7.1 and that of 9.8.0. All entries are signed with the names of the
author.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="971-980-incompatible"></a>
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
  Changed: GridIn::read_partitioned_msh(), the overload of GridIn::read_msh()
  which uses the GMSH API, and the overload of GridOut::write_msh() which uses the
  GMSH API now require that gmsh::initialize() is run prior to calling those
  functions. This change was made to ensure consistent initialization and
  finalization of gmsh as well as its dependencies (such as MPI).
  <br>
  (David Wells, 2026/06/19)
 </li>

 <li>
  Changed: TriangulationDescription::CellData::boundary_ids is now a
  std_cxx26::inplace_vector instead of a std::vector. This change makes this
  class, when compiling with C++20 or newer, trivially copyable.
  <br>
  (David Wells, 2026/06/15)
 </li>

 <li>
  Changed: The `mapping_q.h` header no longer includes `quadrature_lib.h`.
  <br>
  (David Wells, 2026/06/08)
 </li>

 <li>
  Changed: CellData::vertices is now a std_cxx26::inplace_vector instead of a
  `std::vector`. This change was made to improve performance and decrease the
  total number of memory allocations, as the maximum number of vertices per cell
  is known at compile time. Some code, like assigning `vertices` to a
  `std::vector`, will no longer compile, but most constructions (like the
  `assign()` member function, assigning to a `std::initializer_list`, etc.) will
  continue to work. In particular, pre-2021 code (before we switched to
  `std::vector` from a C array) will not require updates.
  <br>
  (David Wells, 2026/05/26)
 </li>

 <li>
  Changed: The maximum number of levels in a Triangulation
  is now capped at 31. This is the same limit used by p4est,
  and was previously implicitly enforced through the CellId
  class.
  <br>
  (David Wells, 2026/05/26)
 </li>

 <li>
  Changed: Increased the maximum number of children after refinement
  (`ReferenceCells::max_n_children()`) from 8 to 10 to accommodate the isotropic
  refinement of pyramids.
  <br>
  (Andreas Steger, 2026/05/01)
 </li>

 <li>
  Changed: FESubfaceValues::reinit() now checks whether the face under
  consideration actually belongs to a cell whose neighbor is more refined. This
  was previously only enforced on cell iterators based on a Triangulation
  object, not those associated to a DoFHandler object. The new behavior is
  now uniform in both cases.
  <br>
  (Martin Kronbichler, 2026/04/17)
 </li>

 <li>
  Changed: There are some functions in deal.II that compute the curl of
  vector fields. The curl is only defined for vector fields in ${\mathbb
  R}^3$, or for two-dimensional vector fields on two-dimensional
  manifolds (perhaps embedded in ${\mathbb R}^3$). Historically, it was
  possible to also call these functions for `dim==1`, which made no
  sense but was syntactically allowed. However, this is now forbidden
  and will lead to syntax errors.
  <br>
  (Wolfgang Bangerth, 2026/04/13)
 </li>

 <li>
  Changed: In three space dimensions, the curl of a vector field is unambiguously
  a three-dimensional vector. In two space dimensions, it is typically defined
  as a scalar, interpreted as the magnitude of a vector aligned with the $z$-axis
  the results from taking the curl of a vector field that lives in the $x-y$ plane.
  For historical reasons, deal.II used to represent this scalar as a Tensor<1,1>
  (i.e., a tensor of rank one -- a tensor with one index -- for which the only
  possible value of the index is zero). This is, of course, entirely equivalent
  to it being a scalar, but it is semantically wrong: The data type should have
  been a scalar, not an array of which one can take the zeroth element. This is
  now fixed: The data type defined by FEValuesViews::Vector::curl_type used for
  curls is not simply a scalar value for `dim==2`. For `dim==3`, it is
  unchanged from what it was before.
  <br>
  (Wolfgang Bangerth, 2026/03/31)
 </li>

 <li>
  Removed: The function TrilinosWrappers::Vector::reinit with the option to copy values
  from a differently partitioned vector has been removed due to inconsistent behavior
  with the other reinit functions. If you require this functionality, use the
  copy-constructor specifying the parallel partitioning, or the equivalent
  assignment operator.
  <br>
  (Rene Gassmoeller, 2026/03/30)
 </li>

 <li>
  Deprecated: Overloads for TrilinosWrappers classes allowing a default number of
  entries per row have been deprecated. The respective overloads for TpetraWrappers
  classes have been removed.
  <br>
  (Daniel Arndt, 2026/03/25)
 </li>

 <li>
  Changed: If Trilinos was built with Tpetra support but without Epetra support,
  the TrilinosWrappers classes are aliases for the respective TpetraWrappers classes
  with NumberType=double and MemorySpace=Host.
  <br>
  (Daniel Arndt, 2026/03/21)
 </li>

 <li>
  Deprecated: All classes in the EpetraWrappers namespace have been deprecated.
  <br>
  (Daniel Arndt, 2026/01/27)
 </li>

 <li>
  Changed: Most of the parallel vector classes (in particular those
  wrapping the functionality in the PETSc, Trilinos Epetra, and Trilinos
  Tpetra packages) observe semantics whereby vectors that have ghost
  elements need to be treated as unchangeable. This is because whenever
  you change a locally-owned element of such a vector, this change is
  not propagated to those processes that store this element among their
  ghost elements. As a consequence, these processes view of the vector
  is now no longer in sync with the view on that process that owns the
  element, and that is a common source of very hard to find bugs.
  <br>
  As a consequence, many functions inside these vector classes that
  modify the vector check that the vector has no ghost elements. But
  some did not, such as TrilinosWrappers::MPI::Vector::operator+=(),
  TrilinosWrappers::MPI::Vector::operator-=(),
  TrilinosWrappers::MPI::Vector::operator*=(),
  TrilinosWrappers::MPI::Vector::operator/=(), and the corresponding
  functions in the Epetra and Tpetra wrappers, along with an assortment
  of other non-`const` functions in these classes that include assigning
  a scalar to all elements of a vector, and writing into (or adding,
  or multiplying/dividing) a single vector entry.
  <br>
  This oversight has now been fixed. We know that some of these
  functions are used in parallel contexts in application codes that will
  now no longer work and will need to be fixed. For this, you will want
  to perform the modifying operations on full-distributed vectors (i.e.,
  vectors of the same type that store no ghost elements) and then copy
  the result into the ghosted vector.
  <br>
  (Wolfgang Bangerth, 2025/12/15)
 </li>

 <li>
  Changed: The class Rol::VectorAdaptor has been renamed TrilinosWrappers::ROLAdaptor.
  Its corresponding header file has been moved from `deal.II/optimization/rol/vector_adaptor.h`
  to `deal.II/trilinos/rol_adaptor.h`.
  <br>
  (Marc Fehling, 2025/11/12)
 </li>

 <li>
  Changed: A template parameter for dim has been added to MGTwoLevelTransferBase.
  <br>
  (Ryan Moulday, 2025/11/10)
 </li>

 <li>
  Changed: The function Timer::stop() no longer returns the CPU time after
  stopping the timer. This functionality is covered by the existing function
  Timer::cpu_time(), which should be called instead.
  <br>
  (Rene Gassmoeller, 2025/10/04)
 </li>

 <li>
  Changed: The TimerOutput class no longer supports producing output while actively timing a subsection.
  This feature was removed to simplify internal data storage and to avoid inconsistencies in parallel computations.
  <br>
  (Rene Gassmoeller, 2025/10/02)
 </li>

 <li>
  Changed: Renamed FETools::cell_to_face_patch() to FETools::cell_to_face_patch_numbering().
  FETools::cell_to_face_patch() was marked deprecated.
  <br>
  (Michał Wichrowski, 2025/09/15)
 </li>

 <li>
  Removed: The classes MGTransferMF, MGTransferGlobalCoarsening, MGTransferBlockMF, and MGTransferBlockGlobalCoarsening
  have been deprecated. Their functionalities are now part of MGTransferMatrixFree and MGTransferBlockMatrixFree.
  <br>
  (Peter Munch, 2025/09/12)
 </li>

 <li>
  Removed: Deprecations from the 9.6 release have been removed.
  See https://github.com/dealii/dealii/pull/18709 for details.
  <br>
  (Daniel Arndt, 2025/07/17)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="971-980-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: Portable::MGTwoLevelTransfer now supports FESystem.
  <br>
  (Ivan Prusak, 2026/07/09)
 </li>

 <li>
  New: step-104 is a new tutorial on solving the Stokes equations matrix-free on GPUs.
  <br>
  (Quang Hoang, Timo Heister, 2026/07/01)
 </li>

 <li>
  New: Added Portable::MGTwoLevelTransfer class supporting geometric multigrid transfer operators on GPUs.
  At the moment, only globally refined meshes and FE_Q are supported.
  <br>
  (Ivan Prusak, Martin Kronbichler, 2026/06/23)
 </li>

 <li>
  Deprecated: deal.II has had a long list of classes in namespace
  MeshWorker that were developed decades ago when we had a different
  programming style and that were intended to hide some of the
  complexity of assembling linear systems as well as evaluating and
  postprocessing solutions. Many parts of the MeshWorker framework
  worked well for a specific purpose, but turned out to be very
  difficult to adapt to other tasks -- say, converting the assembly of a
  term to a term that nonlinearly depends on coefficients that
  themselves depend on the current solution. The framework also dates
  back to a time when deal.II had substantially poorer test coverage and
  documentation, and as a consequence corner cases were often found to
  not actually work as expected, and/or not well documented.
  As a consequence, we have over time developed other, better, better
  documented, and better tested tools for all of these tasks. The only
  remaining hold-out that used the old MeshWorker framework was step-39,
  and with the conversion of step-39 from the old MeshWorker::loop()
  function to the newer MeshWorker::mesh_loop() function, many of the
  old classes have become unused in the tutorial and have been
  deprecated (and will consequently be removed in the second-next
  release). Specifically, these classes and functions are:
  - MeshWorker::DoFInfo
  - MeshWorker::IntegrationInfo
  - MeshWorker::IntegrationInfoBox
  - MeshWorker::loop().
  <br>
  This has also led to many other classes that were previously only used
  within the classes and functions mentioned above, and that are now
  also deprecated:
  - MeshWorker::Assembler::Functional
  - MeshWorker::Assembler::ResidualLocalBlocksToGlobalBlocks
  - MeshWorker::Assembler::MatrixLocalBlocksToGlobalBlocks
  - MeshWorker::Assembler::MGMatrixLocalBlocksToGlobalBlocks
  - MeshWorker::Assembler::ResidualSimple
  - MeshWorker::Assembler::MatrixSimple
  - MeshWorker::Assembler::MGMatrixSimple
  - MeshWorker::Assembler::SystemSimple
  - MeshWorker::Assembler::CellsAndFaces
  - MeshWorker::Assembler::GnuplotPatch
  - MeshWorker::VectorSelector
  - MeshWorker::VectorDataBase
  - MeshWorker::VectorData
  - MeshWorker::MGVectorData
  - MeshWorker::LocalResults
  - MeshWorker::LoopControl
  - MeshWorker::cell_action()
  <br>
  All of these classes will disappear in a future release.
  <br>
  (Wolfgang Bangerth, 2026/06/16)
 </li>

 <li>
  New: The new time-steppers have been added for SUNDIALS::ARKode:
  SUNDIALS::ERKStepper wrapping the ERKStep (Explicit Runge-Kutta)
  module of ARKODE, SUNDIALS::LSKRStepperSSP and SUNDIALS::LSKRStepperSTS
  providing access to the SSP (Strong-Stability-Preserving) and STS
  (Super Time-Stepping) time integration methods within the LSRKStep
  (Low-Storage Runge-Kutta) module of ARKODE, respectively.
  <br>
  (Vladimir Ivannikov, 2026/04/21)
 </li>

 <li>
  New: Portable::MatrixFree now supports systems with more than DoFHandler.
  <br>
  (Timo Heister, 2026/06/01)
 </li>

 <li>
  New: deal.II has a C++17-compliant version of most of C++26's new inplace_vector
  class template.
  <br>
  (David Wells, 2026/05/16)
 </li>

 <li>
  Reworked `execute_refinement_isotropic` to handle cells producing children of
  different reference cells and added support for isotropic pyramid refinement.
  <br>
  (Andreas Steger, 2026/05/01)
 </li>

 <li>
  New/Changed: The SUNDIALS ARKODE wrapper has been substantially
  extended and refactored. The time-stepper functionality previously
  embedded in SUNDIALS::ARKode has been extracted into the new
  SUNDIALS::ARKodeStepper interface. A concrete implementation is
  currently provided for the SUNDIALS implicit/explicit ARKStep
  integration module (SUNDIALS::ARKStepper), but the design is open to
  further ARKODE time-steppers. A stepper object is constructed
  independently and passed to the SUNDIALS::ARKode driver, which handles
  the time-marching loop. Method-specific parameters are configured
  through per-class AdditionalData structures.
  <br>
  The previous interface of SUNDIALS::ARKode that exposed the
  right-hand side functions directly on the ARKode object is still
  available but is now deprecated and will be removed in a future
  version. Users are encouraged to migrate to the new stepper-based
  interface.
  <br>
  (Vladimir Ivannikov, 2026/04/21)
 </li>

 <li>
  New: The step-98 tutorial program illustrates application of the
  FE_Q, FE_Nedelec, FE_RaviartThomas, and FE_DGQ finite elements to
  two-dimensional problems in magnetostatics.
  <br>
  (Siarhei Uzunbajakau, 2026/04/10)
 </li>

 <li>
  Reworked `execute_refinement_isotropic` to support cells with multiple face
  types. Additionally, support for the isotropic refinement of wedges was added.
  <br>
  (Andreas Steger, 2026/01/20)
 </li>

 <li>
  DataOut has been generalized to correctly handle parallel vectors with
  arbitrary ownership layouts, improving robustness and correctness for
  distributed output.
  <br>
  (Mathias Anselmann, 2026/01/11)
 </li>

 <li>
  New: Matrix-free infrastructure for evaluation/integration of Nedelec elements has been developed.
  <br>
  (Natalia Nebulishvili, Martin Kronbichler, 2025/12/30)
 </li>

 <li>
  New: PETSc linear algebra now supports multi-threaded assembly and reading from ghost vectors.
  <br>
  (Timo Heister, 2025/11/23)
 </li>

 <li>
  Improved: The hp adaptivity features for DG of MatrixFree have
  been extended. For face-centric loops, neighboring cells can have
  quadrature rules with different number of quadrature points. On
  each face, the quadrature rule with the highest number of quadrature
  points is used. For element-centric loops, p-adaptivity is enabled. It is
  experimental, since categorization and hp-enabled DoFHandler are not
  supported at the same time by MatrixFree, limiting, in most cases, the
  usablity to VectorizedArray<double, 1>.
  <br>
  (Peter Munch, 2025/10/21)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="971-980-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: TriangulationDescription::Description now correctly sets material_ids
  of cells on finer levels.
  <br>
  (David Wells, 2026/07/13)
 </li>

 <li>
  Fixed: SolverGMRES and SolverFGMRES with orthogonalization strategy
  LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt (the
  default case) could run into a division by zero in the rare case when residual
  vectors became exactly zero. This invalid operation is now avoided, and the
  lucky breakdown is handled gracefully.
  <br>
  (Martin Kronbichler, 2026/06/30)
 </li>

 <li>
  Added: get_divergence() and submit_divergence() have been added to Portable::FEEvaluation.
  <br>
  (Ryan Moulday, 2026/06/29)
 </li>

 <li>
  Fixed: DoFTools::extract_dofs_with_support_contained_within() previously returned
  DoFs as contained within when running in parallel that are actually not contained
  within.
  <br>
  (Torsten Schmid, 2026/06/23)
 </li>

 <li>
  New: parallel::CellWeights can now work with pre-computed weights.
  This requires that the weights are determined only by information
  about finite elements.
  <br>
  (Marc Fehling, 2026/06/19)
 </li>

 <li>
  Improved: The matrix-free framework now stores Dirichlet-constrained indices
  as index numbers::invalid_unsigned_int, rather than as a separate constraint
  that needs to be resolved from a list. This leads to slightly better machine
  code on most machines, thus speeding up matrix-free codes.
  <br>
  (Martin Kronbichler, 2026/06/18)
 </li>

 <li>
  Improved: The functions VectorizedArray::gather() and
  VectorizedArray::scatter() are now able to filter out indices with value
  numbers::invalid_unsigned_int, indicating a mask for the vector access.
  This is useful for indirect index access.
  <br>
  (Martin Kronbichler, 2026/06/18)
 </li>

 <li>
  Improved: Functions in Utilities::MPI no longer compress and decompress their input
  arguments. This makes, e.g., Utilities::MPI::gather() about twice as fast.
  <br>
  (David Wells, 2026/06/12)
 </li>

 <li>
  New: diagonal_operator() function returns a LinearOperator
  representing the action of a diagonal matrix
  with entries given by a vector.
  <br>
  (Quang Hoang, 2026/06/12)
 </li>

 <li>
  New: Added a `DoFRenumbering::component_wise()` overload that accepts FEValuesExtractors, where the ordering is
  defined by the order of the passed extractors. Convenience alias `AnyExtractor` is introduced for a variant type holding
  any of the supported `FEValuesExtractors` (scalar, vector, tensor, symmetric tensor).
  <br>
  (Gordan Segon, 2026/06/11)
 </li>

 <li>
  Changed: The step-39 tutorial program has been rewritten to use
  MeshWorker::mesh_loop() instead of the removed MeshWorker::loop(),
  MeshWorker::DoFInfo, MeshWorker::IntegrationInfo, and
  MeshWorker::IntegrationInfoBox interfaces.
  <br>
  (Wolfgang Bangerth, 2026/06/11)
 </li>

 <li>
  Deprecated: The class MeshWorker::Assembler::Functional has been marked as
  deprecated and will be removed in a future release.
  <br>
  (Wolfgang Bangerth, 2026/06/10)
 </li>

 <li>
  Deprecated: The classes Algorithms::Event, Algorithms::OperatorBase,
  Algorithms::OutputOperator, Algorithms::Newton,
  Algorithms::ThetaTimestepping, Algorithms::TimestepControl, and
  Algorithms::DoFOutputOperator have been marked as deprecated.
  <br>
  (Wolfgang Bangerth, 2026/06/10)
 </li>

 <li>
  New: Added output in NetCDF using CF/UGRID conventions to DataOut.
  <br>
  (Daniel Abele, 2026/06/09)
 </li>

 <li>
  Improved: The function vectorized_load_and_transpose() now produces more
  efficient remainder code for sizes 1, 2, 3 on x86.
  <br>
  (Martin Kronbichler, 2026/06/09)
 </li>

 <li>
  New: Add NetCDF as an optional dependency (to be used for output in an upcoming feature).
  <br>
  (Daniel Abele, 2026/06/08)
 </li>

 <li>
  New: Add new preconditioner class: PreconditionFROSch,
  which adds a wrapper for the FROSch preconditioner, that is part of Trilinos.
  <br>
  (Sebastian Kinnewig, 2026/06/07)
 </li>

 <li>
  New: Extended the TimeStepping namespace by the adaptive Tsitouras 5(4) explicit Runge-Kutta method (TimeStepping::TSITOURAS5).
  <br>
  (Jalil Khan, 2026/06/06)
 </li>

 <li>
  Improved: FilteredIterator now re-uses the provided predicate object
  between different iterators instead of copying it (which required
  an allocation). This lowers the number of allocations in step-32 by
  about 15\%.
  <br>
  (David Wells, 2026/05/25)
 </li>

 <li>
  Fixed: FEEvaluation::distribute_local_to_global() would produce a segmentation
  fault if initialized with DoFHandler::active_cell_iterator and used with
  LinearAlgebra::distributed::Vector. This is now fixed.
  <br>
  (Martin Kronbichler, 2026/05/22)
 </li>

 <li>
  New: Improve VTK dataset import/export utilities by supporting conversion of
  generic VTK datasets to unstructured grids while preserving attached data
  arrays, adding XML `.vtu` reading support, and propagating optional material,
  boundary, and manifold id field names through `read_tria()` and `read_vtk()`.
  <br>
  (Luca Heltai, 2026/05/21)
 </li>

 <li>
  Fixed: The SolverFIRE class assumed that the matrix object has a
  `residual()` function, but this function is not part of the required
  interface for all matrix types that can be passed to linear
  solvers. This is now fixed by instead relying on `vmult()`
  functionality.
  <br>
  (Wolfgang Bangerth, 2026/05/21)
 </li>

 <li>
  Fixed: Fix multigrid W- and F-cycles. Previously, intermediate
  solutions were overwritten between cycles, effectively resulting in a V-cycle.
  <br>
  (Peter Munch, 2026/05/18)
 </li>

 <li>
  Improved: The arrays inside a Triangulation now use strides which depend on the
  present ReferenceCell objects, rather than the maximum values over all possible
  ReferenceCell objects. For example, a Triangulation containing only tetrahedra
  now only stores four face orientations per cell instead of six. This lowers
  memory consumption for tetrahedral Triangulations by about 20\%.
  <br>
  (David Wells, 2026/05/16)
 </li>

 <li>
  New: Extended the PETScWrappers::MPI namespace by a full matrix
  (PETScWrappers::MPI::FullMatrix).
  <br>
  (Jonas Plank, 2026/05/15)
 </li>

 <li>
  Moved: The function ReferenceCell::n_vertices_to_type() has been deprecated.
  Its replacement is ReferenceCells::n_vertices_to_reference_cell().
  <br>
  (Wolfgang Bangerth, 2026/05/12)
 </li>

 <li>
  Fixed: The default mapping provided to FEFieldFunction only worked for
  hypercube meshes. This is now fixed by introducing a second
  constructor that determines what the appropriate linear mapping for a
  given mesh is, and then using this mapping.
  <br>
  (Wolfgang Bangerth, Sam Scheuerman, David Wells, 2026/05/12)
 </li>

 <li>
  Improved: DynamicSparsityPattern's internal allocation scheme was rewritten to
  use about 50\% fewer allocations and be about 10\% faster for scalar problems.
  In some cases, such as FESystem<3>(FE_Q<3>(2), 3), the new allocation scheme is
  about 50\% faster. This scheme uses about 50\% more memory, but since the
  creation of a DynamicSparsityPattern is usually followed by creation of a
  SparsityPattern and a SparseMatrix, each of which is about the same size as the
  original DynamicSparsityPattern, the maximum memory usage of a process is the
  the same and hence the added memory cost of the DynamicSparsityPattern is
  negligible.
  <br>
  (David Wells, 2026/05/06)
 </li>

 <li>
  Improved: TimeStepping::ExplicitRungeKutta now also supports
  TimeStepping::HEUN_EULER.
  <br>
  (Peter Munch, 2026/05/05)
 </li>

 <li>
  Improved: VTU output is now about 5\% faster on large meshes
  and no longer allocates memory for every processed cell.
  <br>
  (David Wells, 2026/05/03)
 </li>

 <li>
  New: Add a new signed distance function for an infinite cylinder
  (Functions::SignedDistance::InfiniteCylinder). It supports calculation of its
  values and gradients.
  <br>
  (Magdalena Schreter-Fleischhacker, Julian Brotz, Bruno Blais, 2026/04/30)
 </li>

 <li>
  New: There is now a version of GridGenerator::cheese() that allows the
  creation of meshes based on true/false masks for pixelized/voxelized
  data.
  <br>
  (Wolfgang Bangerth, 2026/04/28)
 </li>

 <li>
  Fixed: Despite its name, the GridGenerator::cheese() function in 3d
  did not actually create any holes. It just produced a solid
  volume. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2026/04/28)
 </li>

 <li>
  Fixed: TimerOutput::print_wall_time_statistics() no longer requires all processes
  to enter the same sections when the TimerOutput is constructed without an MPI
  communicator. Sections that are not entered on a process are now reported with
  zero timing. Previously, this could trigger an assertion or lead to a deadlock.
  <br>
  (Torsten Schmid, 2026/04/20)
 </li>

 <li>
  Fixed: When configuring for using both MPI and VTK, then VTK's CMake scripts
  unset information about the MPI version. As a consequence, the
  variables `DEAL_II_MPI_VERSION_MAJOR` and `DEAL_II_MPI_VERSION_MINOR`
  were left blank, and we always fell back to using MPI 3 functionality
  even if the MPI installation would have supported newer MPI functions.
  <br>
  (Wolfgang Bangerth, Matthias Maier, 2026/04/13)
 </li>

 <li>
  Improved: The order by which a MatrixFree object passes through cells is now
  using the order implied by the categorization specified by
  MatrixFree::AdditionalData::cell_vectorization_category, if the number of
  categories equals the number of locally owned active cells in case the flag
  MatrixFree::AdditionalData::cell_vectorization_categories_strict is set to
  false, or if at least a quarter of cells is given for the flag set to true. By
  setting a one-to-one numbering of cells, one can now inject a fixed ordering
  of cells.
  <br>
  (Martin Kronbichler, 2026/04/07)
 </li>

 <li>
  Fixed: The function FEValuesViews::Vector::curl() computed the wrong
  value for `dim==2` and `spacedim==3`, i.e., for the case of a two-dimensional
  manifold embedded in 3d. In that case, it assumed that we are dealing with
  a two-dimensional vector field, for which the curl is a scalar, when what is
  meant is the curl of a three-dimensional vector field defined on a surface.
  <br>
  (Wolfgang Bangerth, 2026/03/18)
 </li>

 <li>
  Fixed: Initializing FEEvaluation with a Mapping, FiniteElement, and Quadrature
  object would leave some data pointers unset, leading to problems e.g. when
  evaluating elements of type FE_RaviartThomasNodal. This is now fixed.
  <br>
  (Martin Kronbichler, 2026/03/30)
 </li>

 <li>
  Extended: parallel::distributed::Triangulation::copy_triangulation() now allows to
  control whether to copy also the triangulation settings or to use a developer provided
  setting.
  This can be useful for multigrid computations.
  (Vivienne Ehlert, Martin Kronbichler, 2026/03/27)
 </li>

 <li>
  Changed: Serialization of Table and AlignedVector objects storing
  "trivial" data types (things such as `double` or `int`) now do this by
  serializing a binary blob of bits, rather than individual array
  members. This has no effect for binary archives (as one should use for
  checkpointing, for example), but makes text-based archives
  substantially smaller.
  <br>
  (Wolfgang Bangerth, 2026/03/20)
 </li>

 <li>
  Fixed: If FEFaceValues::reinit() and FESubfaceValues::reinit() were used on
  objects that remained valid after a change in the triangulation and used the
  same cell, they could sometimes produce wrong results for the metric
  terms. This is now fixed.
  <br>
  (Martin Kronbichler, 2026/03/19)
 </li>

 <li>
  Fixed: The function AlignedVector::load() function did not work
  correctly if the AlignedVector object on which it was called had
  non-zero size and the object being deserialized had zero size.
  <br>
  (Wolfgang Bangerth, 2026/03/18)
 </li>

 <li>
  Improved: Triangulation::pack_data_serial() previously had quadratic complexity
  in the preparation of cell relations, which had a negative impact on runtime
  when refining a serial Triangulation in combination with SolutionTransfer.
  This has been improved.
  <br>
  (Peter Munch, Martin Kronbichler, Lukáš Lejdar, 2026/03/17)
 </li>

 <li>
  Fixed: Deserialization of an AlignedVector now works correctly with types
  which are not trivially copyable.
  <br>
  (David Wells, 2206/03/15)
 </li>

 <li>
  Fixed: When using the PETSc solver through PETScWrappers,
  the original code will use the l2 norm of the preconditioned residual.
  This is not consistent with deal.II expectations and example usage.
  This PR fixes the issue by using the unpreconditioned residual norm.
  <br>
  (Yiliang Wang, 2026/03/05)
 </li>

 <li>
  Fixed: The function AlignedVector::replicate_across_communicator() could not deal with tables that took up more than about 2GB because of an integer overflow. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2026/02/25)
 </li>

 <li>
  Fixed: NonMatching::MappingInfo::reinit() would leak a small amount of memory in each call. This is now fixed. This indirectly affects the usage of FEPointEvaluation.
  <br>
  (Timo Heister, 2026/02/24)
 </li>

 <li>
  Improved: `CellAccessor::child_iterators()` return vector size now depends on
  the `constexpr` function `ReferenceCells::max_n_children<dim>()`.
  <br>
  (Andreas Steger, 2026/02/23)
 </li>

 <li>
  New: VTKWrappers::unstructured_grid_to_dealii_triangulation() and
  VTKWrappers::dealii_triangulation_to_unstructured_grid() allow translations of
  a Triangulation from and to vtkUnstructuredGrid.
  <br>
  (Luca Heltai, 2026/02/19)
 </li>

 <li>
  Fixed: FEPointEvaluation::get_gradient() would return wrong results if a mesh
  consisted of both non-affine and affine elements and several points were
  evaluated on the cells of affine shape. This also affected the step-19
  tutorial program. This is now fixed.
  <br>
  (Martin Kronbichler, Zdeněk Bonaventura, 2026/02/13)
 </li>

 <li>
  New: A new class DataOutPoints to produce graphical
  output in a set of points specified by the user.
  <br>
  (Timo Heister, 2026/02/04)
 </li>

 <li>
  New: RemotePointEvaluation::CellData::get_data_view() now
  supports cases with more than one component via an optional
  third argument.
  <br>
  (Timo Heister, 2026/01/28)
 </li>

 <li>
  Fixed: Various minor issues with algorithms operating
  on large problems with DEAL_II_WITH_64BIT_INDICES=ON have
  been addressed.
  <br>
  (Timo Heister, 2026/01/27)
 </li>

 <li>
  Fixed: The function VectorTools::project_boundary_values_div_conforming()
  would previously allocate an array of length equal to the global number of
  degrees of freedom. This is now fixed, making sure that also large-scale
  parallel programs can fit into memory.
  <br>
  (Martin Kronbichler, 2026/01/26)
 </li>

 <li>
  New: step-7 now has extensive sections on comparing numerical methods via
  cost-accuracy diagrams.
  <br>
  (Wolfgang Bangerth, 2026/01/26)
 </li>

 <li>
  Improved: Added an overload of DoFAccessor::get_interpolated_dof_values() which
  takes an ArrayView, instead of a Vector, as an argument. This function enables
  callers like FEValues to avoid allocating memory for many low-order elements.
  <br>
  (David Wells, 2026/01/25)
 </li>

 <li>
  Fixed: MappingFEField::get_vertices() correctly works with threads again.
  <br>
  (David Wells, 2026/01/25)
 </li>

 <li>
  New: The GridGenerator::uniform_channel_with_sphere() function
  creates the three-dimensional mesh of a channel with a sphere inside. The mesh
  created is conformal.
  <br>
  (Quentin Viville, 2026/01/24)
 </li>

 <li>
  New: LinearAlgebra::TpetraWrappers::PreconditionAMGMueLu has been added as
  preconditioner for TpetraWrappers::SparseMatrix.
  <br>
  (Jan Philipp Thiele, Daniel Arndt, 2026/01/08)
 </li>

 <li>
  New: VTKWrappers::read_vtk() allows to read all content of a vtk file (grid and
  data) into a Triangulation, DoFHandler, and Vector.
  <br>
  (Luca Heltai, 2025/12/30)
 </li>

 <li>
  New: Make NonMatching::MeshClassifier handle the case where the zero contour of
  a level set function aligns perfectly with a face in the mesh, by adding a new
  entry, "aligned", in the enum LocationToLevelSet.
  <br>
  (Simon Sticko, Magdalena Schreter-Fleischhacker, 2025/12/25)
 </li>

 <li>
  New: The new function VTKWrappers::read_all_data() loads a VTK unstructured
  grid file (optionally cleaning overlapping points) and concatenates every point
  and cell data array it finds into the provided output vector.
  <br>
  (Luca Heltai, 2025/12/24)
 </li>

 <li>
  New: The function GridTools::parallel_to_serial_vertex_indices() returns a
  vector mapping each vertex of a parallel Triangulation to the corresponding
  vertex index in the originating serial Triangulation.
  <br>
  (Luca Heltai, 2025/12/23)
 </li>

 <li>
  New: Implemented a new RTreeFunctionalVisitor class and corresponding
  visit_rtree() function that allow to visit (read only) a boost RTree data
  structure, executing a user provided function for each node, leaf, or element
  of the tree.
  <br>
  (Luca Heltai, 2025/12/19)
 </li>

 <li>
  Changed: step-42 used to create output of the solid object in the
  deformed configuration by deforming the mesh with the solution before
  output, creating output in the deformed configuration, and then
  undoing the mesh deformation. But we can do all of this without
  deforming the mesh at all: The MappingQEulerian class can do that more
  easily, and step-42 has been changed to use this approach now.
  <br>
  (Wolfgang Bangerth, 2025/12/17)
 </li>

 <li>
  Add: Add support for LinearAlgebra::distributed::Vector to the SparseDirectMUMPS
  class.
  <br>
  (Marco Feder, 2025/12/14)
 </li>

 <li>
  Added: FE_Hermite has now a new function
  FE_Hermite::get_dofs_corresponding_to_outward_normal_derivatives() that
  returns those cell unknowns that have a particular derivative index nonzero on
  all faces, which helps with the implementation of boundary conditions.
  <br>
  (Ivy Weber, Martin Kronbichler, 2025/12/11)
 </li>

 <li>
  Enable Taskflow-based parallelization of parallel::apply_to_subranges().
  This functionality is missed in the current code, which only supports
  TBB and sequential execution.
  <br>
  (Qingyuan Shi, 2025/12/10)
 </li>

 <li>
  New: The new function VTKWrappers::vtk_to_finite_element() generates a
  FiniteElement object from a VTK supported file. The generated FiniteElement can
  be used to create a field that is compatible with the field data contained in
  the file.
  <br>
  (Luca Heltai, 2025/12/09)
 </li>

 <li>
  Improved: The Taskflow based WorkStream::run using graph coloring has been
  improved to be significantly more performant in the cases where work units are
  very cheap.
  <br>
  (Ryan Moulday, Timo Heister, 2025/12/08)
 </li>

 <li>
  New: A nodal variant of Nedelec elements, FE_NedelecNodal, has been added.
  <br>
  (Natalia Nebulishvili, Martin Kronbichler, Katharina Kormann, 2025/12/05)
 </li>

 <li>
  Add: Introduce a fully vectorizable function to approximate the exponential
  function. The function offers faster computation than std::exp while maintaining
  high accuracy and supports SIMD/vectorized operations.
  <br>
  (Julian Brotz, 2025/12/02)
 </li>

 <li>
  New: A new overload of MatrixFree::initialize_dof_vector() for BlockVector creates a block
  for each DoFHandler passed to MatrixFree::reinit().
  <br>
  (Timo Heister, 2025/11/29)
 </li>

 <li>
  New: The function VTKWrappers::read_cell_data() allows to read cell data contained in a VTK file using the VTK library.
  <br>
  (Luca Heltai, 2025/11/25)
 </li>

 <li>
  New: The function VTKWrappers::read_vertex_data() allows to read vertex data
  from vtk files using the VTK library.
  <br>
  (Luca Heltai, 2025/11/25)
 </li>

 <li>
  New: VTKWrappers::read_vtk() can now read any unstructured vtk file using the
  native VTK library. This includes binary vtk files, and any format supported by
  vtkUnstructuredGridReader.
  <br>
  (Luca Heltai, 2025/11/25)
 </li>

 <li>
  Fix: Shorten documentation of ParsedFunction that is printed to parameter files.
  <br>
  (Luca Heltai, 2025/11/21)
 </li>

 <li>
  New: You can now mask vector elements in TrilinosWrappers::ROLAdaptor to
  reduce the size of the optimization space. This is useful for excluding
  constrained degrees of freedom from the optimization process.
  <br>
  (Marc Fehling, 2025/11/20)
 </li>

 <li>
  Fixed: TrilinosWrappers::ROLAdaptor::reduce() now also works on parallel.
  <br>
  (Marc Fehling, 2025/11/19)
 </li>

 <li>
  New: In Portable::MatrixFree the method get_cell_iterator() has been added.
  <br>
  (Michał Wichrowski, 2025/11/11)
 </li>

 <li>
  Changed: Time step information is saved for pvtu output. In addition the default
  value of VtkFlags::cycle is changed to numbers::invalid_unsigned_int and
  VtkFlags::time to std::numeric_limits<double>::lowest() to ensure it is added to
  the output at the first time step (cycle=0).
  <br>
  (Florian Schulze, 2025/11/08)
 </li>

 <li>
  Fixed: Move constructor and assignment operator of Triangulation<dim> did not
  move the data member cell_data_attached_data. This led to failures when trying
  to use a moved to Triangulation e.g. for refinement.
  <br>
  (Alexander Reinhold, 2025/11/07)
 </li>

 <li>
  New: Portable::MatrixFree is now able to support actions on multigrid levels to support local smoothing multigrid.
  <br>
  (Michał Wichrowski, 2025/11/05)
 </li>

 <li>
  New: FEPointEvaluation can now evaluate from a list of values
  that represent multiple components of a scalar FiniteElement.
  <br>
  (Martin Kronbichler, 2025/10/29)
 </li>

 <li>
  New: SegmentDataOut — adds a DataOut-style interface for exporting
  data on 1D segments.
  <br>
  (Michał Wichrowski, 2025/10/27)
 </li>

 <li>
  Changed: The ParticleHandler class now provides references to the underlying
  triangulation and mapping through member functions get_triangulation() and
  get_mapping(). In addition, the member functions get_particle_positions(
  VectorType &) and get_particle_positions(std::vector<Point<spacedim>> &) are
  marked as const member functions, so that they can be called from const
  ParticleHandler objects.
  <br>
  (Yimin Jin, 2025/10/27)
 </li>

 <li>
  New: The GridGenerator::convert_simplex_to_hypercube_mesh() function
  splits all cells of a simplex mesh into hypercubes, producing a
  mesh with cells of only the latter kind.
  <br>
  (Wolfgang Bangerth, Dominik Still, 2025/09/25)
 </li>

 <li>
  New: Added variadic template overloads for FullMatrix::kronecker_product()
  to compute the Kronecker product of multiple matrices in a single call.
  <br>
  (Michał Wichrowski, 2025/10/21)
 </li>

 <li>
  New: Added FullMatrix::permute() method to permute rows and columns of a matrix.
  <br>
  (Michał Wichrowski, 2025/10/20)
 </li>

 <li>
  Improved: The various all_zero() functions (e.g., Vector::all_zero())
  now consistently use `std::all_of()` to detect nonzero entries instead of
  norms or ad-hoc loops.
  <br>
  (David Wells, 2025/10/20)
 </li>

 <li>
  New: FEEvaluationData::read_cell_data now is also able to fetch data from Table
  <br>
  (Michał Wichrowski, 2025/09/27)
 </li>

 <li>
  Fixed: `MGTransferMatrixFree::initialize_dof_vector()` optionally forces the target vector to have the same partitioner as stored in the transfer operator. This fixes a bug where the two partitioners were assumed to be the same in `MGTwoLevelTransfer::restrict_and_add()`, triggering an assertion in debug mode for block-based systems.
  <br>
  (Michele Bucelli, 2025/10/01)
 </li>

 <li>
  New: A safety factor that multiplies the maximum eigenvalue estimate for
  PreconditionChebyshev and PreconditionRelaxation
  can now be provided via AdditionalData.
  <br>
  (Michał Wichrowski, 2025/09/27)
 </li>

 <li>
  Fixed: Resolved an access violation when using a dealii::DataOut object that has been copied or moved.
  <br>
  (Daniel Abele, 2025/09/26)
 </li>

 <li>
  New: The GridIn::read_ugrid() function now reads mesh files written
  in ugrid format.
  <br>
  (Wolfgang Bangerth, Jerett Cherry, 2025/09/25)
 </li>

 <li>
  Changed: The TimerOutput class now reports CPU time spent in sections
  as the sum over all processes in the given MPI communicator.
  Previously, it would only report for each section the CPU time on the current
  process. This meant the time spent in sections would not add up to the total
  CPU time, which was always reported as the sum over all MPI processes.
  <br>
  (Rene Gassmoeller, 2025/09/24)
 </li>

 <li>
  Added the MatrixScaling class. The class implements various known matrix scaling algorithms.
  <br>
  (Davide Polverino, 2025/09/17)
 </li>

 <li>
  New: Added bundled version of magic_enum.
  <br>
  (Marc Fehling, 2025/09/17)
 </li>

 <li>
  Improved: TimeStepping::LowStorageRungeKutta now also supports
  TimeStepping::FORWARD_EULER and TimeStepping::HEUN_EULER.
  <br>
  (Peter Munch, 2025/09/10)
 </li>

 <li>
  Changed: We now require Kokkos 4.0 or newer when device support is
  enabled (CUDA, HIP, SYCL).
  <br>
  (Timo Heister, 2025/09/02)
 </li>

 <li>
  Fixed: Polynomials derived from the class PolynomialsVectorAnisotropic would
  previously not overwrite all previous content in the `values`, `grads` and
  higher derivative fields that are designated as output. The behavior has been
  fixed, making sure that all content is explicitly written in the function.
  <br>
  (David Wells, Martin Kronbichler, 2025/08/19)
 </li>

 <li>
  New: GridIn used to support vtk meshes with version 3.0 only.
  Now it supports vtk meshes with version 5.1.
  <br>
  (Mohamad Ghadban, 2025/08/09)
 </li>

 <li>
  New: The ClosestSurfacePoint class in the Tools Shifted Boundary Method
  computes the closest point on a surface defined by a level set function.
  <br>
  (Michał Wichrowski, 2025/07/25)
 </li>

 <li>
  New: Integration of Gmsh api utilities for parallel distributed mesh handling: the
  function GridIn::read_partitioned_msh() can now be used to read partitioned mesh files in parallel.
  <br>
  (Raksha Devi, Luca Heltai, 2025/07/25)
 </li>

 <li>
  Improved: MappingFEField now stores which DoFs are nonzero in each component,
  which makes subsequent calculations about 30-50\% faster.
  <br>
  (David Wells, 2025/07/14)
 </li>

 <li>
  Fixed: PetscWrappers::SparseMatrix::memory_consumption() previously
  always returned 0. This is now fixed.
  <br>
  (Anna Long, 2025/07/11)
 </li>

</ol>

*/
