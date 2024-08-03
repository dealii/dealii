// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
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
@page changes_between_9_5_2_and_9_6_0 Changes between Version 9.5.2 and 9.6.0

<p>
This is the list of changes made between the release of deal.II version
9.5.2 and that of 9.6.0. All entries are signed with the names of the
author.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="952-960-incompatible"></a>
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
  Changed: ParameterHandler::read_in() now throws a ExcFileNotOpen exception
  instead of PathSearch::ExcFileNotFound.
  <br>
  (Matthias Maier, 2024/01/08)
 </li>

 <li>
  Changed: The GridIn::read() function variant taking a filename string no
  longer uses PathSearch to try to find an appropriate file. Instead, the
  function requires now that the filename string is a valid relative or
  absolute path to mesh file in question. Additionally, GridIn::read() now
  throws ExcFileNotOpen instead of PathSearch::FileNotFound.
  <br>
  (Matthias Maier, 2024/01/08)
 </li>

 <li>
  Changed: GridTools::orthogonal_equality() now returns a
  std::optional<unsigned char> instead of a
  std::optional<std::bitset<3>>. As a consequence,
  PeriodicFacePair now also stores an unsigned char instead of
  a std::bitset<3> to represent the relative orientation of the first face to the
  second face.
  <br>
  (David Wells, 2024/05/25)
 </li>

 <li>
  Changed: The Quadrature class used to have a constructor that took an
  integer argument. This was error-prone because it was easy to accidentally
  write
  @code
    Quadrature<dim> quadrature(3);
  @endcode
  where
  @code
    QGauss<dim> quadrature(3);
  @endcode
  was meant. As a consequence, this constructor has been removed from the `public`
  interface of the class.
  <br>
  (Wolfgang Bangerth, 2024/05/13)
 </li>

 <li>
  Removed: The Tensor<0,dim> template had functions `begin_raw()` and
  `end_raw()`, which had been deprecated for a couple of releases
  already. They have been removed now.
  <br>
  Furthermore, the corresponding functions for the general
  Tensor<rank,dim> template have been deprecated. They can still be used
  for the time being, but only for the specific case `rank==1`. This is
  because the underlying assumption of these functions is that tensors
  store their elements in memory in a contiguous fashion, but that is
  only true for the case `rank==1`. Similarly, the `make_array_view()`
  function that takes tensors as input has been marked as deprecated,
  and can only be used for rank-1 tensors from now on.
  <br>
  (Wolfgang Bangerth, 2024/03/26)
 </li>

 <li>
  Changed: SolverGMRES::AdditionalData now controls the size of the Arnoldi
  basis through the new member variable
  SolverGMRES::AdditionalData::max_basis_size, rather than the number of
  temporary vectors (which is the basis size plus 2). As a result, the default
  value of the basis size is now 30, compared to 28 used before. The old
  variable SolverGMRES::AdditionalData::max_n_tmp_vectors is still available,
  but whenever SolverGMRES::AdditionalData::max_basis_size is set to a non-zero
  value (including the value set by the default constructor), the latter takes
  precedence. Furthermore, the default algorithm has been changed from
  the modified Gram-Schmidt algorithm to the classical Gram-Schmidt algorithm.
  The latter uses unconditional reorthogonalization delayed by one step,
  following the algorithm described in @cite Bielich2022.
  <br>
  (Martin Kronbichler, 2024/03/12)
 </li>

 <li>
  Changed:`ParameterHandler::parse_input_from_json()` reads now a json-file
  without mangled parameter names.
  <br>
  (Magdalena Schreter, Peter Munch, 2024/01/26)
 </li>

 <li>
  Removed: The rarely used long double variants `std_cxx17::legendrel` and
  `std_cxx17::cyl_bessel_jl` have been removed due to a compatibility issue
  with boost-1.83. If you happen to use these functions in your project
  you can change to using `std::legendrel` and `std::cyl_bessel_jl` directly,
  which are provided by the C++ standard template library.
  <br>
  (Matthias Maier, 2024/01/08)
 </li>

 <li>
  Removed: The function FEPointEvaluation::real_point() has been
  renamed to FEPointEvaluation::quadrature_point(). The old function
  has been deprecated.
  <br>
  (Peter Munch, 2024/01/01)
 </li>

 <li>
  Changed: We have always considered PETSc- and Trilinos-based vectors
  that have ghost elements as immutable, i.e., it is possible to read
  elements (including the ghost elements) but not to write into them. On
  the other hand, the `compress()` function available in all vector
  types is meant to communicate values written into non-locally-owned
  vector elements to their proper owners, for example during assembly of
  a right hand side vector. This `compress()` operation clearly only
  makes sense if a vector does not have ghost elements, because only
  then is it possible to write into the vector at all, but this
  restriction was not enforced -- the `compress()` function simply did
  not do anything at all in these cases. This has now changed: Because
  it does not make sense to call `compress()` on vectors that have ghost
  elements, it is now forbidden to call it and this case will be caught
  by an assertion.
  <br>
  (Wolfgang Bangerth, 2023/12/27)
 </li>

 <li>
  Changed: FiniteElement::adjust_quad_dof_index_for_face_orientation(),
  FiniteElement::face_to_cell_index(), and
  FiniteElement::adjust_line_dof_index_for_line_orientation() now use the
  standardized <tt>combined_orientation</tt> orientation encoding as input
  arguments rather than one or three booleans.
  <br>
  (David Wells, 2023/12/21)
 </li>

 <li>
  Changed: Meshes in 2D now use the face_orientation boolean (instead of the
  face_flip boolean) to represent their orientation.
  <br>
  (David Wells, 2023/12/20)
 </li>

 <li>
  Removed: The file `include/deal.II/multigrid/mg_constrained_dofs.h`
  used to include `deal.II/dofs/dof_tools.h` and
  `deal.II/multigrid/mg_tools.h`. This is not the case any more.
  <br>
  (Peter Munch, 2023/10/31)
 </li>

 <li>
  Removed: The functions evaluate() and integrate() of
  CUDAWrappers::FEEvaluation that take Bools have been
  deprecated and replaced by versions that take
  EvaluationFlags.
  <br>
  (Peter Munch, 2023/10/27)
 </li>

 <li>
  Changed: The class FEEvaluation now uses a different internal data layout for
  the gradients, exposed via FEEvaluation::begin_gradients(). Now, the entries
  of the partial derivatives in the space directions are placed adjacent in
  memory. The entries of different components are still separated by the entries
  on all points. This change has been made to simplify the access in the
  FEEvaluation::get_gradient() and FEEvaluation::submit_gradient() functions,
  which is especially useful for the case with many FE components. For the
  regular use of FEEvaluation apart from the plain pointers mentioned above,
  there is no change in behavior.
  <br>
  (Martin Kronbichler, 2023/09/11)
 </li>

 <li>
  Changed: For Triangulations created with
  MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(),
  the MeshSmoothing flag Triangulation::eliminate_unrefined_islands will
  be removed. It caused unintentional refinement during coarsening which
  led to problems in GlobalCoarseningFineDoFHandlerView. See also
  <a href="https://github.com/dealii/dealii/issues/15541">#15541</a>.
  <br>
  (Marc Fehling, 2023/08/22)
 </li>

 <li>
  Changed: MGTransferMatrixFree::interpolate_to_mg()
  returns now non-ghosted vectors.
  <br>
  (Peter Munch, 2023/07/25)
 </li>

 <li>
  Removed: Support for PETSc without MPI has been removed.
  <br>
  (Daniel Arndt, 2023/07/24)
 </li>

 <li>
  Removed: The methods copy_to_mg, copy_from_mg,
  and interpolate_to_mg used of
  MGTransferGlobalCoarsening and MGTransferMatrixFree
  used to have spacedim as template argument. Since it
  is not used and the underlying functionality did not work,
  it has been removed.
  <br>
  (Peter Munch, 2023/07/21)
 </li>

 <li>
  Removed: The deprecated header files
  - deal.II/base/std_cxx14/algorithm.h
  - deal.II/base/std_cxx14/memory.h
  - deal.II/base/std_cxx14/utility.h
  have been removed.
  <br>
  (Daniel Arndt, 2023/07/13)
 </li>

 <li>
  Changed: The CMake configuration does no longer export git revision strings
  and dates in `deal.IIConfig.cmake`. This avoids unnecessary
  reconfigurations of client projects (such as the testsuite). Instead the
  git revision strings and dates are now recorded in
  `deal.IIConfigGit.cmake` whose path is recorded in `${DEAL_II_GIT_CONFIG}`.
  This is similar to how we store the git revision in the header
  `deal.II/base/revision.h` instead of `deal.II/base/config.h` to avoid
  unnecessary recompilations of the library and user programs. User project
  who need the information in CMake must do an
  `include(${DEAL_II_GIT_CONFIG})`.
  <br>
  (Matthias Maier, 2023/07/03)
 </li>

 <li>
  Removed: The LinearAlgebra::Vector class has been removed without deprecation.
  Users should use the standard Vector class instead.
  <br>
  (David Wells, 2023/07/03)
 </li>

 <li>
  Removed: The deprecated functions
  Triangulation::create_triangulation_compatibility(),
  GridTools::cell_measure(),
  MappingQCache::initialize(),
  TriaAccessor::number_of_children(),
  the deprecated class LinearAlgebra::CommunicationPatternBase
  and the PETScWrappers::PreconditionerBase alias
  have been removed.
  <br>
  (Daniel Arndt, 2023/07/03)
 </li>

 <li>
  Removed: The deprecated functions GridReordering::reorder_cells()
  and GridReordering::invert_all_cells_of_negative_grid()
  have been removed.
  <br>
  (Daniel Arndt, 2023/07/03)
 </li>

 <li>
  Changed: The overloads of KellyErrorEstimator::estimate() which compute
  multiple vector error estimates at once now require that vectors be passed
  in ArrayView objects of pointers to ReadVector instead of std::vector objects
  of pointers to VectorType objects.
  <br>
  (David Wells, 2023/07/02)
 </li>

 <li>
  Removed: The deprecated MatrixFree::FEEvaluation
  and MatrixFreeFaceEvaluation member function taking
  bools instead of EvaluationFlags have been removed.
  <br>
  (Daniel Arndt, 2023/07/02)
 </li>

 <li>
  Removed: The deprecated member variables
  SUNDIALS::KINSOL::solve_jacobian_system,
  and SUNDIALS::IDA::solve_jacobian_system
  have been removed.
  <br>
  (Daniel Arndt, 2023/07/02)
 </li>

 <li>
  Removed: The deprecated overloads for DoFTools::extract_boundary_dofs
  returning void have been removed.
  <br>
  (Daniel Arndt, 2023/07/02)
 </li>

 <li>
  Changed: Most of the member functions of FEValues now have different template
  parameters. As a result, some function calls which relied on implicit
  conversions to ArrayView now require explicit conversions instead.
  <br>
  (David Wells, 2023/07/01)
 </li>

 <li>
  Updated: deal.II requires a compiler supporting C++17.
  <br>
  (Daniel Arndt, 2023/07/01)
 </li>

 <li>
  Removed: The deprecated classes Threads::Thread,
  and Threads::ThreadGroup and the function
  Threads::new_thread() has been removed.
  <br>
  (Daniel Arndt, 2023/07/01)
 </li>

 <li>
  Removed: The deprecated classes FEValuesViews::Scalar::OutputType,
  FEValuesViews::Vector::OutputType,
  FEValuesViews::SymmetricTensor::OutputType,
  have been removed.
  <br>
  (Daniel Arndt, 2023/06/30)
 </li>

 <li>
  Removed: Deprecated functions in the FEInterfaceValues nampespace have been removed.
  <br>
  (Daniel Arndt, 2023/06/29)
 </li>

 <li>
  Removed: Deprecated constructors for ConsensusAlgorithms::Interface,
  ConsensusAlgorithms::Process, ConsensusAlgorithms::NBX, ConsensusAlgorithms::PEX,
  ConsensusAlgorithms::Serial, ConsensusAlgorithms::Selector have been removed
  as well as ConsensusAlgorithms::AnonymousProcess.
  <br>
  (Daniel Arndt, 2023/06/29)
 </li>

 <li>
  Removed: The deprecated GridTools::CellDataTransferBuffer class has been removed.
  <br>
  (Daniel Arndt, 2023/06/29)
 </li>

 <li>
  Removed: The deprecated signal Triangulation::Signals::cell_weight
  has been removed along with the deprecated class LegacySignal. Use
  Triangulation::Signals::weight instead.
  <br>
  (Marc Fehling, 2023/06/27)
 </li>

 <li>
  Removed: The deprecated Functions::LevelSet namespace
  has been removed
  <br>
  (Daniel Arndt, 2023/06/30)
 </li>

 <li>
  MatrixFree::reinit() would always set up the data structures for inner faces,
  also in case only
  MatrixFree::AdditionalData::mapping_updates_flags_boundary_faces was set. As
  this can lead to considerably higher memory consumption, the inner faces are
  now only set up when requested, increasing efficiency. When inner faces are
  desired, make sure to set
  MatrixFree::AdditionalData::mapping_updates_flags_inner_faces.
  <br>
  (Martin Kronbichler, 2023/06/19)
 </li>

 <li>
  Changed: When calling Triangulation::get_manifold(manifold_id) for
  a non-existing manifold_id, the function returned a flat manifold.
  This behaviour has been changed and an assert is thrown in this case.
  Furthermore, the functions Triangulation::reset_manifold() and
  Triangulation::reset_all_manifolds() do not actually remove the manifolds
  but make them flat.
  <br>
  (Peter Munch, 2022/12/09)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="952-960-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: LinearAlgebra::TpetraWrappers preconditioner interfaces
  for Trilinos Ifpack2, as listed below:
  <ul>
      <li>Preconditioners that were also available in the old interface.
          <ul>
              <li>PreconditionIdentity</li>
              <li>PreconditionJacobi</li>
              <li>PreconditionSOR</li>
              <li>PreconditionSSOR</li>
              <li>PreconditionChebyshev</li>
              <li>PreconditionILU</li>
              <li>PreconditionILUT</li>
              <li>PreconditionBlockJacobi</li>
              <li>PreconditionBlockSOR</li>
              <li>PreconditionBlockSSOR</li>
          </ul>
      </li>
      <li>Newly introduced preconditioners
          <ul>
              <li>PreconditionL1Jacobi<br>
                  - A new variant to improve coupling in MPI-parallel applications</li>
              <li>PreconditionL1GaussSeidel<br>
                  - A new variant to improve coupling in MPI-parallel applications</li>
              <li>PreconditionIfpack<br>
                  - A generic Ifpack2 preconditioner for any
                  classes or parameters not interfaced (yet)</li>
          </ul>
      </li>
  </ul>
  (Jan Philipp Thiele, 2024/07/18)
 </li>

 <li>
  New: step-1 and step-3 now discuss how to use triangular meshes in their "Results" sections.
  <br>
  (Wolfgang Bangerth, 2024/05/13)
 </li>

 <li>
  Updated: Tetrahedral meshes are now refined quasi-uniformly, preventing degradation
  with higher refinement levels. The changes also include the relevant adjustments in the
  geometric multigrid infrastructure.
  <br>
  (Dominik Still, Niklas Fehn, Peter Munch, Martin Kronbichler, Wolfgang Bangerth, 2024/05/10)
 </li>

 <li>
  New: LinearAlgebra::TpetraWrappers direct solver interfaces
  for Trilinos Amesos2, as listed below
  <dl>
      <dt>SolverDirect</dt>
      <dd>- Same interface as TrilinosWrappers, but with the option to pass a parameter list</dd>
      <dt>SolverDirectKLU2</dt>
      <dd>- An interface to the KLU2 solver, providing solver-specific AdditionalData</dd>
  </dl>
  (Jan Philipp Thiele, 2024/02/12)
 </li>

 <li>
  Updated: The version of BOOST bundled with deal.II has been updated from 1.70
  to 1.84.
  <br>
  (Wolfgang Bangerth, 2024/04/04)
 </li>

 <li>
  New: There is now a function-like macro `DEAL_II_NOT_IMPLEMENTED();`
  that can be used to indicate places where something is not
  implemented. If a code runs into a place where it appears, the program
  is aborted (or, as documented, an exception is thrown). In contrast to
  the old way of indicating such a thing (namely, writing
  `Assert(false, ExcNotImplemented());`, the error is raised even in
  release mode.
  <br>
  Similarly, there is now also `DEAL_II_ASSERT_UNREACHABLE()` that is used
  in places that control flow really should not reach, but where it is
  best to abort a program if it does.
  <br>
  (Wolfgang Bangerth, 2024/01/24)
 </li>

 <li>
  New: Support for hanging nodes in FE_NedelecSZ.
  <br>
  (Sebastian Kinnewig, 2023/11/28)
 </li>

 <li>
  New: Introduce FERemoteEvaluation. This new class provides similar
  interfaces as FEEvaluation, FEFaceEvaluation, and FEPointEvaluation.
  It can be used to access values and/or gradients in quadrature points
  at remote triangulations. In particular the class is designed such that
  it can also be used for Nitsche-type mortaring or L2-Projection with
  CGAL in the future.
  <br>
  (Johannes Heinz, Peter Munch, 2023/11/27)
 </li>

 <li>
  New: The new tutorial step-89 presents the use of FERemoteEvaluation
  during matrix-free operator evaluation for non-matching and Chimera
  methods. The acoustic conservation equations are solved using
  Nitsche-type mortaring and point-to-point interpolation to demonstrate
  that a simple point-to-point interpolation approach is sometimes not
  sufficient.
  <br>
  (Johannes Heinz, Maximilian Bergbauer, Marco Feder, Peter Munch, 2023/11/27)
 </li>

 <li>
  New: LinearAlgebra::TpetraWrappers::SparseMatrix class
  that implements a wrapper for Tpetra::CrsMatrix.
  <br>
  (Sebastian Kinnewig, 2023/11/22)
 </li>

 <li>
  New: FE_SimplexP and FE_SimplexDGP now support cubic elements.
  <br>
  (David Wells and Peter Munch, 2023/10/11)
 </li>

 <li>
  New: The new tutorial step-87 presents the advanced point evaluation
  functionalities of deal.II, specifically useful for evaluating
  finite element solutions at arbitrary points on distributed meshes.
  <br>
  (Magdalena Schreter-Fleischhacker, Peter Munch, 2023/09/05)
 </li>

 <li>
  Removed: Some, but not all, of the vector classes were derived from a
  base class VectorSpaceVector. This class had been intended to provide
  an abstract interface (via `virtual` functions) to vector-vector
  operations such as dot products or norms. But it turns out that that
  is not practical in many cases: Functions still need to either have
  access to individual elements of the vector, or they need to be able
  to do matrix-vector products. As a consequence, it is rarely useful to
  only have a reference to the base class VectorSpaceVector: One
  actually needs a reference to the derived class. Because of this lack
  of use, we have removed the VectorSpaceVector base class from the
  library.
  <br>
  (Wolfgang Bangerth, 2023/07/06)
 </li>

 <li>
  Deprecated: The entire LocalIntegrators namespace has been
  deprecated. The namespace represents an effort at making integration
  simpler by providing a library of terms that appear frequently in
  equations. But this only made the simple case simple -- everything
  that exceeded the simplest case still needed to be done by hand, and
  using these pre-built integrators did not help teaching how to do
  that. As a consequence, these pre-built integrators will be removed in
  a future version of the library.
  <br>
  (Wolfgang Bangerth, 2023/07/03)
 </li>

 <li>
  New: All vector classes in deal.II now inherit from ReadVector, which
  provides some common read operations.
  <br>
  (David Wells, 2023/07/01)
 </li>

 <li>
  New: Added classes for conforming Hermite interpolation polynomials,
  allowing for higher levels of regularity between elements to be directly
  enforced.
  <br>
  (Ivy Weber, 2023/06/21)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="952-960-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Compatibility: deal.II's interface to CGAL has been updated to ensure
  compatibility with the upcoming CGAL 6.0 release.
  <br>
  (S
  bastien Loriot, 2024/07/22)
 </li>

 <li>
  Fixed: The identification of periodic boundaries for
  parallel::distributed::Triangulation did not work properly for some meshes
  that had different face orientations on the two sides of a periodic
  boundary. This is now fixed.
  <br>
  (Martin Kronbichler, Paras Kumar, 2024/07/16)
 </li>

 <li>
  New: There is a new initialization function MGTwoLevelTransfer::reinit() that
  takes two MatrixFree objects to describe a p-coarsening transfer operation,
  besides the previous option of computing the information via two separate
  DoFHandler objects. The new option is faster to set up and more memory
  efficient in case the MatrixFree objects already are available, such as when
  using those to define the level operators.
  <br>
  (Martin Kronbichler, 2024/07/22)
 </li>

 <li>
  New: The functions VectorTools::compute_nonzero_normal_flux_constraints(),
  VectorTools::compute_nonzero_normal_flux_constraints_on_level(),
  VectorTools::compute_no_normal_flux_constraints(),
  VectorTools::compute_no_normal_flux_constraints_on_level(),
  VectorTools::compute_nonzero_tangential_flux_constraints(), and
  VectorTools::compute_normal_flux_constraints() now also can use
  the mapping to compute the normal vectors.
  <br>
  (Peter Munch, 2024/07/19)
 </li>

 <li>
  New: The functions GridGenerator::hyper_shell() and
  GridGenerator::hyper_cube_with_cylindrical_hole() now support
  dim=2 and spacedim=3.
  <br>
  (Peter Munch, 2024/07/16)
 </li>

 <li>
  New: Class DiscreteTime supports now serialization.
  <br>
  (Pasquale Claudio Africa, 2024/07/11)
 </li>

 <li>
  Deprecated: The `interpolate` callback of the
  PETScWrappers::TimeStepper class has been renamed to
  PETScWrappers::TimeStepper::transfer_solution_vectors_to_new_mesh. The former
  name is still available for backward compatibility, but is now
  deprecated.
  <br>
  (Wolfgang Bangerth, 2024/07/09)
 </li>

 <li>
  Fixed: ReadWriteVector's copy constructor now correctly sets up its IndexSet
  when it has zero locally owned elements.
  <br>
  (David Wells, Laryssa Abdala, 2023/07/08)
 </li>

 <li>
  New: All the basic Manifold objects describing a coordinate system now make
  their essential geometric properties (e.g., the center) available through member
  functions.
  <br>
  (David Wells, 2024/07/06)
 </li>

 <li>
  New: Added a function FETools::compute_nodal_quadrature()
  which computes the nodal Quadrature rule for an interpolatory
  FiniteElement.
  <br>
  (David Wells, 2023/07/05)
 </li>

 <li>
  Deprecated: The `decide_for_coarsening_and_refinement` callback of the
  PETScWrappers::TimeStepper class has been renamed to
  PETScWrappers::TimeStepper::decide_and_prepare_for_remeshing. The former
  name is still available for backward compatibility, but is now
  deprecated.
  <br>
  (Wolfgang Bangerth, 2024/07/01)
 </li>

 <li>
  Deprecated: The `distribute` callback of the
  PETScWrappers::TimeStepper class has been renamed to
  PETScWrappers::TimeStepper::update_constrained_components. The former
  name is still available for backward compatibility, but is now
  deprecated.
  <br>
  (Wolfgang Bangerth, 2024/06/30)
 </li>

 <li>
  New: MatrixFree::loop() and MatrixFree::cell_loop() as well as
  FEEvaluation/FEFaceEvaluation can now also be used with ArrayView
  arguments. This allows to wrap externally managed vectors into the matrix-free
  facilities of deal.II.
  <br>
  (Martin Kronbichler, 2024/06/24)
 </li>

 <li>
  Fixed: AffineConstraints::make_consistent_in_parallel() now iteratively
  updates the constraints on user-specified DoFs until the set of constraints
  globally converged on all subdomains.
  <br>
  The function now also updates the locally stored constraints of the
  underlying AffineConstraints object. After using this function, it might be
  necessary to use the IndexSet retrieved by
  AffineConstraints::get_local_lines() when initializing data structures,
  for example:
  <code>
  DynamicSparsityPattern dsp(constraints.get_local_lines());
  LinearAlgebra::distributed::Vector<double> solution(
    locally_owned_dofs, constraints.get_local_lines(), mpi_communicator);
  </code>
  This used to be an issue in parallel hp-adaptive applications, when
  finite elements of different types have constraints on faces between
  locally relevant and artificial cells.
  <br>
  (Wolfgang Bangerth, Marc Fehling, Martin Kronbichler, Peter Munch, 2024/06/20)
 </li>

 <li>
  Fixed: The GridTools::Cache class accessed its members in ways that
  were not thread-safe. This should now be fixed.
  <br>
  (Wolfgang Bangerth, 2024/06/19)
 </li>

 <li>
  Fixed: The ParticleHandler could cause a division by zero
  for certain mappings when sorting particles that are
  located almost exactly on mesh vertices. This is fixed now.
  <br>
  (Rene Gassmoeller, 2024/06/18)
 </li>

 <li>
  Fixed: InitFinalize would initially Kokkos without taking into account the final
  number of threads that are used in deal.II.
  <br>
  (Daniel Arndt, 2024/06/04)
 </li>

 <li>
  Fixed: The SUNDIALS::ARKODE::solve_ode() function returned the number
  of time steps at which the intermediate output function was called,
  rather than the number of time steps performed as promised in the
  documentation. This has now been fixed.
  <br>
  The same bug also existed in the SUNDIALS::IDA::solve_dae() function,
  where it has also been fixed.
  <br>
  (Wolfgang Bangerth, 2024/06/02)
 </li>

 <li>
  New: There is a new class Threads::TaskResult that represents the
  result of a background task.
  <br>
  (Wolfgang Bangerth, 2024/05/21)
 </li>

 <li>
  Deprecated: Namespace FETools::Compositing had three functions
  FETools::Compositing::multiply_dof_numbers(),
  FETools::Compositing::compute_restriction_is_additive(), and
  FETools::Compositing::compute_nonzero_components() that took pointers to five
  finite elements and five multiplicities. These have now been
  deprecated in favor of the versions of these functions that take a
  vector of elements and a vector of multiplicities.
  <br>
  (Wolfgang Bangerth, 2024/05/19)
 </li>

 <li>
  New: MatrixFreeTools::compute_diagonal() and MatrixFreeTools::compute_matrix()
  can now handle face and boundary integrals needed for DG.
  <br>
  (Peter Munch, Magdalena Schreter-Fleischhacker, 2024/05/16)
 </li>

 <li>
  New: Added overloads of missing trigonometric functions (acos, asin, atan) and
  hyperbolic functions (cosh, sinh, tanh, acosh, asinh, atanh) for VectorizedArray.
  <br>
  (Magdalena Schreter-Fleischhacker, Johannes Resch, Martin Kronbichler 2024/05/15)
 </li>

 <li>
  Deprecated: The member variables `input_vector_names` and
  `output_names` of MeshWorker::LocalIntegrator are not used in the
  library. They are now deprecated. If you use them, simply introduce
  variables of the same name and type in your own derived class.
  <br>
  (Wolfgang Bangerth, 2024/05/14)
 </li>

 <li>
  Fixed: FE_SimplexP::get_restriction_matrix() now returns the correct restriction matrix
  for continuous elements, which differs compared to discontinuous elements.
  <br>
  (Dominik Still, Martin Kronbichler, 2024/05/15)
 </li>

 <li>
  Changed: step-39 now implements the local integration routines itself,
  rather than relying on deprecated concepts from namespace
  LocalIntegrator.
  <br>
  (Wolfgang Bangerth, 2024/05/13)
 </li>

 <li>
  Changed: Several tutorials (step-15, step-72, step-77) dealing with
  nonlinear problems have been simplified by using AffineConstraints
  for boundary conditions.
  <br>
  (Timo Heister, 2024/05/11)
 </li>

 <li>
  New: There is a new function ThreadLocalStorage::get_for_thread() that
  returns the thread-local object of another thread.
  <br>
  (Wolfgang Bangerth, 2024/04/24)
 </li>

 <li>
  New: deal.II is now compatible with SUNDIALS v7.0.0.
  <br>
  (Marc Fehling, 2024/04/17)
 </li>

 <li>
  Fixed: When accessing elements of Trilinos vectors that do not exist
  on the current MPI process (for example because a function *expects* a
  vector with ghost elements but got a completely distributed vector
  instead), deal.II tries to produce a reasonable error message pointing
  out the cause. But when using 64-bit indices, the logic that creates
  the error message was faulty, resulting in errors that were very hard
  to debug and not helpful at all. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2024/04/17)
 </li>

 <li>
  New: Added get_shape_value_component(), get_shape_grad_component()
  and get_shape_grad_grad_component() to FE_NedelecSZ.
  <br>
  (Sebastian Kinnewig, 2024/04/16)
 </li>

 <li>
  New: Added get_prolongation_matrix() to FE_NedelecSZ.
  <br>
  (Sebastian Kinnewig, 2024/04/16)
 </li>

 <li>
  New: Add GridGenerator::cylinder, GridGenerator::subdivided_cylinder,
  and GridGenerator::truncated_cone to the python wrappers
  <br>
  (Bruno Turcksin, 2024/04/11)
 </li>

 <li>
  New: The new function DoFTools::extract_level_constant_modes()
  allows to extract constant modes on multigrid levels.
  <br>
  (Peter Munch, 2024/04/07)
 </li>

 <li>
  Improved: VectorTools::interpolate is now faster for the case of a DoFHandler
  based on an FESystem with a single base element that has classical support
  points (like FE_Q or FE_DGQ).
  <br>
  (Martin Kronbichler, 2024/04/03)
 </li>

 <li>
  Improved: Several functions in MappingCartesian have been carefully rewritten
  for better performance, e.g., by avoiding divisions and making use of recurring
  operations.
  <br>
  (Martin Kronbichler, 2024/04/01)
 </li>

 <li>
  New: Add GridGenerator::plate_with_a_hole, GridGenerator::channel_with_cylinder,
  and GridGenerator::hyper_ball_balanced to the python wrappers
  <br>
  (Bruno Turcksin, 2024/03/27)
 </li>

 <li>
  Fixed: Utilities::unpack() could result in a segmentation fault for
  objects that have alignment greater than the default. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2024/03/24)
 </li>

 <li>
  Changed: The class FEEvaluationAccess, which was previously used to control
  access to values, gradients and Hessians at quadrature points in specialized
  form for scalar, vector-valued and general n-component system, has been merged
  into the class FEEvaluationBase, using conditional types deduced from
  std::conditional. This should not change any application code, but reduced the
  amount of code by 1,500 lines, thus making the library simpler to maintain.
  <br>
  (Martin Kronbichler, Magdalena Schreter-Fleischhacker, 2024/03/22)
 </li>

 <li>
  New: SolverGMRES and SolverFGMRES can now use an additional orthogonalization
  strategy, controlled by
  LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt. This
  implements the classical Gram-Schmidt method with delayed reorthogonalization,
  a low-synchronization algorithm (performing one global reduction per GMRES
  iteration for deal.II's own vectors) that has excellent stability properties.
  <br>
  (Martin Kronbichler, 2024/03/19)
 </li>

 <li>
  Changed: hp::Refinement::choose_p_over_h() now communicates
  refinement flags and future FE indices on ghost cells for
  all types of parallel Triangulation objects to decide between
  p- and h-refinement.
  <br>
  (Marc Fehling, 2024/03/18)
 </li>

 <li>
  Improved: The implementations of SolverGMRES and SolverFGMRES have been
  overhauled and made more similar. In particular, SolverFGMRES now uses the
  same internal algorithm to solve the minimization problem in the Arnoldi
  basis, employing Givens rotations in analogy to the setting used by
  SolverGMRES. Since the Arnoldi process is sensitive to roundoff errors, this
  change might slightly affect iteration counts (often giving slightly better
  results).
  <br>
  (Martin Kronbichler, 2024/03/12)
 </li>

 <li>
  Fixed: The memory of particles which are deleted during refinement was not
  correctly freed in the memory pool.  This could lead to a minor memory leak and
  error messages about an inconsistent state during destruction of the property
  pool. This is fixed now.
  <br>
  (Rene Gassmoeller, 2024/03/07)
 </li>

 <li>
  Improved: The parse functions of ParameterHandler used to add subsections
  if subsections/parameters have not been defined. This led to an output
  of ParameterHandler::print_parameters() that also contained these subsections.
  This is not the case anymore.
  <br>
  (Peter Munch, Magdalena Schreter-Fleischhacker, 2024/02/12)
 </li>

 <li>
  Fixed: The function hp::Refinement::limit_p_level_difference() now
  correctly sets future FE indices in case of p-coarsening.
  <br>
  (Marc Fehling, 2024/02/29)
 </li>

 <li>
  Fixed: The colorize option for the cylinder shell could yield inadequate
  boundary ID for some aspect ratios that were very different from the
  default one.
  <br>
  (Bruno Blais, 2024/02/27)
 </li>

 <li>
  Fixed: In a special case where a particle lies on a vertex, the ParticleHandler detected it as lost.
  We introduce a tolerance for GeometryInfo<dim>::is_inside_cell() inside the ParticleHandler to fix this issue.
  <br>
  (Magdalena Schreter-Fleischhacker, Julian Brotz, 2024/02/22)
 </li>

 <li>
  New: There is now functions VectorTools::interpolate_to_coarser_mesh()
  and VectorTools::interpolate_to_finer_mesh() that, different from
  VectorTools::interpolate_to_different_mesh() also work for parallel
  triangulations.
  <br>
  (Wolfgang Bangerth, 2024/02/20)
 </li>

 <li>
  Fixed: Fix make_flux_sparsity_pattern for a FECollection containing
  FE_Nothing elements.
  <br>
  (Magdalena Schreter-Fleischhacker, Andreas Ritthaler, 2024/02/19)
 </li>

 <li>
  Fixed: Added some instantiations to "sparse_matrix.inst.in" such that
  copying real sparse matrices into complex sparse matrices and performing
  matrix vector operations of a complex sparse matrix with a real vector
  into a complex vector are possible now.
  <br>
  (Malik Scheifinger, 2024/02/16)
 </li>

 <li>
  New: Add potentially more efficient vector operation to SUNDIALS wrappers. These operations are used automatically
  for deal.II's distributed vectors.
  <br>
  (Sebastian Proell, 2024/02/15)
 </li>

 <li>
  Fixed: The function GridTools::compute_active_cell_halo_layer() now
  also supports perodic meshes.
  <br>
  (Peter Munch, 2024/02/12)
 </li>

 <li>
  New: Add capacity to colorize boundary in GridGenerator::cylinder_shell
  <br>
  (Bruno Blais, 2023/02/12)
 </li>

 <li>
  New: If the relaxation parameter is set to 0, PreconditionRelaxation now
  determines it based on estimated eigenvalues, similarly as
  PreconditionChebyshev does.
  <br>
  (Peter Munch, Laura Prieto Saavedra, 2024/03/09)
 </li>

 <li>
  Deprecated: Matrix Free classes previously in the CUDAWrappers namespace have
  been moved to the Portable namespace.
  <br>
  (Bruno Turcksin, 2024/02/09)
 </li>

 <li>
  Fixed: Fix initialization of FEFaceEvaluation for hp in 3D.
  <br>
  (Maximilian Bergbauer, Peter Munch, Andreas Ritthaler, Magdalena Schreter-Fleischhacker, 2024/02/09)
 </li>

 <li>
  Fixed: The function DoFCellAccessor::distribute_local_to_global(local_matrix,
  local_vector, global_matrix, global_vector) mistook the argument type
  when calling OutputMatrix::add(row, n_cols, col_indices, values), which
  causes compile errors. This is now fixed.
  <br>
  (Yimin Jin, 2024/01/23)
 </li>

 <li>
  Changed:
  FE_DGQ::get_interpolation_matrix() now also takes as a source element FE_Nothing.
  <br>
  (Magdalena Schreter, Peter Munch, 2024/01/18)
 </li>

 <li>
  Changed:
  parallel::distributed::Triangulation::prepare_coarsening_and_refinement()
  is now a collective call that exchanges coarsen/refinement flags on
  ghost cells. This makes adaptive refinement independent of the number
  of MPI ranks involved.
  <br>
  (Timo Heister, Quang Hoang, Vladimir Yushutin, 2024/01/15)
 </li>

 <li>
  Fixed: In parallel computations, it was not possible to take the end
  iterator of a Trilinos sparsity pattern for the last locally-owned
  row. This is now fixed.
  <br>
  (Simon Wiesheier, Peter Munch, 2024/01/07)
 </li>

 <li>
  New: Add utility functions that setup FERemoteEvaluationCommunicator for
  point-to-point interpolation and Nitsche-type mortaring.
  <br>
  (Johannes Heinz, 2023/12/31)
 </li>

 <li>
  Improved: convert_to_distributed_compute_point_locations_internal() takes
  a pointer to a vector of quadratures as an additional argument. In case
  this pointer is not a null pointer, the vector is populated by the mapped
  quadratures that are internally used during the conversion.
  <br>
  (Johannes Heinz, 2023/12/29)
 </li>

 <li>
  Fixed: Exceptions of type RecoverableUserCallbackError, raised in
  callbacks `solve_with_jacobian` and
  `solve_with_jacobian_and_track_n_linear_iterations` of the
  TrilinosWrappers::NOXSolver class, are now treated as "recoverable", if
  the NOX parameter "Newton/Rescue Bad Newton Solve" is set to `true`,
  which is, in fact, its default value. Exceptions of all other types and
  also all exceptions raised in other callbacks are still treated as
  "irrecoverable".
  <br>
  (Vladimir Ivannikov, 2023/12/26)
 </li>

 <li>
  New: There is now a function Utilities::MPI::partial_and_total_sum().
  <br>
  (Wolfgang Bangerth, 2023/12/19)
 </li>

 <li>
  Fixed: When creating triangulations in dim=1 and spacedim=3, the
  Triangulation class would sometimes crash because of a code path that
  should not have been taken. This is now fixed.
  <br>
  (Vinayak Vijay, Wolfgang Bangerth, 2023/11/28)
 </li>

 <li>
  New: Add AdditionalData struct for configuring the RemotePointEvaluation
  class. The latter is passed into a new constructor. Deprecated: The old
  constructor of RemotePointEvaluation class taking the parameters
  value-by-value is marked as deprecated.
  <br>
  (Magdalena Schreter-Fleischhacker, Peter Munch, 2023/12/06)
 </li>

 <li>
  New: The FEValuesBase class now lazily computes the FEValuesViews
  objects it returns when you do something like `fe_values[velocities]`,
  whenever they are needed first, rather than *always* creating these
  objects.
  <br>
  (Wolfgang Bangerth, 2023/11/28)
 </li>

 <li>
  Fixed: Meshes created with GridGenerator::torus() would previously run
  into an assertion in p4est when using
  parallel::distributed::Triangulation. This is now fixed.
  <br>
  (Martin Kronbichler, Mathias Anselmann, Peter Munch, 2023/11/22)
 </li>

 <li>
  New: Introduce a new function,
  Utilities::Trilinos::teuchos_comm_to_mpi_comm(), to convert a
  Teuchos::RCP<Teuchos::Comm<int>> communicator into an MPI_Comm
  communicator.
  <br>
  (Sebastian Kinnewig, 2023/11/20)
 </li>

 <li>
  New: Introduce a new function, Utilities::Trilinos::internal::make_rcp(),
  to create Teuchos::RCP objects, avoiding the usage of 'operator new'.
  <br>
  (Sebastian Kinnewig, 2023/11/15)
 </li>

 <li>
  Fixed: MatrixFree::get_cell_active_fe_index did not
  return the correct value when used with categorization
  of cells. This affects also MFTools::compute_diagonal().
  This is now fixed.
  <br>
  (David Schneider, 2023/11/13)
 </li>

 <li>
  Fixed: In step-41, the mass matrix B is missing in the expressions
  of the residuals of the nonlinear system. See the equation under the
  line "This suggest a semismooth Newton step of the form" in the tutorial.
  This is fixed.
  <br>
  (Tao Jin, 2023/11/12)
 </li>

 <li>
  Improved: Point now supports constexpr construction and evaluation,
  just like Tensor.
  <br>
  (Chengjiang Yin, 2023/11/09)
 </li>

 <li>
  Fixed: A mistake in marking a `private` member function of the
  SphericalManifold class made it impossible to derive from this class
  in any meaningful way. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2023/11/09)
 </li>

 <li>
  Fixed: It was not possible to compare non-const iterators into block
  sparse matrices because a `friend` declaration referred to the wrong
  class. This is now fixed.
  <br>
  (Wolfgang Bangerth, Simon Wiesheier, 2023/10/31)
 </li>

 <li>
  New: There is a new overload of Threads::new_task() that takes lambda
  functions taking arguments, followed by the arguments the lambda
  function is to be called with on a new task.
  <br>
  (Wolfgang Bangerth, 2023/10/30)
 </li>

 <li>
  New: The new function RefinementCase::all_refinement_cases() returns
  an array of all possible refinement options of a cell.
  <br>
  (Wolfgang Bangerth, 2023/10/26)
 </li>

 <li>
  New: A new wrapper class Lazy<T> provides an efficient, thread-safe
  mechanism for lazy initialization of objects.
  <br>
  (Matthias Maier, 2023/10/25)
 </li>

 <li>
  Fixed: The periodic faces information was not being copied when using the
  `Triangulation::copy_triangulation()` function. This is now fixed.
  <br>
  (Laura Prieto Saavedra, Peter Munch, 2023/10/19)
 </li>

 <li>
  Fixed: There was an inconsistency between the way
  DoFTools::make_periodicity_constraints() and MGConstrainedDoFs
  set up the constraints for periodic constraints. This is fixed
  now.
  <br>
  (Peter Munch, Laura Prieto Saavedra, 2023/10/19)
 </li>

 <li>
  New: LinearAlgebra::TpetraWrappers::VectorReference for easier access
  to LinearAlgebra::TpetraWrappers::Vector. Moreover, several
  new member functions are introduced to LinearAlgebra::TpetraWrappers::Vector
  to be used as the "VectorType" template argument in AffineConstraints.
  <br>
  (Sebastian Kinnewig, 2023/10/18)
 </li>

 <li>
  New: New constructor for IndexSet that takes a
  Teuchos::RCP<Tpetra::Map> as input. And the
  function IndexSet::make_tpetra_map_rcp()
  returning a Teuchos::RCP<Tpetra::Map>.
  <br>
  (Sebastian Kinnewig, 2023/10/17)
 </li>

 <li>
  New: The print_formatted() functions of the following matrix classes
  now offer the option to specify a string that separates row entries:
  <br>
  SparseMatrixEZ, SparseMatrix, LAPACKFullMatrix, FullMatrix,
  ChunkSparseMatrix, BlockSparseMatrix, CUDAWrappers::SparseMatrix
  <br>
  (Marc Fehling, 2023/10/15)
 </li>

 <li>
  New: The ArrayView class now has a constructor that allows creation of
  an ArrayView object from a `boost::container::small_vector`.
  <br>
  (Wolfgang Bangerth, 2023/10/24)
 </li>

 <li>
  Fixed: The querying of line orientations in the case
  of tetrahedra has been fixed.
  <br>
  (Peter Munch, David Wells, 2023/10/10)
 </li>

 <li>
  New: The ArrayView class now has a constructor that allows creation of
  an ArrayView object from a `std::initializer_list`.
  <br>
  (Wolfgang Bangerth, 2023/10/10)
 </li>

 <li>
  New: The function AffineConstraints::get_view() can be used to obtain
  the constraints that correspond to a specific subset of degrees of
  freedom.
  <br>
  (Wolfgang Bangerth, 2023/10/05)
 </li>

 <li>
  Added: An option to sum into the solution values vector was added
  to FEPointEvaluation integrate() and test_and_sum().
  <br>
  (Maximilian Bergbauer, 2023/10/04)
 </li>

 <li>
  Added: An output operator was added to DerivativeForm.
  <br>
  (Maximilian Bergbauer, 2023/10/04)
 </li>

 <li>
  Improved: MGTransferMF can now permute DoFs during
  MGTransferMF::copy_to_mg(), MGTransferMF::copy_to_mg(),
  and MGTransferMF::interpolate_to_mg() if the outer
  solver and the multigrid preconditioner have been set up
  with different DoFHandler objects.
  <br>
  (Peter Munch, Laura Prieto Saavedra, 2023/10/03)
 </li>

 <li>
  New: The new function AffineConstraints::add_constraint() adds an
  entire constraint all at once, rather than splitting the task across
  AffineConstraint::add_line(), multiple calls to
  AffineConstraint::add_entry(), and
  AffineConstraint::set_inhomogeneity(). The new function
  AffineConstraints::constrain_dof_to_zero() is a short-cut that adds a
  constraint that requires a specific degree of freedom to be zero.
  <br>
  (Wolfgang Bangerth, 2023/10/05)
 </li>

 <li>
  New: The function IndexSet::get_view() now has an overload that
  computes the view with regard to an arbitrary mask.
  <br>
  (Wolfgang Bangerth, 2023/10/02)
 </li>

 <li>
  Improved: SparsityTools::reorder_Cuthill_McKee now runs considerably faster,
  especially for the case with many couplings between matrix rows.
  <br>
  (Martin Kronbichler, 2023/09/20)
 </li>

 <li>
  New: Add alternative interfaces to RemotePointEvaluation::evaluate_and_process and
  RemotePointEvaluation::process_and_evaluate.
  <br>
  (Magdalena Schreter, Peter Munch, 2023/09/17)
 </li>

 <li>
  New: Added helper class RemotePointEvaluation::CellData to store
  and access data of cell-specific points.
  <br>
  (Magdalena Schreter, Peter Munch, 2023/09/17)
 </li>

 <li>
  Changed: The interface to QGaussRadauChebyshev<dim>
  now matches that of the new quadrature QGaussRadau<dim>.
  <br>
  (Jan Philipp Thiele, 2023/09/16)
 </li>

 <li>
  New: Added QGaussRadau quadrature up to and including 8 quadrature points.
  <br>
  (Jan Philipp Thiele, 2023/09/15)
 </li>

 <li>
  Added: A new function shink_to_fit() was added to AlignedVector
  in analogy to std::vector.
  <br>
  (Maximilian Bergbauer, 2023/09/14)
 </li>

 <li>
  Improved: FEEvaluation::get_gradient() and FEEvaluation::submit_gradient() are
  now considerably faster for Raviart-Thomas elements on non-Cartesian meshes.
  The previous algorithm used non-optimal loop layouts for the local tensor
  contractions in the derivative of the Piola transform, which have been
  transformed to optimal-complexity variants.
  <br>
  (Martin Kronbichler, 2023/09/12)
 </li>

 <li>
  Improved: The internal implementation of the tensor-product evaluators used
  for FEEvaluation and FEFaceEvaluation has been cleaned up. This reduces the
  compile times for both the deal.II library and code using FEEvaluation with
  template parameters on the polynomial degrees. Also, the code is now simpler
  to maintain, especially for the evaluators for Raviart-Thomas elements.
  <br>
  (Martin Kronbichler, 2023/09/11)
 </li>

 <li>
  Fixed: FEValues requesting only mapping information but initialized with
  elements derived from FE_PolyTensor would previously run into an assertion on
  unstructured 2d meshes. This is now fixed.
  <br>
  (Martin Kronbichler, 2023/09/07)
 </li>

 <li>
  Changed: IndexSet objects could only be compared for equality or
  inequality against other IndexSet objects that had the same size. This
  did not allow for comparison against default-constructed objects, for
  example to test whether an object had been initialized. The
  restriction is therefore relaxed: IndexSet objects can be compared for
  equality and inequality against objects of the same size, or or size
  zero.
  <br>
  (Wolfgang Bangerth, 2023/09/07)
 </li>

 <li>
  Fixed: MGTwoLevelTransferBase now preserves the
  ghost state of the source vector in all call scenarios.
  <br>
  (Richard Schussnig, Martin Kronbichler, Peter Munch, 2023/09/04)
 </li>

 <li>
  Improved: VectorizedArray now also supports ARM Neon intrinsics.
  <br>
  (Maximilian Bergbauer, Peter Munch, 2023/08/24)
 </li>

 <li>
  New: The function IndexSet::is_subset_of() does as its name suggests.
  <br>
  (Wolfgang Bangerth, 2023/08/24)
 </li>

 <li>
  New: Add two tensor functions for the split of a 2nd-order symmetric tensor
  into a positive part and a negative part based on the signs of the eigenvalues
  obtained from the spectrum decomposition. The function positive_negative_split()
  performs the positive-negative split of the 2nd-order symmetric tensor given as
  the input. The function positive_negative_projectors() not only performs the split,
  but also provides the derivatives (two fourth-order tensors) of the positive/negative
  part of the tensor with respect to the input tensor.
  <br>
  (Tao Jin, 2023/08/22)
 </li>

 <li>
  Changed: You can now set MeshSmoothing flags in non-empty triangulations
  with Triangulation::set_mesh_smoothing().
  <br>
  (Marc Fehling, 2023/08/22)
 </li>

 <li>
  New: The Rayleigh--Kothe vortex has been extracted
  from step-68 and is now available as the new class
  Functions::RayleighKotheVortex.
  <br>
  (Bruno Blais, Peter Munch, 2023/08/20)
 </li>

 <li>
  Improved: The classes MGSmootherRelaxation and mg::SmootherRelaxation can now
  also handle matrices of types `MGLevelObject<std::unique_ptr<...>>` or
  `MGLevelObject<std::shared_ptr<...>>` in their `initialize()` functions by the
  use of Utilities::get_underlying_value(), rather than MGLevelObject of the
  actual matrix type only.
  <br>
  (Martin Kronbichler, 2023/08/16)
 </li>

 <li>
  Improved: deal.II now has a flag DEAL_II_USE_VECTORIZATION_GATHER to control
  the use of gather/scatter instructions on the x86 architecture. On a wide
  range of Intel hardware with microcode mitigation for the Intel Gather Data
  Speculation (GDS, aka Downfall) side channel vulnerability, in particular,
  server processors of the Broadwell, Skylake, Cascade Lake, and Ice Lake
  families released between 2015 and 2021, these instructions can be much slower
  than scalar loads. While the default behavior of deal.II is to aggressively
  enable these instructions in the intrinsics-class VectorizedArray, the new
  variable can be used to disable their use if deemed to give better
  performance.
  <br>
  (Martin Kronbichler, Matthias Maier, 2023/08/14)
 </li>

 <li>
  New: Added a function extract_children_of_level() that returns the
  bounding boxes associated to the children of a given level of an Rtree
  and stores them in a vector.
  <br>
  (Marco Feder, 2023/08/10)
 </li>

 <li>
  New: MGTwoLevelTransferNonNested now also supports FE_SimplexP.
  <br>
  (Peter Munch, Marco Feder, 2023/08/02)
 </li>

 <li>
  New: The new class FEValuesViews::RenumberedView allows one to filter an existing
  FEValuesViews object via two renumbering vectors, one acting on the degrees of
  freedom, and the other acting on the quadrature points.
  <br> (Luca Heltai, 2023/08/02)
 </li>

 <li>
  Improved: GridTools::find_active_cell_around_point() and
  GridTools::find_all_active_cells_around_point() now also work
  for simplices.
  <br>
  (Peter Munch, David Wells, 2023/07/25)
 </li>

 <li>
  New: The new function Triangulation::as_dof_handler_level_iterator()
  allows to create level iterators based on other cell iterators.
  <br>
  (Peter Munch, 2023/08/01)
 </li>

 <li>
  New: The function Triangulation::contains_cell()
  allows to check if Tringulation::create_cell_iterator()
  can be called for a specific cell id.
  <br>
  (Peter Munch, 2023/08/01)
 </li>

 <li>
  Improvement: Added gradient() implementation to VectorFunctionFromTensorFunction in
  function.h file.
  <br>
  (Abbas Ballout, 2023/07/30)
 </li>

 <li>
  New: Similar to the DataOut::add_data_vector() case,
  DataOut::add_mg_data_vector() now also copies the vector
  and performs the ghost-vector update internally.
  <br>
  (Peter Munch, 2023/07/25)
 </li>

 <li>
  New: Created a new class FECouplingValues, that helps with assembling of coupled
  finite element methods across different dimensions or different grids.
  <br>
  (Luca Heltai, 2023/07/21)
 </li>

 <li>
  New: Added InitFinalize class. This class is similar to MPI_InitFinalize but it
  allows users to decide which libraries should be initialized/finalized by the
  class.
  <br>
  (Bruno Turcksin, 2023/07/13)
 </li>

 <li>
  Fixed: FETools::get_fe_by_name() now works also for simplex finite elements.
  <br>
  (Luca Heltai, 2023/07/08)
 </li>

 <li>
  Improvement: Step-68 now uses the FEPointEvaluation to calculate the particle velocity from the velocity solution instead of manually interpolating at the particle location.
  <br>
  (Bruno Blais, 2023/07/0)
 </li>

 <li>
  Deprecated: The CellStatus enumeration has been deprecated in the Triangulation<dim, spacedim> class and moved to the global deal.II namespace. As a consequence, all references to Triangulation<dim, spacedim>, parallel::distributed::Triangulation<dim, spacedim>::CellStatus, and similar should be updated to use CellStatus directly. Also, the enumeration values have been renamed: CELL_PERSIST -> cell_will_persist, CELL_REFINE -> cell_will_be_refined, CELL_COARSEN -> children_will_be_coarsened, CELL_INVALID -> cell_invalid).
  <br>
  (Pasquale Claudio Africa, 2023/07/03)
 </li>

 <li>
  Fixed: The SolutionTransfer class writes into output vectors, but does
  not call compress() on them. This is of no consequence for deal.II
  vectors for which this class is mostly used (in contrast to the
  parallel::distributed::SolutionTransfer class), but leads to awkward
  downstream failures with, for example, PETSc vectors. This is now
  fixed.
  <br>
  (Wolfgang Bangerth, 2023/07/01)
 </li>

 <li>
  Changed: The constructor of the
  PETScWrappers::PreconditionSSOR::AdditionalData class is now
  `explicit`, thereby disallowing the implicit conversion of a number
  (the relaxation factor) to an object of this type.
  <br>
  (Wolfgang Bangerth, 2023/07/01)
 </li>

 <li>
  Improved: The `deal_ii_pickup_tests()` macro now prints a status line at
  the end summarizing how many tests have been configured for the given test
  category. Similarly, the top level target `setup_tests` concatenates these
  status lines and prints a summary after invocation.
  <br>
  (Matthias Maier, 2023/06/30)
 </li>

 <li>
  Fixed: IndexSet::add_index() would not recognize single indices already added,
  leading to an unnecessary quadratic complexity in case the same entry is added
  many times. This is now fixed.
  <br>
  (Martin Kronbichler, 2023/06/30)
 </li>

 <li>
  Fixed: It was previously possible to access the return value of a
  Threads::Task object if the underlying task had ended with an
  exception. But that return value was not initialized. This is now
  checked: You can no longer call Threads::Task::return_value() after an
  exception.
  <br>
  (Wolfgang Bangerth, 2023/06/30)
 </li>

 <li>
  Fixed: MGTwoLevelTransfer used within the global coarsening multigrid
  framework did not work when deal.II was compiled without MPI or when MPI_Init
  was not called. This is now fixed.
  <br>
  (Martin Kronbichler, 2023/06/28)
 </li>

 <li>
  Fixed: The function IndexSet::add_indices() was not efficient when
  adding sets of indices that are sorted but contain duplicates. This is
  now fixed.
  <br>
  (Wolfgang Bangerth, 2023/06/20)
 </li>

 <li>
  New: Added a function DoFTools::map_boundary_to_bulk_dof_iterators() that
  generates a mapping of codimension-1 active DoFHandler cell iterators to
  codimension-0 cells and face indices, to couple DoFHandler objects of different
  co-dimensions, initialized on grids generated with
  GridTools::extract_boundary_mesh()
  <br>
  (Luca Heltai, 2023/04/12)
 </li>

</ol>

*/
