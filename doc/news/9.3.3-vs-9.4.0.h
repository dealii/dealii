// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
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
@page changes_between_9_3_3_and_9_4_0 Changes between Version 9.3.3 and 9.4.0

<p>
This is the list of changes made between the release of deal.II version
9.3.3 and that of 9.4.0. All entries are signed with the names of the
author.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="933-940-incompatible"></a>
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
  Changed: For ParameterHandler::OutputStyle::JSON or
  ParameterHandler::OutputStyle::ShortJSON, the function
  ParameterHandler::print_parameters() demangles now the
  parameters before printing the output. Users do not
  need to demangle the output themselves anymore.
  <br>
  (Magdalena Schreter, Peter Munch, 2022/05/25)
 </li>

 <li>
  Improved: GridTools::get_coarse_mesh_description() no longer returns subcell data
  corresponding to default values (i.e., when the manifold_id is
  numbers::flat_manifold_id and when the boundary_id is
  numbers::internal_face_boundary_id). By skipping this data, which is not required by
  Triangulation::create_triangulation(), the total amount of data returned by this
  function is significantly lowered (and, therefore, so is the cost of creating a new
  Triangulation from the output).
  <br>
  (David Wells, 2022/05/09)
 </li>

 <li>
  Changed: All functions `n_nonzero_elements()`, which are members of the
  SparseMatrix and SparsityPattern family of classes from PETSc and
  Trilinos, now always return `std::uint64_t` instead of
  types::global_dof_index.
  <br>
  (Marc Fehling, 2022/04/04)
 </li>

 <li>
  Changed: To be able to guarantee correctness of
  the application of constraints in MatrixFree, we
  require that users call AffineConstraints::close()
  before calling MatrixFree::reinit() if AffineConstraints
  objects are passed that are not empty. One can check
  if an AffineConstraints object is closed with the
  new function AffineConstraints::is_closed().
  <br>
  (Peter Munch, 2022/03/08)
 </li>

 <li>
  Fixed: Callback functions attached to any weighted load balancing signal
  of parallel::distributed::Triangulation objects now need to handle cells
  with CELL_INVALID status explicitly.
  <br>
  If a cell gets refined, only the first child has the CELL_REFINE status,
  while all other children are CELL_INVALID.
  <br>
  (Marc Fehling, 2022/02/25)
 </li>

 <li>
  Changed: For weighted load balancing with parallel::distributed::Triangulation
  objects, an initial weight of `1000` will no longer be assigned to each cell.
  <br>
  The old signal Triangulation::Signals::cell_weight has been deprecated.
  Please use the new signal Triangulation::Signals::weight instead.
  <br>
  You can invoke the old behavior by connecting a function to the new
  signal that returns the base weight like this:
  @code{.cc}
  triangulation.signals.weight.connect(
    [](const typename parallel::distributed::Triangulation<dim>::cell_iterator &,
       const typename parallel::distributed::Triangulation<dim>::CellStatus)
      -> unsigned int { return 1000; });
  @endcode
  (Marc Fehling, 2022/02/25)
 </li>

 <li>
  Updated: Spurious inclusions of the header files `deal.II/grid/grid_tools.h`
  and `deal.II/grid/grid_tools_cache.h` from other headers of deal.II have been
  removed. User codes relying on this implicit inclusion will need to add the
  respective include files.
  <br>
  (Martin Kronbichler, 2022/02/22)
 </li>

 <li>
  Deprecated: We have introduced new versions of
  DoFTools::extract_locally_active_dofs(),
  DoFTools::extract_locally_active_level_dofs(),
  DoFTools::extract_locally_relevant_dofs(), and
  DoFTools::extract_locally_relevant_level_dofs().
  These versions return the index sets directly.
  The old versions have been early deprecated.
  <br>
  (Peter Munch, 2022/02/23)
 </li>

 <li>
  Changed: The function MatrixFreeTools::compute_diagonal()
  does not initialize the output vector anymore. Users
  need to call MatrixFree:initialize_dof_vector() themselves.
  This change is required to be able to also support block
  vectors.
  <br>
  (Magdalena Schreter, Peter Munch, 2022/02/01)
 </li>

 <li>
  Changed: The two algorithms in the ConsensusAlgorithms namespace
  (along with their base class) required a callback function that
  correctly sized the buffer into which the answer to a request was to
  be written. These algorithms have been rewritten in ways that make
  this no longer necessary, and they now correctly size these buffers
  themselves. As a consequence, the ConsensusAlgorithms::Process class
  no longer has the `prepare_for_answer()` function, and the
  ConsensusAlgorithms::AnonymousProcess class no longer has a
  corresponding member variable; in the latter case, the constructor of
  the class also no longer takes an argument to this effect.
  <br>
  (Wolfgang Bangerth, Peter Munch, 2022/01/18)
 </li>

 <li>
  Updated: The header inclusions of several files in the `deal.II/matrix_free`
  subfolder have been cleaned up to make the library more maintainable. As a
  result, some headers are no longer implicitly included and a user code might
  need to add the respective include file.
  <br>
  (Martin Kronbichler, 2021/11/30)
 </li>

 <li>
  Renamed: The namespace Functions::LevelSet has been renamed to
  Functions::SignedDistance.
  <br>
  (Peter Munch, 2021/10/19)
 </li>

 <li>
  Changed: For parallel::distributed::Triangulation, vertices will no longer
  be communicated automatically to p4est in debug mode. Please set the
  corresponding flag
  parallel::distributed::Triangulation::Settings::communicate_vertices_to_p4est
  while constructing your triangulation object.
  <br>
  (Marc Fehling, 2021/10/13)
 </li>

 <li>
  Changed: The 3D implementation of GridTools::rotate() with an integer
  for a Cartesian coordinate direction has been superseded by the version
  that accepts unit vectors as rotation axes. Further,
  Physics::Transformations::Rotations::rotation_matrix_3d() now requires
  a Tensor<1,3> object instead of a Point<3> as an axis.
  <br>
  (Marc Fehling, 2021/10/11)
 </li>

 <li>
  Updated: deal.II now requires ginkgo version 1.4.0 or newer
  (Fabian Castelli, 2021/09/01)
 </li>

 <li>
  Deprecated: The Tensor and SymmetricTensor classes had functions
  `begin_raw()` and `end_raw()` that suggested that the elements of
  these objects are stored in arrays and consequently could be accessed
  consecutively through iterators. While they are indeed stored
  contiguously in memory, they are not part of an array and, as a
  consequence, C++ does not actually allow this kind of access. The
  functions have therefore been deprecated and will be removed in due
  time.
  <br>
  As a consequence, the `make_array_view()` functions that allow turning
  a Tensor or SymmetricTensor into an ArrayView object have now also
  been deprecated.
  <br>
  (Wolfgang Bangerth, 2021/08/19)
 </li>

 <li>
  Changed: The DiagonalMatrix objects set up by MatrixFreeOperators::Base and
  its inheriting classes no set ghost entries.
  <br>
  (David Wells, 2021/08/17)
 </li>

 <li>
  Deprecated: Many of the member functions of the FEInterfaceValues class, as well
  as the FEInterfaceViews::Scalar and FEInterfaceViews::Vector classes, that are
  used to compute jumps and averages of shape function-related quantities have
  been renamed in such a way that they are more expressive about what they
  compute.
  <br>
  (Jean-Paul Pelteret, 2021/07/31)
 </li>

 <li>
  Changed: MatrixFreeOperators::MassOperator now computes its true diagonal (and
  not a lumped mass matrix) as its diagonal.
  <br>
  (David Wells, 2021/07/29)
 </li>

 <li>
  Deprecated: The autopartition parameter has been removed from
  parallel::DistributedTriangulationBase::load() and all inheriting
  classes.
  <br>
  (Marc Fehling, 2021/07/24)
 </li>

 <li>
  Changed: Custom repartitioning via Triangulation::Signals::cell_weight
  is no longer part of parallel::distributed::Triangulation::load(). If
  desired, call parallel::distributed::Triangulation::repartition()
  manually after load().
  <br>
  (Marc Fehling, 2021/07/23)
 </li>

 <li>
  Changed: The class MappingQ now applies a high-order mapping to a cells in the
  mesh, not just the ones touching the boundary. This makes the class completely
  equivalent to the previous class MappingQGeneric. All functionality in the
  latter class has been integrated into the MappingQ class. This avoids errors
  with curved manifolds in the interior of the domain, where MappingQ could lead
  to gaps or overlaps in the computational domain. In case the slightly better
  performance of using MappingQ1 in interior cells is desired,
  hp::MappingCollection can be used to replicate the old behavior.
  <br>
  (Martin Kronbichler, 2021/07/10)
 </li>

 <li>
  Updated: FE_RaviartThomasNodal now uses a different polynomial space to allow
  for a simpler use for faces in non-standard orientation. The new polynomials
  are anisotropic tensor products of Lagrange polynomials on the points of a
  Gauss--Lobatto quadrature formula. This change leads to different entries in
  the matrices and constraints, for example, but as the resulting polynomial
  space spans the same polynomials, no change in accuracy should be expected.
  <br>
  (Martin Kronbichler, Katharina Kormann, Konrad Simon, 2021/06/12)
 </li>

 <li>
  Removed: The PETSc Eisenstat preconditioner wrapper, which has never worked
  correctly in parallel, has been removed.
  <br>
  (David Wells, 2021/06/11)
 </li>

 <li>
  Removed: The deprecated exception ParameterHandler::ExcInvalidEntryForPatternXML
  and the deprecated overload of SparsityTools::gather_sparsity_pattern()
  have been removed.
  <br>
  (Daniel Arndt, 2021/06/07)
 </li>

 <li>
  Removed: The deprecated member function
  TriangulationBase::compute_vertices_with_ghost_neighbors() has been removed.
  <br>
  (Daniel Arndt, 2021/06/07)
 </li>

 <li>
  Removed: The deprecated constructors for FESeries::Fourier and
  FESeries::Legendre have been removed.
  <br>
  (Daniel Arndt, 2021/06/07)
 </li>

 <li>
  Removed: The deprecated specializations of VectorTools::integrate_difference()
  and VectorTools::project_boundary_values_curl_conforming() have been removed.
  <br>
  (Daniel Arndt, 2021/06/04)
 </li>

 <li>
  Removed: The deprecated overloads of
  FETools::lexicographic_to_hierarchic_numbering() and
  FETools::hierarchic_to_lexicographic_numbering() have been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member variable VectorizedArray::n_array_elements has
  been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member function CellAccessor::active() has been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member functions DataOut::first_cell() and DataOut::next_cell()
  have been removed.
  <br>
  (Daniel Arndt, 2021/06/02)
 </li>

 <li>
  Changed: Trilinos support in deal.II requires both deal.II and Trilinos to be
  configured with MPI support.
  <br>
  (Daniel Arndt, 2021/06/02)
 </li>

 <li>
  Removed: The deprecated build_matrices() member functions
  in various multigrid transfer classes have been removed.
  <br>
  (Daniel Arndt, 2021/06/01)
 </li>

 <li>
  Removed: The deprecated versions of DoFTools::count_dofs_per_component()
  and DoFTools::count_dofs_per_block() have been removed.
  <br>
  (Daniel Arndt, 2021/06/01)
 </li>

 <li>
  Removed: The deprecated overloads for DoFTools::extract_hanging_node_dofs(),
  and DoFTools::extract_dofs() have been removed.
  <br>
  (Daniel Arndt, 2021/05/28)
 </li>

 <li>
  Removed: Deprecated parallel::CellWeights member functions have been removed.
  <br>
  (Daniel Arndt, 2021/05/27)
 </li>

 <li>
  Removed: The overloads in CUDAWrappers::FEEvaluation that take a local dof index
  or a quadrature point as argument have been removed.
  Use the ones that don't use these arguments in their interface instead.
  <br>
  (Daniel Arndt, 2021/05/26)
 </li>

 <li>
  Removed: The deprecated parallel::Triangulation class has been removed.
  Use parallel::TriangulationBase instead.
  <br>
  (Daniel Arndt, 2021/05/26)
 </li>

 <li>
  Removed: The deprecated member field
  MatrixFree::AdditionalData::level_mg_handler.
  <br>
  (Peter Munch, 2021/05/25)
 </li>

 <li>
  Removed: The deprecated
  parallel::distributed::CellDataTransfer::CoarseningStrategies struct has been removed.
  Use AdaptationStrategies::Coarsening instead.
  <br>
  (Daniel Arndt, 2021/05/25)
 </li>

 <li>
  Removed: The deprecated member functions
  DoFHandler::locally_owned_dofs_per_processor(),
  DoFHandler::n_locally_owned_dofs_per_processor(), and
  DoFHandler::n_locally_owned_mg_dofs_per_processor() have been removed.
  <br>
  (Daniel Arndt, 2021/05/24)
 </li>

 <li>
  Deprecated: The template arguments of the following classes have been
  changed to avoid the legacy argument `DoFHandlerType`:
  <ul>
    <li> `SolutionTransfer<dim, VectorType, DoFHandlerType> -> SolutionTransfer<dim, VectorType, spacedim>`
    <li> `parallel::distributed::SolutionTransfer<dim, VectorType, DoFHandlerType> -> parallel::distributed::SolutionTransfer<dim, VectorType, spacedim>`
    <li> `Functions::FEFieldFunction<dim, DoFHandlerType, VectorType> -> Functions::FEFieldFunction<dim, VectorType, spacedim>`
    <li> `DataOut<dim, DoFHandlerType> -> DataOut<dim, spacedim>`
    <li> `DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim> -> DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>`
    <li> `DataOutFaces<dim, DoFHandlerType> -> DataOutFaces<dim, spacedim>`
    <li> `DataOutRotation<dim, DoFHandlerType> -> DataOutRotation<dim, spacedim>`
  </ul>
  Please change your code accordingly.
  <br>
  If for some reason, you need a code that is compatible with deal.II
  9.3 and the subsequent release, a Legacy namespace has been introduced
  with aliases of these classes with their original interface. You can
  make the following substitutions to your code for each of the affected
  classes:
  <ul>
    <li>X &rarr; Legacy::X
  </ul>
  To perform this substitution automatically, you may use a *search and
  replace* script like the following made for the *bash* shell:
  @code{.sh}
  classes=(SolutionTransfer parallel::distributed::SolutionTransfer Functions::FEFieldFunction DataOut DataOut_DoFData DataOutFaces DataOutRotation)
  for c in \${classes[@]}; do
    find /path/to/your/code -type f -exec sed -i -E "/(\w\${c}|\${c}[^<]|Particles::\${c}|distributed::\${c}|^\s*(\/\/|\*))/! s/\${c}/Legacy::\${c}/g" {} \;
  done
  @endcode
  (Marc Fehling, 2020/11/21)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="933-940-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: The namespace Utilities::MPI::LargeCount provides functionality
  to do MPI operations with more than 2^31 number of objects.
  <br>
  (Timo Heister, Sean Ingimarson, Penfei Jia, Jiaqi Zhang, 2022/05/29)
 </li>

 <li>
  New: Several classes have been added to the NonMatching namespace for using cut
  finite element methods (CutFEM):
  <ul>
    <li>
      MeshClassifier identifies how active cells and faces are located relative
      to the zero contour of a level set function, by associating one of the
      values of the enum LocationToLevelSet: {inside, outside, intersected}, with
      each cell/face.
    </li>
    <li>
      QuadratureGenerator creates immersed quadrature rules over the 3 regions of
      a BoundingBox defined by the sign of a level set function.
      FaceQuadratureGenerator creates the same type of quadratures but over a
      face of a BoundingBox. DiscreteQuadratureGenerator and
      DiscreteFaceQuadratureGenerator create the same type of quadratures over a
      cell/face when the level set function lies in a finite element space.
    </li>
    <li>
      FEImmersedSurfaceValues is a FEValues-like class for evaluating real space
      values based on an ImmersedSurfaceQuadrature, i.e., for integrating over
      the zero contour in a cell.
    </li>
    <li>
      NonMatching::FEValues simplifies assembling by creating dealii::FEValues
      and FEImmersedSurfaceValues objects initialized with the 3 different
      immersed quadratures of a cell. Analogously, NonMatching::FEInterfaceValues
      creates dealii::FEInterfaceValues objects initialized with immersed
      quadrature rules over a face.
    </li>
  </ul>
  <br>
  (Simon Sticko, Maximilian Bergbauer, Martin Kronbichler, 2022/05/29)
 </li>

 <li>
  New: Added support for the CGAL library (www.cgal.org). The following features
  are supported:
  <ul>
    <li>Conversion from deal.II to CGAL point types and vice-versa</li>
    <li>Conversion from deal.II cell types to CGAL::Surface_mesh and vice-versa</li>
    <li>Conversion from deal.II Triangulation to CGAL::Surface_mesh and vice-versa</li>
    <li>Insertion of deal.II points in CGAL triangulation types</li>
    <li>Conversion from CGAL::Triangulation_2 to Triangulation<dim, 2></li>
    <li>Conversion from CGAL::Triangulation_3 to Triangulation<dim, 3></li>
    <li>Conversion from CGAL::Surface_mesh or CGAL::Polyhedron_3 to Triangulation<2, 3></li>
  </ul>
  <br>
  (Marco Feder, Luca Heltai, 2022/04/26)
 </li>

 <li>
  New: Added a new function DoFRenumbering::support_point_wise()
  which renumbers DoFs so that DoFs with identical support points
  are now numbered adjacently.
  <br>
  (David Wells, 2021/03/11)
 </li>

 <li>
  New: A tutorial, Step-85, which shows how to use functionality in the
  NonMatching namespace to solve a PDE with the cut finite element method,
  when the domain is described by a level set function.
  <br>
  (Simon Sticko, 2022/01/07)
 </li>

 <li>
  New: Geometry GridGenerator::pipe_junction() creates an intersection of
  arbitrary truncated cones in 3D. It can be used to generate T-pipes,
  Y-pipes, corner pieces and other fittings.
  <br>
  (Marc Fehling, Maximilian Bergbauer, Peter Munch, Martin Kronbichler, 2021/12/01)
 </li>

 <li>
  New: The step-49 tutorial now has an extended section about generating meshes from Gmsh with a more interesting example.
  <br>
  (Timo Heister, Sean Ingimarson, 2021/10/08)
 </li>

 <li>
  New: Added an interface to p4est_search.h (>=v2.2) to find the MPI
  rank of cells owning remote points without communication. This
  works in the special (but not uncommon) case that no nontrivial
  manifold is attached to the triangulation, i.e., the triangulation
  has only FlatManifold attached to it.
  <br>
  (Konrad Simon, 2021/09/24)
 </li>

 <li>
  New: The tutorial program step-82 (local discontinuous Galerkin method
  for the bi-Laplacian problem) has been added.
  <br>
  (Diane Guignard, 2021/09/17)
 </li>

 <li>
  Added: Evaluation and integration of hessians both on cells and on faces is now
  available in the matrix-free framework.
  <br>
  (Maximilian Bergbauer, Martin Kronbichler, Peter Munch, 2021/09/06)
 </li>

 <li>
  New: VTU files have traditionally been used for visualization
  purposes, but more recently they have also been used for more general
  postprocessing tasks since they can easily be imported into Python,
  for example. By setting the DataOutBase::VtkFlags::physical_units
  variable and attaching such a flags object to a DataOut,
  Particles::DataOut, or MatrixOut object, one can now specify the
  physical units associated with variables that are written into VTU
  files, so that they can be parsed by scripts that do further
  postprocessing.
  <br>
  (Wolfgang Bangerth, 2021/07/09)
 </li>

 <li>
  Improved: DataOut_DoFData and the derived classes create internally
  a copy of vectors attached via add_data_vector() and are now
  responsible for updating ghost values so that users do not
  need to call update_ghost_values() and potentially zero_out_ghost_values().
  <br>
  (Peter Munch, Magdalena Schreter, 2021/05/27)
 </li>

 <li>
  New: The new postprocessing class DataOutResample is similar to
  DataOut but allows to interpolate the result onto a second
  triangulation.
  <br>
  (Peter Munch, Martin Kronbichler, Magdalena Schreter, 2021/06/08)
 </li>

 <li>
  Improved: The internal particle handler storage structure has been reworked to
  allow for significantly faster particle sorting and iteration as well as for
  future improvements and extensions.
  <br>
  (Rene Gassmoeller, Bruno Blais, Martin Kronbichler, Peter Munch, 2021/06/06)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="933-940-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: VectorTools::project can now be instantiated for complex numbers, bypassing undefined operations when compiling the parallel code path with MatrixFree evaluators.
  <br>
  (Pascal Kraft, Martin Kronbichler, 2022/05/30)
 </li>

 <li>
  Fixed: Several I/O routines are now handling large files/requests
  correctly. This includes Triangulation::save(), Triangulation::load(),
  DataOutInterface::write_vtu_in_parallel().
  <br>
  +(Timo Heister, Sean Ingimarson, Penfei Jia, Jiaqi Zhang, 2022/05/29)
 </li>

 <li>
  New: PreconditionChebyshev can now also estimate eigenvalues with the power
  method, rather than the Lanczos algorithm underlying the SolverCG
  estimation. This is useful for matrices that are slightly non-symmetric such
  that the eigenvalue estimate of SolverCG breaks down, but not too much so that
  the Chebyshev iteration is still a useful iteration.
  <br>
  (Martin Kronbichler, 2022/05/23)
 </li>

 <li>
  New: The Utilities::pack() and Utilities::unpack() functions are now
  optimized for the rather common case where the argument is of type
  `std::vector<T>` where `T` is a type that can be copied bit by bit
  (using `memcpy()`).
  <br>
  (Wolfgang Bangerth, 2022/05/05)
 </li>

 <li>
  Improved: The thread-local storage in FunctionParser has been reworked to avoid
  extra locking and unlocking calls. This makes evaluations of simple expressions
  about three times faster.
  <br>
  (David Wells, 2022/05/16)
 </li>

 <li>
  New: We added the possibility to average the contributions to a DoF from
  different cells in parallel::distributed::SolutionTransfer. For this purpose,
  set the corresponding flag in the constructor.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/05/14)
 </li>

 <li>
  New: The hdf5 DataOut now correctly supports very
  large output files with more than 2 billion nodes/vertices.
  <br>
  (Sean Ingimarson, Timo Heister, Jiaqi Zhang, 2022/05/10)
 </li>

 <li>
  New: The class SolverCG now supports the interleaving of vector operations
  with the matrix-vector product. The prerequisite is an associated `MatrixType`
  class to provide a `vmult` class with two `std::function` objects to specify
  the operation before and after the matrix-vector product, and a
  `PreconditionerType` class that provides either a function `apply` to apply
  the preconditioner action on a single element or and `apply_to_subrange(const
  unsigned int, const unsigned int) const` that can selectively apply the
  precondition on a part of a vector. For optimal performance, the matrix and
  preconditioner types need to agree on suitable sizes for the sub-ranges.
  <br>
  (Dmytro Sashko, Martin Kronbichler, Peter Munch, 2022/05/09)
 </li>

 <li>
  New: The class SolverFlexibleCG implements a flexible variant of the conjugate
  gradient method that runs the so-called Polak-Ribiere update formula instead
  of the default Fletcher-Reeves formula of SolverCG.
  <br>
  (Martin Kronbichler, 2022/05/06)
 </li>

 <li>
  New: There is now a function ReferenceCell::get_midpoint_quadrature()
  that generalizes QMidpoint.
  <br>
  (Wolfgang Bangerth, 2022/05/05)
 </li>

 <li>
  New: There are now functions ReferenceCell::volume() and
  ReferenceCells::barycenter() that compute what their names suggest for
  the various kinds of reference cells deal.II supports.
  <br>
  (Wolfgang Bangerth, 2022/05/04)
 </li>

 <li>
  New: MGConstrainedDoFs now has a function add_boundary_indices()
  to add boundary indices on levels.
  <br>
  (Jiaqi Zhang, 2021/05/02)
 </li>

 <li>
  Changed: Improved: The internal indexing in FESystem::fill_fe_values() has been
  improved and is now about twice as fast.
  <br>
  (David Wells, 2022/04/30)
 </li>

 <li>
  New: The class DataPostprocessors::BoundaryIds can be used to
  visualize the boundary indicators used in a mesh.
  <br>
  (Wolfgang Bangerth, 2022/04/27)
 </li>

 <li>
  Improved: The step-44 tutorial now uses FE_DGP elements to discretize the
  pressure and dilatation fields, rather than FE_DGPMonomial. This improves
  the condition number of the linear system, necessitating fewer iterations
  to compute its iterative inverse.
  <br>
  (Jean-Paul Pelteret, 2022/04/26)
 </li>

 <li>
  Fixed: The build process now works well with MPI option on Windows systems.
  If configured properly, this allows exploring the capability of distributed computing
  on Windows system with deal.II.
  <br>
  (Daniel Sun, 2022/04/26)
 </li>

 <li>
  New: There is a new function DoFRenumbering::matrix_free_data_locality() that
  re-orders the indices of the unknowns in a DoFHandler to give optimal data
  locality of pre- and post-operations in matrix-free loops.
  <br>
  (Martin Kronbichler, Peter Munch, Alexander Roschlaub, Dmytro Sashko, 2022/04/21)
 </li>

 <li>
  Added: New function enable_abort_on_exception.
  <br>
  (Pasquale Claudio Africa, 2022/04/20)
 </li>

 <li>
  New: FEInterfaceValues now has functions get_*_function_values() to evaluate solution
  values across cell faces, such as the jump in (or average of) values,
  gradients, and Hessians.
  <br>
  (Jiaqi Zhang, 2022/04/11)
 </li>

 <li>
  Fixed: Renumbering of the support points of `FE_RaviartThomasNodal`, which previously had x and z coordinates switched for the second component with `dim=3`.
  <br>
  (Niklas Wik, 2022/04/04)
 </li>

 <li>
  Changed: The bundled version of muParser has been updated to version 2.3.3.
  <br>
  (David Wells, 2022/03/28)
 </li>

 <li>
  Fixed: deal.II, compiled with the bundled copy of muParser, can now be
  safely linked alongside other projects which use a separate copy of muParser.
  <br>
  (David Wells, 2022/03/25)
 </li>

 <li>
  Fixed: The functions in DoFRenumbering now work correctly (they do not deadlock)
  when a processor has zero locally owned degrees of freedom.
  <br>
  (David Wells, 2022/03/23)
 </li>

 <li>
  Improved: MeshWorker::CopyData gains a new template argument for the scalar
  type, which allows the collected local matrix or vector contributions to be
  stored as types other than `double`.
  <br>
  (Jean-Paul Pelteret, 2022/03/14)
 </li>

 <li>
  New: Add new wrappers ArborXWrappers::SphereIntersectPredicate and
  ArborXWrappers::SphereNearestPredicate to perform geometrical search on spheres
  using ArborX.
  <br>
  (Bruno Turcksin, 2022/03/11)
 </li>

 <li>
  New: There is now a function ReferenceCell::contains_point() that
  provides the generalized equivalent of a function in the GeometryInfo
  class that will eventually be deprecated.
  type.
  <br>
  (Wolfgang Bangerth, 2022/03/05)
 </li>

 <li>
  Improved: The Utilities::MPI::broadcast() has been optimized for
  types that are natively supported by MPI.
  <br>
  (Francesco Andreuzzi, Wolfgang Bangerth, 2022/03/04)
 </li>

 <li>
  New: There is now a template variable Utilities::MPI::is_mpi_type that
  can be used to query whether a data type is a natively supported MPI
  type.
  <br>
  (Wolfgang Bangerth, 2022/03/01)
 </li>

 <li>
  Changed: We now require the MPI library to support at least the MPI 3.0 standard.
  <br>
  (Pengfei Jia, Matthias Maier, Timo Heister, 2022/02/28)
 </li>

 <li>
  New: The template variable Utilities::MPI::mpi_type_id_for_type
  provides information about the `MPI_Datatype` that corresponds to a
  given type. This information is necessary in all MPI calls where one
  needs to provide an `MPI_Datatype` to identify the kinds of objects
  stored in a send or receive buffer.
  <br>
  (Wolfgang Bangerth, 2022/02/27)
 </li>

 <li>
  Fixed: We can now load certain msh files generated by Gmsh that contain some cells with inverted volume.
  <br>
  (Sean Ingimarson, Timo Heister, 2022/02/26)
 </li>

 <li>
  Improved: The asserts in parallel::distributed::SolutionTransfer have
  been improved. Before this change, an assert has been only called
  in the parallel case for distributed vectors for cells neighboring
  interprocess boundaries. Now, we check for all DoFs that the
  contributions from all cells have the same values. This condition
  might be not given if the solution has high gradients at hanging
  nodes.
  <br>
  (Peter Munch, Magdalena Schreter, 2022/02/25)
 </li>

 <li>
  Fixed: ParsedConvergenceTable::difference() was not working. This is now fixed.
  <br>
  (Luca Heltai, 2022/02/22)
 </li>

 <li>
  Fixed: The function parallel::distributed::copy_triangulation() used to
  copy the reference to the attached SolutionTransfer instance, which resulted
  in calling the same SolutionTransfer instance multiple times if the
  triangulation and the new triangulation are refined/coarsened independently.
  This has been fixed.
  <br>
  (Ivo Dravins, Peter Munch, 2022/02/22)
 </li>

 <li>
  Fixed: GridOut::write_msh() with gmsh api now works also with non standard combinations of boundary and manifold ids.
  <br>
  (Luca Heltai, 2022/02/20)
 </li>

 <li>
  Fixed: NonMatching::create_coupling_sparsity_pattern() and
  NonMatching::create_coupling_mass_matrix() were not working properly with
  hanging node constraints.
  <br>
  (Luca Heltai, 2022/02/18)
 </li>

 <li>
  New: AffineConstraints::add_entries_local_to_global() now supports also the case
  where columns don't use the same constraints as rows.
  <br>
  (Luca Heltai, 2022/02/18)
 </li>

 <li>
  Added: The DEAL_II_QUERY_GIT_INFORMATION() CMake macro now also records the
  date and time of the commit in the CMake variable GIT_TIMESTAMP.
  <br>
  (Matthias Maier, 2022/02/06)
 </li>

 <li>
  New: GridGenerator::generate_from_name_and_argument() now supports also GridGenerator::hyper_ball_balanced()
  <br>
  (Luca Heltai, 2022/02/17)
 </li>

 <li>
  Fixed: Patterns::Tools::Convert now supports also std::array types.
  <br>
  (Luca Heltai, 2022/02/16)
 </li>

 <li>
  New: Add new DistributedTree class based on ArborX. This class performs multiple
  kinds of geometric search on distributed objects: intersection of bounding boxes
  and intersection of bounding boxes with points. This class is the distributed
  equivalent of ArborXWrappers::BVH.
  <br>
  (Bruno Turcksin, 2022/02/15)
 </li>

 <li>
  Fixed: BlockDynamicSparsityPattern in step-55 will be initialized with
  the correct IndexSet for locally relevant dofs.
  <br>
  (Marc Fehling, 2022/02/10)
 </li>

 <li>
  Fixed: GridGenerator::hyper_ball_balanced() now works correctly with nonzero
  centers.
  <br>
  (David Wells, 2022/02/09)
 </li>

 <li>
  New: The FE_Nothing class can now return its (empty) constant modes.
  DoFTools::extract_constant_modes can now work with DoFHandlers that
  contain FE_Nothing.
  <br>
  (Sebastian Proell, 2022/02/07)
 </li>

 <li>
  Improved: To be consistent with the FEValues classes, one can now
  loop over all DoF and quadrature-point indices in a
  range-based loop with the help of the new functions dof_indices() and
  quadrature_point_indices().
  <br>
  (Peter Munch, 2022/02/06)
 </li>

 <li>
  Changed: The minimum CMake version required to configure and build deal.II
  has been raised to 3.3.0.
  <br>
  (Matthias Maier, 2022/02/06)
 </li>

 <li>
  Improved: In MatrixFreeTools::compute_diagonal a template argument VectorType
  has been introduced to be applicable for arbitrary vector types.
  <br>
  (Magdalena Schreter, Peter Munch, 2022/02/02)
 </li>

 <li>
  New: Added support for hessians to FE_NedelecSZ
  <br>
  (Sebastian Kinnewig, 2022/01/28)
 </li>

 <li>
  New: FE_Nedelec<2> now includes the function get_embedding_dofs(),
  which permits direct identification of DoF indices between two
  realizations of the FE_Nedelec cell.
  <br>
  (Jake Harmon, 2022/01/25)
 </li>

 <li>
  Fixed: MappingQ::transform_points_real_to_unit_cell has been made more robust,
  correcting an unsuitable initial guess for the Newton iteration to compute the
  inverse transformation in the case of bilinear/trilinear cells.
  <br>
  (Bruno Blais, Martin Kronbichler, 2022/01/22)
 </li>

 <li>
  New: MeshWorker::ScratchData and MeshWorker::CopyData now have hp-capabilities.
  Creation of a ScratchData object with an hp::FECollection and hp::QCollection
  enables those capabilities. A few specialized MeshWorker::ScratchData::reinit()
  were added to support face integration where the integration rule and mapping
  differs on either side of an interface. Additionally, some
  MeshWorker::CopyData::reinit() methods have been implemented to cater for the
  change in number of degrees-of-freedom between cells iterated over during
  assembly.
  <br>
  (Jean-Paul Pelteret, 2022/01/09)
 </li>

 <li>
  Fixed: GridGenerator::convert_hypercube_to_simplex_mesh() now correctly copies
  line manifold ids in 3D.
  <br>
  (David Wells, 2022/01/07)
 </li>

 <li>
  Changed: The step-43 tutorial program used `Amg_preconditioner` and
  `Mp_preconditioner` to denote the preconditioners for the top left and
  bottom right block of a mixed-Laplace matrix. These variable names
  were taken from step-31/step-32 where they made sense, but in the
  current context, the names do not reflect what is in the corresponding
  matrix blocks. They have been renamed to `top_left_preconditioner` and
  `bottom_right_preconditioner`, respectively.
  <br>
  (Wolfgang Bangerth, 2022/01/07)
 </li>

 <li>
  Fixed: Serializing empty sparsity patterns ran into a segmentation
  fault. This is now fixed.
  <br>
  (Lucas Myers, Wolfgang Bangerth, 2022/01/05)
 </li>

 <li>
  Fixed: deal.II has been updated to support the new OneAPI api interface
  introduced by the Intel Threading Building Blocks Library. Note that the
  threading support of the matrix-free backend will be disabled in this case:
  For now, the MatrixFree::AdditionalData::tasks_parallel_scheme control has
  no effect and simply defaults to the serial loop.
  <br>
  (Wolfgang Bangerth, Matthias Maier, 2022/01/06)
 </li>

 <li>
  Fixed: The FEInterfaceValues class is now made compatible with filtered
  iterators. Specifically, it is now possible to call FEInterfaceValues::reinit()
  with the active `cell` being a filtered iterator and the type returned by
  `cell->neighbor()` being some other non-filtered iterator type.
  <br>
  (Simon Sticko, Jean-Paul Pelteret, 2022/01/05)
 </li>

 <li>
  Improved: The class PreconditionRelaxation has been refactored.
  Similarly to PreconditionChebyshev, it can be now used
  stand-alone and users can attach their own preconditioners.
  Furthermore, users can specify multiple iteration steps,
  which is particularly useful if PreconditionRelaxation is
  used as smoother in a multigrid algorithm.
  <br>
  (Peter Munch, 2022/01/01)
 </li>

 <li>
  New: More "getter" functions have been added to MeshWorker::ScratchData,
  facilitating the construction of other FEValues-type objects with some of
  the same objects as those used in a ScratchData instance.
  <br>
  (Jean-Paul Pelteret, 2021/12/31)
 </li>

 <li>
  Fixed: The SparseDirectUMFPACK class had a bug where it did not
  correctly set up the data passed to the UMFPACK package if (i) the
  matrix is a block matrix, and (ii) the elements are
  complex-valued. This is now fixed.
  <br>
  (Kuljit S. Virk, 2021/12/30)
 </li>

 <li>
  New: The IteratorFilters::BoundaryIdEqualTo and IteratorFilters::ManifoldIdEqualTo
  filtered iterator predicates have been added to help select cell faces that
  have, respectively, one or more specified boundary ID(s) and manifold ID(s).
  <br>
  (Jean-Paul Pelteret, 2021/12/30)
 </li>

 <li>
  Fixed: GridGenerator::flatten_triangulation() did not copy information
  on the manifold ids of interior faces and edges. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2021/12/27)
 </li>

 <li>
  Fixed: GridGenerator::flatten_triangulation() did not copy information
  of edges of three-dimensional triangulations, even though it copied
  information from faces. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2021/12/26)
 </li>

 <li>
  Changed: The step-60 tutorial now works also with non-zero Dirichlet BCs and non-zero right-hand side
  <br>
  (Marco Feder, 2021/12/24)
 </li>

 <li>
  New: Build also a docker image with root user, to be used in github actions directly.
  <br>
  (Luca Heltai, 2021/12/24)
 </li>

 <li>
  Deprecated: The Threads::new_thread() functions have been
  deprecated. These functions return a Threads::Thread object, which is
  a class that had previously been deprecated, and should no longer be
  used: use the `std::thread` or `std::jthread` classes and related
  functionality instead.
  <br>
  (Wolfgang Bangerth, 2021/12/23)
 </li>

 <li>
  Fixed: The matrix-free FEFaceEvaluation would produce wrong results in case
  hanging nodes were present on faces with non-standard orientation. This is now
  fixed.
  <br>
  (Martin Kronbichler, 2021/12/20)
 </li>

 <li>
  Fixed: IndexSet::subtract_set would previously run in quadratic complexity in
  the number of ranges. This is now fixed.
  <br>
  (Martin Kronbichler, 2021/12/13)
 </li>

 <li>
  New: GridOut::write_msh() now also works for triangular and tetrahedral meshes.
  <br>
  (Wolfgang Bangerth, 2021/12/06)
 </li>

 <li>
  New: The DataPostprocessorInputs::CommonInputs class now also
  stores the face number that is currently being worked on when
  using a DataOutFaces object. This allows accessing face-related
  data like the boundary id from a DataPostprocessor.
  <br>
  (Wolfgang Bangerth, 2021/12/03)
 </li>

 <li>
  Changed: The function TrilinosWrappers::Vector::is_non_negative() now also works in parallel.
  <br>
  (Justin O'Connor, 2021/12/02)
 </li>

 <li>
  Fixed: MatrixFreeOperators::CellwiseInverseMassMatrix with template argument
  `-1` for the polynomial degree now chooses the pre-compiled code with fast
  tensor-product algorithms if available, rather than a slower general-purpose
  code.
  <br>
  (Martin Kronbichler, 2021/11/29)
 </li>

 <li>
  New: Namespace Physics::VectorRelations features functions to compute
  angles between (spatial) vectors.
  <br>
  (Marc Fehling, 2021/11/29)
 </li>

 <li>
  New: The SD::BatchOptimizer::copy_from() function can be used to duplicate the
  data from one batch optimizer instance to another.
  <br>
  (Jean-Paul Pelteret, 2021/11/27)
 </li>

 <li>
  Fixed: The step-22 tutorial now uses FEValuesExtractors also for the right-hand side
  <br>
  (Marco Feder, 2021/11/23)
 </li>

 <li>
  New: Utilities::MPI::create_mpi_data_type_n_bytes() helps with large MPI communication and IO.
  <br>
  (Timo Heister, 2021/11/17)
 </li>

 <li>
  New: In the same spirit as provided by the C++20
  [range adaptors](https://en.cppreference.com/w/cpp/ranges) feature, it is now
  possible to "filter" ranges over the active cells of Triangulation
  and DoFHandler objects.
  <br>
  This allows to replace:
  @code
    DoFHandler<dim> dof_handler;
    ...
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            ...do the local integration on 'cell'...;
          }
      }
  @endcode
  by:
  @code
    DoFHandler<dim> dof_handler;
    ...
    for (const auto &cell :
            dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      {
        fe_values.reinit (cell);
        ...do the local integration on 'cell'...;
      }
  @endcode
  Here, the `operator|` is to be interpreted in the same way as is done in
  the [range adaptors](https://en.cppreference.com/w/cpp/ranges) feature
  that is part of [C++20](https://en.wikipedia.org/wiki/C%2B%2B20). It has
  the same meaning as the `|` symbol on the command line: It takes what is
  on its left as its inputs, and filters and transforms to produce some
  output. In the example above, it "filters" all of the active cell iterators
  and removes those that do not satisfy the predicate -- that is, it produces
  a range of iterators that only contains those cells that are both active
  and locally owned.
  <br>
  (Wolfgang Bangerth, 2021/11/03)
 </li>

 <li>
  New: The GridOut::write_gnuplot() function can now deal with
  tetrahedral meshes, as well as with wedges and pyramids.
  <br>
  (Wolfgang Bangerth, 2021/11/15)
 </li>

 <li>
  New: The GridIn::read_comsol_mphtxt() function adds the ability to
  read `.mphtxt` files generated by COMSOL.
  <br>
  (Wolfgang Bangerth, 2021/11/15)
 </li>

 <li>
  New: The new class Functions::SignedDistance::Ellipsoid computes
  a signed distance function for ellipsoids. It is implemented for 1D and
  2D ellipsoids (ellipses).
  <br>
  (Nils Much, 2021/11/05)
 </li>

 <li>
  New: One can now output simplex, pyramid, and wedge cells in GNUPLOT
  format via DataOut and other classes that are based on DataOutBase and
  call DataOutBase::write_gnuplot().
  <br>
  (Wolfgang Bangerth, 2021/11/03)
 </li>

 <li>
  Changed: The DataOutBase::write_deal_II_intermediate() function had a
  bug when outputting data that corresponds to cells other than
  hypercubes. This is now fixed. In the process, the format version for
  this function's output has been incremented from 3 to 4, which also
  includes changing the order in which vertices are outputted.
  <br>
  The DataOutReader class was also updated to match the new output
  format.
  <br>
  (Wolfgang Bangerth, 2021/10/25)
 </li>

 <li>
  New: The ReferenceCell class now has input and output operators (`>>`
  and `<<`).
  <br>
  (Wolfgang Bangerth, 2021/10/25)
 </li>

 <li>
  Fixed: GridIn::read_unv() can now ignore any sections before 2411. In other
  words, the unv file can start with any number sections (before 2411) as long as
  they have a proper format (must begin and end with a "-1").
  <br>
  (Vachan Potluri, 2021/10/22)
 </li>

 <li>
  New: There is now a function ReferenceCell::vertex() that returns the
  coordinates of a vertex of a reference cell.
  <br>
  (Wolfgang Bangerth, 2021/10/13)
 </li>

 <li>
  Fixed: The DataOutFaces class can now also be used with simplex
  meshes. It previously aborted with exceptions.
  <br>
  (Wolfgang Bangerth, 2021/10/13)
 </li>

 <li>
  Improved: PreconditionChebyshev can now run the vector updates on a part of
  the vector for better data locality. The prerequisite is an associated
  `MatrixType` class to provide a `vmult` class with two `std::function` objects
  to specify the operation before and after the matrix-vector product.
  <br>
  (Martin Kronbichler, 2021/10/12)
 </li>

 <li>
  Fixed: The KellyErrorEstimator for 1D can now handle the hp case.
  <br>
  (Fabian Castelli, 2021/10/12)
 </li>

 <li>
  Fixed: Deleted duplicate occurrence of VectorTools::interpolate_boundary_values()
  in setup_dofs() of step-45.
  <br>
  (Raghunandan Pratoori, 2021/10/11)
 </li>

 <li>
  Changed: The 3D implementation of GridTools::rotate() now accepts unit
  vectors as rotation axes.
  <br>
  (Marc Fehling, 2021/10/11)
 </li>

 <li>
  New: VectorizedArray can now be constructed from an initializer list.
  <br>
  (Peter Munch, 2021/10/08)
 </li>

 <li>
  Fixed: The FEPointEvaluation class always evaluated a solution starting at the component 0 of the selected vector-valued base element, even if a different component was selected by the caller. This has been fixed to correctly start evaluating at the selected component.
  <br>
  (Rene Gassmoeller, 2021/10/01)
 </li>

 <li>
  Fixed: The payload for Trilinos-based linear operators is now able to use its
  own type as an exemplar operator. This means that one can now use linear
  operators based on Trilinos matrices as exampler operators when initializing
  a LinearOperator for a Trilinos preconditioner.
  <br>
  (Wenyu Lei, Jean-Paul Pelteret, 2021/09/12)
 </li>

 <li>
  Fixed: ParticleHandler::insert_particles() runs into a deadlock in parallel
  when a processor is not inserting any particles. This is now fixed.
  <br>
  (Sebastian Fuchs, 2021/08/31)
 </li>

 <li>
  Fixed: There was a bug in the default constructor of
  TimeStepping::EmbeddedExplicitRungeKutta(). Some parameters were left
  uninitialized which triggered a segmentation fault when calling
  `evolve_one_time_step()`.
  <br>
  (Praveen Chandrashekar, Bruno Turcksin, 2021/08/30)
 </li>

 <li>
  Changed: For 1d cells, the following code did not compile:
  @code
      for (const auto &face : cell->face_iterators())
        if (face->at_boundary())
          {
            face->set_boundary_id(1);
          }
  @endcode
  That was because, unlike in the 2d and 3d case, the class that
  describes the faces of 1d cells had a `set_boundary_id()` function
  that was not marked as `const`. This has now been fixed by adding the
  `const` specification for this function, along with the one on the
  `set_all_boundary_ids()` function.
  <br>
  (Wolfgang Bangerth, Tyler Anderson, 2021/08/21)
 </li>

 <li>
  New: The hp::FECollection::hp_vertex_dof_identities() function (and
  related functions for lines and quads) allows for the computation of
  multiway computations of DoF identities in the hp-context.
  <br>
  (Wolfgang Bangerth, 2021/08/20)
 </li>

 <li>
  New: There is now an overload for
  Particles::ParticleAccessor::set_properties() that takes a Tensor
  as argument. It is used in step-19.
  <br>
  (Wolfgang Bangerth, 2021/08/19)
 </li>

 <li>
  Fixed: Exceptions thrown via `AssertThrow` did not have a stacktrace
  attached to them, unlike exceptions thrown with `Assert`. This is now
  fixed.
  <br>
  (Wolfgang Bangerth, Paras Kumar, 2021/08/20)
 </li>

 <li>
  Fixed: It was not possible to default-construct objects of type
  `DoFHandler::face_iterator` for `dim==1`. As a consequence, one
  could also not call `cell->face_iterators()` for the cells of 1d
  DoFHandler objects. This is now fixed.
  <br>
  (Wolfgang Bangerth, Tyler Anderson, 2021/08/19)
 </li>

 <li>
  New: MeshWorker::ScratchData is now able to evaluate the jumps in, and averages
  of, finite element functions and their derivatives across an interface.
  <br>
  (Jean-Paul Pelteret, 2021/08/14)
 </li>

 <li>
  New: FEInterfaceViews now has functions get_*_function_values() and
  get_*_function_values_from_local_dof_values() to evaluate solution
  values across cell faces, such as the jump in (or average of) values,
  gradients, and Hessians. With these functions, FEInterfaceValues
  can now be used with automatic differentiation.
  <br>
  (Jiaqi Zhang, 2021/08/11)
 </li>

 <li>
  New: The functions FEFaceValuesBase::face_number(), FEInterfaceValues::get_cell()
  and FEInterfaceValues::face_number() have been implemented to provide
  introspection as to which entities these objects have been initialized for.
  <br>
  (Jean-Paul Pelteret, 2021/08/06)
 </li>

 <li>
  New: The function FEInterfaceValues::dof_indices() has been added to make
  writing range-based for loops over all local interface degrees of freedom much
  simpler.
  <br>
  (Jean-Paul Pelteret, 2021/08/05)
 </li>

 <li>
  New: Diagnostic functions AffineConstraints::n_identities() and
  AffineConstraints::n_inhomogeneities() to analyse your constraints.
  <br>
  (Marc Fehling, 2021/08/04)
 </li>

 <li>
  Fixed: Use tolerances consistently in RemotePointEvaluation,
  VectorTools::point_values(), and GridTools::distributed_compute_point_locations().
  <br>
  (Peter Munch, Magdalena Schreter, 2021/07/31)
 </li>

 <li>
  New: QWitherdenVincentSimplex now implements even-order rules
  in addition to the standard Gauss-like odd-order rules.
  <br>
  (David Wells, 2021/07/29)
 </li>

 <li>
  Fixed: AlignedVector objects for which
  AlignedVector::replicate_across_communicator() had previously been
  called, could not be move-constructed or move-assigned: The resulting
  object was usable, but destruction of the target object led to memory
  errors that aborted the program. This is now fixed. As a side effect
  of the patch, memory consumption of AlignedVector objects has also
  been shrunk for the most common case where
  AlignedVector::replicate_across_communicator() has not actually been called.
  <br>
  (Wolfgang Bangerth, 2021/07/28)
 </li>

 <li>
  New: The new function TableBase::clear() allows to
  empty a Table object.
  <br>
  (Jiaqi Zhang, 2021/07/21)
 </li>

 <li>
  Fixed: MappingFEField now works with simplex faces correctly.
  <br>
  (David Wells, 2021/07/19)
 </li>

 <li>
  New: Added a class QIteratedSimplex for building composite simplex quadrature rules.
  <br>
  (David Wells, 2021/07/15)
 </li>

 <li>
  New: ParticleAccessor::get_local_index() is a new function returning particle
  indices local to each MPI process, which is useful for direct access into
  arrays with temporary results associated to particles.
  <br>
  (Martin Kronbichler, 2021/07/06)
 </li>

 <li>
  Fixed: step-34 had a minor memory leak whereby the function
  `get_singular_quadrature()` allocated some memory the
  first time it was called, and then never released it. This is now
  fixed. Separately, this function was not thread-safe, and now it is.
  <br>
  (Wolfgang Bangerth, 2021/07/03)
 </li>

 <li>
  New: When calling `cell->vertex_dof_index(...)` on a DoFHandler that
  uses the hp capability, then one needed to provide the last argument
  indicating what the active fe index is. But for cells, there can only
  be one active fe index, so leaving the last argument at its default
  value is now interpreted as using the active fe index of the cell. The
  argument still needs to be provided if the `vertex_dof_index()`
  function is called on faces, edges, or vertices.
  <br>
  (Wolfgang Bangerth, 2021/07/02)
 </li>

 <li>
  New: The new function Triangulation::n_global_coarse_cells() allows to
  query the global number of coarse cells. This function is particularly useful
  for parallel::fullydistributed::Triangulation, since in this case the coarse
  cells might differ on each process.
  <br>
  (Peter Munch, 2021/06/26)
 </li>

 <li>
  New: The new unified function MGTwoLevelTransfer::reinit() selects automatically
  if MGTwoLevelTransfer::reinit_geometric_transfer() or
  MGTwoLevelTransfer::reinit_polynomial_transfer() is needed.
  <br>
  (Peter Munch, 2021/06/23)
 </li>

 <li>
  New: The new function CellAccessor::is_ghost_on_level() allows to
  check if a cell is a ghost cell on a multigrid level.
  <br>
  (Peter Munch, 2021/06/23)
 </li>

 <li>
  Improved: i) Standardized the format of the references in the file
               doc/doxygen/references.bib using in particular the convention
               "M. Mo" or "R.B. Kellogg" for names and the full name for journals.
           ii) Fixed the display of some titles (e.g., upper case for proper noun).
          iii) Added an URL entry for the references with a DOI.
  <br>
  (Diane Guignard, 2021/06/22)
 </li>

 <li>
  Fixed: The Gauss-like rules implemented by QGaussSimplex now integrate
  2 * n_points_1D - 1 degree polynomials exactly.
  <br>
  (David Wells, 2021/06/22)
 </li>

 <li>
  Fixed: The step-21 program called a function that was disabled on
  Microsoft Windows to work around a bug in the Microsoft Visual Studio
  compiler. This is now fixed and the program should run successfully
  again on Windows.
  <br>
  (Wolfgang Bangerth, 2021/06/08)
 </li>

 <li>
  Fixed: MappingQCache::get_vertices() would previously always return the actual
  vertices, ignoring a possible displacement that is possible with this
  class. This is now fixed.
  <br>
  (Martin Kronbichler, 2021/06/07)
 </li>

 <li>
  Improved: An assert has been added in the constructor of
  FEPointEvaluation to prevent accessing non-existing components.
  <br>
  (Judith Pauen, 2021/06/05)
 </li>

 <li>
  New: The new function Utilities::MPI::reduce() allows to reduce
  arbitrary types with a user-specified binary operation.
  <br>
  (Peter Munch, 2021/05/27)
 </li>

 <li>
  New: Fixed a bug in clear user data for standard Triangulation.
  <br>
  (Nicola Giuliani, 2021/05/25)
 </li>

 <li>
  Added conformity tests for several vector elements
  in 2D/3D as well as minor fixes for 2d Nedelec and
  Raviart-Thomas to make them work on non-standard meshes.
  <br>
  (Konrad Simon, 2021/05/08)
 </li>

 <li>
  New: The class MGTransferMatrixFree::build() now also
  accepts an optional function for initializing the internal level vectors.
  This is useful if one uses the transfer operators in the context of
  smoothers that are built around MatrixFree objects.
  <br>
  (Peter Munch, 2021/03/09)
 </li>

 <li>
  New: You can now perform range-based iterations on hp::FECollection, hp::QCollection,
  and hp::MappingCollection objects.
  <br>
  (Peter Munch, 2021/10/27)
 </li>

</ol>

*/
