// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
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
@page changes_between_9_6_0_and_9_7_0 Changes between Version 9.6.0 and 9.7.0

<p>
This is the list of changes made since the last release of deal.II.
All entries are signed with the names of the authors.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="960-970-incompatible"></a>
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
  Changed: The hp-variants of KellyErrorEstimator::estimate now take
  hp::MappingCollection objects. The versions that require Mapping objects
  are now deprecated.
  <br>
  (Marc Fehling, 2025/06/20)
 </li>

 <li>
  Changed: Template parameters for ArborXWrappers::BVH and
  ArborXWrappers::DistributedTree vary between the ArborX 1.X and 2.X series
  <br>
  (Bruno Turcksin, 2025/06/17)
 </li>

 <li>
  Changed: Given that the MatrixOut class now allows to create "sparse"
  representations (see the "Minor changes" section of this page), the
  default for what kind of output is being produced has now also been
  changed to create this "sparse" output, rather than the traditional
  "dense" output. You can switch back to the traditional way in which
  zero matrix entries are also represented in the visualization output
  using the MatrixOut::Options::create_sparse_plot variable.
  <br>
  (Wolfgang Bangerth, 2025/04/21)
 </li>

 <li>
  Changed: deal.II now requires Trilinos version 13.2.
  <br>
  (Timo Heister, 2025/03/21)
 </li>

 <li>
  Changed: Every deal.II header file includes the deal.II header
  `deal.II/base/config.h` that contains some basic configuration
  information such as what deal.II version is being used, whether
  deal.II was configured with Trilinos/PETSc/... support, etc. In turn,
  `config.h` included `deal.II/base/types.h` and
  `deal.II/base/numbers.h`, making the contents of these latter two
  files available throughout deal.II.
  For conceptual reasons, for example to support C++20 modules in the
  long run, the `deal.II/base/config.h` file must stand separate from
  the rest of the deal.II header files. As a consequence, it can no
  longer have includes for other deal.II header files (namely,
  `deal.II/base/types.h` and `deal.II/base/numbers.h` as mentioned
  above), and no longer does so now. If your code requires knowledge of
  the contents of `deal.II/base/types.h` and `deal.II/base/numbers.h`,
  you will have to explicitly include these files like you do with any
  other deal.II header file already.
  <br>
  (Wolfgang Bangerth, 2024/01/16)
 </li>

 <li>
  Updated: We now require Boost version 1.74 or newer.
  <br>
  (Matthias Maier, 2025/01/30)
 </li>

 <li>
  Changed: The function FETools::compute_face_embedding_matrices() used
  to take its second argument as a fixed-size (C-style) array of
  matrices of size `GeometryInfo<dim>::max_children_per_face`. This has
  been changed to be an object of type `ArrayView<FullMatrix<number>>`
  to allow for storage of that array in other ways, say as an object of
  type `std::vector<FullMatrix<double>>`. To convert old code, you may
  have to explicitly call the function with its template arguments, say
  as `FETools::compute_face_embedding_matrices<dim,double,spacedim(...)`
  or to convert the second argument from its actual type to the
  ArrayView object by wrapping it in a call to make_array_view().
  <br>
  (Wolfgang Bangerth, 2024/01/16)
 </li>

 <li>
  Fixed: Previously, FiniteElement::has_generalized_support_points() and
  FiniteElement::has_face_support_points() returned false for the FE_Nothing
  element. Now, FE_Nothing::has_generalized_support_points() and
  FE_Nothing::has_face_support_points() correctly return true, as the empty
  arrays returned by FE_Nothing::get_generalized_support_points() and
  FE_Nothing::get_unit_face_support_points() accurately describe the support
  points of the element (i.e., they don't have any, as there are no degrees of
  freedom).
  <br>
  (Oreste Marquis, 2024/08/21)
 </li>

 <li>
  Fixed: The FiniteElement::has_support_points() function is poorly
  named because it does not return *whether* an element has support
  points, but whether it *implements* the
  FiniteElement::get_unit_support_points() function correctly. Because
  of this misunderstanding, it returned `false` for the FE_Nothing
  element. This is now fixed: FE_Nothing::has_support_points() now
  returns `true`, because the empty array returned by
  FE_Nothing::get_unit_support_points() correctly describes the support
  points the element has (namely: it does not have any, as there are no
  degrees of freedom).
  <br>
  (Wolfgang Bangerth, 2024/08/20)
 </li>

 <li>
  Changed: The header file `deal.II/grid/tria.h` used to automatically
  include the header file `deal.II/grid/tria_description.h` that
  declares classes such as CellData and SubCellData. But this led to a
  circular set of inclusions because `deal.II/grid/tria_description.h`
  also included `deal.II/grid/tria.h`, something that C++ allows but
  that leads to other problems. As a consequence, `deal.II/grid/tria.h`
  now no longer includes `deal.II/grid/tria_description.h`. If your
  program uses the CellData, SubCellData, or similar classes previously
  declared in `deal.II/grid/tria_description.h`, then you may want to
  explicitly include `deal.II/grid/cell_data.h` (for the
  CellData and SubCellData classes) or `deal.II/grid/tria_description.h`
  (for the classes in namespace TriaDescription) in your program.
  The file `deal.II/grid/cell_data.h` is new. As a consequence, if
  you need to ensure that your code compiles with both deal.II 9.6
  and 9.7, you can't directly include it. However,
  `deal.II/grid/tria_description.h` now includes it, so a backward
  compatible solution is to only include that file, even if you
  only need the CellData and SubCellData classes but not the ones
  in namespace TriaDescription.
  <br>
  (Wolfgang Bangerth, 2024/05/16)
 </li>

 <li>
  Changed: The SolverControl, ReductionControl, and
  IterationNumberControl classes have constructors that take two boolean
  flags `log_history` and `log_result`. The first of these represents
  whether each time these classes are asked to evaluate progress of an
  iterative solver, it should output information about the current
  number of the iteration and the current residual to the `deallog`
  object. The latter flag determines whether this kind of output should
  be generated once the solver has been determined to have either
  succeeded or failed, i.e., at the end of the iteration. By default,
  these flags were set to `false` and `true`, respectively.
  Creating and formatting this kind of output turns out to be
  surprisingly expensive. This is a nuisance because by default, the
  `deallog` variable (via the constructor arguments of the LogStream
  class) is instructed to simply ignore whatever information is sent to
  it, rather than putting the output onto the screen or into a log file.
  As a consequence, the defaults of the constructor arguments of the
  three classes mentioned above have been changed from `false` and
  `true` to `false` and `false`. If you have instructed `deallog` to
  pass information given to it to the screen or into an output file, you
  can always explicitly also set the constructor arguments of
  SolverControl, ReductionControl, or IterationNumberControl objects to
  retain the previous behavior.
  <br>
  (Wolfgang Bangerth, 2025/01/21)
 </li>

 <li>
  Changed: Mapping::transform_points_real_to_unit_cell() used to
  designate points for which it could not successfully find reference
  coordinates by setting the first vector component to
  `std::numeric_limits<double>::infinity()`. Unfortunately, on some
  platforms, the use of infinities leads to dramatically slower code
  execution. As a consequence, we now use the (finite) value
  `std::numeric_limits<double>::lowest()` (somewhere around `-1e308`) to
  denote invalid values.
  <br>
  (Wolfgang Bangerth, 2024/12/17)
 </li>

 <li>
  Updated: The interface of Rol::VectorAdaptor class now uses `ROL::Ptr`.
  It is a shared pointer wrapper for either `Teuchos::RCP` or
  `std::shared_ptr`, and can be specified while configuring Trilinos.
  See `Trilinos/packages/rol/cmake/BuildOptions.cmake` for details.
  <br>
  (Marc Fehling, 2024/11/28)
 </li>

 <li>
  Changed: The minimum version for Trilinos has been bumped to 12.14.1.
  <br>
  (Marc Fehling, 2024/11/28)
 </li>

 <li>
  Changed: FEInterfaceValues::normal() has been renamed FEInterfaceValues::normal_vector().
  <br>
  (Peter Munch, 2024/11/23)
 </li>

 <li>
  Changed: The deprecated entries Text and ShortText in OutputStyle
  have been removed.
  <br>
  (Peter Munch, 2024/08/21)
 </li>

 <li>
  Changed: The minimum version of p4est compatible with deal.II has been
  increased to 2.2.
  <br>
  (David Wells, 2024/08/15)
 </li>

 <li>
  Removed: The CMake option `DEAL_II_COMPILE_EXAMPLES` has been removed.
  As a consequence, examples will no longer be compiled if
  `DEAL_II_COMPONENT_EXAMPLES` is set to true.
  <br>
  You can use the testsuite to compile all examples at once. First, set up
  the `setup_tests_examples` target, i.e., with `make setup_tests_examples`.
  Then compile and run the examples with `ctest -R examples`.
  <br>
  (Marc Fehling, 2024/08/12)
 </li>

 <li>
  Changed: The configuration option DEAL_II_WITH_COMPLEX_VALUES that
  enables support for linear algebra classes to be used with std::complex
  is now disabled by default.
  <br>
  (Timo Heister, 2024/08/11)
 </li>

 <li>
  Removed: All support for CUDAWrappers and CUDA-related macros have been removed.
  GPU support is provided through Kokkos.
  <br>
  (Daniel Arndt, 2024/08/11)
 </li>

 <li>
  Removed: Deprecations from the 9.5 release have been removed.
  See https://github.com/dealii/dealii/pull/17444 for details.
  <br>
  (Daniel Arndt, 2024/08/07)
 </li>

 <li>
  Removed: step 52 has been removed.
  step-86 provides an alternative tutorial on time-stepping approaches.
  <br>
  (Daniel Arndt, 2024/07/28)
 </li>

 <li>
  Deprecated: The classes SolutionTransfer and parallel::distributed::SolutionTransfer
  have been unified to SolutionTransfer. The class now supports both serial and parallel
  meshes. The class has lost some functions: the function SolutionTransfer::interpolate()
  that takes the input vector as well as the less frequently used functions
  SolutionTransfer::prepare_for_pure_refinement() and
  SolutionTransfer::refine_interpolate(). Please use the other functions to accomplish the
  same functionality. For the time being, the old implementation has been moved to the
  `Legacy` namespace. The old parallel::distributed::SolutionTransfer has been now early deprecated.
  <br>
  (Pasquale Claudio Africa, Bruno Blais, Peter Munch, 2024/01/03)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="960-970-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: Added a new class MappingP1, which is now the default linear simplex mapping
  (i.e., it is the mapping returned by ReferenceCell::get_default_linear_mapping()).
  This mapping is significantly more efficient than MappingFE for linear simplex
  elements: e.g., Mapping::fill_fe_values() is about 6x faster in a simple load
  vector benchmark.
  <br>
  (David Wells, 2025/06/26)
 </li>

 <li>
  New: Re-introduced an interface for the MUMPS (Multifrontal Massively Parallel
  Sparse Direct Solver) library. MUMPS support had been removed in 2015, since we
  opted for its use through Trilinos and PETSc. As part of the dealii-X project,
  we plan to integrate the new features of MUMPS (i.e., CUDA support, currently
  under development) within deal.II, enabling its use with native deal.II objects
  and classes.
  <br> (Luca Heltai, Martin Kronbichler, Alfredo Buttari, 2025/04/02)
 </li>

 <li>
  Added support for the PSBLAS library: Parallel Sparse BLAS. A library of
  distributed sparse linear algebra with support for GPU and Multithread
  acceleration. Part of the PSCToolkit: Parallel Sparse Computation Toolkit.
  <br>
  (Luca Heltai, 2025/03/12)
 </li>

 <li>
  New: The step-97 tutorial program illustrates application of the
  FE_Nedelec and FE_RaviartThomas finite elements to problems in
  electromagnetics.
  <br>
  (Siarhei Uzunbajakau, Wolfgang Bangerth, 2025/01/28)
 </li>

 <li>
  Updated: deal.II now bundles TaskFlow 3.10 to support parallel computations.
  <br>
  (Wolfgang Bangerth, 2025/01/18)
 </li>

 <li>
  New: GridIn::read_vtk() was extended to read vtk meshes with field data defined at the cells, to enable importing and use of cell related data beyond a material ID or Manifold ID. The imported field data can be accessed using a new function get_cell_data(), which returns a map with a key containing the field data identifier and a vector containing the field data values associated with the cell indices of the imported mesh.
  <br>
  (Vaishnavi Kale, Marc Secanell, Mohamad Ghadban and Mayank Sabharwal, 2024/12/11)
 </li>

 <li>
  New: The step-93 tutorial program shows how to use nonlocal dofs in the
  deal.II framework, in the context of a simple optimization problem.
  <br>
  (Sam Scheuerman, Wolfgang Bangerth, 2024/11/06)
 </li>

 <li>
  Changed: The `SmartPointer` class has been renamed to ObserverPointer,
  since this name is a much better reflection of the purpose of this
  class. The documentation of the class has also been updated
  significantly to better explain what the class does.
  Correspondingly, the `Subscriptor` class has been renamed to
  EnableObserverPointer, since that is what the class actually does.
  <br>
  (Wolfgang Bangerth, 2024/09/09)
 </li>

 <li>
  Fixed: The AffineConstraints::make_consistent_in_parallel() function
  made sure that constraints for a degree of freedom stored on different
  processes were the same, but they were occasionally wrong. This has
  now been fixed.
  <br>
  (Wolfgang Bangerth, 2024/08/30)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="960-970-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Changed: The Intel Threading Building Blocks (TBB) library is no
  longer bundled as part of deal.II. Instead, the taskflow library is
  bundled and used for multithreaded tasks by default. TBB functionality
  is still being used if taskflow is explicitly disabled, but this will
  likely be removed in the future.
  <br>
  (Timo Heister, 2025/07/08)
 </li>

 <li>
  New: GridTools::extract_ordered_boundary_vertices.
  The function returns a vector with all closed
  boundaries of a triangulation. Each closed boundary
  is represented as a vector with vertex indices and
  coordinates (as a pair).
  The outer boundary is ordered counter clockwise,
  inner boundaries are ordered clockwise.
  <br>
  (Sascha Hofstetter, 2025/07/06)
 </li>

 <li>
  New: The ordering strategy DoFRenumbering::lexicographic() has been added.
  <br>
  (Michał Wichrowski, 2025/07/01)
 </li>

 <li>
  New: Implementation of CGAL functionalities for 2D.
  They allow to convert a dealii cell or a
  dealii triangulation to CGAL::Polygon_with_holes_2.
  Moreover, boolean operations can be performed
  on two CGAL::Polygon_with_holes_2.
  <br>
  (Sascha Hofstetter, 2025/07/01)
 </li>

 <li>
  Fixed: TriaAccessor::measure() now returns 1 for all `structdim == 0` objects.
  <br>
  (David Wells, 2025/06/25)
 </li>

 <li>
  Fixed: The function GridGenerator::half_hyper_shell would colorize
  the boundary ids of the triangulation incorrectly if the created shell
  was too thin. In particular the top and bottom boundary ids would
  be set incorrectly. This is fixed now.
  <br>
  (Rene Gassmoeller, 2025/06/19)
 </li>

 <li>
  Changed: ParameterHandler::print_parameters now pretty prints XML output.
  <br>
  (Marc Fehling, 2025/06/19)
 </li>

 <li>
  Fixed: The function GridGenerator::quarter_hyper_shell would colorize
  the boundary ids of the triangulation incorrectly if the created shell
  was too thin. In particular the top and bottom boundary ids would
  not be set. This did not matter for the bottom boundary since it is 0,
  which is also the default value, but the top boundary would not be
  marked correctly for these geometries. This is fixed now.
  <br>
  (Rene Gassmoeller, 2025/06/18)
 </li>

 <li>
  New: There is now a new class VectorFunctionFromTensorFunctionObject
  that can be used to convert between function objects returning a
  Tensor<1,dim>, and a vector-valued Function object.
  <br>
  (Wolfgang Bangerth, 2025/06/18)
 </li>

 <li>
  Augmented: deal.II provides functions deviator(), deviator_tensor(),
  and Physics::Elasticity::StandardTensors::dev_P() that all relate to
  the computation of the "deviator" of a tensor. These functions use a
  factor of $\frac{1}{\text{dim}}$ in their definition. This factor is
  unquestionably correct for `dim==3`, but for `dim==2` it depends on
  whether the model represents a truly two-dimensional situation, or is
  thought of as a cross-section through a three-dimensional body. This
  is, in other words, a modeling assumption. The documentation of these
  functions now explicitly describes these sorts of considerations.
  <br>
  (Wolfgang Bangerth, 2025/06/18)
 </li>

 <li>
  Added: A new mesh generator (uniform_channel_with_cylinder) for a channel
  in which there is a cylindrical obstacle. All dimensions of the
  channel can be customized by the user at the moment of the grid
  generation.
  <br>
  (Bruno Blais, 2025/06/10)
 </li>

 <li>
  New: Enable high-order VTK output in 1D.
  <br>
  (Peter Munch, 2025/04/28)
 </li>

 <li>
  Improved: The SparseDirectMUMPS class now takes an AdditionalData to
  control the MUMPS execution.
  <br>
  (Davide Polverino, 2025/05/20)
 </li>

 <li>
  Fixed: CGALWrappers dealii_cell_to_cgal_surface_mesh
  and dealii_tria_to_cgal_surface_mesh now consistently
  return closed and oriented surface meshes for tets
  and hex. This avoids errors in boolean operations.
  <br>
  (Sascha Hofstetter, 2025/05/20)
 </li>

 <li>
  Improved: GridGenerator::convert_hypercube_to_simplex_mesh() can now split
  quadrilaterals into two triangles and hexahedra into six tetrahedra. These
  anisotropic splits are useful in contexts where it is important to add as few
  new elements as possible.
  <br>
  (Kyle Schwiebert, David Wells, 2025/05/11)
 </li>

 <li>
  Fix: Ensure that undefined parameters within nested JSON subnodes
  trigger proper exception handling when using
  ParameterHandler::parse_input_from_json() with skip_undefined set
  to false. Previously, the exception ExcEntryUndeclared was only
  thrown when undefined parameters appeared at the top level of the
  parameter file. This fix extends the check to nested subnode
  entries as well.
  <br>
  (Magdalena Schreter-Fleischhacker, 2025/05/08)
 </li>

 <li>
  Fixed: FETools::extrapolate previously failed to correctly compress ghost values on
  MPI simulations with distributed triangulation and LinearAlgebra::distributed
  vectors, due to the different behaviours of the compress(insert) operation on
  Petsc/Trilinos and LinearAlgebra::distributed vectors. This is now fixed.
  <br>
  (Guilhem Poy, 2025/05/08)
 </li>

 <li>
  Changed: When the SUNDIALS::KINSOL solver fails to converge, it
  returns an error. In the past, deal.II then aborted the program (in
  debug mode), but there are cases where KINSOL's failure can make
  legitimate sense, and where a user program could catch the error and
  re-start, perhaps from a better chosen starting point. As a
  consequence, the behavior has been changed: Instead of aborting the
  program, SUNDIALS::KINSOL now throws an exception (in both debug and
  release mode) that can be caught and processed by user code.
  <br>
  (Simon Wiesheier, Wolfgang Bangerth, 2025/05/06)
 </li>

 <li>
  New: A new function FETools::cell_to_face_lexicographic() to generate a
  mapping from cell-local DoFs to lexicographic ordering of DoFs on two adjacent cells
  has been added.
  <br>
  (Michał Wichrowski, 2025/05/03)
 </li>

 <li>
  Added: 1D matrices: mass, laplace, ghost penalty.
  FullMatrix::kronecker_product.
  <br>
  (Michał Wichrowski, 2025/05/01)
 </li>

 <li>
  Improved: The MappingQEulerian class now supports
  the use of higher order mappings for the underlying
  undeformed geometry. This for example allows
  for higher accuracy in models with curved boundaries.
  <br>
  (Rene Gassmoeller, 2025/04/29)
 </li>

 <li>
  New: Fifth- and sixth order Runge--Kutta schemes have been
  added.
  <br>
  (Peter Munch, 2025/04/28)
 </li>

 <li>
  Fix: Add explicit instantiation for TpetraWrappers::BlockSparseMatrix::reinit
  (std::vector<IndexSet>, BlockDynamicSparsityPattern, MPI_Comm, bool) to fix
  missing symbol error when initializing Tpetra block matrices from
  distributed sparsity patterns.
  <br>
  (Qingyuan Shi, 2025/04/27)
 </li>

 <li>
  New: MatrixOut now allows to create "sparse" representations of
  matrices by setting the MatrixOut::Options::create_sparse_plot
  variable.
  <br>
  (Wolfgang Bangerth, 2025/04/18)
 </li>

 <li>
  Changed: Flexible GMRES can be viewed as a special case
  of the MPGMRES algorithm. Consequently, the SolverFGMRES class
  is now derived from SolverMPGMRES and its solve method now
  supports the function signature
  <code>
  solve(A, x, b, preconditioner_1, preconditioner_2, ...)
  </code>
  <br>
  (Wyatt Smith, Matthias Maier, 2025/04/10)
 </li>

 <li>
  New: Added a new class SolverMPGMRES that implements a variant of
  the GMRES algorithm allowing for the use of multiple preconditioners
  at once. For more details, see @cite Greif2017.
  <br>
  (Wyatt Smith, Matthias Maier, 2025/04/10)
 </li>

 <li>
  Improved: Portable::MatrixFree now uses a more beneficial indexing into
  multi-dimensional arrays for the cell degrees of freedom and quadrature
  data, leading to more coalesced access on GPUs.
  <br>
  (Martin Kronbichler, Urvij Saroliya, 2025/04/04)
 </li>

 <li>
  New: FE_SimplexPoly now implements FiniteElement::face_to_cell_index(),
  enabling periodicity for derived classes, with restrictions similar to the
  implementation in FE_Q_Base.
  <br>
  (Kyle Schwiebert, 2025/04/04)
 </li>

 <li>
  Fixed: DoFTools::map_dofs_to_support_points() did not work when there
  were cells in the mesh that had no degrees of freedom (e.g., because
  the cell's finite element was FE_Nothing). This is now fixed.
  <br>
  (Davit Gyulamiryan, Wolfgang Bangerth, 2025/03/23)
 </li>

 <li>
  Fixed: GridGenerator::subdivided_hyper_L() may now be used via
  GridGenerator::generate_from_name_and_arguments().
  (Bruna Campos, 2025/03/20)
 </li>

 <li>
  New: Portable::FEEvaluation::get_current_cell_index() and
  Portable::FEEvaluation::get_matrix_free_data() allow querying
  the current Portable::MatrixFree::Data object and the current cell index respectively.
  <br>
  (Daniel Arndt, Timo Heister, 2025/03/19)
 </li>

 <li>
  New: The bundled version of Kokkos has been upgraded to 4.5.1.
  <br>
  (Wolfgang Bangerth, 2025/03/10)
 </li>

 <li>
  Fixed: The particle generator function 'probabilistic_locations' could
  under certain conditions cause a division by zero if the user-provided
  probability density function evaluated to 0 everywhere in a local domain,
  but not globally. This has been fixed.
  <br>
  (Rene Gassmoeller, 2025/02/28)
 </li>

 <li>
  Improved: A default expression can now be use when
  declaring a ParsedFunction. This change allows for
  a default value different than zero which wasn't
  previously possible.
  <br>
  (Olivier Gaboriault, 2025/02/28)
 </li>

 <li>
  New: Matrix-Free now can build data structures for ghosted cells. The
  functionality can be accessed  by AdditionalData::store_ghost_cells.
  <br>
  (Michał Wichrowski, 2025/01/30)
 </li>

 <li>
  Improved: Local operation classes of the portable matrix-free
  infrastructure do not need n_dofs_1d anymore.
  <br>
  (Peter Munch, 2024/01/28)
 </li>

 <li>
  Fixed: When checking whether a particle is in a specific cell,
  ParticleHandler always assumed that the cell is a quadrilateral or
  hexahedron. This is of course correct if the mesh only contains such
  cells, but is wrong if you are working on simplex or mixed meshes and
  sometimes led to putting particles into the wrong cells. This is now
  fixed.
  <br>
  (Wolfgang Bangerth, Bruno Blais, 2025/01/27)
 </li>

 <li>
  New: There function GridGenerator::hyper_ball() now also works for
  `dim=2`, `spacedim=3`, where it creates a two-dimensional disk
  embedded in three-dimensional space.
  <br>
  (Wolfgang Bangerth, 2025/01/23)
 </li>

 <li>
  Deprecated: The functions IndexSet::pop_front() and
  IndexSet::pop_back() are now deprecated. An IndexSet should be seen as
  a set, not an ordered collection, and so speaking of "front" and
  "back" is not well defined.
  <br>
  (Wolfgang Bangerth, 2025/01/23)
 </li>

 <li>
  Deprecated: The functions ReferenceCell::unit_tangential_vectors() and
  ReferenceCell::unit_normal_vectors() were poorly renamed, given the
  plural (when they only returned a single vector) and the use of the
  word "unit" when they return information about faces. As a
  consequence, these functions have been renamed to
  ReferenceCell::face_tangent_vector() and
  ReferenceCell::face_normal_vector(), and the old names are now
  deprecated.
  <br>
  (Wolfgang Bangerth, 2025/01/23)
 </li>

 <li>
  New: There are now functions
  FiniteElement::shape_function_belongs_to() that allow testing whether
  a given shape function corresponds to a vector component described by
  a FEValuesExtractors object.
  <br>
  (Wolfgang Bangerth, 2024/12/18)
 </li>

 <li>
  Deprecated: The function parallel::transform() has been deprecated.
  <br>
  (Wolfgang Bangerth, 2024/12/18)
 </li>

 <li>
  Fixed: For `dim = 1`, the setup of MatrixFree did not correctly
  identify faces between cells of different refinement level. This is
  now fixed.
  <br>
  (Sean Johnson, 2024/12/12)
 </li>

 <li>
  New: Configuration and CI job for pre-commit.
  The tool pre-commit can install hooks into the local git repo.
  These will be run on every call to git commit and perform some quick checks
  on the changeset.
  The CI job runs these checks on the whole codebase.
  <br>
  (Jan Philipp Thiele, 2024/11/27)
 </li>

 <li>
  Improved: The function FEInterfaceValues::quadrature_point() has
  been added.
  <br>
  (Peter Munch, 2024/11/23)
 </li>

 <li>
  New: FiniteElement::get_local_dof_sparsity_pattern() can be used to
  provide the coupling between DoFs within a cell and is used in
  make_sparsity_pattern() to generate matrices with fewer nonzero
  entries depending on the element. FE_Q_iso_Q1 now provides this
  coupling information.
  <br>
  (Timo Heister, Luca Heltai, 2024/11/20)
 </li>

 <li>
  Deprecated: The Utilities::MPI::ConsensusAlgorithms::Payload base
  class has been deprecated. It used to serve as an interface class for
  the algorithms in namespace Utilities::MPI::ConsensusAlgorithms, but
  we have since come up with ways of formulating these algorithms in
  terms of function objects that are more flexible than the base
  class/virtual function interface now deprecated.
  <br>
  (Wolfgang Bangerth, 2024/11/12)
 </li>

 <li>
  New: The cmake configuration system now allows enabling
  interprocedural and link-time optimization by the compiler. These
  kinds of optimization often result in substantially faster executables
  because they allow the compiler to see a bigger part of the whole
  program when optimizing.
  To enable interprocedural and link-time optimization, pass the
  `-DDEAL_II_USE_LTO=ON` flag to cmake when configuring the library.
  <br>
  (Wolfgang Bangerth, 2024/09/24)
 </li>

 <li>
  Improved: The orthogonalization done within SolverGMRES and SolverFGMRES for
  the deal.II vectors would previously lead to data access pattern that are
  unfriendly to data prefetchers on modern CPUs. This has been addressed by
  implementing a suitable loop blocking.
  <br>
  (Martin Kronbichler, 2024/11/01)
 </li>

 <li>
  Improved: DoFRenumbering::matrix_free_data_locality() now also works for
  elements that are not FE_Q, simply using the ordering of the degrees of
  freedom on the cells.
  <br>
  (Martin Kronbichler, 2024/11/01)
 </li>

 <li>
  Improved: The new function
  Utilities::MPI::compute_index_owner_and_requesters()
  allows to compute index owners but also returns the requesters
  for locally owned indices.
  <br>
  (Peter Munch, 2024/10/15)
 </li>

 <li>
  Improved: The new function LAPACKFullMatrix::get_state() allows to
  query the current state of LAPACKFullMatrix, which allows to decide
  which method (e.g., LAPACKFullMatrix::solve() vs.
  LAPACKFullMatrix::vmult()) to use.
  <br>
  (Peter Munch, 2024/10/15)
 </li>

 <li>
  Improved: FunctionFromFunctionObjects can now accept
  a single function within which individual components need
  to be handled.
  <br>
  (Peter Munch, 2024/10/15)
 </li>

 <li>
  Fixed: FESystem now works correctly with cubic FE_SimplexP elements.
  <br>
  (David Wells, 2024/10/15)
 </li>

 <li>
  Improved: The deal_ii_invoke_autopilot CMake macro gained support for
  multi-configuration generators.
  <br>
  (Matthias Maier, 2024/10/12)
 </li>

 <li>
  Improved: TableHandler now correctly outputs Integers if
  TableHandler::set_scientific() is set.
  <br>
  (Peter Munch, 2024/10/05)
 </li>

 <li>
  Improved: The class ScalarFunctionFromFunctionObject can now
  also handle time-dependent functions.
  <br>
  (Peter Munch, 2024/09/26)
 </li>

 <li>
  Improved: Added AlignedVector::insert(), which works the same way as
  std::vector::insert().
  <br>
  (David Wells, 2024/09/26)
 </li>

 <li>
  New: Utilities::MPI::NoncontiguousPartitioner::export_to_ghosted_array()
  can now handle multiple components.
  <br>
  (Peter Munch, 2024/09/25)
 </li>

 <li>
  Improved: MGTransferMF now also supports FE_DGP.
  <br>
  (Peter Munch, Nils Margenberg, 2024/09/24)
 </li>

 <li>
  New: The constructors of the classes within namespace
  FEValuesExtractors are now all marked as `constexpr`. As a
  consequence, FEValuesExtractors objects can now also be `constexpr`.
  <br>
  (Wolfgang Bangerth, 2024/09/24)
 </li>

 <li>
  New: The new functions get_normal_hessian() and submit_normal_hessian() for
  FEFaceEvaluation implement second derivatives in normal direction for the
  matrix-free framework.
  <br>
  (Maximilian Bergbauer, Andreas Koch, 2024/09/18)
 </li>

 <li>
  Patched: GridTools::transform() now supports 3D
  anisotropically refined grids. Note that FEValues
  still do not support hanging node constraints for
  anisotropically refined grids.
  <br>
  (Robin Hiniborch, 2024/09/10)
 </li>

 <li>
  Fixed: Triangulation::get_boundary_ids() does not return duplicates when
  `dim == 1` any more.
  <br>
  (David Wells, 2024/09/06)
 </li>

 <li>
  Fix: Add the missing `compressed = false` in TpetraWrappers::Vector::add()
  function to ensure distributed Tpetra vector with writable nonlocal
  entries will always be non-compressed after any entry addition.
  <br>
  (Qingyuan Shi, 2024/08/28)
 </li>

 <li>
  New: The new setting OutputStyle::KeepOnlyChanged allows
  to print only changed parameters with
  ParameterHandler::print_parameters().
  <br>
  (Peter Munch, 2024/08/22)
 </li>

 <li>
  New: deal.II now contains a python script to indent ParameterHandler .prm
  files. The script is located in contrib/utilities/prm_tools.py and can be
  used to update .prm file in-place or write the output into a new file.
  <br>
  (Rene Gassmoeller, 2024/08/21)
 </li>

 <li>
  New: The TaskResult::emplace_object() function sets the value of a
  TaskResult object to a specific value without requiring running a
  function on a separate task.
  <br>
  (Wolfgang Bangerth, 2024/08/21)
 </li>

 <li>
  Improved: DataOut now also supports cubic output of simplices.
  <br>
  (Peter Munch, 2024/08/14)
 </li>

 <li>
  New: FESystem::compare_for_domination() now accepts two FESystem objects
  with different number of base elements.
  <br>
  (Mohamad Ghadban,  2024/08/20)
 </li>

 <li>
  Fixed: Objects of type FEFaceEvaluation can now be collected in an std::vector.
  <br>
  (Martin Kronbichler, John Coughlin, 2024/08/19)
 </li>

 <li>
  Changed: The concepts::is_vector_space_vector now also requires the vector
  class to provide a function `VectorType::get_mpi_communicator() -> MPI_Comm`,
  which returns the underlying MPI communicator or, if a sequential vector,
  @p MPI_COMM_SELF.
  <br>
  (Martin Kronbichler, David Wells, 2024/08/19)
 </li>

 <li>
  Improved: MappingFE now supports the case that not
  all geometric objects of a cell have the same
  manifold id.
  <br>
  (Peter Munch, 2024/08/16)
 </li>

 <li>
  Improved: The classical Gram--Schmidt orthonormalization in SolverGMRES now
  uses SIMD-optimized routines also for dealii::Vector and dealii::BlockVector.
  <br>
  (Martin Kronbichler, 2024/08/16)
 </li>

 <li>
  New: Added SparseVanka::Tvmult() and SparseVanka::clear()
  SparseVanka can now be passed to MGSmootherPrecondition to be used as a multigrid smoother.
  <br>
  (Chayapol Chaoveeraprasit, 2024/08/14)
 </li>

 <li>
  Improved: Portable::MatrixFree now supports FESystem(FE_Q, n_components)
  with arbitrary number of components.
  <br>
  (Peter Munch, Daniel Arndt, 2024/08/14)
 </li>

 <li>
  New: Added support for magic_enum.hpp. Now PatternsTools::Convert<T> works also when T is an enum type.
  <br>
  (Luca Heltai, 2024/08/12)
 </li>

 <li>
  New: deal.II contributions are now automatically checked
  for typos using the 'typos' software.
  <br>
  (Rene Gassmoeller, 2024/08/12)
 </li>

 <li>
  New: Parameter entries in the ParameterHandler class can now be marked as
  deprecated using the function 'mark_as_deprecated'. If deprecated parameters
  are found in a parameter file an exception of type
  'ExcEncounteredDeprecatedEntries' is thrown that can be caught in
  the calling code, or if ignored will emit an error message.
  <br>
  (Rene Gassmoeller, 2024/08/11)
 </li>

 <li>
  New: The new functions DoFTools::extract_rigid_body_modes
  and DoFTools::extract_level_rigid_body_modes
  allow you to extract translational and rotational modes,
  need to set up AMG for solving elasticity problems.
  <br>
  (Marco Feder, Peter Munch, Richard Schussnig, 2024/08/10)
 </li>

 <li>
  New: Apply the 'dealii' namespace import change to step-49 and step-50,
  ensuring consistency with the other tutorial programs.
  <br>
  (L&oacute;r&aacute;nt Hadnagy, 2024/08/08)
 </li>

 <li>
  Fixed: MatrixFreeTools::compute_matrix() now
  correctly handles the case that individual
  components are requested.
  <br>
  (Peter Munch, 2024/08/04)
 </li>

 <li>
  Fixed: MatrixFreeTools::compute_diagonal() used
  to use first_selected_component also for accessing
  (block) vectors. This has been fixed.
  <br>
  (Peter Munch, Nils Much, 2024/07/30)
 </li>

 <li>
  Fix: Corrected pointer iteration range in VectorBase::mean_value()
  and VectorBase::lp_norm() in PETScWrappers namespace to proper handle
  distributed PETSc vectors. Now pointers obtained via VecGetArrayRead()
  iterate only over locally_owned_size(), and results are correctly
  aggregated across processes.
  <br>
  (Qingyuan Shi, 2024/08/01)
 </li>

 <li>
  New: Fourth order tensors can now be symmetrized based
  on the required symmetry type (major or minor).
  <br>
  (Vinayak, 2024/07/31)
 </li>

 <li>
  New: TrilinosWrappers::SolverDirect can now also
  be used for a matrix of type Epetra_Operator and
  vectors of type Epetra_MultiVector. In addition,
  FullMatrix objects can be used as multi vectors.
  <br>
  (Peter Munch, 2024/07/31)
 </li>

 <li>
  New: Uniformly apply the practice of importing the 'dealii' namespace
  only within the StepXX namespace across all tutorial programs,
  except for step-49 and step-50. These steps will be handled separately
  to ensure proper integration and testing.
  <br>
  (L&oacute;r&aacute;nt Hadnagy, 2024/07/30)
 </li>

 <li>
  New: VectorTools::interpolate() can now also used on multigrid levels.
  <br>
  (Peter Munch, Richard Schussnig, 2024/07/26)
 </li>

 <li>
  Deprecated: ObserverPointer is clearly documented as not owning the
  object it points to. Yet, ObserverPointer::clear() deletes the object
  pointed to, thereby assuming that the ObserverPointer object actually
  owns the object. This must surely be confusing. As a consequence the
  function is now deprecated.
  <br>
  (Wolfgang Bangerth, 2024/07/25)
 </li>

 <li>
  New: Variables that represent a mapping (and are generally called `mapping`)
  are now all `const` in tutorial programs.
  <br>
  (L&oacute;r&aacute;nt Hadnagy, 2024/07/25)
 </li>

</ol>

*/
