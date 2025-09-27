// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
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
@page changes_between_9_0_1_and_9_1_0 Changes between Version 9.0.1 and 9.1.0

<p>
This is the list of changes made between the release of deal.II version
9.0.1 and that of 9.1.0. All entries are signed with the names of the
author.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="901-910-incompatible"></a>
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
  Changed: Tasks of member function
  hp::FECollection::find_least_face_dominating_fe() are now divided
  into two functions hp::FECollection::find_common_fes()
  and hp::FECollection::find_dominated_fe().
  <br>
  (Marc Fehling, 2019/04/08)
 </li>

 <li>
  Removed: The IndexSet instantiations for FEValues::get_function_values()
  and DataOut_DoFData::add_data_vector() have been removed.
  <br>
  (Daniel Arndt, 2019/03/20)
 </li>

 <li>
  Changed: The "debug" and "release" targets that are created by the
  DEAL_II_INVOKE_AUTOPILOT() macro no longer automatically rebuild the
  project.
  <br>
  (Matthias Maier, 2019/02/05)
 </li>

 <li>
  Changed: VectorSlice has been deprecated in favor of the more general ArrayView
  class. All public interfaces now take ArrayView arguments instead of VectorSlice
  arguments. Since VectorSlice now inherits from ArrayView this should be
  compatible with the overwhelming majority of use cases.
  <br>
  (David Wells, 2019/01/21)
 </li>

 <li>
  Changed: The class VectorView has been removed. The suggested replacements are
  to either use an ArrayView, a BlockVector, or to copy the relevant subset into a
  Vector.
  <br>
  (David Wells, 2019/01/11)
 </li>

 <li>
  Changed: The particle MPI datatype macro <code>PARTICLE_INDEX_MPI_TYPE</code> is
  now properly namespaced: its new name is
  <code>DEAL_II_PARTICLE_INDEX_MPI_TYPE</code>.
  <br>
  (David Wells, 2018/11/30)
 </li>

 <li>
  Changed: Subscriptor::subscribe() requires providing a pointer to a boolean
  that can be used to signal validity of the object pointed to by the subscribing
  object.
  <br>
  (Daniel Arndt, 2018/11/02)
 </li>

 <li>
  Deprecated: deal.II's NetCDF bindings have been deprecated. The present set of
  bindings are incompatible with recent releases of NetCDF and, at the present
  time, no one has volunteered to either maintain or upgrade our interface. The
  bindings will be removed in a future version if no one steps forward to update
  the bindings.
  <br>
  (David Wells, 2018/10/27)
 </li>

 <li>
  Changed: Time is described by a real-valued scalar in Function, TensorFunction,
  ConstantTensorFunction and ZeroTensorFunction, i.e. all of them derive from
  FunctionTime with a real-valued scalar template parameter.
  <br>
  (Daniel Arndt, 2018/10/26)
 </li>

 <li>
  Deprecated: FiniteElement::compare_for_face_domination() has been
  deprecated and will be replaced by the more versatile function
  FiniteElement::compare_for_domination().
  <br>
  (Marc Fehling, 2018/10/11)
 </li>

 <li>
  Changed: Triangulation::create_triangulation() doesn't allow to set a boundary
  id different from numbers::internal_face_boundary_id for internal faces anymore.
  <br>
  (Daniel Arndt, 2018/08/07)
 </li>

 <li>
  Changes: The data transfer interface of the class
  parallel::distributed::Triangulation now requires different kinds of
  callback functions. parallel::distributed::Triangulation::register_data_attach()
  now takes `std::function<std::vector<char>(cell_iterator &, CellStatus)>`,
  which returns the buffer of the packed data.
  parallel::distributed::Triangulation::notify_ready_to_unpack()
  requires `std::function<void(const cell_iterator &, const CellStatus,
  const boost::iterator_range<std::vector<char>::%const_iterator &>`, where
  the last argument describes an iterator range, from which the callback
  function is allowed to read.
  <br>
  Further, parallel::distributed::Triangulation::register_data_attach()
  now requires a boolean argument `returns_variable_size_data`, which
  denotes if the registered pack_callback function interacts with the fixed
  size (`=false`) or variable size (`=true`) buffer for data transfer.
  <br>
  (Marc Fehling, 2018/07/20)
 </li>

 <li>
  Changed: MatrixCreator::create_mass_matrix and
  MatrixCreator::create_boundary_mass_matrix that had mixed number type (real
  valued for the matrix and complex valued for vectors) have been changed to
  support complex number types uniformly. This implies that now all
  underlying number types of matrix, vectors and AffineConstraints objects
  have to match; mixed variants are no longer supported.
  <br>
  (Matthias Maier, 2018/05/25)
 </li>

 <li>
  Changed: The AffineConstraints object (formerly ConstraintMatrix) gained a
  template parameter for the underlying storage type. In order to facilitate
  this change, the interface for AffineConstraints has become more strict: In
  particular, the AffineConstraints::distribute_local_to_global(),
  AffineConstraints::add_entries_local_to_global(), and other functions now
  require that all matrix and vector arguments have the same matching number
  type as the AffineConstraints object.
  <br>
  (Matthias Maier, 2018/05/25)
 </li>

 <li>
  Changed: The OpenCASCADE Manifold classes with names ending in Boundary
  (i.e., OpenCASCADE::NormalProjectionBoundary,
  OpenCASCADE::DirectionalProjectionBoundary, and
  OpenCASCADE::NormalToMeshProjectionBoundary) have been deprecated in favor
  of renamed classes ending in Manifold (i.e.,
  OpenCASCADE::NormalProjectionManifold,
  OpenCASCADE::DirectionalProjectionManifold, and
  OpenCASCADE::NormalToMeshProjectionManifold).
  <br>
  (David Wells, 2018/05/16)
 </li>

 <li>
  Changed: The PolynomialSpace::compute_index() and
  PolynomialsP::directional_degrees() functions used to return their
  information through an array that they received as a reference
  argument. Instead, they now return a `std::array<unsigned int,dim>`
  by value.
  <br>
  (Wolfgang Bangerth, 2018/05/09)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="901-910-general"></a>
<h3>General</h3>

<ol>

 <li>
  New: The Step-63 tutorial has been added which discusses block
  smoothers vs. point smoothers inside a geometric multigrid
  v-cycle. This tutorial also implements a GMG solver for a
  non-symmetric PDE, and demonstrates the effect of DoF/cell renumbering
  for multiplicative smoothers for advection-dominated problems.
  <br>
  (Thomas C. Clevenger, Timo Heister, 2019/05/10)
 </li>

 <li>
  Improved: All tutorial programs have been reviewed and updated with new features
  in the library and C++11.
  <br>
  (Daniel Arndt, David Wells, Matthias Maier, Denis Davydov, Wolfgang Bangerth, Jean-Paul Pelteret, Martin Kronbichler, 2019/05/13)
 </li>

 <li>
  Deprecated: deal.II's nanoflann bindings (i.e., the KDTree class) have been
  deprecated in favor of using <code>boost::geometry::index::rtree</code>.
  <br>
  (David Wells, 2019/05/10)
 </li>

 <li>
  New: The step-64 tutorial demonstrates how to use the CUDAWrappers::MatrixFree
  framework (possibly with MPI) and discusses the peculiarities of using CUDA
  inside deal.II in general.
  <br>
  (Daniel Arndt, Bruno Turcksin, 2019/05/10)
 </li>

 <li>
  New: The new tutorial program step-61 shows how to use "weak Galerkin" finite element method to solve the Poisson equation.
  <br>
  (Zhuoran Wang, 2019/05/02)
 </li>

 <li>
  New: The new tutorial program step-62 shows how to calculate the resonance
  frequency and the bandgap of a phononic crystal. The program solves the elastic
  wave equation with Perfectly Matches Layers in the frequency domain using
  complex algebra and parallelization via MUMPS and MPI.
  <br>
  (Daniel Garcia-Sanchez, 2019/05/01)
 </li>

 <li>
  New: The tutorial examples now all use a quadrature formula with a degree depending on the degree of the
  finite element, so that one changing the finite element degree does not also require changing by hand the number
  of quadrature points.
  <br>
  (Roland Richter, 2019/04/25)
 </li>

 <li>
  New: The ParsedConvergenceTable class allows convenient construction of convergence tables, exploiting
  parameter files.
  <br>
  (Luca Heltai, 2019/04/10)
 </li>

 <li>
  New: A new class Differentiation::SD::Expression has been added which allows
  the creation and manipulation of scalar symbolic functions using the SymEngine
  library. It will later be used as the basis of symbolic tensor calculus, amongst
  other things. It incorporates numerous features, which currently include:
  <ul>
    <li>expression parsing,</li>
    <li>comparison operations,</li>
    <li>logical operations,</li>
    <li>basic math operations,</li>
    <li>conditional expression construction,</li>
    <li>differentiation,</li>
    <li>substitution (partial and complete), and</li>
    <li>serialization.</li>
  </ul>
  <br>
  (Jean-Paul Pelteret, 2019/03/29)
 </li>

 <li>
  New: A new class Differentiation::AD::VectorFunction has been
  added to help implement point-wise vector functions using automatic
  differentiation. In particular, this class is designed to compute the
  Jacobian of a vector function that is parameterized in
  terms of scalar, vector, and tensor arguments. One example of its use
  would be to compute the linearization of a multi-field constitutive
  law that is expressed in terms of some kinetic variables.
  <br>
  (Jean-Paul Pelteret, 2019/03/21)
 </li>

 <li>
  New: Triangulation::execute_coarsening_and_refinement() now also performs
  p-coarsening and p-refinement on all attached hp::DoFHandler objects.
  <br>
  (Marc Fehling, 2019/03/18)
 </li>

 <li>
  New: A new class Differentiation::AD::ScalarFunction has been
  added to help implement point-wise scalar functions using automatic
  differentiation. In particular, this class is designed to compute the
  gradient and Hessian of a scalar function that is parameterized in
  terms of scalar, vector, and tensor arguments. One example of its use
  would be to compute the derivatives of a multi-field constitutive law
  that is expressed in terms of an energy function.
  <br>
  (Jean-Paul Pelteret, 2019/02/20)
 </li>

 <li>
  Improved: The preconditioner and solver setup in step-20 has been rewritten
  to highlight the new LinearOperator and PackagedOperation classes.
  <br>
  (Matthias Maier, 2019/02/08)
 </li>

 <li>
  New: Support for Ginkgo, a high-performance numerical linear algebra library
  has been added with classes inheriting from GinkgoWrappers::SolverBase. Ginkgo
  provides advanced highly optimized linear solvers, matrix formats and an
  abstraction to easily create linear operators on both the cpu and the gpu. The
  deal.II's Ginkgo interface can currently be used to solve linear systems using
  Krylov solvers with the cuda and OpenMP paradigms.
  <br>
  (Ginkgo developers, 2019/01/30)
 </li>

 <li>
  Changed: Class hp::DoFHandler now transfers the active_fe_index of each
  cell automatically when refining/coarsening a Triangulation,
  parallel::shared::Triangulation, or
  parallel::distributed::Triangulation. However, serialization of a
  parallel::distributed::Triangulation still requires a user to
  explicitly call the functions
  hp::DoFHandler::prepare_for_serialization_of_active_fe_indices() and
  hp::DoFHandler::deserialize_active_fe_indices().
  <br>
  (Marc Fehling, 2019/01/27)
 </li>

 <li>
  New: A new HDF5 interface has been added. The Hierarchical Data Format (HDF) is
  a cross platform and a high I/O performance format designed to store large
  amounts of data. It supports serial and MPI I/O access. The deal.II's HDF5
  interface can be used to write and read data, such the results of a simulation,
  in a HDF5 file.
  <br>
  (Daniel Garcia-Sanchez, 2019/01/09)
 </li>

 <li>
  Changed: hp::DoFHandler now automatically sets active_fe_indices
  on cells that will be refined and coarsened whenever a Triangulation or
  parallel::shared::Triangulation is used. Upon refinement, the
  active_fe_index will be passed on to the children. For coarsening, the
  parent's active_fe_index will be determined from its former children
  using the FiniteElementDomination logic.
  <br>
  (Marc Fehling, 2018/10/15)
 </li>

 <li>
  New: %Function FiniteElement::compare_for_domination() inspects two
  FiniteElement objects upon FiniteElementDomination, and uses a codim
  parameter that determines in which subspace we actually compare them.
  <br>
  (Marc Fehling, 2018/10/11)
 </li>

 <li>
  New: A new class Differentiation::AD::ResidualLinearization has been
  added to help compute the linearization of a residual vector defined on the
  level of a cell (a finite element residual), or for local nonlinear equations.
  <br>
  (Jean-Paul Pelteret, 2018/09/30)
 </li>

 <li>
  Changed: hp::DoFHandler::distribute_dofs() called on
  parallel::distributed::Triangulation objects now assigns exactly
  one index for each degree of freedom located on interfaces with
  ghost cells. Previously, those degrees of freedom were assigned
  indices by each adjacent subdomain.
  <br>
  (Marc Fehling, 2018/08/31)
 </li>

 <li>
  New: A new class Differentiation::AD::EnergyFunctional has been
  added to help implement (incremental) variational formulations using automatic
  differentiation. In particular, this class is designed to compute the finite
  element residuals and their linearizations.
  <br>
  (Jean-Paul Pelteret, 2018/08/27)
 </li>

 <li>
  New: The parallel::distributed::SolutionTransfer class is now capable
  of transferring data that has been associated with a hp::DoFHandler
  object. See the class documentation of the former on how to use this
  new functionality.
  <br>
  (Marc Fehling, 2018/08/15)
 </li>

 <li>
  Improved: GridGenerator::extrude_triangulation can now optionally extrude manifold
  ids.
  <br>
  (David Wells, 2018/08/11)
 </li>

 <li>
  New: FE_BernardiRaugel is the implementation of an element developed by Christine Bernardi and
  Genevieve Raugel in the 1985 Mathematics of Computation paper. The Bernardi-Raugel element is an
  enrichment of the P1-P0 element for Stokes problems by the addition of edge-based quadratic
  bubble functions. The original description is for simplicial meshes, but this implementation is
  for quadrilateral and hexahedral meshes. The BR1-P0 pair is LBB stable for Stokes problems. Due
  to the edge-based bubble functions, features such as hanging nodes are not supported.
  <br>
  (Graham Harper, 2018/07/09)
 </li>

 <li>
  New: A new function GridGenerator::channel_with_cylinder has been added. This
  function generates a grid corresponding to the classic flow past a cylinder test
  problem, and uses manifolds to blend between the description of the cylinder and
  the Cartesian bulk geometry.
  <br>
  (David Wells, 2018/07/04)
 </li>

 <li>
  Changed: The ConstraintMatrix class has been renamed to AffineConstraints.
  Further, AffineConstraints gained a template parameter for the underlying
  storage type. This in particular enables the definition of complex-valued
  constraints (for example, complex-valued Dirichlet boundary conditions, or
  periodic boundary conditions with an additional phase shift).
  <br>
  (Matthias Maier, Daniel Arndt, David Wells, 2018/05/31)
 </li>

 <li>
  New: The introduction to step-6 has been completely rewritten. It
  previously mentioned only the bare minimum about what this important
  tutorial program actually does. It now gives a much broader overview
  of why we use adaptive meshes and how this is achieved in practice.
  <br>
  (Wolfgang Bangerth, 2018/05/25)
 </li>

 <li>
  Changed: clang-format replaces astyle as source code formatting tool.
  <br>
  (Matthias Maier, Daniel Arndt, 2018/05/24)
 </li>

 <li>
  New: FE_NedelecSZ is a new H(curl)-conforming element which overcomes the sign conflict problem
  which can be encountered for Nedelec elements on hexahedral meshes. This is an implementation of
  the method described in the PhD thesis of Sabine Zaglmayr and should be considered over the
  FE_Nedelec element where meshes with non-standard orientation are used. Note that not all
  functionality (e.g. Hessians), has been implemented for this element.
  <br>
  (Ross Kynch, 2018/04/28)
 </li>

 <li>
  New: The new class ColorEnriched::Helper constructs an hp::FECollection, a
  collection of finite element objects used by hp::DoFHandler, in a domain
  with multiple, possibly overlapping, sub-domains with individual
  enrichment functions. Note that the overlapping regions may have
  multiple enrichment functions associated with them. This is implemented
  using a general constructor of FE_Enriched which allows different
  enrichment functions.
  <br>
  (Nivesh Dommaraju, Denis Davydov, 2018/04/20)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="901-910-specific"></a>
<h3>Specific improvements</h3>

<ol>

 <li>
  Improved: GridIn::read_msh has been extended to allow for gmsh mesh format 4.1.
  <br>
  (Daniel Arndt, 2019/05/21)
 </li>

 <li>
  Improved: The point identification of TransfiniteInterpolationManifold has
  been made more robust. It would previously sometimes fail for strongly curved,
  long and skinny cells.
  <br>
  (Martin Kronbichler, 2019/05/16)
 </li>

 <li>
  Improved: DoFHandler::renumber_dofs() on parallel::distributed::Triangulation
  would scale linearly with the number of MPI ranks on more than 50k ranks,
  reaching times of more than 10 seconds on 150k ranks. This is now fixed for
  the case the renumbering does not lead to a change in the underlying IndexSet
  field of locally owned unknowns, i.e., in case the numbers merely change
  within each MPI rank but not across them.
  <br>
  (Martin Kronbichler, 2019/05/10)
 </li>

 <li>
  New: Introduced indices for future finite elements that cells will be assigned
  to after the triangulation changes. The future_fe_index of each cell can be
  accessed and changed via member functions of the DoFCellAccessor class,
  namely DoFCellAccessor::future_fe_index(), DoFCellAccessor::set_future_fe_index(),
  DoFCellAccessor::future_fe_index_set(), and DoFCellAccessor::clear_future_fe_index().
  <br>
  (Marc Fehling, 2019/05/10)
 </li>

 <li>
  New: Add Utilities::MPI::mean_and_standard_deviation() to calculate mean and standard deviation.
  <br>
  (Denis Davydov, 2019/05/09)
 </li>

 <li>
  New: Add DynamicSparsityPattern::get_view() to create a subset of a sparsity pattern based on a selection of some rows.
  <br>
  (Denis Davydov, 2019/05/08)
 </li>

 <li>
  New: Some utility functions that perform symbolic substitution on, and
  evaluation of, symbolic scalar and tensor expressions have been added to the
  Differentiation::SD namespace.
  <br>
  (Jean-Paul Pelteret, 2019/05/07)
 </li>

 <li>
  Improved: Extend Utilities::inverse_Hilbert_space_filling_curve to handle degenerate cases.
  <br>
  (Denis Davydov, 2019/05/07)
 </li>

 <li>
  New: Some utility functions that resolve explicit (and cyclic) dependencies
  between entries in symbolic substitution maps.
  <br>
  (Jean-Paul Pelteret, 2019/05/05)
 </li>

 <li>
  New: Some utility functions that facilitate the creation of symbolic maps and
  symbolic substitution maps have been added to the Differentiation::SD namespace.
  <br>
  (Jean-Paul Pelteret, 2019/05/05)
 </li>

 <li>
  New: Mapping::get_bounding_box() returns the correct bounding box for mappings that do not preserve vertex
  locations.
  <br>
  (Luca Heltai, 2019/05/03)
 </li>

 <li>
  Improved: The `Meshworker::mesh_loop()` function is now capable of working with an
  IteratorRange, and also supports iterator ranges constructed from
  FilteredIterators.
  <br>
  (Jean-Paul Pelteret, 2019/05/02)
 </li>

 <li>
  Improved: The WorkStream::run() function is now capable of working with iterator
  ranges, or any general iterable object that defines the `begin()` and `end()`
  of a range of elements to iterate over.
  <br>
  (Jean-Paul Pelteret, 2019/05/02)
 </li>

 <li>
  New: Some utility functions that perform symbolic differentiation of scalar
  symbolic expressions, as well as tensorial expressions, have been added to the
  Differentiation::SD namespace.
  <br>
  (Jean-Paul Pelteret, 2019/04/23)
 </li>

 <li>
  Changed: The chunk_size used in CUDA kernels has been reduced from 8 to 1 to
  improve performance on newer architectures.
  <br>
  (Bruno Turcksin, 2019/04/23)
 </li>

 <li>
  New: Some utility functions that facilitate the creation of both scalar and
  tensor symbolic variables and symbolic functions have been added to the
  Differentiation::SD namespace.
  <br>
  (Jean-Paul Pelteret, 2019/04/22)
 </li>

 <li>
  New: Functions::CutOffFunctionBase now supports rescaling the function to keep its integral equal to one.
  <br>
  (Luca Heltai, 2019/04/18)
 </li>

 <li>
  New: The new FunctionFromFunctionObjects class allows one to wrap a vector of
  std::function objects into a Function object, to allow fast prototyping of user codes.
  <br>
  (Luca Heltai, 2019/04/12)
 </li>

 <li>
  Fixed: Make sure that GridTools::Cache::get_cell_bounding_boxes_rtree() honors mappings.
  <br>
  (Luca Heltai, 2019/04/11)
 </li>

 <li>
  New: Added BoundingBox::extend() that allows extending and shrinking of BoundingBox objects.
  <br>
  (Luca Heltai, 2019/04/11)
 </li>

 <li>
  New: The method Mapping::get_center() allows one to retrieve the cell center of a cell when the
  mapping object does not preserve vertex locations.
  <br>
  (Luca Heltai, 2019/04/11)
 </li>

 <li>
  New: Member function hp::FECollection::find_dominating_fe_extended()
  returns the index of the most dominating finite element out of a given
  set of indices. If none was found, the search will be extended on the
  complete collection.
  <br>
  (Marc Fehling, 2019/04/08)
 </li>

 <li>
  New: Member function hp::FECollection::find_dominating_fe()
  returns the index of the most dominating finite element out of a given
  set of indices.
  <br>
  (Marc Fehling, 2019/04/08)
 </li>

 <li>
  New: Member function hp::FECollection::find_enclosing_fes()
  returns a set of indices from the full FECollection whose associated
  finite elements are dominated by all elements of a given set of indices.
  <br>
  (Marc Fehling, 2019/04/08)
 </li>

 <li>
  New: Member function hp::FECollection::find_common_fes()
  returns a set of indices from the full FECollection whose associated
  finite elements dominate all elements of a given set of indices.
  <br>
  (Marc Fehling, 2019/04/08)
 </li>

 <li>
  Changed: SmartPointer and Subscriptor use a `std::string`
  instead of a `const char *` for subscriber identification. As a result,
  subscriber strings are no longer compared by their memory address but instead
  by their content.
  <br>
  (Daniel Arndt, 2019/04/08)
 </li>

 <li>
  New: The GeneralDataStorage class facilitates the storage of any general data.
  It offers the ability to store any amount of data, of any type, which is then
  made accessible by an identifier string.
  <br>
  (Luca Heltai, Jean-Paul Pelteret, 2019/04/02)
 </li>

 <li>
  New: Add ArrayView::ArrayView() and ArrayView::reinit(value_type *, const std::size_t).
  <br>
  (Denis Davydov, 2019/04/02)
 </li>

 <li>
  New: The Differentiation::SD::Expression class can now be used to perform
  symbolic math calculations.
  <br>
  (Jean-Paul Pelteret, 2019/03/30)
 </li>

 <li>
  New: The Vector class can now be initialized using an object of type
  std::initializer_list. Such objects are, in particular, created when
  the compiler sees a brace-enclosed list of numbers.
  <br>
  (Wolfgang Bangerth, 2019/03/28)
 </li>

 <li>
  Changed: parallel::distributed::SolutionTransfer::prepare_serialization() has
  been deprecated in favor of
  parallel::distributed::SolutionTransfer::prepare_for_serialization().
  <br>
  (Daniel Arndt, 2019/03/26)
 </li>

 <li>
  New: Added a variant of MeshWorker::mesh_loop() that takes a class and its member functions as workers
  and copiers.
  <br>
  (Luca Heltai, 2019/03/23)
 </li>

 <li>
  New: The type alias FEValuesViews::View now allows one to infer the correct FEValuesViews type
  that would be returned by FEValuesBase::operator[]() when called with a specified extractor type from the
  FEValuesExtractors namespace.
  <br>
  (Luca Heltai, 2019/03/20)
 </li>

 <li>
  New: Add DynamicSparsityPattern::nonempty_cols() and
  DynamicSparsityPattern::nonempty_rows() to return columns/rows stored in the
  sparsity pattern.
  <br>
  (Denis Davydov, 2019/03/20)
 </li>

 <li>
  New: The new method GridTools::assign_co_dimensional_manifold_indicators() propagate manifold ids from
  cells to faces and edges, according to a disambiguation function that takes the set of manifold ids of
  the cells that share the given face or edge.
  <br>
  (Luca Heltai, 2019/03/19)
 </li>

 <li>
  New: FEValuesExtractors classes now have a new method get_name() that returns a unique string identifier
  for each extractor.
  <br>
  (Luca Heltai, 2019/03/19)
 </li>

 <li>
  New: Add SparsityTools::gather_sparsity_pattern().
  <br>
  (Denis Davydov, 2019/03/19)
 </li>

 <li>
  New: Add Utilities::MPI::create_ascending_partitioning() to create a one-to-one
  ascending partitioning from locally owned sizes.
  <br>
  (Denis Davydov, 2019/03/19)
 </li>

 <li>
  New: VectorTools::interpolate_to_different_mesh() now works with
  hp::DoFHandler objects.
  <br>
  (Sebastian Stark, 2019/03/16)
 </li>

 <li>
  Improved: DoFRenumbering::cell_wise() now support parallel::Triangulation .
  <br>
  (Denis Davydov, 2019/03/16)
 </li>

 <li>
  Improved: The apply_all_indicators_to_manifold flag in GridIn::read_ucd()
  lets the indicators be used for cells as manifold id as well.
  <br>
  (Daniel Arndt, 2019/03/16)
 </li>

 <li>
  New: Added two new classes to the MeshWorker namespace: MeshWorker::ScratchData and
  MeshWorker::CopyData, that can be used as good default classes with the
  WorkStream::run() and MeshWorker::mesh_loop() functions.
  <br>
  (Luca Heltai, 2019/03/13)
 </li>

 <li>
  New: If present, we will detect the LLD linker ld.lld at configuration time.
  <br>
  (Timo Heister, 2019/03/11)
 </li>

 <li>
  Changed: The TrilinosWrappers::PreconditionAMG::AdditionalData data structure
  is now able to return a parameter list, which can be adjusted and fine-tuned by
  the user and later used to initialize the AMG preconditioner. It can also
  initialize the null space settings of an existing parameter list.
  <br>
  (Jean-Paul Pelteret, 2019/03/10)
 </li>

 <li>
  Fixed: DoFHandler::renumber_dofs() called on levels would previously run into
  an assertion on distributed triangulation when used with continuous
  elements. This is now fixed.
  <br>
  (Martin Kronbichler, 2019/03/08)
 </li>

 <li>
  New: Added a new function Utilities::type_to_string() to demangle type names.
  <br>
  (Luca Heltai, 2019/03/07)
 </li>

 <li>
  Fixed: LinearAlgebra::CUDAWrappers::Vector::add_and_dot was giving a wrong result
  if the size of vector was greater than 4096.
  <br>
  (Bruno Turcksin, 2018/03/06)
 </li>

 <li>
  Fixed: EllipticalManifold::push_forward_gradient() did previously not take the
  rotation into account. This is now fixed.
  <br>
  (Daniel Appel, Martin Kronbichler, 2019/03/05)
 </li>

 <li>
  Fixed: CylindricalManifold::get_new_point would previously return wrong
  results if the resulting point is on the axis, and the cylinder axis does not
  pass through the origin.
  <br>
  (Daniel Appel, Martin Kronbichler, 2019/03/05)
 </li>

 <li>
  Fixed: the DoFTools::make_sparsity_pattern variant that takes, as arguments, two
  different DoFHandler objects now works correctly with
  parallel::shared::Triangulation and parallel::distributed::Triangulation.
  <br>
  (David Wells, 2019/03/04)
 </li>

 <li>
  New: Signals Triangulation::pre_distributed_repartition and
  Triangulation::post_distributed_repartition which will be triggered in
  parallel::distributed::Triangulation::repartition().
  <br>
  (Marc Fehling, 2019/03/03)
 </li>

 <li>
  Fixed: CellAccessor::material_id interfered with the manifold object
  for refining the mesh in case its dimension was less than the space dimension.
  <br>
  (Daniel Arndt, Sebastian Stark, 2019/03/02)
 </li>

 <li>
  New: Hierarchy of finite elements in hp::FECollection objects. Get succeeding
  and preceding indices via hp::FECollection::next_in_hierarchy() and
  hp::FECollection::previous_in_hierarchy(). By default, a hierarchy corresponding
  to indices is established. Hierarchy can be overridden via
  hp::FECollection::set_hierarchy().
  <br>
  (Marc Fehling, 2019/02/28)
 </li>

 <li>
  Fixed: Tensor-valued pvtu output.
  <br>
  (Daniel Jodlbauer, 2019/02/27)
 </li>

 <li>
  New: %Function GridTools::guess_point_owner(), which uses a covering rtree
  to guess which processes own the points of the given vector of points.
  <br>
  (Giovanni Alzetta, 2019/02/25)
 </li>

 <li>
  Fixed: FEFaceEvaluation::evaluate() and FEFaceEvaluation::integrate() would
  provide wrong results for faces in non-standard orientation layout, namely
  with `orientation=true, face_flip={false,true}, face_rotation=true`. This is
  now fixed.
  <br>
  (Martin Kronbichler, 2019/02/21)
 </li>

 <li>
  New: A new function
  DoFRenumbering::random(DoFHandlerType &, const unsigned int) which allows for a
  random renumbering of degrees of freedom on a level in a mutilevel hierarchy.
  <br>
  (Conrad Clevenger, 2019/02/19)
 </li>

 <li>
  Fixed: Handle the case properly that a local matrix with all diagonal elements equal to zero is distributed with AffineConstraints::distribute_local_to_global() to the global matrix while constraints apply.
  <br>
  (Sebastian Stark, 2019/02/13)
 </li>

 <li>
  Fixed: inverse_operator() now handles a (composite) LinearOperator, the
  temporary PreconditionIdentity(), or no argument as preconditioner argument
  correctly.
  <br>
  (Matthias Maier, 2019/02/13)
 </li>

 <li>
  New: %Function MGConstrainedDoFs::make_no_normal_flux_constraints() which adds
  functionality for no normal flux constraints during geometric mutigrid computations.
  Currently, this function is limited to meshes with no normal flux boundaries
  normal to the x-, y-, or z-axis.
  <br>
  (Conrad Clevenger, 2019/02/11)
 </li>

 <li>
  New: Member functions DoFHandler::set_fe() and hp::DoFHandler::set_fe()
  that will register a finite element object without enumerating all
  degrees of freedom. In the hp case, the active_fe_indices will be
  initialized and communicated amongst all processors, additionally.
  <br>
  (Marc Fehling, 2019/02/10)
 </li>

 <li>
  New: make_const_array_view() creates a constant view from a non-const object.
  <br>
  (Daniel Arndt, 2019/02/09)
 </li>

 <li>
  Fixed: PreconditionChebyshev would previously run one iteration more than
  requested, i.e., perform `degree+1` matrix-vector products rather than
  `degree` which is the convention in literature. This is now fixed. Note that
  the quality of PreconditionChebyshev is obviously slightly reduced by this
  change, and some solvers might need more outer iterations due to a lighter
  Chebyshev iteration.
  <br>
  (Martin Kronbichler, 2019/02/08)
 </li>

 <li>
  New: Add default constructor to SparsityPatternIterators::Accessor().
  <br>
  (Denis Davydov, 2019/02/08)
 </li>

 <li>
  Improved: The iteration of PreconditionChebyshev has been overhauled to reduce
  the number of vectors written per iteration to one, leading to slightly faster
  execution of the vector updates.
  <br>
  (Martin Kronbichler, 2019/02/07)
 </li>

 <li>
  Improved: The FEEvaluation::evaluate() and FEEvaluation::integrate() routines
  would previously unconditionally use an evaluation routine with a
  transformation to a collocation basis and subsequent collocation derivative
  for `n_q_points_1d > fe_degree` because it is faster for the usual case of
  `n_q_points_1d` close to the polynomial degree. Now, we switch back to the
  standard evaluation if `n_q_points_1d > 3 * (fe_degree+1) / 2` where that
  variant is faster.
  <br>
  (Martin Kronbichler, 2019/02/06)
 </li>

 <li>
  Improved: Use binary search in BlockIndices::global_to_local() to improve performance
  for large number of blocks.
  <br>
  (Denis Davydov, 2019/02/06)
 </li>

 <li>
  Fixed: MGLevelGlobalTransfer::copy_to_mg would initialize the level vectors
  without ghosts. When combined with Chebyshev smoother and a DiagonalMatrix as
  well as specific matrix-free loops, this could lead to ghosted vectors in
  places where they were not supposed to be ghosted, leading to exceptions. This
  is now fixed.
  <br>
  (Martin Kronbichler, 2019/02/05)
 </li>

 <li>
  Documented: The minimal supported MPI version is MPI-2.0.
  <br>
  (Daniel Arndt, 2019/02/02)
 </li>

 <li>
  New: added constructor to the class Point, which converts
  a boost::geometry::model::point to a dealii::Point
  <br>
  (Giovanni Alzetta, Daniel Arndt, Wolfgang Bangerth, 2019/01/31)
 </li>

 <li>
  Changed: %Function GridGenerator::general_cell() can now generate
  a cell of dimension `dim` inside a space of dimension `spacedim`
  with `dim <= spacedim`
  <br>
  (Giovanni Alzetta, 2019/01/30)
 </li>

 <li>
  New: Add internal::MatrixFreeFunctions::DoFInfo::get_dof_indices_on_cell_batch() to
  return locally owned DoFs used by matrix-free framework on a given cell.
  <br>
  (Denis Davydov, 2019/01/30)
 </li>

 <li>
  Improved: A new dummy enumeration Differentiation::AD::NumberTypes::none has
  been added. It exists to represent number types that are scalar arithmetic
  types, i.e those that hold no derivatives. They are implemented primarily to
  facilitate the use of template meta-programming techniques to switch between
  different AD types. This covers the case when the user does not want an AD
  type at all, but rather a primitive type.
  <br>
  (Jean-Paul Pelteret, 2019/01/25)
 </li>

 <li>
  New: Add default constructor to DynamicSparsityPatternIterators::Accessor() and
  DynamicSparsityPatternIterators::Iterator() to be able to store such objects
  in STL containers.
  <br>
  (Denis Davydov, 2019/01/18)
 </li>

 <li>
  New: Indentation script <code>deal.II/contrib/utilities/indent.py</code> is added
  to format code in external projects that use deal.II
  <br>
  (Vishal Boddu, 2019/01/13)
 </li>

 <li>
  New: GridTools::Cache::get_covering_rtree() automatically stores and generates
  a rtree of bounding boxes covering the whole triangulation. In a distributed setting
  it can identify which processes have locally owned cells in any part of the mesh.
  <br>
  (Giovanni Alzetta, 2019/01/11)
 </li>

 <li>
  Fixed: Tensor indices of set_k_vectors() in FESeries::Fourier.
  Coefficients will now be distributed in every
  coordinate direction correctly.
  <br>
  (Marc Fehling, 2019/01/08)
 </li>

 <li>
  New: LinearAlgebra::distributed::Vector with CUDA memory space uses CUDA-aware
  MPI in case deal.II is configured with DEAL_II_MPI_WITH_CUDA_SUPPORT=ON.
  <br>
  (Daniel Arndt, Bruno Turcksin, 2019/01/04)
 </li>

 <li>
  New: Add LAPACKFullMatrix::apply_givens_rotation() and Vector::apply_givens_rotation() to
  apply a Givens rotation matrix.
  <br>
  (Denis Davydov, 2018/12/30)
 </li>

 <li>
  New: Add SolverBFGS to minimize a function using the limited memory BFGS approach.
  <br>
  (Denis Davydov, 2018/12/26)
 </li>

 <li>
  Replaced: TrilinosWrappers::MPI::Vector::vector_partitioner() has been
  deprecated in favor of TrilinosWrappers::MPI::Vector::trilinos_partitioner().
  The latter one returns an Epetra_BlockMap instead of an Epetra_Map avoiding an
  invalid downcast.
  The IndexSet constructor taking an Epetra_Map object has been replaced to
  accept a base class object of type Epetra_BlockMap instead.
  <br>
  (Daniel Arndt, 2018/12/22)
 </li>

 <li>
  Fixed: Reinit blockvector only when necessary in MGTransferBlockMatrixFree::copy_to_mg.
  <br>
  (Daniel Jodlbauer, 2018/12/19)
 </li>

 <li>
  New: Class parallel::distributed::CellDataTransfer has been introduced to transfer
  cell by cell data across distributed meshes in case of refinement/serialization.
  <br>
  (Marc Fehling, 2018/12/18)
 </li>

 <li>
  New: %Function GridTools::build_global_description_tree which exchanges a given
  vector of bounding boxes on all processes and uses the result to build an Rtree
  with packing algorithm.
  <br>
  (Giovanni Alzetta, 2018/12/17)
 </li>

 <li>
  New: Add helpers functions to use varying coefficients when using
  CUDAWrappers::MatrixFree.
  <br>
  (Bruno Turcksin, 2018/12/10)
 </li>

 <li>
  New: GridTools::find_active_cell_around_point now allows you to specify an (optional) rtree, constructed
  from the used vertices of the triangulation. Once you have built a tree, querying for a nearest
  vertex is an O(log(N)) operation, where N is the number of used vertices. You can ask a GridTools::Cache
  object to return a tree that is compatible with the new function signature.
  The previous version of this function had a cost that was O(N^2) when the point was not in the cell_hint
  object.
  <br>
  (Luca Heltai, 2018/12/05)
 </li>

 <li>
  New: Add LAPACKFullMatrix::transpose() to perform out-of-place transposition.
  <br>
  (Denis Davydov, 2018/12/05)
 </li>

 <li>
  New: Add DEAL_II_LAPACK_WITH_MKL macro which will be set if Intel MKL is being used.
  <br>
  (Denis Davydov, 2018/12/05)
 </li>

 <li>
  Changed: FESeries::Fourier::calculate and FESeries::Legendre::calculate
  accept both Vector<float> and Vector<double> as input parameters for
  local dof values.
  <br>
  (Marc Fehling, 2018/12/04)
 </li>

 <li>
  New: Added entries to GridTools::Cache to build a vertices RTree, and a bounding boxes RTree.
  <br>
  (Luca Heltai, 2018/12/03)
 </li>

 <li>
  New: Class EllipticalManifold derived from ChartManifold,
  valid only for dim=2 and spacedim=2. It maps points from a system of
  cartesian coordinates to a system of elliptical coordinates and vice-versa.
  <br>
  (Stefano Dominici, 2018/11/30)
 </li>

 <li>
  New: A new header <code>deal.II/base/undefine_macros.h</code> has been
  added. This header, as the name suggests, undefines all deal.II macros that are
  not prefixed with either <code>DEAL</code> or <code>deal</code>.
  <br>
  (David Wells, 2018/11/30)
 </li>

 <li>
  Fixed: The internal::MatrixFreeFunctions::FaceSetup::initialize() function
  would sometimes run into an assertion when ghost faces occur on subdomain
  interfaces. This is now fixed.
  <br>
  (Martin Kronbichler, 2018/11/28)
 </li>

 <li>
  New: Triangulation::n_faces, Triangulation::n_active_faces, and Triangulation::n_raw_faces now work also
  in one dimension.
  <br>
  (Luca Heltai, 2018/11/22)
 </li>

 <li>
  New: Add support for MPI when using CUDAWrappers::MatrixFree. This only works
  with LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
  <br>
  (Bruno Turcksin, 2018/11/21)
 </li>

 <li>
  Fixed: The TriaAccessor<0,dim,spacedim> and TriaAccessor<0,1,spacedim>
  classes were missing member functions `get_triangulation()` and
  `operator<()`. This prevented these objects from being stored in
  `std::map`s, among other issues.
  <br>
  (Wolfgang Bangerth, Luca Heltai, 2018/11/20)
 </li>

 <li>
  Fix bug with FEFaceValues initialization when used with the FE_NedelecSZ
  <br>
  (Alexander Grayver, 2018/11/19)
 </li>

 <li>
  New: A new function GridTools::get_coarse_mesh_description has been added. This
  function extracts the vertices, CellData, and SubCellData objects from a
  Triangulation (i.e., it is the inverse of Triangulation::create_triangulation).
  <br>
  (David Wells, 2018/11/11)
 </li>

 <li>
  Fixed: The operation FEEvaluation::read_dof_values performed invalid array
  accesses with AVX512 vectorization in rare circumstances of particular 1d
  meshes. This is now fixed.
  <br>
  (Martin Kronbichler, 2018/11/10)
 </li>

 <li>
  Deprecated: The binary Tecplot output using the DataOutBase::tecplot_binary
  format has been deprecated.
  <br>
  (Daniel Arndt, 2018/11/10)
 </li>

 <li>
  New: GridOut::write_vtk and GridIn::read_vtk are now compatible with each other, and can be used to save
  and restore a Triangulation with all manifold_ids and boundary_ids correctly saved.
  This is now the only format that supports both boundary_id and manifold_id both in reading and writing.
  <br>
  (Luca Heltai, 2018/11/07)
 </li>

 <li>
  Improved: The class ArrayView checks the memory location of the data it presents
  and gained a template parameter indicating whether the data is stored in CPU or
  CUDA memory.
  <br>
  (Daniel Arndt, 2018/11/04)
 </li>

 <li>
  Improved: GridIn::read_msh also support version 4 of MSH file format
  extracting the same information as for the older versions.
  <br>
  (Daniel Arndt, 2018/11/04)
 </li>

 <li>
  Fixed: The grid generated by GridGenerator::truncated_cone now uses more level 0
  cells to improve conditioning of cells when they are refined.
  <br>
  (David Wells, 2018/11/03)
 </li>

 <li>
  Changed: The class Subscriptor and SmartPointer check for dangling pointers and
  use-after-move. The destructor of Subscriptor doesn't signal an error when there
  is still an object subscribed to it. Instead, validity is checked when
  dereferencing the SmartPointer object.
  <br>
  (Daniel Arndt, 2018/11/02)
 </li>

 <li>
  Fixed: SparsityPattern::begin() and DynamicSparsityPattern::begin()
  resulted in an invalid iterator if the sparsity pattern was actually
  empty.
  <br>
  (Wolfgang Bangerth, Denis Davydov, 2018/11/01)
 </li>

 <li>
  New: Add Functions::IncrementalFunction which returns the difference in values of
  a given function at two time steps.
  <br>
  (Denis Davydov, Jean-Paul Pelteret, 2018/10/05)
 </li>

 <li>
  New: Classes FESeries::Fourier and FESeries::Legendre now take two
  template arguments dim and spacedim, which allow them to use corresponding
  hp::FECollection objects.
  <br>
  (Marc Fehling, 2018/10/17)
 </li>

 <li>
  Fixed: The project_boundary_values_div_conforming() function now works correctly
  with parallel::distributed::Triangulation.
  <br>
  (Jonathan Matthews, 2018/10/13)
 </li>

 <li>
  New: Functions DoFAccessor::get_active_fe_indices() and
  internal::DoFAccessorImplementation::get_active_vertex_fe_indices()
  return all active finite element indices on the corresponding object.
  <br>
  (Marc Fehling, 2018/10/08)
 </li>

 <li>
  New: Added RTree class and pack_rtree() function to create a boost::geometry::index::rtree from
  containers and/or iterators of BoundingBox, Point, or Segment objects.
  <br>
  (Luca Heltai, 2018/10/06)
 </li>

 <li>
  The computation of normal vectors on very small cells could fail. This is now fixed.
  <br>
  (Timo Heister, 2018/10/05)
 </li>

 <li>
  The class VectorizedArray has gained support for the intrinsics of Altivec
  of IBM Power processors.
  <br>
  (Martin Kronbichler, Bruno Turcksin, Sambit Das, 2018/10/05)
 </li>

 <li>
  New: Add Utilities::inverse_Hilbert_space_filling_curve() to map points in
  dim-dimensional space to line index of the Hilbert space filling curve.
  <br>
  (Denis Davydov, 2018/10/05)
 </li>

 <li>
  Deprecated: The Threads::n_existing_threads() and
  Threads::this_thread_id() functions have been deprecated. There are
  equivalent C++11 functions that provide the same kind of
  functionality.
  <br>
  (Wolfgang Bangerth, 2018/10/05)
 </li>

 <li>
  New: There is now a constructor of hp::QCollection that, like the
  existing one of hp::FECollection, takes a variable number of
  quadrature objects and creates a collection from it.
  <br>
  (Wolfgang Bangerth, 2018/10/04)
 </li>

 <li>
  New: Class parallel::CellWeights allows to set the
  weight of a cell depending on the number of degrees of freedom
  associated with it. This guarantees an improved load balancing
  in case a hp::DoFHandler is used on a
  parallel::distributed::Triangulation.
  <br>
  (Marc Fehling, 2018/10/01)
 </li>

 <li>
  Removed: The Threads::DummyMutex, Threads::DummyConditionVariable, and
  Threads::DummyThreadBarrier classes have been removed. These were only
  used if deal.II was compiled without thread support, in which case
  certain functionality from the TBB was not available. But the
  implementations of the underlying functionality is now unconditionally
  available via the C++11 standard library, and so we don't need the
  dummy implementations any more.
  <br>
  (Wolfgang Bangerth, 2018/09/27)
 </li>

 <li>
  Deprecated: The Threads::PosixThreadBarrier class and the
  corresponding Threads::Barrier `using` declaration are deprecated. The
  functionality implemented by these classes is easily replaced by C++11
  functionality.
  <br>
  (Wolfgang Bangerth, 2018/09/27)
 </li>

 <li>
  Deprecated: The Threads::ConditionVariable class is now
  deprecated. Use `std::condition_variable` instead.
  <br>
  (Wolfgang Bangerth, 2018/09/27)
 </li>

 <li>
  Deprecated: The Threads::Mutex class is now derived from
  std::mutex. As a consequence, the Threads::Mutex::ScopedLock class is
  now no longer necessary and is deprecated in favor of std::lock_guard
  or any of the other similar classes provided by C++11 and
  later. Similarly, the Threads::Mutex::acquire() and
  Threads::Mutex::release() functions are deprecated in favor of the
  corresponding functions in the std::mutex base class.
  <br>
  (Wolfgang Bangerth, 2018/09/27)
 </li>

 <li>
  New: Support for complex values can be controlled using the CMake flag
  DEAL_II_WITH_COMPLEX_VALUES.
  <br>
  (Daniel Arndt, 2019/09/26)
 </li>

 <li>
  Fixed: The GridIn::read_abaqus() function now correctly reads in grids designed
  for the co-dimension 1 case (specifically with dim == 2 and spacedim == 3).
  <br>
  (Yuxiang Wang, Jean-Paul Pelteret, 2018/09/20)
 </li>

 <li>
  Fixed: MatrixFree::loop() would not zero a vector if all degrees of freedom
  are constrained. This is now fixed.
  <br>
  (Martin Kronbichler, 2018/09/20)
 </li>

 <li>
  New: Add optional flag to ParameterHandler to skip undefined
  sections and parameters.
  <br>
  (Denis Davydov, 2018/09/20)
 </li>

 <li>
  New: extend Patterns::Tools::Convert to work with ComponentMask
  and std::unique_ptr<FunctionParser>.
  <br>
  (Denis Davydov, 2018/09/19)
 </li>

 <li>
  New: SolverSelector gained a new constructor initializing the object completely.
  <br>
  (Daniel Arndt, 2018/09/17)
 </li>

 <li>
  New: added the function TrilinosWrappers::MPI::Vector::import() to get data from a ReadWriteVector.
  <br>
  (Timo Heister, 2018/09/16)
 </li>

 <li>
  New: Add line minimization functions.
  <br>
  (Denis Davydov, 2018/09/13)
 </li>

 <li>
  Fixed: <code>make test</code> now runs in parallel automatically (i.e., it uses
  a CMake command to determine the number of cores available). Since <code>make
  test</code> just immediately invokes <code>ctest</code> as a subprocess, it is
  currently not possible to link the number of globally available
  <code>make</code> processes with the number of concurrently run tests; i.e.,
  the number of cores <code>make test</code> uses is not changeable by users.
  <br>
  (David Wells, 2018/09/13)
 </li>

 <li>
  Fixed: VectorTools::project_boundary_values_div_conforming works correctly
  in case the DoFHandler object contains more than one FE_RaviartThomas element.
  <br>
  (Daniel Arndt, 2018/09/13)
 </li>

 <li>
  New: The function DoFRenumbering::hierarchical() is now also available
  for hp::DoFHandler objects, as well as for DoFHandlers build on
  triangulations with `dim != spacedim`.
  <br>
  (Wolfgang Bangerth, 2018/09/06)
 </li>

 <li>
  Added: New function
  Utilities::MPI::compute_n_point_to_point_communications()
  is implemented to be used for computing number of processes
  in an MPI universe to expect communication from.
  <br>
  (Eldar Khattatov, 2018/09/04)
 </li>

 <li>
  New: Utilities::CUDA::allocate_device_data() and
  Utilities::CUDA::delete_device_data() help to manage data on CUDA devices
  using std::unique_ptr.
  <br>
  (Daniel Arndt, 2018/09/01)
 </li>

 <li>
  Fixed: Member function AffineConstraints::add_entries() now skips already
  existing entries.
  <br>
  (Dustin Kumor, 2018/08/24)
 </li>

 <li>
  Fixed: TrilinosWrappers::SparseMatrix::mmult() sometimes gave the wrong result
  in parallel.
  <br>
  (Damien Lebrun-Grandie and Bruno Turcksin, 2018/08/22)
 </li>

 <li>
  Fixed: The function Particles::ParticleAccessor::serialize() was
  declared but not implemented. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2018/08/21)
 </li>

 <li>
  New: Add a mask argument to filter out some vector lanes in FEEvaluationBase.
  This is required for local time stepping using the MatrixFree framework,
  because some cells of batches must be excluded from read/write operations
  when operating on different time steps compared to their neighbor cells.
  <br>
  (Svenja Schoeder, 2018/08/21)
 </li>

 <li>
  New: Table and FullMatrix now have complete (both mutable and const)
  random-access iterators.
  <br>
  (David Wells, 2018/08/20)
 </li>

 <li>
  New: CUDAWrappers::PreconditionILU and CUDAWrappers::PreconditionIC can be used
  for preconditioning CUDAWrappers::SparseMatrix objects.
  <br>
  (Daniel Arndt, 2018/08/20)
 </li>

 <li>
  New: There are nothrow versions of AssertCUDA and AssertCusparse, namely
  AssertNothrowCUDA and AssertNothrowCusparse.
  <br>
  (Daniel Arndt, 2018/08/20)
 </li>

 <li>
  New: There are new methods CUDAWrappers::SparseMatrix::print() and
  CUDAWrappers::SparseMatrix::print_formatted(). Additionally,
  CUDAWrappers::SparseMatrix gained a move assignment operator.
  <br>
  (Daniel Arndt, 2018/08/20)
 </li>

 <li>
  Fixed: Throwing an exception in Assert doesn't slice objects derived from
  ExceptionBase anymore.
  <br>
  (Daniel Arndt, 2018/08/17)
 </li>

 <li>
  New: There is LinearAlgebra::CUDAWrappers::Vector::norm_sqr() for computing
  the squared l2-norm and LinearAlgebra::CUDAWrappers::Vector::swap()
  for swapping two vectors.
  <br>
  (Daniel Arndt, 2018/08/17)
 </li>

 <li>
  Fixed: LinearAlgebra::CUDAWrappers::Vector performs now a deep copy instead of
  a shallow copy in its copy assignment operator.
  <br>
  (Daniel Arndt, 2018/08/16)
 </li>

 <li>
  Fixed: MGConstrainedDoFs::initialize() now handles refinement edges across periodic boundaries correctly.
  <br>
  (Alexander Knieps, 2018/08/13)
 </li>

 <li>
  Fixed: Bug using ScaLAPACKMatrix::save() and load() for configurations with pHDF5 in combination with serial HDF5 matrices from user code.
  <br>
  (Benjamin Brands, 2018/08/10)
 </li>

 <li>
  Improved: Triangulation::prepare_coarsening_and_refinement also takes care of
  limiting the level difference across periodic boundary faces.
  <br>
  (Daniel Arndt, Stefan Kaessmair, 2018/08/09)
 </li>

 <li>
  New: The new function
  MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis()
  enables to change the basis from a representation in quadrature points into
  the original basis of FEEvaluation.
  <br>
  (Martin Kronbichler, 2018/08/06)
 </li>

 <li>
  Improved: FEEvaluation and FEFaceEvaluation have gained support for
  renumbering of degrees of freedom for fully discontinuous finite element
  spaces which can improve performance by around 10 percent.
  <br>
  (Martin Kronbichler, 2018/08/06)
 </li>

 <li>
  New: Add support for constraint when using matrix-free on GPU. Note you cannot
  set constraints in 2D if the degree of the finite element is even.
  <br>
  (Bruno Turcksin, 2018/08/06)
 </li>

 <li>
  New: A convenience function Triangulation::active_face_iterators() (which is like
  Triangulation::active_cell_iterators(), but for faces) has been added.
  <br>
  (David Wells, 2018/08/05)
 </li>

 <li>
  Fixed: Minor issues when using matrix-free cell categorization
  with multi-grid operators.
  <br>
  (Denis Davydov, 2018/08/02)
 </li>

 <li>
  Fixed: The import function for LinearAlgebra::distributed::Vector would fail if
  the IndexSet of the ReadWriteVector has no locally owned element
  <br>
  (Bruno Turcksin, 2018/08/02)
 </li>

 <li>
  New: GridGenerator::merge_triangulations can merge multiple
  Triangulation objects at once.
  <br>
  (Daniel Arndt, 2018/08/01)
 </li>

 <li>
  New: LinearAlgebra::distributed::Vector has an extra template parameter
  MemorySpace. This new template parameter can be used if the data should be
  allocated on the host or on the device.
  <br>
  (Bruno Turcksin, 2018/08/01)
 </li>

 <li>
  Added support for high-order VTK output by using newly
  introduced Lagrange VTK cells.
  <br>
  (Alexander Grayver, 2018/07/27)
 </li>

 <li>
  New: Add FiniteSizeHistory class to store a finite-size history of objects.
  <br>
  (Denis Davydov, 2018/07/26)
 </li>

 <li>
  Deprecate VectorTools::project_boundary_values_curl_conforming in favor
  of VectorTools::project_boundary_values_curl_conforming_l2. The former
  has a number of known limitations and bugs.
  <br>
  (Alexander Grayver, 2018/07/26)
 </li>

 <li>
  New: Add boost adaptor for dealii::BoundingBox.
  <br>
  (Luca Heltai, 2018/07/24)
 </li>

 <li>
  New: Added geometry adaptors for boost.
  <br>
  (Luca Heltai, Bruno Turcksin, 2018/07/22)
 </li>

 <li>
  New: Add new constructor type for FunctionParser class.
  <br>
  (Luca Heltai, 2018/07/21)
 </li>

 <li>
  Fixed: hp::DoFHandler::distribute_dofs() now works on
  parallel::distributed::Triangulation<3> objects that
  contain artificial cells.
  <br>
  (Marc Fehling, 2018/07/20)
 </li>

 <li>
  Changed: The Subscriptor class has been made thread-safe, and its behavior in
  debug and release modes are now identical. Although this makes subscription and
  un-subscription slightly more expensive in release mode, debugging should be
  easier overall.
  <br>
  (Jean-Paul Pelteret, 2018/07/15)
 </li>

 <li>
  New: A new function GridGenerator::concentric_hyper_shells has been added. This
  function generates a grid comprised of shells of varying widths, where the
  shells near the center are much thinner than those near the outer boundary.
  <br>
  (David Wells, 2018/07/15)
 </li>

 <li>
  Fixed: Iterating over the entries of a TrilinosWrappers::SparsityPattern object
  in parallel and after calling TrilinosWrappers::SparsityPattern::compress()
  works now.
  <br>
  (Mathias Anselmann, Daniel Arndt, 2018/07/11)
 </li>

 <li>
  Fixed: GridTools::delete_duplicated_vertices also treats the SubCellData
  object correctly when renumbering the vertices.
  <br>
  (Daniel Arndt, 2018/07/08)
 </li>

 <li>
  New: GridGenerator::merge_triangulations gained an option to copy manifold ids.
  <br>
  (Daniel Arndt, 2018/07/08)
 </li>

 <li>
  New: Add GridGenerator::plate_with_a_hole() that generates a
  rectangular plate with an (offset) cylindrical hole.
  <br>
  (Denis Davydov, 2018/06/21)
 </li>

 <li>
  Fixed: The functions FEFaceEvaluation::evaluate and
  FEFaceEvaluation::integrate would access invalid memory when only evaluating
  function values and using multiple components in 1D. This is now fixed.
  <br>
  (Martin Kronbichler, 2018/07/04)
 </li>

 <li>
  New: The class MatrixFreeOperators::CellwiseInverseMassMatrix can now also be
  used in 1D.
  <br>
  (Martin Kronbichler, 2018/07/04)
 </li>

 <li>
  New: Add support for output of tensor-valued data in
  DataOut::write_vtu().
  <br>
  (Denis Davydov, 2018/06/21)
 </li>

 <li>
  New: Extend DataOutInterface::get_vector_data_ranges() to return extra element
  in the tuple to represent DataComponentInterpretation. The functions were
  renamed to get_nonscalar_data_ranges().
  <br>
  (Denis Davydov, 2018/06/21)
 </li>

 <li>
  Fixed: FETools::back_interpolate can be used
  for TrilinosWrappers::MPI::BlockVector.
  <br>
  (Daniel Arndt, 2018/06/20)
 </li>

 <li>
  Changed: Calling GridRefinement::refine_and_coarsen_fixed_fraction()
  and GridRefinement::refine_and_coarsen_fixed_fraction() with values for
  the coarsen and refine fractions that added to 1.0 in exact
  arithmetic, but something slightly larger in floating point
  arithmetic, resulted in an error. These functions are now slightly
  more relaxed in that they allow arguments that add to 1.0 up to a
  small floating point round off error.
  <br>
  (Wolfgang Bangerth, 2018/06/18)
 </li>

 <li>
  New: The LinearAlgebra::distributed::Vector::compress() function has
  gained two new combine operations VectorOperation::min and
  VectorOperation::max to define the entry at the owning processor as
  the minimum/maximum of its own value and the one coming from a ghost.
  <br>
  (Svenja Schoeder, 2018/06/18)
 </li>

 <li>
  Fixed: TrilinosWrappers::MPI::Vector::reinit now checks if
  IndexSet::is_ascending_and_one_to_one is `true` for
  TrilinosWrappers::MPI::Vector::parallel partitioner before calling
  IndexSet::make_trilinos_map because IndexSet::make_trilinos_map may be able
  to make a linear map if the property holds.
  <br>
  (Joshua Hanophy, 2018/06/16)
 </li>

 <li>
  Fixed: The eigenvalue approximation in SolverGMRES failed if no iterations
  were performed. Now, there is simply no output in this case.
  <br>
  (Daniel Arndt, 2018/06/15)
 </li>

 <li>
  Bugfix: Fixed bugs for TrilinosWrappers::SparseMatrix::add()
  and TrilinosWrappers::SparseMatrix::copy_from() for the case of
  non-contiguous rows.
  <br>
  (Uwe Koecher, 2018/06/14)
 </li>

 <li>
  Fixed: Mismatched numbers of components between the exact solution and
  underlying FE in Step-51.
  <br>
  (Pi-Yueh Chuang, 2018/06/14)
 </li>

 <li>
  Bugfix: Fixed a bug where internal::VectorDataExchange::ghosts_were_set was
  not set inside MatrixFree::cell_loop() for block vectors
  with many blocks.
  <br>
  (Denis Davydov, 2018/06/12)
 </li>

 <li>
  Fixed: CMake now supports again in-source builds (This functionality was
  accidentally broken in the previous release).
  <br>
  (Matthias Maier, 2018/06/11)
 </li>

 <li>
  Deprecated: FunctionMap is deprecated. It is recommended to
  use the aliased type directly.
  <br>
  (Daniel Arndt, 2018/06/11)
 </li>

 <li>
  Fixed: BlockLinearOperator::vmult(), BlockLinearOperator::vmult_add(),
  BlockLinearOperator::Tvmult() and BlockLinearOperator::Tvmult_add()
  can be used with identical source and destination vectors.
  <br>
  (Daniel Arndt, 2018/06/10)
 </li>

 <li>
  Improved: Extend ScaLAPACKMatrix::invert() to use pXtrtri for inversion of triangular matrices.
  <br>
  (Sambit Das, 2018/06/06)
 </li>

 <li>
  New: TimerOutput::cpu_and_wall_times_grouped will print CPU and
  wallclock time in a single table.
  <br>
  (Denis Davydov, 2018/05/29)
 </li>

 <li>
  Improved: TimerOutput now dynamically adjust to the width of Section column.
  <br>
  (Denis Davydov, 2018/05/29)
 </li>

 <li>
  New: The function GridTools::find_all_active_cells_around_point() returns a
  list of all cells around a point, rather than only a single one for the old
  GridTools::find_active_cell_around_point.
  <br>
  (Martin Kronbichler, Niklas Fehn, 2018/05/28)
 </li>

 <li>
  Improved: MatrixFree::cell_loop() and MatrixFree::loop() directly use
  LinearAlgebra::distributed::BlockVector::update_ghost_values() and
  LinearAlgebra::distributed::BlockVector::compress() calls on block vectors
  with many blocks, rather than splitting each method into two parts for
  overlapping communication and computation. The latter is inefficient once
  too many MPI requests are in flight.
  <br>
  (Martin Kronbichler, Denis Davydov, 2018/05/24)
 </li>

 <li>
  Fixed: Step-16 was broken caused by uninitialized variables, it is now fixed.
  <br>
  (Ce Qin, 2018/05/21)
 </li>

 <li>
  Improved: MGTransferPrebuilt now supports
  PETScWrappers::MPI::Vector and
  PETScWrappers::MPI::SparseMatrix.
  <br>
  (Alexander Knieps, 2018/05/21)
 </li>

 <li>
  Fixed: DoFTools::make_flux_sparsity_pattern was accessing invalid dof indices
  in case neighboring cells were not using FiniteElement objects with the same
  number of degrees of freedom. This is fixed now.
  <br>
  (Daniel Arndt, 2018/05/19)
 </li>

 <li>
  Fixed: Two invalid off-by-one data accesses in the initialization of
  MatrixFree with face data enabled, that appeared in some rare cases on certain
  configurations of processors, have been fixed.
  <br>
  (Daniel Arndt, Martin Kronbichler, 2018/05/18)
 </li>

 <li>
  New: Classes that act as an interface to the drivers of internally supported automatic
  differentiation libraries (ADOL-C and Sacado) have been implemented. A brief summary
  of what these classes are used for is given in the @ref auto_symb_diff module.
  <br>
  (Jean-Paul Pelteret, 2017/11/18)
 </li>

</ol>

*/
