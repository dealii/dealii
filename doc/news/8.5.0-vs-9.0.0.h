// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
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
@page changes_between_8_5_0_and_9_0_0 Changes between Version 8.5.0 and 9.0.0

<p>
This is the list of changes made between the release of deal.II version
8.5.0 and that of 9.0.0. All entries are signed with the names of the
author.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="850-900-incompatible"></a>
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
  New: Every function in the GridGenerator namespace now attaches a default manifold to the
  curved parts of the boundary of the domain, and sets reasonable defaults for manifold indicators in the
  domain, independently on the boundary indicators.
  <br>
  (Luca Heltai, Daniel Arndt, Jean-Paul Pelteret, David Wells, 2018/04/17)
 </li>

 <li>
  New: The implementation of the divergence in FEValuesExtractors::Tensor
  was changed so that the gradient operator is contracted from the right.
  This is done in order to make it consistent with gradient of the
  second order tensor, namely $Grad(T) : I = Div(T)$.
  <br>
  (Denis Davydov, 2018/04/12)
 </li>

 <li>
  New: The Manifold class now requires the implementation of a Manifold::clone() method in
  all derived classes.
  <br>
  (Luca Heltai, 2018/04/08)
 </li>

 <li>
  New: The class FE_DGQHermite has been switched to the more stable polynomial
  basis Polynomials::HermiteLikeInterpolation. This changes the representation
  of the numbers by the FE_DGQHermite class. The change is yet highly beneficial
  because it significantly improves the accuracy (in terms of roundoff) for this
  basis and also reduces iteration counts for some iterative solvers with simple
  preconditioner ingredients such as point-Jacobi or point-Jacobi smoothing.
  <br>
  (Martin Kronbichler, 2018/02/09)
 </li>

 <li>
  Changed: FETools::get_fe_by_name returns a std::unique_ptr instead of a
  owning raw pointer to prevent creating memory leaks.
  <br>
  (Daniel Arndt, 2018/02/09)
 </li>

 <li>
  Changed: The function Mapping::clone (and all inheriting classes' implementations, such as MappingQ1Eulerian::clone) now return a <code>unique_ptr<Mapping<dim,spacedim>></code>.
  <br>
  (David Wells, 2018/02/08)
 </li>

 <li>
  Changed: The GridRefinement functions only allow Vector objects
  to be passed in as criteria.
  <br>
  (Daniel Arndt, 2018/02/04)
 </li>

 <li>
  Changed: TableIndices uses std::size_t for index access,
  LAPACKFullMatrix::size_type is types::blas_int and
  FullMatrix::size_type is std::size_t.
  <br>
  (Daniel Arndt, 2018/02/02)
 </li>

 <li>
  Changed: The classes OpenCASCADE::NormalProjectionBoundary,
  OpenCASCADE::DirectionalProjectionBoundary
  and OpenCASCADE::NormalToMeshProjectionBoundary are rebased
  from Boundary to FlatManifold.
  <br>
  (Daniel Arndt, 2018/01/05)
 </li>

 <li>
  Changed: It is now necessary to explicitly enable MPI support when
  configuring deal.II via -DDEAL_II_WITH_MPI=ON, or -DWITH_MPI=ON on the
  command line. The old behavior of just specifying an MPI compiler wrapper
  via environment, or CMAKE_(C|CXX)_COMPILER does not work any more.
  <br>
  (Matthias Maier, 2018/01/01)
 </li>

 <li>
  Changed: The type <code>LAPACKFullMatrix::size_type</code> is now <code>unsigned int</code>, which matches <code>FullMatrix::size_type</code>.
  <br>
  (David Wells, 2017/12/20)
 </li>

 <li>
  Changed: Change the logic of DoFTools::get_subdomain_association to have the
  degrees of freedom along a refinement edge be now all given to the processor
  with the smallest subdomain_id.
  <br>
  (Alexander Grayver, 2017/12/11)
 </li>

 <li>
  Changed: The parameters contained in struct SolverQMRS::AdditionalData have
  been changed in order to meet the rewritten algorithm. See the class documentation
  for more details.
  <br>
  (Ingo Kligge, 2017/12/06)
 </li>

 <li>
  Changed: The default partitioner for parallel::shared::Triangulation
  is Zoltan.
  <br>
  (Daniel Arndt, 2017/11/28)
 </li>

 <li>
  Changed: The minimal supported Trilinos version is 12.4.
  <br>
  (Daniel Arndt, 2017/11/18)
 </li>

 <li>
  Removed: The default constructor for ArrayView was removed since
  an object constructed in this way can't be used in any sensible way.
  <br>
  (Daniel Arndt, 2017/11/17)
 </li>

 <li>
  Changed: The many KellyErrorEstimator::estimate() overloads all have
  an argument that represents the Neumann boundary values (if the
  solution of the PDE satisfies Neumann boundary conditions on parts of
  the boundary). They used to just be of type Function<dim>, but this
  has now been changed to Function<dim,number> (or, more precisely,
  FunctionMap<dim,number>, which stores a map to Function<dim,number>
  objects) where @p number is the underlying data type of the input
  vector. In other words, if you are using a regular Vector<double> (or
  PETSc or Trilinos equivalents), then nothing has changed since @p
  number is @p double, which is the default for the Function class. On
  the other hand, if you are using a Vector<std::complex<double>>, for
  example, then you need to pass a Function<dim,std::complex<double>> --
  which of course makes sense because the Neumann boundary values of a
  complex-valued solution should also be complex-valued.
  <br>
  (Wolfgang Bangerth, 2017/11/17)
 </li>

 <li>
  Changed: The distribution of degrees of freedom along a refinement edge for a parallel::shared::Triangulation has been changed to mimic that of a parallel::distributed::Triangulation, that is, they are now all given to the processor with the smallest subdomain_id.
  <br>
  (Conrad Clevenger, 2017/11/13)
 </li>

 <li>
  Changed: ParameterHandler::print_parameters() in LaTeX format now correctly
  escapes special characters. As a consequence, the names of labels changed
  because they are now mangled. See the documentation for more details.
  <br>
  (Timo Heister, 2017/11/03)
 </li>

 <li>
  Removed: The BZIP2 dependency has been removed as it was not used
  in the library anymore.
  <br>
  (Daniel Arndt, 2017/10/25)
 </li>

 <li>
  Changed: The direction argument of the DoFRenumbering::downstream()
  and related functions have been changed from a Point<dim> to a
  Tensor<1,dim> type, since that is the correct data type for a
  direction vector.
  <br>
  (Wolfgang Bangerth, 2017/10/23)
 </li>

 <li>
  Changed: LAPACKSupport::Properties is renamed to LAPACKSupport::Property to be
  consistent with LAPACKSupport::State.
  <br>
  (Denis Davydov, 2017/08/24)
 </li>

 <li>
  Changed: The virtual function SolverRichardson::criterion() now
  receives the residual and preconditioned residual vectors as
  arguments, rather than accessing it through the member variables of
  the class. It has also been made @p const.
  <br>
  (Wolfgang Bangerth, 2017/09/12)
 </li>

 <li>
  Changed: The methods
  <ol>
    <li>Manifold::get_new_point</li>
    <li>Manifold::project_to_manifold</li>
  </ol>
  have had their declarations changed in an incompatible manner: these three
  methods now take ArrayView arguments instead of <code>std::vector</code>s. This
  change was done for performance reasons: before this change about a third of the
  time spent generating a curved grid (in this particular benchmark, a circle
  described with polar coordinates and a transfinite interpolation) was used in
  allocating and freeing memory used by the <code>std::vector</code> arguments to
  Manifold::get_new_point() even though the sizes of the arrays are usually known
  at compilation time. This change completely eliminates this allocation cost.
  The method <code>Manifold::add_new_points()</code> has been removed in favor of
  Manifold::get_new_points(), which also uses ArrayView arguments instead of
  <code>std::vector</code>s.
  <br>
  (David Wells, 2017/09/10)
 </li>

 <li>
  Changed: deal.II now requires BOOST version 1.59 or newer.
  <br>
  (David Wells, 2017/09/04)
 </li>

 <li>
  Changed: LAPACKSupport::upper_triangle and LAPACKSupport::lower_triangle are
  renamed to LAPACKSupport::upper_triangular and LAPACKSupport::lower_triangular,
  respectively.
  <br>
  (Denis Davydov, 2017/08/24)
 </li>

 <li>
  Changed: SolverGMRES used to notify via deallog when re-orthogonalization
  of the Arnoldi vectors kicks in. Now, there is a slot one can connect to
  retrieve this information. By default nothing is printed.
  <br>
  (Daniel Arndt, 2017/09/02)
 </li>

 <li>
  Fixed: FiniteElement::get_generalized_support_points() now always returns a
  list of unique points. This is in contrast to the old behavior where
  get_generalized_support_points() returned a repeated list of (nodal)
  support points for an FESystem consisting of Lagrangian elements.
  <br>
  (Matthias Maier, 2017/08/31)
 </li>

 <li>
  Changed: Specialization of the ProductType class are now implemented through
  specialization of the class internal::ProductTypeImpl . This was done in order
  to ensure that product operations performed with qualified number types do not
  result in the intended specializations being overlooked by the compiler.
  <br>
  (Jean-Paul Pelteret, Wolfgang Bangerth, 2017/08/24)
 </li>

 <li>
  Changed: The VectorMemory::Pointer class used to have an automatic
  conversion operator to the underlying pointer-to-vector object. As
  part of a rewrite of this class, this conversion (rarely used in
  practice) has been removed.
  <br>
  (Wolfgang Bangerth, 2017/08/22)
 </li>

 <li>
  Changed: The constructors SymmetricTensor (const Tensor &) and
  Tensor (const array_type &) have been marked as 'explicit'.
  <br>
  (Daniel Arndt, 2017/08/20)
 </li>

 <li>
  Changed: The configuration file <code>deal.II/base/config.h</code> has been
  thoroughly cleaned up. As a result, the following preprocessor symbols (which
  were either never used or were workarounds for older, now unsupported compilers)
  are no longer defined:
  <ul>
    <li> <code>DEAL_II_EXPLICIT_CONSTRUCTOR_BUG</code>
    <li> <code>DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG</code>
    <li> <code>DEAL_II_BOOST_BIND_COMPILER_BUG</code>
    <li> <code>DEAL_II_HAVE_ISNAN</code>
    <li> <code>DEAL_II_STD_ISNAN</code>
    <li> <code>DEAL_II_HAVE_UNDERSCORE_ISNAN</code>
    <li> <code>DEAL_II_HAVE_ISFINITE</code>
    <li> <code>DEAL_II_HAVE_SYS_TIMES_H</code>
    <li> <code>DEAL_II_HAVE_TIMES</code>
    <li> <code>DEAL_II_HAVE_SYS_TYPES_H</code>
    <li> <code>DEAL_II_HAVE_SYS_TIMES_H</code>
  </ul>
  <br> (Matthias Maier and David Wells, 2017/08/08 - 2017/08/24)
 </li>

 <li>
  Changed: The member function Boundary::project_to_surface is now ignored in the
  function GridTools::fix_up_distorted_child_cells. Instead,
  GridTools::fix_up_distorted_child_cells uses an internal function that
  duplicates the algorithm used by StraightBoundary::project_to_surface, which was
  the only place project_to_surface was implemented in the library.
  <br>
  (Luca Heltai, David Wells, 2017/08/12)
 </li>

 <li>
  Changed: Triangulation::get_manifold() will now return a FlatManifold instead of
  a StraightBoundary in the case that no manifold description has been attached to
  a provided manifold ID.
  <br> (Luca Heltai, David Wells, 2017/08/12)
 </li>

 <li>
  Removed: Some old functionality has been removed from LogStream:
  <ul>
    <li>The reproducible test functionality has been removed.
      This includes <code>test_mode()</code>, <code>threshold_double()</code>,
      and <code>threshold_float()</code>. Use numdiff, or other tools that can
      do output comparison with small numerical differences instead.</li>
    <li>The <code>log_cerr()</code> function has been removed.
    <li>The timing mechanism has been removed. This includes
      <code>timestamp()</code>, <code>log_execution_time()</code> and
      <code>log_time_differences()</code>. Use the Timer, TimerOutput and
      TimerOutput::Scope classes instead.</li>
    <li>LogStream now defaults to write to std::cout instead of
      std::cerr</li>
  </ul>
  <br>
  (Matthias Maier, 2017/08/10)
 </li>

 <li>
  Changed: DataOutBase (and derived classes) no longer sets floating point
  output precision of ostreams.
  <br>
  (Matthias Maier, 2017/08/09)
 </li>

 <li>
  Changed: The types::boundary_id and types::material_id types have been changed
  from unsigned char to unsigned int.
  <br>
  (Jean-Paul Pelteret, 2017/08/07)
 </li>

 <li>
  Deprecated: A number of classes have been deprecated in favor of the
  LinearOperator concept:
  <ul>
    <li>PointerMatrix, PointerMatrixAux
  </ul>
  Use the LinearOperator class instead, see the module on
  @ref LAOperators "linear operators" for more details.
  <br>
  (Matthias Maier, 2017/07/14)
 </li>

 <li>
  Changed: In 1D, MappingQGeneric now considers a cell to have multiple manifolds if the manifold attached to the cell does not match one of the manifolds attached to a face (i.e., vertex).
  <br>
  (David Wells, 2017/07/10)
 </li>

 <li>
  Changed: The SynchronousIteartors::iterators member variable has been
  made private, as mentioned in the changelog of the previous
  release. It can be accessed via SynchronousIterators::operator*(),
  however.
  <br>
  (Wolfgang Bangerth, 2017/07/11)
 </li>

 <li>
  Changed: The TableIndices class had a single constructor that just
  took as many indices as one wanted, and padded any unassigned indices
  with numbers::invalid_unsigned_int. This invited mistakes and indeed
  led to difficult to track down errors.
  <br>
  This constructor is therefore now deprecated, and has been replaced by
  a series of constructors with between one and five arguments that can
  be used for objects of type TableIndices<1> to TableIndices<5>.
  <br>
  (Wolfgang Bangerth, 2017/07/10)
 </li>

 <li>
  Changed: In 1D, the face values calculated by <code>MappingManifold</code> now use the manifold attached to the relevant vertex instead of the manifold on the current face.
  <br>
  (David Wells, 2017/07/08)
 </li>

 <li>
  Changed: Timer::print_data doesn't restrict its output to the first MPI process
  anymore but leaves such a choice to the stream given.
  <br>
  (Daniel Arndt, 2017/06/27)
 </li>

 <li>
  Changed: The restart_parameter member variable of
  TrilinosWrappers::SolverGMRES::AdditionalData has been removed.
  This value is now stored as gmres_restart_parameter in the base class
  TrilinosWrappers::SolverBase::AdditionalData.
  <br>
  (Jean-Paul Pelteret, 2017/06/26)
 </li>

 <li>
  Changed: The deprecated data in SolverCG::AdditionalData
  and SolverGMRES::AdditionalData::compute_eigenvalues have been removed.
  Use the respective connect_* member functions instead.
  <br>
  (Daniel Arndt, 2017/06/12)
 </li>

 <li>
  Changed: The deprecated member functions add(), normalize(), conjugate(),
  abs(), sadd(), equ() and mult() in the vector classes have been removed.
  <br>
  (Daniel Arndt, 2017/06/12)
 </li>

 <li>
  Deprecated: The ParameterHandler::print_parameters_section() function
  has been deprecated.
  <br>
  (Wolfgang Bangerth, 2017/06/07)
 </li>

 <li>
  Changed: Versions of PETSc prior to 3.3.0 are no longer supported.
  <br>
  (David Wells, 2017/06/03)
 </li>

 <li>
  Changed: The FiniteElement::clone() function and all of its
  implementations in derived classes now return a std::unique_ptr rather
  than a plain pointer to a finite element object. User-implemented
  finite element classes need to be adjusted accordingly.
  <br>
  (Wolfgang Bangerth, 2017/06/01)
 </li>

 <li>
  Changed: The deprecated function SparsityTools::reorder_Cuthill_McKee()
  acting on a SparsityPattern object has been removed.
  Use the one acting on a DynamicSparsityPattern instead.
  <br>
  (Daniel Arndt, 2017/05/15)
 </li>

 <li>
  Changed: The deprecated functions BlockSparseMatrixEZ::n_rows()
  and BlockSparseMatrixEZ::n_cols() have been removed.
  Use BlockSparseMatrixEZ::m() and BlockSparseMatrixEZ::n() instead.
  <br>
  (Daniel Arndt, 2017/05/14)
 </li>

 <li>
  Changed: The deprecated function Utilities::System::job_supports_mpi()
  has been removed. Use Utilities::MPI::job_supports_mpi() instead.
  <br>
  (Daniel Arndt, 2017/05/14)
 </li>

 <li>
  Removed: The deprecated serial Trilinos vector classes have been removed.
  <br>
  (Daniel Arndt, 2017/05/08)
 </li>

 <li>
  Removed: The deprecated serial PETSc vector classes have been removed.
  <br>
  (David Wells, 2017/05/06)
 </li>

 <li>
  Changed: The deprecated member functions in the classes
  TrilinosWrappers::BlockVector, TrilinosWrappers::MPI::BlockVector,
  TrilinosWrappers::Vector and TrilinosWrappers::MPI::Vector
  have been removed.
  <br>
  (Daniel Arndt, 2017/05/05)
 </li>

 <li>
  Changed: The deprecated version of Manifold::get_new_point() that took
  an argument of type Quadrature has been removed. Use the other variant
  of that function instead. The same has been done to all
  implementations of that interface in derived classes.
  <br>
  (Wolfgang Bangerth, 2017/05/05)
 </li>

 <li>
  Changed: The deprecated typedef
  FEFieldFunction::ExcPointNotAvailableHere has been removed. Use
  VectorTools::ExcPointNotAvailableHere instead.
  <br>
  (Wolfgang Bangerth, 2017/05/04)
 </li>

 <li>
  Changed: The deprecated variant of
  GridTools::get_face_connectivity_of_cells() has been removed. Use the
  other variant instead.
  <br>
  (Wolfgang Bangerth, 2017/05/04)
 </li>

 <li>
  Changed: The deprecated constructors of MappingQEulerian
  and MappingQ1Eulerian have been removed. Use the other constructor
  of each class instead.
  <br>
  (Wolfgang Bangerth, 2017/05/04)
 </li>

 <li>
  Changed: The deprecated function FEValuesBase::transform()
  has been removed. It only forwarded the call to the Mapping object
  used by the FEValuesBase object. Use the corresponding function of the
  mapping instead.
  <br>
  (Wolfgang Bangerth, 2017/05/04)
 </li>

 <li>
  Changed: The deprecated member variables @p supports_distributed_data
  that was present in all vector classes has been removed. If you needed
  this functionality, use the type trait is_serial_vector instead.
  <br>
  (Wolfgang Bangerth, 2017/04/30)
 </li>

 <li>
  Changed: The deprecated functions DoFHandler::get_tria() and
  hp::DoFHandler::get_tria() have been removed. Use
  DoFHandler::get_triangulation() and
  hp::DoFHandler::get_triangulation() instead.
  <br>
  (Wolfgang Bangerth, 2017/04/25)
 </li>

 <li>
  Changed: The deprecated functions DataOutInterface::write_pvd_record()
  and DataOutBase::write_visit_record() have been removed. Use the
  corresponding functions in namespace DataOutBase.
  <br>
  (Wolfgang Bangerth, 2017/04/25)
 </li>

 <li>
  Changed: The deprecated function FEValuesBase::get_normal_vectors()
  that returned a vector of Point objects has been removed. Its
  replacement, FEValuesBase::get_all_normal_vectors() has now itself
  been deprecated, and we have created a new function
  FEValuesBase::get_normal_vectors() that returns a vector of
  Tensor<1,dim> objects. The net effect is that the function with the
  old name has simply gotten a new return type.
  <br>
  (Wolfgang Bangerth, 2017/04/25)
 </li>

 <li>
  Removed: DataPostprocessor had member functions
  <code>compute_derived_quantities_*()</code> that had previously
  already been deprecated. These have now been removed.
  <br>
  (Wolfgang Bangerth, 2017/04/24)
 </li>

 <li>
  Removed: The ParameterHandler::read_input() function and friends have
  been removed. They were already deprecated in the previous release.
  <br>
  (Wolfgang Bangerth, 2017/04/21)
 </li>

 <li>
  Changed: It was previously allowed to copy one FESystem object to
  another via the copy constructor. There are probably few reasons to do
  so, but they complicated the design of data structures. Consequently,
  it is now no longer allowed to copy such objects.
  <br>
  (Wolfgang Bangerth, 2017/04/17)
 </li>

 <li>
  Removed: The FETools::compute_node_matrix() variant that takes two
  arguments has been removed. It was previously already deprecated. Use
  the variant with just one argument instead.
  <br>
  (Wolfgang Bangerth, 2017/04/16)
 </li>

 <li>
  Changed: The <code>AssertGlobalIndexRange</code> macro has been removed: the expansion of this macro involved an undeclared template and always generates compiler errors if used.
  <br>
  (David Wells, 2017/04/13)
 </li>

 <li>
  Removed: A number of deprecated classes have been removed in favor of the
  new LinearOperator concept:
  <ul>
    <li>BlockDiagonalMatrix
    <li>InverseMatrixRichardson
    <li>IterativeInverse
    <li>MeanValueFilter
    <li>ProductMatrix
    <li>ProductSparseMatrix
    <li>ScaledMatrix
    <li>SchurMatrix
    <li>ShiftedMatrix
    <li>ShiftedMatrixGeneralized
    <li>TransposeMatrix
  </ul>
  Use the LinearOperator class instead, see the module on
  @ref LAOperators "linear operators" for more details.
  <br>
  (Matthias Maier, 2017/04/06 - 2018/05/03)
 </li>

 <li>
  Removed: The FiniteElement::interpolate() function and all of its
  implementations in derived classes has been removed. It was
  previously already deprecated. Use
  FiniteElement::convert_generalized_support_point_values_to_nodal_values()
  instead.
  <br>
  (Wolfgang Bangerth, 2017/04/05)
 </li>

 <li>
  Removed: The deprecated MPI communicator and the constructor which used it were
  removed from the MatrixFree::AdditionalData class.
  <br>
  (Denis Davydov, 2017/03/26)
 </li>

 <li>
  Removed: The deprecated CMake flag <code>DEAL_II_CXX11_FLAG</code> has been
  removed.
  <br>
  (David Wells, 2017/03/23)
 </li>

 <li>
  Changed: deal.II now requires a compiler supporting (very nearly) the entire C++11 standard. The minimal version of GCC supported is now 4.8.
  <br>
  (David Wells, 2017/03/23)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="850-900-general"></a>
<h3>General</h3>

<ol>

 <li>
  New: The tutorial step-59 shows a matrix-free solver for the Poisson
  equation discretized with the symmetric interior penalty discontinuous
  Galerkin method.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/05/04)
 </li>

 <li>
  New: The MatrixFree framework has been reworked to also support face integrals
  in DG. The new function MatrixFree::loop() takes three function pointers for
  workers on cells, interior faces and boundary faces. The loop can perform the
  data exchange with MPI and is also parallelized with threads taking into
  account the face dependencies. The new class FEFaceEvaluation implements the
  respective evaluation and access routines for face terms.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/04/27)
 </li>

 <li>
  New: Triangulation now uses std::unique_ptr to store Manifolds, and duplicates the Manifolds when
  Triangulation::set_manifold() is called. This allows the construction of the Manifold objects
  to be independent of the Triangulation. In particular, it will be possible to associate
  to the Triangulation reasonable default Manifold descriptors, that will be copied
  around together with the Triangulation itself.
  <br>
  (Luca Heltai, 2018/04/06)
 </li>

 <li>
  Changed: The bundled version of Intel Threading Building Blocks has been updated to 2018 U2.
  <br> (David Wells, 2018/03/02)
 </li>

 <li>
  Improved: The manifold smoothing algorithms applied in the Triangulation class
  and MappingQGeneric have been changed from the old Laplace-style smoothing to
  a transfinite interpolation that linearly blends between the descriptions on
  the faces around a cell. The old transformation introduced boundary layers
  inside cells that prevented convergence rates from exceeding 3.5 in the global
  L2 errors on typical settings. This change also considerably improves mesh
  quality on settings where curved descriptions are only applied to the boundary
  rather than the whole volume.
  <br>
  (Martin Kronbichler, 2017/12/04)
 </li>

 <li>
  New: Add ScaLAPACKMatrix and ProcessGrid wrappers for high-performance dense linear
  algebra routines for parallel distributed memory machines.
  <br>
  (Denis Davydov, Benjamin Brands, 2017/11/20)
 </li>

 <li>
  New: The new ParticleHandler class can store and organize a collection of particles,
  and provide information about their properties. It is currently limited to parallel
  distributed computations and specific applications, but will be extended over time.
  <br>
  (Rene Gassmoeller, 2017/11/10)
 </li>

 <li>
  New: The majority of deal.II classes that deal with solution vectors
  are now also instantiated for complex-valued vectors. In other words,
  complex-valued vectors should now be supported at the same level as
  real-valued vectors.
  <br>
  (Wolfgang Bangerth, 2017/11/08)
 </li>

 <li>
  New: The new namespace Rol contains an adaptor class that provides
  the implementation of the ROL::Vector interface.
  The Trilinos package, Rapid Optimization Library (ROL), can solve unconstrained
  and constrained optimization problems, and optimization problems under
  uncertainty.
  <br>
  (Vishal Boddu 2017/11/02)
 </li>

 <li>
  New: Using the Physics::Notation::Kelvin class, it is possible to store and
  retrieve tensors and symmetric tensors in or from a compressed format.
  <br>
  (Jean-Paul Pelteret, 2017/10/16)
 </li>

 <li>
  New: The function distribute_mg_dofs has been written for a parallel::shared::Triangulation. This allows for geometric multigrid computations on an adaptively refined mesh using a shared triangulation with the possibility of a user defined partitioning of the active and level cells.
  <br>
  (Conrad Clevenger, 2017/10/10)
 </li>

 <li>
  New: Added support for the KINSOL solver of the SUNDIALS
  library. KINSOL is a solver for nonlinear algebraic systems.
  It includes a Newton-Krylov solver as well as Picard and
  fixed point solvers, both of which can be accelerated with
  Anderson acceleration.
  <br>
  (Luca Heltai, 2017/09/28)
 </li>

 <li>
  New: Added support for the ARKode solver of the SUNDIALS
  library. ARKode is a solver library that provides
  adaptive-step time integration of the initial value problem
  for systems of stiff, nonstiff, and multi-rate systems of
  ordinary differential equations (ODEs) given in linearly
  implicit form.
  <br>
  (Luca Heltai, 2017/09/27)
 </li>

 <li>
  New: There is now support for the storage of Particles and their
  properties in the new namespace Particles.
  <br>
  (Rene Gassmoeller, 2017/09/22)
 </li>

 <li>
  New: A new ParameterAcceptor class has been added to the library.
  The class is intended to be used as a base class for any class
  that wants to handle parameters using the ParameterHandler class.
  If you derive all your classes from ParameterAcceptor, and declare
  your parameters either with parse/declare_parameters methods or
  via the ParameterAcceptor::add_parameter() method, then the declaring
  and parsing of your parameter files will be automatically managed
  by the ParameterAcceptor::initialize() function.
  <br>
  (Luca Heltai, 2017/09/20)
 </li>

 <li>
  New: Added support for the Open Asset Import Library (Assimp)
  (https://assimp.sourceforge.net/). This library can be used
  to read about 40 different 3D graphics formats, used in 3D
  modelers (such as Blender, Maya, etc.). Some of these formats
  contain mesh information, that in turn can be read
  into deal.II Triangulation<2,3> objects.
  <br>
  (Luca Heltai, 2017/09/16)
 </li>

 <li>
  New: A new GridTools::Cache class
  allows caching of some expensive data of the
  Triangulation class, computed generally using
  functions in GridTools.
  <br>
  (Luca Heltai, 2017/09/14)
 </li>

 <li>
  New: A new SUNDIALS::IDA class has been added that interfaces
  the SUNDIALS IDA library (https://computation.llnl.gov/projects/sundials)
  This class can be used to solve differential algebraic equations
  using IDA (an implicit differential algebraic equation solver, that
  supports variable step and variable order BDF schemes).
  <br>
  (Luca Heltai, 2017/09/12)
 </li>

 <li>
  Improved: The support for non-Lagrangian elements with generalized support
  points has been vastly improved. FESystem::get_generalized_support_point()
  now returns a list of unique generalized support points for the finite
  element system. In order to do interpolation a function
  FiniteElement::convert_generalized_support_points_to_dof_values() can be
  used. This interface is implement for a wide variety of base classes and
  FESystem. TODO
  <br> (Matthias Maier, Luca Heltai, 2017/08/31)
 </li>

 <li>
  New: A new KDTreeDistance class has been added that interfaces
  the nanoflann library (https://github.com/jlblancoc/nanoflann).
  This class can be used to extract nearest neighbour information
  on collection of points, query for the closest points to a target
  point or all points contained within a given distance.
  <br>
  (Luca Heltai, 2017/08/13)
 </li>

 <li>
  New: A new MeshWorker::mesh_loop() function has been added
  that performs the same tasks of the MeshWorker::loop() function
  without forcing the users to adhere to a specific interface.
  <br>
  (Luca Heltai, 2017/08/12)
 </li>

 <li>
  New: The eigenvectors of a rank-2 symmetric tensor can now be computed using one
  of three approaches through the eigenvectors() function. The three algorithms
  that have been implemented are:
  1. The QL algorithm with implicit shifting.
  2. A hybrid algorithm that preferentially uses an analytical algorithm
  and falls back to the QL algorithm if the calculations are deemed
  inaccurate.
  3. The Jacobi algorithm.
  <br>
  (Joachim Kopp, Jean-Paul Pelteret, Ester Comellas, 2017/07/27)
 </li>

 <li>
  New: The eigenvalues of a rank-2 symmetric tensor can now be computed using an
  analytical approach via the eigenvalues() function.
  <br>
  (Jean-Paul Pelteret, Ester Comellas, 2017/07/27)
 </li>

 <li>
  New: The new namespace Patterns::Tools contains
  utilities that can be used to convert from complicated types
  to strings and vice versa. These tools have been put to use
  in the method ParameterHandler::add_parameter() that allows
  users to perform in one single call the operations
  ParameterHandler::declare_parameter(),
  ParameterHandler::get() and to convert the string to a valid
  value stored in the variable that is given as a parameter to
  ParameterHandler::add_parameter().
  <br>
  (Luca Heltai, 2017/07/20)
 </li>

 <li>
  New: The cmake configuration now supports unity builds with the option <code>-DDEAL_II_UNITY_BUILD=ON</code>. This option speeds up the build by about 10 to 25%.
  <br> (David Wells, 2017/07/20)
 </li>

 <li>
  New: The hp::DoFHandler class can now work with shared triangulations
  of type parallel::shared::Triangulation.
  <br>
  (Wolfgang Bangerth, 2017/07/17)
 </li>

 <li>
  New: SolverFIRE implements FIRE (Fast Inertial Relaxation Engine) for solving
  the problem of minimization of a given objective function.
  <br>
  (Vishal Boddu, Denis Davydov 2017/07/11)
 </li>

 <li>
  New: A new manifold class TransfiniteInterpolationManifold implementing an
  interpolation from a curved boundary description to the interior has been
  added. This class enables high-order convergence rates of more than three in
  the power of the mesh size for situations where a curved manifold can only be
  prescribed to the boundary but not in a whole volume.
  <br>
  (Martin Kronbichler, Luca Heltai, 2017/06/01)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="850-900-specific"></a>
<h3>Specific improvements</h3>

<ol>

 <li>
  New: A new tutorial program, step-60, shows how to deal with non-matching discretizations using
  Distributed Lagrange Multipliers.
  <br>
  (Luca Heltai, Giovanni Alzetta, 2018/05/07)
 </li>

 <li>
  New: Gmsh::create_triangulation_from_boundary_curve is a new method that constructs a grid from an
  OpenCASCADE TopoDS_Edge object. It offers the reverse functionality of the method
  OpenCASCADE::create_curves_from_triangulation_boundary().
  <br>
  (Luca Heltai, 2018/05/03)
 </li>

 <li>
  New: NonMatching::create_coupling_sparsity_pattern and NonMatching::create_coupling_mass_matrix now have an
  overloaded version that takes an additional GridTools::Cache, instead of
  building it inside the functions. In the process, NonMatching::create_coupling_sparsity_pattern
  gained also an additional ConstraintMatrix argument to reflect the same argument list of
  its companion function NonMatching::create_coupling_mass_matrix.
  <br>
  (Luca Heltai, 2018/05/02)
 </li>

 <li>
  New: Utilities::MPI::create_group allows using the functionality provided by
  MPI_Comm_create_group to be used also for a non-MPI-3.0-conforming MPI
  dependency.
  <br>
  (Daniel Arndt, 2018/04/28)
 </li>

 <li>
  New: There are new functions FEEvaluation::gather_evaluate() and
  FEEvaluation::integrate_scatter() that combine the vector access through
  FEEvaluation::read_dof_values() and FEEvaluation::evaluate() as well as
  FEEvaluation::integrate() and
  FEEvaluation::distribute_local_to_global(), respectively. This allows to
  write somewhat more compact code and is also more efficient in the case
  of FEFaceEvaluation where equivalent functions exists.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/04/27)
 </li>

 <li>
  New: The function MatrixFree::cell_loop takes an optional boolean
  argument that enables the class to zero vectors within the loop for
  better cache locality, clearing only that part of the vector that will
  be accessed soon.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/04/27)
 </li>

 <li>
  New: The new field
  MatrixFree::AdditionalData::cell_vectorization_categories allows to
  control the way batches of cells are formed in MatrixFree by assigning
  different numbers to different parts in the domain.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/04/27)
 </li>

 <li>
  New: The class FEEvaluation can now be constructed for only a subset of
  the components inside an FESystem, picking e.g. the pressure part
  outside of a Taylor-Hood element.
  <br>
  (Katharina Kormann, Martin Kronbichler, 2018/04/27)
 </li>

 <li>
  New: FECollection gained an equality comparison operator. Both FiniteElement
  and FECollection have a non-equality comparison operator now.
  <br>
  (Daniel Arndt, 2018/04/25)
 </li>

 <li>
  New: Added missing implementation of MappingFEField::get_vertices.
  <br>
  (Luca Heltai, 2018/04/24)
 </li>

 <li>
  Updated: deal.II is now compatible with PETSc version 3.9.0, SLEPC version
  3.9.0
  <br>
  (Matthias Maier, 2018/04/23)
 </li>

 <li>
  New: GridTools::minimal_cell_diameter and GridTools::maximal_cell_diameter now take an optional Mapping
  argument, that allows one to compute minimal and maximal cell diameters of deformed grids.
  <br>
  (Luca Heltai, 2018/04/23)
 </li>

 <li>
  Fixed: GridTools::distributed_compute_point_locations now returns the correct
  maps values (third component of the output tuple). Added a test to check it.
  <br>
  (Giovanni Alzetta, 2018/04/23)
 </li>

 <li>
  New: Add direct solvers (Cholesky and LU factorization) to solve problems on the
  GPU
  <br>
  (Bruno Turcksin, 2018/04/23)
 </li>

 <li>
  New: Utilities::System::get_current_vectorization_level() returns the currently used vectorization support in string format.
  <br>
  (Timo Heister, 2018/04/21)
 </li>

 <li>
  Changed: In accordance with dealii::DoFHandler, hp::DoFHandler stores a copy
  of the FECollection instead of a pointer.
  <br>
  (Daniel Arndt, 2018/04/21)
 </li>

 <li>
  Add default constructor and member function `initialize` to hp::DoFHandler,
  making the interface consistent with DoFHandler.
  <br>
  (Ce Qin, 2018/04/20)
 </li>

 <li>
  New: Cell weighting can now be taken into account by the Metis partitioner. This can be done directly
  using the new GridTools::partition_triangulation() function that accepts a vector of cell weights,
  or by adding the appropriate signal to the triangulation itself (Triangulation::Signals::cell_weight).
  <br>
  (Jean-Paul Pelteret, 2018/04/17)
 </li>

 <li>
  Changed: GridTools::distributed_compute_point_locations now takes as input the
  global description of the manifold using bounding boxes.
  <br>
  (Giovanni Alzetta, 2018/04/16)
 </li>

 <li>
  New: Introduced Triangulation::reset_manifold(), as a substitute for Triangulation::set_manifold()
  with only one argument. The set_manifold() method with one argument only was used to reset the manifold
  object. This method is now deprecated, in favor of the new, more explicative, method.
  An additional Triangulation::reset_all_manifolds() method has been added to remove all manifold
  objects from the triangulation.
  <br>
  (Luca Heltai, 2018/04/14)
 </li>

 <li>
  New: GridTools::map_boundary_to_manifold_ids() allows you to set manifold ids on the
  boundary according to a given map of boundary to manifold ids.
  <br>
  (Luca Heltai, 2018/04/13)
 </li>

 <li>
  Fixed: Ensure that numbers::NumberTraits and its corresponding operations are
  well defined for all supported auto-differentiable numbers.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  Fixed: The SymmetricTensor class was previously not instantiated for
  auto-differentiable numbers.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  Fixed: When using diagonal SymmetricTensors with automatic-differentiable numbers, computations
  using the eigenvalue() and eigenvector() functions would return correct values of the
  eigenvalues/vectors. However, the derivatives of these values/vectors were incorrect as
  the sensitivities of the eigenvalues/vectors with respect to one another was not encoded
  in the returned result. Therefore, under these specific conditions the returned result is now
  a more coarse approximation for the eigenvalues/vectors but with the trade-off that a
  meaningful approximation of the derivative of the result can now be computed.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  Fixed: The Tensor class was previously not instantiated for
  auto-differentiable numbers.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  Fixed: The Tensor::invert() function would previously not work with some Sacado
  number types. This has now been fixed.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  Fixed: Some compile-times errors would previously appear for some arithmetic functions
  when using Sacado::Fad::DFad types (e.g. Physics::Elasticity::Kinematics::F_iso() ).
  By defining the ProductType Sacado expression templates, this issue is now avoided.
  <br>
  (Jean-Paul Pelteret, 2018/04/12)
 </li>

 <li>
  New: FEValuesExtractors::Tensor now supports calculation of gradients.
  <br>
  (Denis Davydov, 2018/04/12)
 </li>

 <li>
  New: The two new functions Polynomials::jacobi_polynomial_value() and
  Polynomials::jacobi_polynomial_roots() provide a user-visible
  implementation of Jacobi polynomials. This functionality has previously only
  been available internally in the library.
  <br>
  (Martin Kronbichler, 2018/04/06)
 </li>

 <li>
  New: Patterns::Tools::to_string() and Patterns::Tools::to_value() simplify the conversion to and from
  strings of arbitrarily complex types.
  <br>
  (Luca Heltai, 2018/04/06)
 </li>

 <li>
  New: Extend GridGenerator::extrude_triangulation to using exact slicing z-axis
  values and reimplement the existing GridGenerator::extrude_triangulation using
  the newly developed overload.
  <br>
  (Weixiong Zheng, 2018/04/05)
 </li>

 <li>
  New: There is now a function Utilities::dynamic_unique_cast() that does for
  `std::unique_ptr` objects what `dynamic_cast` does for regular pointers.
  <br>
  (Wolfgang Bangerth, 2018/04/03)
 </li>

 <li>
  Improved: VectorTools::get_position_vector() now supports parallel triangulations.
  <br>
  (Denis Davydov, 2018/04/01)
 </li>

 <li>
  New: Multigrid classes obtained the new signals
  mg::Signals::transfer_to_mg,
  mg::Signals::transfer_to_global,
  mg::Signals::coarse_solve,
  mg::Signals::restriction,
  mg::Signals::prolongation,
  mg::Signals::pre_smoother_step,
  mg::Signals::post_smoother_step
  that functions can be connected to.
  <br>
  (Daniel Arndt, Timo Heister, 2018/03/31)
 </li>

 <li>
  Fixed: The VectorTools::integrate_difference() function allows users
  to provide a weight function that can also serve as a component mask
  to select individual components of the solution vector for error
  computation. For components not selected, such a mask would then
  simply be zero.
  <br>
  In some cases, the solution vector contains NaN numbers, for example
  when one uses the FE_FaceQ element for certain components of the
  solution vector and uses a quadrature formula for error evaluation
  that has quadrature points in the interior of the cell. For any
  "regular" solution component for which the component mask has a zero
  weight, the value of that component will be multiplied by zero and
  consequently does not add anything to the error computation. However,
  if the NaNs of a FE_FaceQ are multiplied with zero weights, the result
  is still a NaN, and adding it to the values times weights of the other
  components results in NaNs -- in effect rendering it impossible to get
  any information out of the VectorTools::integrate_difference()
  function if one of the finite elements involved is FE_FaceQ.
  <br>
  This is now fixed by simply skipping vector components for which the
  weight vector is zero. This has the same result as before for all
  "normal" situations, but also properly skips the NaN case outlined
  above.
  <br>
  (Wolfgang Bangerth, 2018/03/29)
 </li>

 <li>
  Fixed: Threads::Thread<T> and Threads::Task<T> objects can now also be
  used when @p T is a type that is only move-constructable or movable,
  but not necessarily copy-constructible or copyable. Consequently, one
  can now also call functions that return such objects on tasks and
  threads. In particular, this is relevant for functions that return
  std::unique_ptr objects that have this property.
  <br>
  (Wolfgang Bangerth, 2018/03/28)
 </li>

 <li>
  Fixed: LinearAlgebra::distributed::Vector::compress(VectorOperation::insert) now flashes
  ghost part of the vector in Release mode.
  <br>
  (Denis Davydov, 2018/03/17)
 </li>

 <li>
  Fixed: GridTools::minimal_cell_diameter() and GridTools::maximal_cell_diameter()
  return the maximal respectively minimal cell diameter over the whole mesh for
  parallel::distributed::Triangulation object, too.
  <br>
  (Daniel Arndt, 2018/03/25)
 </li>

 <li>
  New: Added a third test variant to the testsuite driver: For a test
  consisting of a test.prm.in (and a test.output) file the test.prm.in file
  will be configured/preprocessed to a test.prm file. This is done with the
  CMake macro CONFIGURE_FILE that replaces all strings \@VARIABLE\@ with the
  contents of the corresponding CMake variable. This is useful in particular
  to conveniently substitute \@SOURCE_DIR\@ with the full source directory path
  of the test.
  <br>
  (Matthias Maier, 2018/03/21)
 </li>

 <li>
  Fixed: Previously, it was not possible to iterate over the local range
  of rows of PETSc matrix objects because one would have to call
  `matrix.end(row)` where `row` is the last locally owned row, and that
  triggered an assertion because this end iterator is also the begin
  iterator of the next row -- which is not locally owned any more.
  <br>
  This is now fixed.
  <br>
  (Feimi Yu, Wolfgang Bangerth, 2018/03/20)
 </li>

 <li>
  Fixed: The SparsityPattern::copy_from() variant that takes a
  FullMatrix argument was previously of quadratic complexity in the
  number of nonzero entries per row. This is now fixed, and the function
  is now linear.
  <br>
  (Ben Shields, Wolfgang Bangerth, 2018/03/20)
 </li>

 <li>
  Fixed: FETools::get_interpolation_difference_matrix() used to not
  clear content of the matrix passed as argument, but instead just add
  to it. This resulted in wrong output, but is now fixed.
  <br>
  (Wolfgang Bangerth, 2018/03/18)
 </li>

 <li>
  New: Add LAPACKFullMatrix<number>::set(const size_type, const size_type, const number)
  to set an element of the matrix.
  <br>
  (Denis Davydov, 2018/03/17)
 </li>

 <li>
  Fixed: FullMatrix::residual(), FullMatrix::add_col(), and FullMatrix::add_row()
  now work correctly with rectangular matrices.
  <br> (David Wells, 2018/03/17)
 </li>

 <li>
  New: JSON files which are written out by the
  parameter handler can be read in again with the
  function ParameterHandler::parse_input_from_json().
  <br>
  (Menno Fraters, 2018/03/12)
 </li>

 <li>
  New: The testsuite can be run using valgrind via
  'ctest -S ../tests/run_memorycheck.cmake'.
  <br>
  (Daniel Arndt, 2018/03/11)
 </li>

 <li>
  New: The Utilities::dealii_version_string() function returns a string
  representation of the deal.II version being used.
  <br>
  (Wolfgang Bangerth, 2018/03/09)
 </li>

 <li>
  Fixed: Copying a Trilinos::MPI::Vector to a local deal.II Vector using
  operator=() resulted in only a partial copy of data to the local vector. In
  addition, the copied  elements  were offset incorrectly on each process.
  Similar held for copying a  Trilinos::MPI::BlockVector to a local deal.II
  BlockVector. This has been fixed and now works as expected.
  <br>
  (Jean-Paul Pelteret, 2018/03/07)
 </li>

 <li>
  New: LAPACKFullMatrix::compute_inverse_svd_with_kernel allows to
  set an expected kernel size for an inverse singular value decomposition.
  <br>
  (Daniel Arndt, 2018/03/07)
 </li>

 <li>
  New: The ParameterAcceptorProxy class allows you to wrap a class that provides a static member
  declare_parameters and a member parse_parameters into a ParameterAcceptor subclass.
  <br>
  (Luca Heltai, 2018/03/06)
 </li>

 <li>
  New: NonMatching::create_coupling_sparsity_pattern and NonMatching::create_coupling_mass_matrix functions
  to handle L2 projections between arbitrary non-aligned grids.
  <br>
  (Luca Heltai, 2018/03/05)
 </li>

 <li>
  Fixed: GridTools::distort_random works with
  parallel::shared::Triangulation objects.
  <br>
  (Daniel Arndt, 2018/03/05)
 </li>

 <li>
  New: A template class LinearIndexIterator has been added with the intent of
  using it to generalize iterators over contiguously stored data.
  <br> (David Wells, 2018/03/01)
 </li>

 <li>
  Changed: The GridRefinement::hyper_sphere() function used to have two
  template arguments (`dim` and `spacedim`), but it really only existed
  if <code>dim == spacedim-1</code>. Consequently, it now has lost its
  `dim` template argument and has only retained `spacedim`.
  <br>
  (Wolfgang Bangerth, 2018/01/25)
 </li>

 <li>
  New: TransposeTable (the base class of LAPACKFullMatrix and ScaLAPACKMatrix) now
  has a random access iterator implementation similar to the one provided by
  SparseMatrix.
  <br> (David Wells, 2018/02/25)
 </li>

 <li>
  New: A top level target 'expand_all_instantiations' generates all .inst files.
  <br>
  (Timo Heister, Daniel Arndt, 2018/02/20)
 </li>

 <li>
  New: LAPACKFullMatrix: make operators *= and /= use Lapack function,
  LAPACKFullMatrix::add() uses BLAS 1 routine instead of handwritten loops
  <br>
  (Benjamin Brands, 2018/02/20)
 </li>

 <li>
  Changed: The GridOut::write_mesh_per_processor_as_vtu function now
  only includes processor info in the .vtu filename and writes a .pvtu file
  when using a <code>parallel::Triangulation</code>.
  <br>
  (Conrad Clevenger, 2018/02/19)
 </li>

 <li>
  New: Add a new class CUDAWrappers::SparseMatrix, i.e., wrappers for cuSPARSE csr
  sparse matrix. The matrix is copied from  deal.II own's SparseMatrix and copied
  to the device.
  <br>
  (Bruno Turcksin, 2018/02/18)
 </li>

 <li>
  New: Add LAPACKFullMatrix<number>::Tmmult(LAPACKFullMatrix<number> &, const LAPACKFullMatrix<number> &, const Vector<number> &, const bool) const
  to do a triple matrix product with a diagonal matrix in the middle. Add
  LAPACKFullMatrix::scale_rows() to scale rows via a given vector. Make
  LAPACKFullMatrix::Tmmult() and LAPACKFullMatrix::mTmult() use Xsyrk if A==B.
  <br>
  (Denis Davydov, 2018/02/14)
 </li>

 <li>
  New: add parallel::distributed::BlockVector::mmult(BlockVector &, const FullMatrixType &, const NumberType, const NumberType) const
  and parallel::distributed::BlockVector::multivector_inner_product_with_metric(const FullMatrixType &, const BlockVector &V, const bool) const
  to operate on multivectors with a metric tensor.
  <br>
  (Denis Davydov, 2018/02/10)
 </li>

 <li>
  New: The new polynomial class Polynomials::HermiteLikeInterpolation is a
  modification of the Hermite polynomials with good conditioning of
  interpolation for all degrees, as opposed to the
  Polynomials::HermiteInterpolation class.
  <br>
  (Anian Fuchs, Martin Kronbichler, 2018/02/09)
 </li>

 <li>
  Extend VectorTools::project() function with MatrixFree quadrature data to
  optionally take the finite element component index.
  <br>
  (Denis Davydov, 2018/02/07)
 </li>

 <li>
  New: Added routines to perform addition, multiplication and scaling for ScaLAPACKMatrix
  <br>
  (Benjamin Brands, 2018/02/07)
 </li>

 <li>
  Fixed: MGTransferMatrixFree::restrict() and MGTransferMatrixFree::prolongate()
  would produce a segmentation fault when used with multi-component systems for
  polynomial degrees larger than 10. This is now fixed.
  <br>
  (Martin Kronbichler, 2018/01/30)
 </li>

 <li>
  New: Add AssertCusparse macro to assert the error code of cuSPARSE function
  <br>
  (Bruno Turcksin, 2018/01/29)
 </li>

 <li>
  New: There are new overloads of make_array_view for Tensor, SymmetricTensor,
  LAPACKFullMatrix, C-style array and Vector.
  <br>
  (Daniel Arndt, 2018/01/25)
 </li>

 <li>
  New: Tensor::begin_raw, Tensor::end_raw, SymmetricTensor::begin_raw
  and SymmetricTensor::end_raw provide access to the underlying storage for Tensor
  and SymmetricTensor.
  <br>
  (Daniel Arndt, 2018/01/25)
 </li>

 <li>
  New: Add save/load functions for ScaLAPACKMatrix to save/load distributed matrix to/from disc using HDF5. If HDF is configured with MPI, parallel I/O is used to save/load the matrix.
  <br>
  (Benjamin Brands, 2018/01/25)
 </li>

 <li>
  New: Added Utilities::MPI::gather wrapper and tests to gather objects from all to one MPI process.
  <br>
  (Benjamin Brands, 2018/01/25)
 </li>

 <li>
  Improved: The implementation of the SolverQMRS iteration has been renewed to
  perform much faster and to be adaptable to left and right side preconditioning
  of the system matrix.
  <br>
  (Ingo Kligge, 2018/01/23)
 </li>

 <li>
  New: New function GridTools::distributed_compute_point_locations ; similarly to GridTools::compute_point_locations , given
  a vector of points, it returns vectors containing them, their reference position and the process owning them as it works
  with shared and distributed meshes.
  <br>
  (Giovanni Alzetta, 2018/01/19)
 </li>

 <li>
  Fixed: In some situations,
  DoFTools::locally_owned_dofs_per_subdomain() returned a vector of the
  wrong size. Specifically, this happened if the processor with the
  highest subdomain id owned no degrees of freedom. This is now fixed.
  <br>
  (Jean-Paul Pelteret, Wolfgang Bangerth, 2018/01/18)
 </li>

 <li>
  Fixed: In FE_Enriched element, avoid evaluating quadrature points if no dofs are
  assigned. This happens when FE_Nothing is used together with other FE
  (i.e. FE_Q) as enrichments in FE_Enriched constructor.
  <br>
  (Nivesh Dommaraju, 2018/01/15)
 </li>

 <li>
  New: Enhanced Raviart-Thomas finite element FE_RT_Bubbles,
  allows for local elimination of a vector variable in
  multipoint flux mixed finite element methods and similar.
  <br>
  (Eldar Khattatov, 2018/01/11)
 </li>

 <li>
  New: Enhanced Raviart-Thomas polynomial space PolynomialsRT_Bubbles
  consisting of RT_k + curls bubbles.
  <br>
  (Eldar Khattatov, 2018/01/09)
 </li>

 <li>
  Fixed: The DataOut::write_vtu output pretended to write double precision data,
  while in reality the data was first converted to float and then written as
  double. This was fixed by writing all data (including vertex positions) as
  single precision float and adjusting the output types accordingly.
  <br>
  (Rene Gassmoeller, 2018/01/04)
 </li>

 <li>
  New: Implementing CylindricalManifold::push_forward_gradient allows to compute
  normal vectors for boundaries described by CylindricalManifold objects.
  <br>
  (Daniel Arndt, 2018/01/04)
 </li>

 <li>
  Fixed: In parallel computations, the DoFRenumbering::hierarchical()
  function created DoF indices that were dependent on the previous DoF
  indices owned by each processor. This was not intended: the new DoF
  indices were supposed to only depend on the order of cells, not any
  previous numbering. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/12/31)
 </li>

 <li>
  New: There is now a function VectorizedArray::streaming_store() that
  stores the content of a vectorized array into another array which bypasses the
  CPU's caches on supported architectures using _mm_stream_pd() intrinsics. This
  operation is useful because it can avoid the read-for-ownership memory access
  on large arrays and directly write to the destination for streaming
  stores. However, it can also be detrimental to performance in case the store
  would fit into caches. This new function is not used within the library
  because only the final user code knows the properties of hardware and whether
  a particular store pattern is so large that it will not fit into caches and
  thus benefit from this optimization.
  <br>
  (Martin Kronbichler, 2017/12/27)
 </li>

 <li>
  Fixed: implement missing instantiation of TrilinosWrappers::SparsityPattern::copy_from().
  <br>
  (Timo Heister, 2017/12/22)
 </li>

 <li>
  Improved: SparseMatrix::mmult(..) is now using the
  DynamicSparsityPattern::compute_mmult_pattern(..) function to create
  the sparsity pattern of the final matrix C.
  <br> (Christoph Goering, 2017/12/22)
 </li>

 <li>
  New: Added QSimplex, QTrianglePolar, QDuffy, and QSplit classes to perform quadratures
  on reference simplices, on their affine transformations, and on hyper cubes
  split by a given point.
  <br>
  (Luca Heltai, 2017/12/21)
 </li>

 <li>
  New: add LAPACK_WITH_64BIT_BLAS_INDICES configure option and introduce
  types::blas_int to support 64bit BLAS indices.
  <br>
  (Denis Davydov, 2017/12/21)
 </li>

 <li>
  New: add LinearAlgebra::distributed::BlockVector::norm_sqr() and
  LinearAlgebra::distributed::Vector::norm_sqr() that return square of the l2 norm.
  <br>
  (Denis Davydov, 2017/12/21)
 </li>

 <li>
  New: LAPACKFullMatrix::remove_row_and_column() to remove certain row and column
  of the matrix.
  <br>
  (Denis Davydov, 2017/12/20)
 </li>

 <li>
  Fixed: parallel::distributed::Triangulation::add_periodicity missed to
  update the ghost_owners in the NumberCache member variable.
  <br>
  (Daniel Arndt, Sambit Das, 2017/12/11)
 </li>

 <li>
  New: DynamicSparsityPattern::compute_mmult_pattern(left, right) with two arguments of either
  a DynamicSparsityPattern or a SparsityPattern; or any combination of those.
  The result is the pattern which is obtained by multiplying the two sparse matrices on
  the given sparsity patterns.
  <br>
  (Christoph Goering, 2017/12/19)
 </li>

 <li>
  New: LAPACKFullMatrix::rank1_update(const number, const Vector<number> &) performs
  a rank-1 update of a matrix.
  <br>
  (Denis Davydov, 2017/12/17)
 </li>

 <li>
  New: Add move constructor for TrilinosWrappers::SparseMatrix and
  TrilinosWrappers::SparsityPattern
  <br>
  (Bruno Turcksin, 2017/12/10)
 </li>

 <li>
  New: Add Vector::grow_or_shrink() and LAPACKFullMatrix::grow_or_shrink()
  to (partly) keep the previous values upon resizing.
  <br>
  (Denis Davydov, 2017/12/16)
 </li>

 <li>
  Fixed: parallel::distributed::Triangulation::communicate_locally_moved_vertices
  treats periodic faces correctly now.
  <br>
  (Daniel Arndt, Sambit Das, 2017/12/11)
 </li>

 <li>
  New: Add Jupyter Notebook explaining how to use the python wrappers.
  <br>
  (Bruno Turcksin, 2017/12/10)
 </li>

 <li>
  Improved: DoFTools::extract_hanging_node_dofs works for
  parallel::shared::Triangualtion and parallel::distributed::Triangulation
  and reports locally relevant DoFs.
  <br>
  (Daniel Arndt, 2017/12/21)
 </li>

 <li>
  New: Now opencascade also works for spacedim == 2.
  <br>
  (Luca Heltai, 2017/12/08)
 </li>

 <li>
  New: Improve support of triangular matrices in LAPACKFullMatrix wrappers.
  Add LAPACKFullMatrix::solve() to solve a system of equations either when
  the matrix is factorized (Cholesky/LU) or triangular.
  <br>
  (Denis Davydov, 2017/12/07)
 </li>

 <li>
  New: Created a new method OpenCASCADE::Utilities::create_curves_from_triangulation_boundary
  that smoothly interpolates the boundary of two dimensional triangulations into a vector of
  OpenCASCADE TopoDS_Edge closed curves, representing the connected components of the boundary.
  <br>
  (Dirk Peschka, Luca Heltai, 2017/12/05)
 </li>

 <li>
  Fixed: SparsityPattern::print_svg(...) has the right size for the
  white rectangle so there is a small margin between the first and last
  row/column of the printed matrix.
  <br>
  (Christoph Goering, 2017/12/05)
 </li>

 <li>
  New: Add color_sparsity_pattern function to color a graph represented by SparsityPattern Object.
  The function uses coloring function from ZOLTAN library.
  <br>
  (Nivesh Dommaraju, 2017/12/04)
 </li>

 <li>
  New: GridGenerator::extract_boundary_mesh now correctly extracts also Manifold information.
  <br>
  (Luca Heltai, Dirk Peschka, 2017/12/04)
 </li>

 <li>
  Fixed: Manifold::normal_vector() used to pick two directions that are almost
  linearly dependent quite often, which resulted in very bad accuracy of the
  normal vectors. This has been fixed. The downstream use of this function on
  curved geometries in VectorTools::compute_no_normal_flux_constraints() has now
  become much more accurate and robust.
  <br>
  (Martin Kronbichler, 2017/12/02)
 </li>

 <li>
  Improved: SphericalManifold::get_new_points now computes a lot
  of information outside of the loop over all points in get_new_point,
  which saves a significant amount of time for spherical geometries.
  <br>
  (Rene Gassmoeller, 2017/12/01)
 </li>

 <li>
  Fixed: MappingQGeneric::transform_real_to_unit_cell() with degree 1 would
  previously fail for Cartesian meshes with certain combinations of roundoff
  errors. This is now fixed.
  <br>
  (Martin Kronbichler, 2017/12/01)
 </li>

 <li>
  New: Added CellId::serialize function and a test for it using pack/unpack
  <br>
  (Giovanni Alzetta, 2017/11/28)
 </li>

 <li>
  New: Added BoundingBox::serialize function and a test for it using pack/unpack functions
  <br>
  (Giovanni Alzetta, 2017/11/24)
 </li>

 <li>
  Extend: partition_triangulation with options for choosing between METIS or ZOLTAN partitioners.
  METIS (the original partitioner) is the default.
  <br>
  (Nivesh Dommaraju, 2017/11/23)
 </li>

 <li>
  New: Added random_value<T>() and random_point<dim>() in tests.h, to simplify and unify our way to generate random numbers and points.
  <br>
  (Luca Heltai, 2017/11/21)
 </li>

 <li>
  Fixed: Added a tolerance inside GridTools::find_active_cell_around_point to errors with boundary points and a test with such points.
  <br>
  (Giovanni Alzetta, 2017/11/21)
 </li>

 <li>
  Fixed: Added convert_generalized_support_point_values_to_dof_values into
  FE_FaceQ and FE_TraceQ. This is needed by VectorTools::interpolate
  <br>
  (Praveen Chandrashekar, 2017/11/19)
 </li>

 <li>
  Fixed: ConvergenceTable::evaluate_convergence_rates now also supports values
  and keys of type 'unsigned long long int' as used for types::global_dof_index
  with 64 bit integers enabled.
  <br>
  (Martin Kronbichler, 2017/11/19)
 </li>

 <li>
  New: Add Utilities::pack and Utilities::unpack, to serialize and unserialize arbitrary
  objects that support boost::serialization/deserialization.
  <br>
  (Luca Heltai, 2017/11/18)
 </li>

 <li>
  New: Utilities::MPI::all_gather and Utilities::MPI::some_to_some functions have been
  added to perform general collective communications between processors.
  <br>
  (Giovanni Alzetta, Luca Heltai, 2017/11/18)
 </li>

 <li>
  Fixed: The function GridRefinement::refine_and_coarsen_fixed_number
  no longer produces an EXC_BAD_ACCESS exception when either of the
  arguments top_fraction or bottom_fraction are set to 1.0.
  <br>
  (Oliver Sutton, 2017/11/17)
 </li>

 <li>
  Changed: The specializations of the scalar_product function for Sacado::Fad::DFad
  numbers were redundant and have been removed.
  <br>
  (Jean-Paul Pelteret, 2017/11/17)
 </li>

 <li>
  New: ArrayView has a data() member function making it more similar to the
  standard containers.
  <br>
  (Daniel Arndt, 2017/11/17)
 </li>

 <li>
  New: Math operations for some automatically differentiable numbers have been
  extended. We also import some math functions for these number types into the
  standard namespace. This gives full compatibility with the Tensor and
  SymmetricTensor classes, and facilitates the integration of these number types
  into existing code.
  <br>
  (Jean-Paul Pelteret, 2017/11/16)
 </li>

 <li>
  Changed: Core functions in the FEValues and FEValuesViews classes have been
  updated to allow the use of certain automatically differentiable numbers to
  represent degree-of-freedom values.
  <br>
  (Jean-Paul Pelteret, 2017/11/16)
 </li>

 <li>
  New: In parallel, the starting indices of
  DoFRenumbering::Cuthill_McKee() can now also be locally active
  indices. They no longer need to be locally owned DoF indices.
  <br>
  (Wolfgang Bangerth, 2017/11/16)
 </li>

 <li>
  Fixed and improved: The GridTools::compute_point_locations() now always return the correct number of cells. The new algorithm is also significantly faster.
  <br>
  (Giovanni Alzetta, 2017/11/15)
 </li>

 <li>
  New: Given a vector of points, the GridTools::guess_point_owner() function uses bounding boxes describing the mesh to guess the processes that own these points. A test for the function has been added.
  <br>
  (Giovanni Alzetta, 2017/11/13)
 </li>

 <li>
  Fixed: The VectorTools::interpolate_boundary_values() function that
  takes a hp::MappingCollection argument was declared, but not
  instantiated, and could consequently not be called. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/11/12)
 </li>

 <li>
  New: TrilinosWrappers::DirectSolver::solve() with two
  LinearAlgebra::distributed::Vector arguments now re-uses the factorization of
  a matrix set in the TrilinosWrappers::DirectSolver::initialize() call.
  <br>
  (Martin Kronbichler, 2017/11/10)
 </li>

 <li>
  Modified: Function Functions::FEFieldFunction::compute_point_locations now calls GridTools::compute_point_locations
  Modified: Function GridTools::compute_point_locations now accepts a cell hint (a test for this hint has been added)
  <br>
  (Giovanni Alzetta, 2017/11/06)
 </li>

 <li>
  Fixed: The algorithm in Manifold::normal_vector could
  cause division-by-zero errors for special cases. This
  was fixed.
  <br>
  (Rene Gassmoeller, 2017/11/03)
 </li>

 <li>
  New: A uniform interface to automatically differentiable numbers that will
  be supported by the deal.II library has been implemented. A corresponding
  set of type traits classes that help query information about these numbers
  has also been implemented.
  <br>
  (Jean-Paul Pelteret, 2017/11/03)
 </li>

 <li>
  Improved: Both Sundials-2.7.0 and Sundials-3.0.0 are supported.
  <br>
  (Daniel Arndt, 2017/11/03)
 </li>

 <li>
  Fixed: ChartManifold now implements the method
  Manifold::get_intermediate_point().
  <br>
  (Martin Kronbichler, 2017/11/02)
 </li>

 <li>
  Improved: SphericalManifold::get_new_point uses a Newton-type
  optimization algorithm instead of a simple averaging.
  <br>
  (Daniel Arndt, 2017/11/02)
 </li>

 <li>
  Fixed: The spherical manifold used to do a division
  by zero inside of SphericalManifold::get_new_point()
  if the new point was exactly at the center of the
  manifold. This is fixed now.
  <br>
  (Rene Gassmoeller, 2017/10/31)
 </li>

 <li>
  New: Added GridTools::exchange_local_bounding_boxes function and tests to exchange
  vectors of bounding boxes and gather them in every process.
  <br>
  (Giovanni Alzetta, 2017/10/30)
 </li>

 <li>
  Fixed: The treatment of inhomogeneous boundary conditions in VectorTools::project
  has been fixed.
  <br>
  (Daniel Arndt, 2017/10/30)
 </li>

 <li>
  Improved: The step-37 tutorial program now describes how to deal with
  inhomogeneous Dirichlet boundary values with the matrix-free framework in the
  results section of the tutorial.
  <br>
  (Martin Kronbichler, 2017/10/29)
 </li>

 <li>
  Fixed: The member variable FEEvaluation::dofs_per_cell and the static variable
  FEEvaluation::tensor_dofs_per_cell only returned the number of degrees of
  freedom per scalar component also for vector-valued evaluators, rather than
  returning the number accumulated over all components in the evaluator. This is
  now fixed.
  <br>
  (Martin Kronbichler, 2017/10/28)
 </li>

 <li>
  New: Created a new Patterns::Tuple class, that allows parsing arbitrary tuples
  of objects in parameter files.
  <br>
  (Luca Heltai, 2017/10/26)
 </li>

 <li>
  New: Removed a typo in QTelles and added a new singular integration test.
  <br>
  (Nicola Giuliani, 25/10/2017
 </li>

 <li>
  New: Add LAPACKFullMatrix<number>::add (const number a, const LAPACKFullMatrix<number> &A)
  <br>
  (Denis Davydov, 2017/10/25)
 </li>

 <li>
  Fixed: The computation of the inverse jacobians in the
  MappingCartesian class was broken. It would trigger an
  assertion at runtime complaining that only an assignment
  with zero is allowed. This is fixed now.
  <br>
  (Rene Gassmoeller, 2017/10/24)
 </li>

 <li>
  New: The function
  LinearAlgebra::distributed::Vector::copy_locally_owned_data_from() can copy
  the locally owned range in a parallel vector without considering the ghost
  entries which might be different between the two vectors.
  <br>
  (Martin Kronbichler, 2017/10/24)
 </li>

 <li>
  Added: Function GridTools::compute_point_locations and a test for it. The function uses the recently added Cache to improve performance.
  <br>
  (Giovanni Alzetta, 2017/10/24)
 </li>

 <li>
  Added: Test for Functions::FEFieldFunction::compute_point_locations
  <br>
  (Giovanni Alzetta, 2017/10/23)
 </li>

 <li>
  New: DoFHandler::renumber_dofs() for multigrid levels is now implemented for
  parallel::distributed::Triangulation for the special case that the numbering
  does not change over processor boundaries.
  <br>
  (Martin Kronbichler, 2017/10/20)
 </li>

 <li>
  New: Added BoundingBox::get_neighbor_type() which returns the neighboring relation between two bounding boxes. Added the function GridTools::compute_mesh_predicate_bounding_box() which returns a vector of BoundingBoxes covering the cells for which a given predicate is true.
  <br>
  (Giovanni Alzetta, 2017/10/20)
 </li>

 <li>
  Changed: The <code>enum</code> VectorOperation is now declared in its own header instead of the header for Vector.
  <br>
  (David Wells, 2017/10/20)
 </li>

 <li>
  New: MatrixFreeOperators::Base now also stores the diagonal of the operator that
  shall be populated in derived class in compute_diagonal() method. The read access
  is provided by MatrixFreeOperators::Base::get_matrix_diagonal().
  <br>
  (Denis Davydov, 2017/10/19)
 </li>

 <li>
  New: A new version of GridTools::find_active_cell_around_point()
  has been added that exploits a GridTools::Cache object.
  <br>
  (Luca Heltai, 2017/10/18)
 </li>

 <li>
  New: The class Utilities::MPI::Partitioner can now handle the MPI import and
  export calls that were previously hidden in the vector class. This allows to
  use MPI communication of integers besides the previously supported floating
  point data.
  <br>
  (Martin Kronbichler, 2017/10/17)
 </li>

 <li>
  New: The class Utilities::MPI::Partitioner now supports handling only a
  selection of ghost indices within a larger ghost set. Use this feature by
  adding a second, larger, IndexSet to Partitioner::set_ghost_indices().
  <br>
  (Martin Kronbichler, 2017/10/17)
 </li>

 <li>
  Fixed: The transformation operation for rank-2 and rank-4 (non-symmetric) tensors
  in the  Physics::Transformations was incorrectly implemented.
  This has now been fixed.
  <br>
  (Jean-Paul Pelteret, 2017/10/17)
 </li>

 <li>
  New: Add LAPACKFullMatrix::trace() to calculate trace of a matrix.
  <br>
  (Denis Davydov, 2017/10/09)
 </li>

 <li>
  New: FiniteElement::get_sub_fe() allows you to return a contained FiniteElement based on a ComponentMask.
  <br>
  (Timo Heister, 2017/10/07)
 </li>

 <li>
  Fixed: It was previously not possible to create objects of type Table
  with elements that are not copyable or copy constructible. An example
  of such a type is `std::unique_ptr`. This has now been fixed.
  <br>
  (Wolfgang Bangerth, 2017/10/02)
 </li>

 <li>
  New: MGTransferBlockMatrixFree allows to use a separate DoFHandler
  for each block. The required interface is supported by PreconditionMG.
  <br>
  (Daniel Arndt, 2017/10/01)
 </li>

 <li>
  New: Add Utilities::MPI::sum() for LAPACKFullMatrix objects.
  <br>
  (Denis Davydov, 2017/09/04)
 </li>

 <li>
  New: Add LinearAlgebra::distributed::BlockVector<Number>::multivector_inner_product()
  to perform inner product between each block of the two block vectors.
  <br>
  (Denis Davydov, 2017/09/30)
 </li>

 <li>
  Fixed: One of the functions GridTools::find_active_cell_around_point()
  supports an additional `mapping` argument. This mapping was used
  to actually transform points from the real to the reference space,
  but somewhere in the algorithm we look for the `closest_vertex` to
  the target point. Such search was done ignoring if the `mapping`
  actually changes the location of the vertices. This is now fixed,
  and the internal function that looks for the closest vertex now
  takes into account the `mapping` transformation.
  <br>
  (Luca Heltai, 2017/09/27)
 </li>

 <li>
  Changed: Split TensorProductMatrixSymmetricSum into base and derived class.
  Added a template specialization for TensorProductMatrixSymmetricSum
  in case of having a VectorizedArray as arithmetic type.
  In addition, TensorProductMatrixSymmetricSum is now able to handle
  distinct 1D matrices in each tensor dimension.
  <br>
  (Julius Witte, 2017/09/27)
 </li>

 <li>
  New: Add MGTransferMatrixFree::interpolate_to_mg() to transfer fine-level
  solution to a multi-grid vector
  <br>
  (Denis Davydov, 2017/09/27)
 </li>

 <li>
  New: Extend MappingQEulerian to work with multi-grid vectors.
  <br>
  (Denis Davydov, 2017/09/27)
 </li>

 <li>
  Fixed: The hp version of DoFTools::make_flux_sparsity_pattern() with masks for
  DoF couplings now works with non-primitive FE spaces, similarly to its non-hp version.
  The hp version of VectorTools::project_boundary_values() is modified to work properly with Hdiv
  conforming finite elements.
  <br>
  (Eldar Khattatov, 2017/09/26)
 </li>

 <li>
  New: Extend LAPACKFullMatrix to provide (i) Cholesky factorization (and use it for inversion);
  (ii) an estimate of the reciprocal condition number; (iii) l1_norm(),
  linfty_norm() and frobenius_norm().
  <br>
  (Denis Davydov, 2017/09/26)
 </li>

 <li>
  New: There are new signals in the Triangulation class
  that signal the beginning and end of a parallel::distributed
  refinement and serialization. This allows to connect functions
  that should be called uniquely once before and after
  refinement and serialization for parallel distributed
  Triangulations (something not possible with the existing
  signals).
  <br>
  (Rene Gassmoeller, 2017/09/25)
 </li>

 <li>
  Changed: The new methods ArrayView::cbegin() and ArrayView::cend() return
  const_iterators and are const themselves. ArrayView::begin() and
  ArrayView::end() return (non-const) iterators and are const themselves.
  <br>
  (Daniel Arndt, Julius Witte, 2017/09/24)
 </li>

 <li>
  New: Added two methods to the class BoundingBox. The first, merge_with(BoundingBox) merges the current BoundingBox with BBox. The second, volume(), returns the volume (dim-dimensional measure) of the current BoundingBox
  <br>
  (Giovanni Alzetta, 2017/09/22)
 </li>

 <li>
  New: A GridTools::find_closest_vertex() has been added
  that takes a map of index-vertex position and a target
  point, and returns the index of the closest point in
  the map.
  <br>
  (Luca Heltai, 2017/09/19)
 </li>

 <li>
  New: A GridTools::extract_used_vertices() has been added
  that returns a map of index-vertex position, using an
  optional Mapping<dim,spacedim> object to compute the
  vertex locations.
  <br>
  (Luca Heltai, 2017/09/18)
 </li>

 <li>
  New: A new GridIn::read_assimp() function allows to
  read in a file in one of the (40) formats supported by Assimp.
  <br>
  (Luca Heltai, 2017/09/16)
 </li>

 <li>
  New: Added the base class BoundingBox. The BoundingBox is constructed from a std::pair of points in real space, ordered following the convention bottom-left, top-right.
  The class has the method BoundingBox::point_inside(Point<spacedim> p) which checks for the presence of a point inside it.
  <br>
  (Giovanni Alzetta, 2017/09/15)
 </li>

 <li>
  Fixed: The destructor of TimerOutput::Scope will now always exit the subsection
  it entered in its constructor, instead of just exiting the last subsection that
  was created in the provided TimerOutput object. This fixes timing output
  measurements for multithreaded applications and cases where subsections are
  nested.
  <br>
  (David Wells, 2017/09/14)
 </li>

 <li>
  New: FiniteElement and FESystem have move constructors.
  <br>
  (Daniel Arndt, 2017/09/11)
 </li>

 <li>
  New: There is now a conversion constructor that converts std::vector
  objects to ArrayView objects.
  <br>
  (Wolfgang Bangerth, 2017/09/11)
 </li>

 <li>
  New: Add more Python wrappers to generate meshes and to flatten a Triangulation.
  <br>
  (Bruno Turcksin, 2017/09/10)
 </li>

 <li>
  New: The newly implemented FiniteElement::operator^ allows in combination
  with the new variadic template constructor or the new std::initializer_list
  constructor for FESystem to construct
  FESystem objects using syntax like
  @code
    FESystem<dim, spacedim> fe_system1 ({FiniteElementType1^n1,
                                         FiniteElementType2^n2});
    FESystem<dim, spacedim> fe_system2 (FiniteElementType1^n3,
                                        FiniteElementType2^n4);
  @endcode
  <br>
  (Daniel Arndt, 2017/09/05)
 </li>

 <li>
  New: Add Utilities::MPI::sum() for FullMatrix objects.
  <br>
  (Denis Davydov, 2017/09/04)
 </li>

 <li>
  Changed: The various overloads for Threads::new_thread() and
  Threads::new_task() taking from 0 to 9 parameters have been replaced
  by two variadic template versions.
  <br>
  (Daniel Arndt, 2017/09/01)
 </li>

 <li>
  New: The new DataPostprocessorTensor can help visualize tensor-valued
  quantities.
  <br>
  (Wolfgang Bangerth, 2017/08/30)
 </li>

 <li>
  New: The DataPostprocessorVector class now has an extensive example of
  its use.
  <br>
  (Wolfgang Bangerth, 2017/08/29)
 </li>

 <li>
  Changed: DoFHandler::get_finite_element(unsigned int) has been removed again in
  favor of equipping DoFHandler::get_fe() with a (defaulted) unsigned int
  parameter. hp::DoFHandler::get_finite_element(unsigned int)  has been renamed to
  hp::DoFHandler::get_fe(unsigned int).
  <br>
  (Daniel Arndt, 2017/08/25)
 </li>

 <li>
  New: Extend Point::distance_square() to work with VectorizedArray numbers.
  <br>
  (Denis Davydov, 2017/08/24)
 </li>

 <li>
  Fixed: On Vector::add(..), a pointer to the first element of a std::vector
  is gotten by std::vector::data() instead of using the "&v[0]" idiom, which
  results in an undefined behaviour.
  <br>
  (Tulio Ligneul, 2017/08/23)
 </li>

 <li>
  Fixed: It turns out that it was not possible to copy an invalid
  operator, i.e., the following code would yield an exception:
  @code
    typename DoFHandler<dim>::active_cell_iterator invalid_1;
    typename DoFHandler<dim>::active_cell_iterator invalid_2;
    invalid_1 = invalid_2;  // resulted in an error
  @endcode
  This made no sense, and has now been fixed.
  <br>
  (Wolfgang Bangerth, 2017/08/23)
 </li>

 <li>
  New: The struct is_base_of_all is a generalization of std::is_base_of
  to template parameter packs and tests if all classes in the parameter pack
  have a given class as base class or are this class.
  <br>
  (Daniel Arndt, 2017/08/21)
 </li>

 <li>
  New: The constructors for hp::FECollection taking a fixed number of
  FiniteElement objects as argument were replaced by a variadic
  template constructor taking an arbitrary number of objects.
  <br> (Daniel Arndt, 2017/08/20)
 </li>

 <li>
  New: DoFHandler::get_fe() and hp::DoFHandler::get_fe()
  have been deprecated in favor of DoFHandler::get_finite_element(),
  DoFHandler::get_fe_collection(), hp::DoFHandler::get_finite_element()
  and hp::DoFHandler::get_fe_collection() that share a common signature each.
  <br>
  (Daniel Arndt, 2017/08/19)
 </li>

 <li>
  New: Add MGTransferBlockMatrixFree which allows using geometric multi-grids with
  LinearAlgebra::distributed::BlockVector<Number> where each block uses the same
  index space of degrees of freedom.
  <br>
  (Denis Davydov, 2017/08/18)
 </li>

 <li>
  Fixed: The non-member invert() function operating on SymmetricTensors has
  been generalized for all number types.
  <br>
  (Wolfgang Bangerth, Jean-Paul Pelteret, 2017/08/17)
 </li>

 <li>
  New: An additional non-member operator*() function, performing a single
  contraction between a symmetric and normal tensor of arbitrary rank, has been
  added to the SymmetricTensor class.
  <br>
  (Matthias Maier, Jean-Paul Pelteret, 2017/08/17)
 </li>

 <li>
  Fixed: The following non-member functions operating on SymmetricTensors have
  been generalized for all number types:
    - scalar_product()
    - double_contract()
    - invert() for rank-4 symmetric tensors
    - operator*() performing a double contraction between two symmetric tensors.
    - operator*() performing a single contraction of a symmetric tensor with a
      rank-1 tensor.
  <br>
  (Jean-Paul Pelteret, 2017/08/17)
 </li>

 <li>
  Fixed: The Tensor::l1_norm() and Tensor::linfty_norm() functions have now been
  generalized for all number types.
  <br>
  (Jean-Paul Pelteret, 2017/08/17)
 </li>

 <li>
  Improved: The CMake CUDA detection and configure has been completely
  rewritten. CUDA support now requires CMake version 3.9 or newer.
  <br>
  (Matthias Maier, 2017/08/08)
 </li>

 <li>
  Fixed: ConsecutiveControl now asserts that check() is being called for zeroth
  iteration.
  <br>
  (Denis Davydov, 2017/08/16)
 </li>

 <li>
  Fixed: Constructing an IndexSet from an Epetra_Map generated the wrong size.
  <br>
  (Timo Heister, 2017/08/15)
 </li>

 <li>
  Deprecated: types_are_equal has been deprecated in favor of std::is_same.
  <br> (David Wells, 2017/08/14)
 </li>

 <li>
  Deprecated: constraint_and_return_value has been deprecated
  in favor of std::enable_if.
  <br> (Daniel Arndt, 2017/08/14)
 </li>

 <li>
  New: A new version of GridTools::find_active_cell_around_point()
  has been added that exploits local maps constructed
  using standard GridTools utilities.
  <br>
  (Rene Gassmoeller, Luca Heltai, 2017/08/13)
 </li>

 <li>
  New: GridTools::exchange_cell_data_to_ghosts() to exchange ghost data.
  <br>
  (Timo Heister, 2017/08/12)
 </li>

 <li>
  New: All of the FEValuesViews classes now have a set of functions
  named get_function_*_from_local_dof_values that return the function
  values, gradients etc. from a local vector of degree-of-freedom
  values. They also now contain a data structure called OutputType
  that that provides the output type for the product of the value
  and derivatives of basis functions of the view and any number type.
  <br>
  (Jean-Paul Pelteret, Luca Heltai, 2017/08/10)
 </li>

 <li>
  Deprecated: internal::bool2type and int2type has been deprecated
  in favor of std::intergral_constant.
  <br>
  (Daniel Arndt, 2017/08/10)
 </li>

 <li>
  Fixed: Creation of coverage information via
  compiling with 'DEAL_II_SETUP_COVERAGE=ON'
  and 'ctest -S ../tests/run_coverage.cmake' works again.
  <br>
  (Matthias Maier, Daniel Arndt, 2017/08/09)
 </li>

 <li>
  Improved: ConstraintMatrix now has a member function
  ConstraintMatrix::get_lines() that returns a range object that allows
  iteration over all constraint entries.
  <br>
  (Matthias Maier, 2017/08/08)
 </li>

 <li>
  Changed: ConstraintMatrix::ConstraintLine is now a public member class
  definition.
  <br>
  (Matthias Maier, 2017/08/08)
 </li>

 <li>
  New: MappingQGeneric uses tensorized evaluation for maybe_update_Jacobians,
  maybe_compute_q_points and maybe_update_jacobian_grads in case
  the Quadrature object represents a tensor product of identical
  one-dimensional quadrature formulas.
  <br>
  (Daniel Arndt, 2017/08/04)
 </li>

 <li>
  New: Quadrature has a (defaulted) move assignment operator.
  <br>
  (Daniel Arndt, 2017/08/04)
 </li>

 <li>
  New: Add Utilities::LinearAlgebra::chebyshev_filter() to apply Chebyshev filter of a
  given degree.
  <br>
  (Denis Davydov, 2017/08/02)
 </li>

 <li>
  New: Add Utilities::LinearAlgebra::lanczos_largest_eigenvalue() to estimate the largest
  eigenvalue of a symmetric linear operator using k-steps of Lanczos algorithm.
  <br>
  (Denis Davydov, 2017/08/01)
 </li>

 <li>
  Fixed: DiagonalMatrix<VectorType> can now be used as a linear operator.
  <br>
  (Denis Davydov, 2017/08/01)
 </li>

 <li>
  New: Extend pArpack solver to cover mode 1 (standard eigenvalue problem) and
  mode 2 (generalized eigenvalue problem without spectral transformation).
  <br>
  (Denis Davydov, 2017/08/01)
 </li>

 <li>
  New: Quadrature::is_tensor_product() returns if the corresponding quadrature
  formula can be represented as a tensor product and the quadrature points are
  sorted lexicographically. The corresponding 1D objects can be queried using
  Quadrature::get_tensor_basis().
  <br>
  (Daniel Arndt, 2017/07/26)
 </li>

 <li>
  Fixed: DataOutBase::write_hdf5_parallel() can also be used in case deal.II
  and HDF5 are not compiled with MPI support.
  <br>
  (Daniel Arndt, 2017/07/26)
 </li>

 <li>
  Fixed: The FunctionParser class demonstrated an incompatibility with
  (very) old versions of the Threading Building Blocks used in
  deal.II. This is now worked around.
  <br>
  (Alberto Sartori, Wolfgang Bangerth, 2017/07/26)
 </li>

 <li>
  Initialization of Chebyshev smoother: Make initial guess robust with
  respect to number of processors by operating on the global index.
  <br>
  (Niklas Fehn, 2017/07/25)
 </li>

 <li>
  Fixed: MappingQGeneric::transform_real_to_unit_cell had a bug where, for certain
  parallelograms, a problem with subtracting nearly equal floating point numbers
  ruined most of the digits of accuracy in the coordinate transformation.
  <br> (Giovanni Di Ilio, David Wells, Jean-Paul Pelteret, 2017/07/25)
 </li>

 <li>
  New: There is now a function ConstraintMatrix::copy_from() that allows
  copying objects of type ConstraintMatrix.
  <br>
  (Wolfgang Bangerth, 2017/07/27)
 </li>

 <li>
  Fixed: The ghost layer computation for the case when periodic boundary
  conditions were set on all boundaries in 3D missed some cells leading to a
  deadlock when used on 43 MPI ranks or more. This is now fixed.
  <br>
  (Martin Kronbichler, 2017/07/20)
 </li>

 <li>
  New: DataOutBase now supports writing zero-dimensional
  data in a higher spatial dimension (e.g. vertices, quadrature points, particles)
  for a number of output formats (vtu,vtk,gnuplot,hdf5). A consequence of the
  change is that HDF5 now in general supports output for dim<spacedim.
  <br>
  (Rene Gassmoeller, 2017/07/18)
 </li>

 <li>
  Improved: linear_operator gained an overload that directly takes a
  LinearOperator as exemplar object to directly copy the
  reinit_(domain|range)_vector functions from.
  <br>
  (Matthias Maier, 2017/07/14)
 </li>

 <li>
  New: Adding the functions FE_FaceQ::hp_vertex_dof_identities(),
  FE_FaceQ::hp_line_dof_identities() and FE_FaceQ::hp_quad_dof_identities()
  allows FE_FaceQ to be used in combination with a hp::DoFHandler.
  <br>
  (Samuel Imfeld, Daniel Arndt, 2017/07/11)
 </li>

 <li>
  New: A new class TensorProductMatrixSymmetricSum that represents a tensor
  product of a mass and derivative matrix has been added. It implements both the
  matrix-vector operation (TensorProductMatrixSymmetricSum::vmult()) as well as
  the inverse operation (TensorProductMatrixSymmetricSum::apply_inverse())
  operation in tensorial form with optimal complexity of
  <tt>size<sup>dim+1</sup></tt> operations. The inverse operation uses the fast
  diagonalization method.
  <br>
  (Martin Kronbichler, 2017/07/11)
 </li>

 <li>
  Fixed: The DoFHandler::n_boundary_dofs() and
  hp::DoFHandler::n_boundary_dofs() functions did not take into account
  that a one-dimensional triangulation may consist of more than one
  segment (and then may have more than two vertices on the
  boundary). This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/07/10)
 </li>

 <li>
  New: MGConstrainedDoFs object now stores information for periodicity
  constraints which is used during the multigrid transfer by
  MGTransferMatrixFree and MGTransferPrebuilt.
  <br>
  (Conrad Clevenger, Martin Kronbichler, 2017/07/08)
 </li>

 <li>
  Changed: The DoFRenumbering::compute_component_wise() function used to
  take two template arguments that denote iterator types, one for the
  begin and one for the end iterator. This has been changed: It now only
  takes an iterator type argument for the begin iterator, and the end
  iterator is automatically casted to that same type.
  <br>
  (Wolfgang Bangerth, 2017/07/08)
 </li>

 <li>
  New: The DoFRenumbering::component_wise() function was available for
  all dimensions and space dimensions for regular DoFHandler objects,
  but only for dimension == space_dimension for hp::DoFHandler. This is
  now fixed.
  <br>
  (Wolfgang Bangerth, 2017/07/07)
 </li>

 <li>
  Fixed: TrilinosWrappers::PreconditionBlock* classes now work correctly if one of the processors does not own any rows.
  <br>
  (Timo Heister, 2017/07/04)
 </li>

 <li>
  Deprecated: The ParpackSolver::Shift class has been deprecated.
  <br>
  (Matthias Maier, 2017/07/03)
 </li>

 <li>
  New: The determinant of a LAPACKMatrix can now be computed via a
  call to LAPACKMatrix::determinant().
  <br>
  (Jean-Paul Pelteret, 2017/07/02)
 </li>

 <li>
  Fixed: FullMatrix::left_invert() and FullMatrix::right_invert() no longer
  fail when the matrix is square.
  <br>
  (Andreas Kergassner, Jean-Paul Pelteret, 2017/07/02)
 </li>

 <li>
  Fixed: FullMatrix::determinant() is now implemented for square matrices
  of size greater than 3. In this case, it is computed through the
  LU decomposition of the matrix using LAPACK.
  <br>
  (Andreas Kergassner, Jean-Paul Pelteret, 2017/07/02)
 </li>

 <li>
  Changed: The complete_index_set() function now returns a compressed
  version of the IndexSet that contains all indices.
  <br>
  (Wolfgang Bangerth, 2017/06/24)
 </li>

 <li>
  New: The Timer class now has new functions Timer::get_total_data(),
  Timer::print_total_data, Timer::last_wall_time, Timer::cpu_time and
  Timer::last_cpu_time() to access information about the total run as
  well as the last lap.
  <br>
  (Daniel Arndt, 2017/06/23)
 </li>

 <li>
  Deprecated: The static member variables DoFHandler::invalid_dof_index
  and hp::DoFHandler::invalid_dof_index are now deprecated. Use
  numbers::invalid_dof_index instead.
  <br>
  (Wolfgang Bangerth, 2017/06/20)
 </li>

 <li>
  New: Patterns::Map now accepts a new optional parameter,
  that allows to specify also the separator between
  key-value pairs.
  <br>
  (Luca Heltai, 2017/06/14)
 </li>

 <li>
  New: Utilities::split_string_list() now has a
  version that allows to pass a string instead of a
  character. With this it is possible to split on
  double chars, for example: "alpha ;; beta ;; gamma".
  <br>
  (Luca Heltai, 2017/06/14)
 </li>

 <li>
  Fixed: The passing of additional settings to
  TrilinosWrappers::SolverBase from derived solvers
  was never completed. This has now been corrected.
  <br>
  (Julian Andrej, Jean-Paul Pelteret, 2017/06/14)
 </li>

 <li>
  Fixed: Patterns derived from PatternBase now return
  std::unique_ptr<PatternBase> in their clone() and
  create() functions.
  <br>
  (Luca Heltai, 2017/06/12)
 </li>

 <li>
  Changed: ParameterHandler::print_parameters() was at times looking at
  individual bits instead of just the declared values
  of ParameterHandler::OutputStyle. This presumably allowed for calling that
  function with a combination of the OutputStyle flags, for reasons that no
  longer seem particularly relevant nor obvious.
  <br>
  This possibility has now been removed from the current implementation of
  the function, but is for the moment retained for the (deprecated) function
  ParameterHandler::print_parameters_section().
  <br>
  (Wolfgang Bangerth, 2017/06/09)
 </li>

 <li>
  Fixed: The ParameterHandler::print_parameters() function
  used the previously set fill character of the output stream, but it
  should have filled with spaces instead. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/06/07)
 </li>

 <li>
  Fixed: FullMatrix::print_formatted() and LAPACKFullMatrix::print_formatted
  would previously print NaN values as zero. This is now fixed.
  <br>
  (Simon Sticko, Martin Kronbichler, 2017/06/08)
 </li>

 <li>
  Initialization of Chebyshev smoother: project vector onto space of vectors
  with zero mean which is necessary in some cases where the matrix is singular.
  Rename variable is_initialized -> eigenvalues_are_initialized
  in order to improve code readability.
  <br>
  (Niklas Fehn, 2017/06/07)
 </li>

 <li>
  Improved: The multigrid smoothers now provide an additional apply() method that
  specializes the smoothing for the case the input vector is set to zero. This
  method is used in the Multigrid class and saves one matrix-vector product,
  which can speed up the multigrid algorithm by up to 10 percent, depending on
  the cost of the smoother.
  <br>
  (Martin Kronbichler, 2017/06/07)
 </li>

 <li>
  Changed: The ParameterHandler::print_parameters() function is now
  @p const, as one would expect it to be.
  <br>
  (Wolfgang Bangerth, 2017/06/07)
 </li>

 <li>
  New: Added a new polynomial class IntegratedLegendreSZ. This implements
  the integrated Legendre polynomials described in the 2006 PhD thesis
  of Sabine Zaglmayr.
  <br>
  (Ross Kynch, 2017/06/06)
 </li>

 <li>
  New: A function TriaAccessor::real_to_unit_cell_affine_approximation that
  computes an approximation to a point in unit coordinates has been added.
  <br>
  (Martin Kronbichler, 2017/06/05)
 </li>

 <li>
  Improved: CylindricalManifold is now really using cylindrical coordinates
  instead of taking an average in space coordinates.
  <br>
  (Daniel Arndt, 2017/06/05)
 </li>

 <li>
  Deprecated: The ParameterHandler::print_parameters_section() function
  has been deprecated.
  <br>
  (Wolfgang Bangerth, 2017/06/02)
 </li>

 <li>
  New: The ConeBoundary class now has a member function ConeBoundary::normal_vector().
  <br> (Anne Glerum, 2017/05/30)
 </li>

 <li>
  Fixed: DataOutInterface::set_flags() was broken when passed an object
  of type DataOutBase::GnuplotFlags. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/05/30)
 </li>

 <li>
  New: deal.IIFeatureConfig.cmake will contain detailed information about
  feature configuration and is installed along with other CMake files.
  These variables can be loaded by calling INCLUDE(\${DEAL_II_FEATURE_CONFIG}).
  This is useful when users want to link in additional parts of dependencies, not
  required by deal.II (i.e. filesystem library from Boost package).
  <br>
  (Denis Davydov and Matthias Maier, 2017/05/29)
 </li>

 <li>
  Removed: The build_tests category of tests has been removed. Build tests
  now compile all configurable example steps as part of the regular build
  stage.
  <br>
  (Matthias Maier, 2017/04/06)
 </li>

 <li>
  New: There is a new header, <tt>lac/trilinos_index_access.h</tt>, that provides
  index and size functions for some common Trilinos objects that work correctly for
  both 32 and 64 bit code.
  <br> (David Wells, 2017/05/10)
 </li>

 <li>
  New: Add method to access timing data from TimerOutput from code.
  <br>
  (Jonathan Robey, 2017/05/22)
 </li>

 <li>
  New: Added a GridTools::regularize_corner_cells function that
  detects if the boundary cells of a mesh at corner positions
  (with dim adjacent faces on the boundary) need to be split into
  cells with smaller angles.
  <br> (Luca Heltai, Martin Kronbichler, 2017/05/18)
 </li>

 <li>
  New: Add a new class ConsecutiveControl which returns SolverControl::State::success
  if and only if a certain positive number of consecutive iterations satisfy the
  specified tolerance.
  <br>
  (Denis Davydov, 2017/05/17)
 </li>

 <li>
  New: The vertex (and face) iterator for dimension 1 triangulations now implements
  the get_manifold method.
  <br> (David Wells, 2017/05/10)
 </li>

 <li>
  New: Add support for matrix-free operator application on GPU using CUDA.
  Constraints (Dirichlet boundary conditions, hanging nodes, etc.) cannot be
  applied.
  <br>
  (Bruno Turcksin and Karl Ljungkvist, 2017/05/10)
 </li>

 <li>
  New: TrilinosWrappers::MPI::BlockMatrix can now return its MPI_Comm via
  get_mpi_communicator() and its range and domain indices via
  locally_owned_range_indices() and locally_owned_domain_indices()
  relying on the information of the TrilinosWrappers::MPI::Matrix
  object it is based on.
  <br>
  (Daniel Arndt, 2017/05/09)
 </li>

 <li>
  New: ConstraintMatrix::distribute_local_to_global can now also assemble to
  rectangular matrices where rows and columns are described by different
  ConstraintMatrix objects.
  <br>
  (Martin Kronbichler, 2017/05/04)
 </li>

 <li>
  Changed: The PETScWrappers::MPI::Vector class always corresponds to a PETSc
  vector with type <code>mpi</code>, even when empty.
  <br> (David Wells, 2017/03/29)
 </li>

 <li>
  New: TriaAccessor::enclosing_ball() computes and return a pair of Point
  and double corresponding to the center and the radius of a reasonably small
  enclosing ball of the TriaAccessor object.
  <br>
  (Vishal Boddu, Denis Davydov 2017/04/26)
 </li>

 <li>
  Fixed: Copying Patterns::List or Patterns::Map objects previously led
  to memory corruption. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2017/04/25)
 </li>

 <li>
  New: Add Point<dim,Number>::distance_square() which calculates the squared
  Euclidean distance.
  <br>
  (Denis Davydov, 2017/04/24)
 </li>

 <li>
  New: The ParameterHandler::parse_input() function and friends now make
  a guarantee that if a parameter value does not match its pattern, or
  an associated action throws an exception, that the enclosing
  ParameterHandler object will be reset to the same subsection level
  it was in before ParameterHandler::parse_input() was called.
  <br>
  (Wolfgang Bangerth, 2017/04/21)
 </li>

 <li>
  Changed: We no longer support Visual Studio 2013 because it
  lacks important c++11 features like constexpr. It is now possible
  to use MSVC 2017 in addition to MSVC 2015.
  <br>
  (Timo Heister, 2017/04/20)
 </li>

 <li>
  New: The new function ParameterHandler::add_action() allows to
  register actions that should be performed when a parameter
  is read from somewhere. This allows, in particular, initializing
  member variables that store parameter values without having to
  explicitly call ParameterHandler::get(),
  ParameterHandler::get_integer(), or similar functions.
  <br>
  (Wolfgang Bangerth, 2017/04/16)
 </li>

 <li>
  Fixed: The MeshWorker framework can also be used in case dim != spacedim.
  <br>
  (Daniel Arndt, 2017/04/16)
 </li>

 <li>
  Changed: The PointValueHistory class still used the old, deprecated
  interface for DataPostprocessor objects in that variant of the
  PointValueHistory::evaluate_field() function that takes such an
  object. This is now fixed: it uses
  DataPostprocessor::evaluate_scalar_field() and
  DataPostprocessor::evaluate_vector_field(), like all other
  users of the DataPostprocessor class.
  <br>
  (Wolfgang Bangerth, 2017/04/16)
 </li>

 <li>
  Fixed: FETools::extrapolate also works for BlockVector types.
  <br>
  (Daniel Arndt, 2017/04/13)
 </li>

 <li>
  Fixed: When initializing a LinearAlgebra::distributed::Vector and omitting
  zeroing the entries, ghost entries were left undefined, leading to use of
  invalid memory when calling compress() without a previous zero_out_ghosts()
  also when the local range was completely valid. This is now fixed.
  <br>
  (Martin Kronbichler, 2017/04/12)
 </li>

 <li>
  New: Augment python interface of Point: add operators to modify point, a
  function to compute the distance between points, and functions to compute norms.
  <br>
  (Bruno Turcksin, 2017/04/11)
 </li>

 <li>
  New: There is now a function TableHandler::declare_entry() that
  creates a column of a table without actually putting a value into it.
  <br>
  (Wolfgang Bangerth, 2017/04/11)
 </li>

 <li>
  New: There is now a function TableHandler::start_new_row() that fills
  all entries of the current row (if any) that hadn't been filled, and
  thus starts a new row of the table.
  interface.
  <br>
  (Wolfgang Bangerth, 2017/04/09)
 </li>

 <li>
  New: The class FESystem now implements
  the
  FiniteElement::convert_generalized_support_point_values_to_dof_values()
  interface.
  <br>
  (Wolfgang Bangerth, Matthias Maier, 2017/08/31)
 </li>

 <li>
  New: The classes FE_Q, FE_DGQ, and FE_DGQArbitraryNodes now implement
  the
  FiniteElement::convert_generalized_support_point_values_to_dof_values()
  interface.
  <br>
  (Wolfgang Bangerth, 2017/04/05)
 </li>

 <li>
  New: The FEEvaluation::read_dof_values and
  FEEvaluation::distribute_local_to_global now use dedicate gather (AVX2,
  AVX-512) and scatter (AVX-512) instructions for faster vector access if those
  are available.
  <br>
  (Martin Kronbichler, 2017/03/24)
 </li>

 <li>
  New: The tensor product kernels used by FEEvaluation have gained a faster
  variant for evaluating gradients derived from spectral identities. This
  improves performance of the Laplace operator evaluation by 10-20% in the
  computation-bound case.
  <br>
  (Martin Kronbichler, 2017/03/24)
 </li>

 <li>
  New: DoFTools::extract_dofs_with_support_contained_within() returns a set of
  degrees of freedom whose support is entirely contained within the cells for
  which the predicate returns true.
  <br>
  (Denis Davydov, 2017/03/24)
 </li>

 <li>
  New: GridTools::compute_bounding_box() computes a bounding box of a subdomain
  whose active cells conform to a given predicate.
  GridTools::compute_active_cell_layer_within_distance() computes a collection
  of active cells that are within a given distance from
  predicate subdomain.
  GridTools::compute_ghost_cell_layer_within_distance() computes a collection
  of active cells that are within a given distance from
  locally owned active cells.
  <br>
  (Vishal Boddu, Denis Davydov 2017/03/23)
 </li>

 <li>
  Extend: The two overloaded functions GridTools::find_active_cell_around_point() now take an optional custom mask for vertices
  to narrow down the search for surrounding cells.
  <br>
  (Vishal Boddu, 2017/03/22)
 </li>

</ol>

*/
