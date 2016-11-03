// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2016 by the deal.II authors
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
@page changes_after_8_4_2 Changes after Version 8.4.2

<p>
This is the list of changes made after the release of deal.II version
8.4.1. All entries are signed with the names of the authors.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

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
 <li> Changed: Some typedefs (both names and definitions) within the
 <tt>generic_linear_algebra.h</tt> header have been altered.
 <br>
 (Jean-Paul Pelteret, 2016/10/21)
 </li>

 <li> Deprecated: ParameterHandler::read_input,
 ParameterHandler::read_input_from_xml, and
 ParameterHandler::read_input_from_string are now deprecated in favor of
 ParameterHandler::parse_input, ParameterHandler::parse_input_from_xml, and
 ParameterHandler::parse_input_from_string. These new functions throw
 exceptions to indicate failure instead of using return codes.
 <br>
 (David Wells, 2016/09/15)
 </li>

 <li> Deprecated: MGCoarseGridLACIteration got deprecated in favor of
 MGCoarseGridIterativeSolver.
 <br>
 (Timo Heister, 2016/09/14)
 </li>

 <li> Changed: The template parameter order in many VectorTools functions is now
 different; this was done so that the order is the same across similar functions.
 This will only effect code that explicitly specifies template parameters for
 overloaded VectorTools functions (no known deal.II-based projects do this).
 <br>
 (David Wells, 2016/09/06)
 </li>

 <li> New: deal.II now requires at least BOOST version 1.56, rather than the
 previous minimal version of 1.54. This is because 1.54 does not support
 serializing objects of type std::unique_ptr if C++11 is used, but we now
 use such objects in a variety of places in classes that can be serialized.
 BOOST 1.56, on the other hand, supports this. deal.II bundles BOOST 1.62
 for cases where no or no sufficiently new version of BOOST is found on
 a system.
 <br>
 (Wolfgang Bangerth, 2016/08/22)
 </li>

 <li> Removed: Deprecated classes CompressedSparsityPattern,
 CompressedSimpleSparsityPattern, CompressedSetSparsityPattern, and their
 block variants got removed.
 <br>
 (Timo Heister, 2016/08/13)
 </li>

 <li> Deprecated: MGLevelObject::clear() deprecated in favor of
 MGLevelObject::clear_elements() due to clear() being inconsistent with
 behavior of other container objects.
 <br>
 (Jonathan Robey, 2016/08/08)
 </li>

 <li> Changed: Several operators from LocalIntegrators::Divergence got moved
 to LocalIntegrators::GradDiv and the never used/tested
 LocalIntegrators::Divergence::grad_div() function was removed.
 <br>
 (Timo Heister, Guido Kanschat, 2016/08/02)
 </li>

 <li> Changed: DoFTools::make_cell_patches() only accepts block lists of type
 SparsityPattern. The reason is that it has to initialize the size of the
 pattern on distributed triangulations by computing the number of locally
 owned cells. Initialization differs between sparsity pattern classes, so no
 generic function would be possible. On the other hand, the block list is an
 object, which only extends over locally owned grid cells and its size can be
 determined efficiently upon initialization. Therefore, SparsityPattern is a
 good choice here.
 <br>
 At the same time, we changed the dof handler template to the type DoFHandler,
 since hp::DoFHandler requires a different setup of the SparsityPattern.
 <br>
 (Guido Kanschat, 2016/08/02)
 </li>

 <li> Changed: The conversion constructors of class Vector from the
 PETScWrappers::Vector, PETScWrappers::MPI::Vector,
 TrilinosWrappers::Vector, and TrilinosWrappers::MPI::Vector classes
 are now marked as <code>explicit</code>, i.e., they will no longer
 allow implicit, silent conversions. Such conversions lead to awkward
 errors that are hard to debug, and in cases where they are necessary,
 are best described in code through explicit casts.
 <br>
 (Wolfgang Bangerth, 2016/06/25)
 </li>

 <li> Changed: The deal.II distributed vector classes do now derive from
 LinearAlgebra::VectorSpaceVector and have been moved to the
 LinearAlgebra::distributed namespace. In the definition of the new
 interfaces, several old vector functions have been marked as deprecated. The
 methods <tt>operator==</tt>, <tt>operator!=</tt>, and
 <tt>is_non_negative</tt> have been removed from the new interface.
 <br>
 (Martin Kronbichler, 2016/06/15)
 </li>

 <li> Changed: The Triangulation::Signals::clear signal is now triggered
 <i>before</i>, not <i>after</i> the internal data structures of the
 triangulation are destroyed. This allows functions attached to the signal to
 save information associated with the triangulation.
 <br>
 (Wolfgang Bangerth, 2016/06/07)
 </li>

 <li> Changed: deal.II used to create template instantiations for scalar
  types <tt>double</tt>, <tt>float</tt>, and <tt>long double</tt>. Since
  <tt>long double</tt> is rarely used and the additional precision does
  usually not pay off because most of the other arithmetics in deal.II are
  only done using <tt>double</tt> variables, it is not instantiated by default
  any more. This reduces the library size by up to 20 percent. In case you
  need instantiations of certain methods using <tt>long double</tt> data
  structures and get linker errors stating undefined symbols involving
  <tt>long double</tt>, include the respective <tt>.templates.h</tt> file(s)
  with the code definitions. See the section on @ref Instantiations in the
  manual for further information.
  <br>
  (Martin Kronbichler, 2016/04/26)
  </li>

  <li> Changed: FlatManifold takes as argument a periodicity option. This
  used to be a Point<dim>, but it should have been a Tensor<1,dim>. This
  is now changed.
  <br>
  (Luca Heltai, 2016/04/09)
  </li>

  <li> Changed: The default nodal point distribution of FE_Q, FE_DGQ,
  FE_Q_DG0, FE_Q_Bubbles, and FE_TraceQ has been changed from equidistant
  points to the node points of the corresponding Gauss-Lobatto quadrature
  formula. For degrees one and two, the Gauss-Lobatto quadrature is
  equidistant and thus the unit support points are as before. However, the
  Gauss-Lobatto points are more dense towards the element boundaries at higher
  degrees. This gives well-conditioned interpolation at arbitrary orders and
  much more stable computations. While these node distribution was available
  before, it was not very visible and often lead to misunderstandings by
  inexperienced users. Most codes will not be affected by this change, even
  those using cubic and higher degree polynomials, apart from slightly
  different (better) interpolation behavior and different entries in solution
  vectors. If you explicitly need equidistant points, use the constructors
  <tt>FE_Q<dim>(QIterated<1>(QTrapez<1>(),degree))</tt> or
  <tt>FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),degree))</tt>.
  <br>
  (Martin Kronbichler, 2016/04/05)
  </li>

  <li> Removed: Support for the legacy <code>Make.global_options</code>
  file has been removed.
  <br>
  (Matthias Maier, 2016/03/17)
  </li>

  <li> Removed: Functions with names containing <code>boundary_indicator</code>
  have been removed. They had previously already been deprecated, and replaced
  by functions containing the string <code>boundary_id</code> instead, to keep
  with the style used for <code>material_id</code>, <code>subdomain_id</code>,
  etc.
  <br>
  (Wolfgang Bangerth, 2016/02/28)
  </li>

  <li> Changed: Many functions in VectorTools and MatrixTools now require
  matching data types between vectors, matrices, and Function arguments.
  <br>
  (Denis Davydov, 2016/02/27)
  </li>

  <li> Changed: ConstraintMatrix::distribute_local_to_global() and numerous
  functions in VectorTools namespace now require matching data types.
  This is done to correctly handle complex-valued case.
  <br>
  (Denis Davydov, 2016/02/22)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>
 <li> Improved: The code in class GridReordering has been rewritten from
 scratch. It now follows the algorithm described in the paper by
 Agelek, Anderson, Bangerth and Barth mentioned in the documentation
 of that class.
 <br>
 (Wolfgang Bangerth, 2016/10/28)
 </li>
 
 <li> Changed: The <tt>generic_linear_algebra.h</tt> has been updated to
 include more type definitions for vectors, sparse matrices, solvers and
 preconditioners for native deal.II, PETSc and Trilinos data types.
 <br>
 (Jean-Paul Pelteret, 2016/10/21)
 </li>

 <li> Improved: deal.II now bundles a subset of BOOST 1.62 instead of a subset
 of BOOST 1.56.
 <br>
 (David Wells, 2016/10/20)
 </li>

 <li> New: Add a new FiniteElement class, FE_P1NC, to implement the scalar
 version of the P1 nonconforming finite element which is a piecewise linear
 element on quadrilaterals in 2d.
 <br>
 (Jaeryun Yim, 2016/10/01)
 </li>

 <li> New: FE_Enriched finite element class implements the partition of unitity
 method which allows to enrich the finite element space based on a priori
 knowledge about solution.
 <br>
 (Denis Davydov, 2016/09/28)
 </li>

 <li> Improved: The doxygen documentation now contains nicely formatted
 boxes containing the text message of each exception. Several messages
 haven been clarified and improved.
 <br>
 (Timo Heister, 2016/09/06)
 </li>

 <li> New: There are 6 new video lectures that explain the
 basics of Linux and the command line, how mesh refinement works, and some
 more complicated time stepping schemes.
 (@dealiiVideoLectureSeeAlso{2.9,2.91,17.25,17.5,17.75,30.25})
 <br>
 (Wolfgang Bangerth, 2016/08/19)
 </li>

 <li> New: deal.II no longer uses features of the C++ language that
 were deprecated with C++11, C++14, or that are scheduled to be
 deprecated for C++17.
 <br>
 (David Wells, Jonathan Robey, Wolfgang Bangerth, 2016/08/11)
 </li>

 <li>
 New: Added Python bindings to generate and manipulate a Triangulation from
 Python. The Triangulation generated in Python can be saved and later, loaded
 inside a C++ code.
 <br>
 (Bruno Turcksin, 2016/08/03)
 </li>

 <li>
 Improved: A few of the introductory examples (steps five through eight) no
 longer use the Function class; they use plain functions instead.
 <br>
 (David Wells, 2016/07/25)
 </li>

 <li> Improved: The build system now checks for usable compiler/linker
 flags during various stages of the configure run. This should catch the
 majority of issues by user supplied flags/libraries and unusable final
 link interfaces before we actually proceed to compile the library.
 <br>
 (Matthias Maier, 2016/07/13)
 </li>

 <li> Improved: The testsuite now supports fine grained feature constraints
 of the form <code>test.with_[feature]_with_[...]=true</code> corresponding
 to variables <code>DEAL_II_<FEATURE>_WITH_[...]</code> exported to
 <code>deal.IIConfig.cmake</code>.
 <br>
 (Matthias Maier, 2016/07/11)
 </li>

 <li> New: The library is now compatible with PETSc 3.7.0. Part of this change
 included adding a new header, <tt>petsc_compatibility.h</tt>, which provides
 some version-independent functions for using common PETSc functions.
 <br>
 (David Wells, 2016/07/07)
 </li>

 <li> Refactored: The contrib/ directory has been cleaned up and the
 Parameter GUI has be reloacted into its own repository:
 https://github.com/dealii/parameter_gui
 <br>
 (Matthias Maier, Timo Heister, 2016/07/06)
 </li>

 <li> New: Add new classes to expand a scalar finite element solution into
 the orthogonal bases FESeries::Fourier and FESeries::Legendre. Also
 provide auxiliary functions to calculate norms of subsets of expansion
 coefficients FESeries::process_coefficients and linear regression
 FESeries::linear_regression. Update step-27 to use this namespace to drive
 the hp-adaptive FEM solution process.
 <br>
 (Denis Davydov, 2016/06/23)
 </li>

 <li> New: The tutorial step-55 shows how to solve the Stokes system
 in parallel with PETSc or Trilinos.
 <br>
 (Timo Heister, 2016/06/17)
 </li>

 <li> New: The tutorial step-56 demonstrates Geometric Multigrid for the
 Stokes equations.
 <br>
 (Ryan Grove, Timo Heister, 2016/06/01)
 </li>

 <li> Improved: The step-44 tutorial now uses the new CellDataStorage class to
 store and retrieve local quadrature point data. An alternative approach to
 solving the linear system using the LinearOperator class has been implemented.
 <br>
 (Jean-Paul Pelteret, 2016/05/20)
 </li>

 <li> New: Add a collection of classes to manage user's quadrature point data:
 CellDataStorage, TransferableQuadraturePointData and
 parallel::distributed::ContinuousQuadratureDataTransfer.
 The implementation of CellDataStorage is flexible to support different types of
 data object at different cells. parallel::distributed::ContinuousQuadratureDataTransfer
 provides a convenient interface to transfer quadrature point data between cells
 of parallel::distributed::Triangulation.
 <br>
 (Denis Davydov, Jean-Paul Pelteret, 2016/04/30)
 </li>

 <li> New: Added an interface to the GNU Scientific Library. Also introduce a
 cubic spline interpolation function Functions::CSpline.
 <br>
 (Denis Davydov, 2016/04/28)
 </li>

 <li> New: Added move operations to BlockIndices, BlockVectorBase and
 BlockVector; Vector move operations nullify old object instead of
 using swap.
 <br>
 (Daniel Shapero, 2016/04/13)
 </li>

 <li> New: Manifold objects were previously only used to compute the
 locations of individual new points on a manifold. Now, they are also
 used to compute tangent vectors (via Manifold::get_tangent_vector()), and this
 functionality provides the basis for computing normal vectors to manifolds
 as well.
 <br>
 In many cases, tangent vectors can be computed quite easily if the
 manifold has a functional description, i.e., if it can be
 represented via the ChartManifold class. In those cases, it is only
 necessary to overload the ChartManifold::push_forward_gradient()
 function that computes the derivatives of the push forward operation.
 <br>
 (Luca Heltai, Wolfgang Bangerth, 2016/04/08)
 </li>

 <li> New: Added a new Mapping class, MappingManifold, to use exact
 geometrical information extracted from the Manifold description instead
 of a polynomial approximation when computing transformations from the
 reference to the real cell. This class allows the computation of
 quadrature points, tangent vectors, and normal vectors which are exact
 with respect to the geometrical description, and it uses the underlying
 Manifold objects of the Triangulation. MappingManifold coincides with
 MappingQ1 for the FlatManifold descriptor.
 <br>
 (Luca Heltai, 2016/04/09)
 </li>

 <li> New: Added indent target to indent all headers and source
 files. Now you can do make (or ninja) indent inside the build
 directory.
 <br>
 (Alberto Sartori, 2016/03/02)
 </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
 <li> New: Add PArpackSolver::reinit(const VectorType &distributed_vector) to
 initialize internal data structures based on a vector. This makes PArpack
 usable with MatrixFree operators.
 <br>
 (Denis Davydov, 2016/10/31)
 </li>

 <li> New: Add MatrixFreeOperators::LaplaceOperator representing a Laplace matrix.
 <br>
 (Denis Davydov, 2016/10/30)
 </li>

 <li> New: Add VectorTools::project() to do L2 projection
 of scalar-valued quadrature point data in parallel.
 <br>
 (Denis Davydov, 2016/10/28)
 <li>

 <li> Fixed: Increased precision of timesteps in DataOutInterface::write_pvd_record().
 <br>
 (Rajat Arora, 2016/10/29)
 </li>
 
<li> New: Add MatrixFreeOperators::MassOperator representing a mass matrix.
 <br>
 (Daniel Arndt, 2016/10/27)
 </li>

 <li> New: There is a new class DiagonalMatrix which represents a diagonal
 matrix via a vector. This is useful for representing Jacobi preconditioners
 with matrix-free methods.
 <br>
 (Martin Kronbichler, 2016/10/27)
 </li>

 <li> New: The class PreconditionChebyshev now offers a third template
 parameter PreconditionerType that is passed to the preconditioner setup via
 AdditionalData::preconditioner. This allows using other preconditioners than
 the default (and previous) selection of a point-Jacobi preconditioner.
 <br>
 (Martin Kronbichler, 2016/10/27)
 </li>

 <li> New: Add ArpackSolver::set_shift() to set the shift value in spectral
 transformation.
 <br>
 (Denis Davydov, 2016/10/25)
 </li>

 <li> New: PreconditionChebyshev now offers a PreconditionChebyshev::step()
 and PreconditionChebyshev::Tstep() methods for usage in relaxation smoothers.
 <br>
 (Martin Kronbichler, 2016/10/21)
 </li>

 <li> Fixed: GridIn::read_vtk() accidentally only read material ids of
 input cells correctly if the file listed them as integers. If they were
 listed them as floating point numbers, then unpredictable numbers were used.
 <br>
 (Wolfgang Bangerth, 2016/10/20)
 </li>

 <li> New: Add a base class for matrix-free operators MatrixFreeOperators::Base.
 <br>
 (Denis Davydov, 2016/10/16)
 </li>

 <li> New: There is now a function FEEvaluation::JxW() to return the Jacobian
 determinant times the quadrature weight in the matrix-free evaluation
 routines similary to FEValues.
 <br>
 (Martin Kronbichler, 2016/10/14)
 </li>

 <li> Fixed: GridGenerator::hyper_cube_slit() with colorized set to
 true is now working correctly.
 <br>
 (Timo Heister, 2016/10/04)
 </li>

 <li> Fixed: SphericalManifold now behaves correctly also when R>>1
 and the center is not the origin.
 <br>
 (Luca Heltai, 2016/10/01)
 </li>

 <li> New: FETools::extrapolate allows for using the
 extrapolate algorithm on parallel::distributed::Triangulations.
 <br>
 (Daniel Arndt, Martin Steigemann, 2016/09/28)
 </li>

 <li> Improved: Some parts of mesh refinement are now parallelized.
 <br>
 (Wolfgang Bangerth, 2016/09/27)
 </li>

 <li> Improved: MGSmootherBlock is now able to use the shared memory pool for
 temporary vector allocation. The constructor requiring an external memory
 allocation has therefore been deprecated.
 <br>
 (Jonathan Robey, 2016/09/21)
 </li>

 <li> New: DoFTools::write_gnuplot_dof_support_point_info outputs
 support point locations and dof indices to a format readable by
 gnuplot.
 <br>
 (Timo Heister, 2016/09/16)
 </li>

 <li> Fixed: EmbeddedRungeKutta methods now correctly increase delta_t_guess
 when the error is below coarsen_tol.
 <br>
 (Vaibhav Palkar, Bruno Turcksin, 2016/09/16)
 </li>

 <li> Fixed: TrilinosWrappers::MPI::Vector::locally_owned_elements()
 now returns the correct IndexSet also if initialized with two
 IndexSets.
 <br>
 (Daniel Arndt, 2016/09/16)
 </li>

 <li> Improved: The multigrid V-cycle has been rewritten for performance on
 large-scale machines. Rather than transferring parts of the defect
 immediately to all coarser levels with a complexity of O(n_levels) global
 communication steps per V-cycle, we now transfer the full defect once to the
 next coarser level only, resulting in crossing all processors only once.
 <br>
 (Martin Kronbichler, 2016/09/16)
 </li>

 <li> Fixed: The Multigrid W-cycle and F-cycle have been fixed (for uniform
 grids).
 <br>
 (Martin Kronbichler, 2016/09/16)
 </li>

 <li> New: LinearAlgebra::Vector is now instantiated for float and double.
 <br>
 (Bruno Turcksin, 2016/09/15)
 </li>

 <li> New: The class MGCoarseGridIterativeSolver is replacing
 MGCoarseGridLACIteration with a simpler interface.
 <br>
 (Timo Heister, 2016/09/14)
 </li>

 <li> Improved: FEValues no longer generates the mapping's internal database if
 the mapping will not be required for the set of update flags specified.
 <br>
 (Jonathan Robey, 2016/09/14)
 </li>

 <li> Fixed: Instantiating class Vector with non-standard template
 arguments did not work because of duplicate function symbols. This
 is now fixed.
 <br>
 (Dragan Nikolic, 2016/09/14)
 </li>

 <li> Fixed: IndexSet::make_trilinos_map now treats non-ascending but
 contiguous IndexSets correctly. It creates a linear EpetraMap only
 if the IndexSets are ascending and 1:1.
 <br>
 (Daniel Arndt, 2016/09/11)
 </li>

 <li> New: IndexSet::is_ascending_and_one_to_one allows to find out
 whether the nth range of indices is stored on the nth process in case
 the IndexSets are contiguous.
 <br>
 (Daniel Arndt, 2016/09/11)
 </li>

 <li> Fixed: The CMake macros <code>DEAL_II_(ADD_TEST|SETUP_TARGET)</code>
 now enforce a stricter <code>CMAKE_BUILD_TYPE</code> handling. This helps
 to avoid situations where targets with different build flavors might
 accidentally get linked against each other.
 <br>
 (Matthias Maier, 2016/09/08)
 </li>

 <li> Fixed: FE_TraceQ now provides unit support points.
 <br>
 (Martin Kronbichler, 2016/09/08)
 </li>

 <li> Fixed: Reimplement copy_triangulation and load in
 dealii::parallel::shared::Triangulation, this avoids the loss of
 partition information which causes parallel::shared::Triangulation to be in an invalid state.
 <br>
 (Ce Qin, 2016/09/05)
 </li>

 <li> Fixed: The build system now uses -fPIC instead of -fpic
 <br>
 (Matthias Maier, 2016/08/31)
 </li>

 <li> Fixed: Fix MPI_InitFinalize by correctly initializing and destroying
   all p4est/libsc related objects by calls to sc_init(), p4est_init(), and
   sc_finalize(); compatibility with p4est versions >1.1.
 <br>
 (Jonathan Perry-Houts, 2016/08/31)
 </li>

 <li> Improved: SparsityPattern::copy_from() copying from a
 DynamicSparsityPattern argument had quadratic complexity in the number of
 rows for sparsity patterns where most of the rows are of length zero. The bad
 algorithm has been replaced by a linear complexity one.
 <br>
 (Dustin Kumor, Martin Kronbichler, 2016/08/31)
 </li>

 <li> New: There is now the possibility to store information about the
 time of an output time step within the .visit file created by
 the DataOutInterface<dim,spacedim>::write_visit_record function.
 <br>
 (Rene Gassmoeller, Juliane Dannberg, 2016/08/24)
 </li>

 <li> New: It is now possible to generate a cell_iterator to a cell
 that is identified by a CellId. CellIds are unique even across
 processes in distributed computations, therefore this change allows
 to identify a particular cell (e.g. a ghost cell of the local process) in
 another domain.
 <br>
 (Rene Gassmoeller, 2016/08/17)
 </li>

 <li> New: Rank-4 symmetric tensors of type SymmetricTensor can now
 be converted to rank-4 tensors of type Tensor.
 <br>
 (Wolfgang Bangerth, 2016/08/11)
 </li>

 <li> New: Implement MGCoarseGridApplySmoother class to do a few steps of a
 smoother at the coarsest level.
 <br>
 (Denis Davydov, 2016/08/09)
 </li>

 <li> New: PreconditionMG can now be used as a LinearOperator.
 <br>
 (Denis Davydov, 2016/08/09)
 </li>

 <li> New: RelaxationBlock classes for geometric multigrid now support parallel
 computations using Trilinos.
 <br>
 (Timo Heister, Guido Kanschat, 2016/08/08)
 </li>

 <li>
 New: Added a new PolarManifold descriptor, that uses a polar coordinate
 system to compute new points, and modified the existing SphericalManifold
 descriptor to use geodesics on the surface of the sphere.
 <br>
 (Luca Heltai, Mauro Bardelloni, 2016/08/04)
 </li>

 <li> Improved: The regular and hp versions of
 DoFTools::make_flux_sparsity_pattern() no longer use the user flags of the
 underlying triangulation to determine if entries along a certain face have been
 added to the sparsity pattern.
 <br>
 (David Wells, 2016/03/02 - 2016/08/02)
 </li>

 <li> New: Introduce operators for residuals and interior penalty terms for
 the Grad-Div operator in LocalIntegrators::GradDiv.
 <br>
 (Timo Heister, Guido Kanschat, 2016/08/02)
 </li>

 <li> Improved: DoFTools::make_cell_patches() can create block lists
 only extending over local cells of distributed triangulations.
 <br>
 (Guido Kanschat, 2016/08/02)
 </li>

 <li> Fixed: (P)ARPACK interface for non-symmetric matrices.
 <br>
 (Joscha Gedicke, 2016/08/01)

 <li> Fixed: The TrilinosWrappers::SparsityPattern::print() and
 TrilinosWrappers::SparsityPattern::print_gnuplot() methods did not produce
 correct output on distributed computations. This is now fixed.
 <br>
 (Martin Kronbichler, 2016/07/30)
 </li>

 <li> Fixed: CMake now tries to pick up the full link interface for gsl.
 This works around an underlinkage issue with libgsl.so not correctly
 stating all shared object dependencies.
 <br>
 (Matthias Maier, 2016/07/28)
 </li>

 <li> Fixed: Level indices for geometric multigrid queried through
 DoFAccessor::get_mg_dof_indices() would return wrong indices on lines
 and faces in non-standard orientation in 3D. This is now fixed.
 <br>
 (Martin Kronbichler, 2016/07/27)
 </li>

 <li> New: There is now a new DoFTools::make_flux_sparsity_pattern()
 which takes a constraint matrix and flux and internal dof masks, in
 parallel. This is useful in the case where some components of a
 finite element are continuous and some discontinuous, allowing
 constraints to be imposed on the continuous part while also building
 building the flux terms needed for the discontinuous part.
 <br>
 (Sam Cox, 2016/07/25)
 </li>

 <li>
 Improved: VectorTools::interpolate() may now be used on FESystems with mixed
 interpolating and non-interpolating FEs, if all of the selected components for
 interpolation originate from interpolating FEs.
 <br>
 (Jonathan Robey, 2016/07/24)
 </li>

 <li>
 Improved: VectorTools::interpolate() now takes a ComponentMask to select the
 components to interpolate.
 <br>
 (Jonathan Robey, 2016/07/21)
 </li>

 <li> Improved: Split out pattern descriptions for LaTeX and Description
 ParameterHandler OutputStyles, and add better description text.
 <br>
 (Jonathan Robey, 2016/07/21)
 </li>

 <li> Improved: Allow for including dofs for individual components on
   boundary in DoFTools::make_vertex_patches().
 <br>
 (Ryan Grove, Daniel Arndt, 2016/07/21)
 </li>

 <li> Improved: Allow for initializing the constrained
 boundary DoFs in MGConstrainedDoFs using a std::set
 instead of a FunctionMap whose function values were not used.
 Allow for non-primitive FiniteElements.
 <br>
 (Daniel Arndt, 2016/07/20)
 </li>

 <li> New: Added GridGenerator::quarter_hyper_ball() to generate the
 intersection of a hyper ball with the positive orthant relative
 to its center.
 <br>
 (Daniel Arndt, 2016/07/19)
 </li>

 <li> Fixed: Work around an issue with the OpenMPI installation on certain
   Ubuntu versions: The build system now automatically drops the
   "-fuse-ld=gold" linker flag if openmpi is incompatible with it.
 <br>
 (Wolfgang Bangerth, Martin Kronbichler, Matthias Maier, 2016/07/13)
 </li>

 <li> Fixed: CMake now handles mixed compiler and linker setup via
 <code>DEAL_II_CXX_FLAGS*</code> / <code>DEAL_II_LINKER_FLAGS*</code> and
 <code>CMAKE_CXX_FLAGS*</code> properly.
 <br>
 (Matthias Maier, 2016/07/13)
 </li>

 <li> Fixed: FEValues::reinit() would sometimes try to be overly
 clever and not re-compute information when called with the same
 cell twice in a row, even if the underlying triangulation had
 been moved, translated, stretched, or otherwise had its vertex
 locations changed between the two calls to FEValues::reinit().
 This is now fixed.
 <br>
 (Wolfgang Bangerth, Jean-Paul Pelteret, Rajat Arora, 2016/07/11)
 </li>

 <li> Fixed: Allow to use FETools::get_fe_by_name for all
 availabale FiniteElements.
 <br>
 (Daniel Arndt, 2016/07/10)
 </li>

 <li> Fixed: SymmetricTensor::access_raw_entry() erroneously produced
 an indexing error for rank-4 symmetric tensors. This is now fixed.
 <br>
 (Wolfgang Bangerth, 2016/07/08)
 </li>

 <li> Fixed: SymmetricTensor::norm() did not work correctly for complex
 underlying scalar types. This is now fixed.
 <br>
 (Wolfgang Bangerth, 2016/07/08)
 </li>

 <li> New: There is now a function DerivativeForm::norm().
 <br>
 (Wolfgang Bangerth, 2016/07/08)
 </li>

 <li> New: A move constructor has been added to Triangulation.
 <br>
 (Daniel Shapero, 2016/07/07)
 </li>

 <li> Fixed: The function DoFTools::dof_couplings_from_component_couplings
 for hp::FECollection arguments was compiled but not exported from the
 object file. This is now fixed.
 <br>
 (Wolfgang Bangerth, 2016/07/01)
 </li>

 <li> New: Added TrilinosWrappers::SolveDirect::initialize and
 TrilinosWrappers::SolverDirect::solve to solve distributed linear systems
 with multiple right hand sides without needing to refactorize the matrix
 everytime.
 <br>
 (Michael Harmon, 2016/06/30)
 </li>

 <li> Changed: Similar to MappingQEulerian, MappingQ1Eulerian has gained
 a second constructor that reverses the order of the arguments to indicate
 which DoFHandler a vector is based on. The old constructor is now
 deprecated and will be removed in a future version.
 <br>
 (Wolfgang Bangerth, 2016/06/25)
 </li>

 <li> New: The MappingFEField class was previously only instantiated
 if the vector type was dealii::Vector. It is now also instantiated
 for PETSc and Trilinos wrapper vector types.
 <br>
 (Wolfgang Bangerth, 2016/06/25)
 </li>

 <li> New: GnuplotFlags now supports user specified space dimension labels
 through a member variable.
 <br>
 (David Wells, 2016/06/24)
 </li>

 <li> New: Added an operator* to SynchronousIterators which returns a
 reference to the stored tuple of iterators iterators. The iterators member
 may be made private in a future release.
 <br>
 (Daniel Shapero, 2016/06/24)
 </li>

 <li> New: IndexSet::at(idx) returns an iterator pointing to the given index
 or the next larger element in the set if idx is not contained.
 <br>
 (Timo Heister, 2016/06/21)
 </li>

 <li> Fixed: Performance of DynamicSparsityPattern::begin(r) and
 DynamicSparsityPattern::end(r) has been improved dramatically in parallel
 computations and if the pattern is empty.
 <br>
 (Timo Heister, 2016/06/21)
 </li>

 <li> Fixed: FEFieldFunction now works correctly in distributed computations,
 where before exceptions of type ExcPointNotAvailableHere could occur for
 evaluation points on or close to a boundary to a ghost cell.
 <br>
 (Timo Heister, 2016/06/06)
 </li>

 <li> Fixed: The Tensor class was not explicitly instantiated. This did
 not matter in almost all contexts because its members are all defined
 as @p inline in the header file. The only cases where it matters if one
 (or the compiler) were to take the address of one of the static member
 variables.
 <br>
 (Wolfgang Bangerth, 2016/06/03)
 </li>

 <li> New: Return value std::vector<unsigned int> vertex_mapping for the
 DoFTools::make_vertex_patches() function, including the optional inversion
 of the vertex mapping.
 <br>
 (Joscha Gedicke, 2016/05/25)
 </li>

 <li> Fixed: Fix a bug where the SparsityPattern could not have more than 4
 billions entries when using 32bit indices.
 <br>
 (Bruno Turcksin, 2016/05/22)
  </li>

 <li> New: There are now additional functions in the FETools::Compositing namespace that build
 finite elements out of simpler finite elements, either by forming tensor
 products or by combining the set of shape functions.
 <br>
 (Denis Davydov, Wolfgang Bangerth, 2016/05/20)
 </li>

 <li> New: Added PArpackSolver::reinit() when dealing with BlockVectors.
 <br>
 (Alberto Sartori, 2016/05/19)
 </li>

 <li> New: Add VectorTools::compute_global_error that computes global
 errors from cellwise errors obtained by VectorTools::integrate_difference()
 and do MPI collectives if necessary.
 <br>
 (Timo Heister, 2016/05/15)
 </li>

 <li> New: Add functions to transform Cartesian coordinates to spherical and back:
 GeometricUtilities::Coordinates::to_spherical and
 GeometricUtilities::Coordinates::from_spherical.
 <br>
 (Denis Davydov, 2016/05/10)
 </li>

 <li> Fixed: Corrected the sign of curl calculated in the functions:
 LocalIntegrators::curl_curl_matrix, LocalIntegrators::curl_matrix,
 LocalIntegrators::nitsche_curl_matrix and LocalIntegrators::ip_curl_matrix in
 integrators/maxwell.h.
 <br>
 (Jihuan Tian, 2016/05/09)
 </li>

 <li> Improved: The method Triangulation::create_triangulation will now throw an
 exception if any cells have negative measure. This check is not run if the
 triangulation keeps track of distorted cells or if the codimension is not zero.
 This check was previously only run in 3D.
 <br>
 (David Wells, 2016/05/07)
 </li>

 <li> Fixed: Bug in the RelaxationBlock class function do_step. Before, the
 corrections were not added together, which leads to a wrong update whenever the
 Jacobi blocks are overlapping. For SOR, SSOR and non-overlapping Jacobi this was
 not an issue.
 <br>
 (Joscha Gedicke, 2016/05/07)
 </li>

 <li> Fixed: The function GridGenerator::subdivided_parallelepiped and its
 variants could generate meshes with cells that had negative Jacobians.
 The function now detects when this will happen and raises a descriptive
 exception instead of going on to produce cells which may have negative measure.
 <br>
 (David Wells, 2016/05/11)
 </li>

 <li> New: Added function GridOut::write_mesh_per_processor_as_vtu. This allows
 the visualization of a parallel finite element mesh that can be separated into each
 processor's owned and ghost cells. It also allows for the visualization of each level
 of a multilevel mesh.
 <br>
 (Conrad Clevenger, 2016/04/28)
 </li>

 <li> Fixed: TrilinosWrappers::SparseMatrix will now exit early if there are no
 entries to add to the matrix. This usually occurs when zero elision is on. This
 fixes a bug where the matrix raises an exception if there are no entries to add
 to a matrix and the provided row and column values are not locally stored.
 <br>
 (David Wells, 2016/04/24)
 </li>

 <li> Fixed: TrilinosWrappers::MPI::Vector and TrilinosWrappers::Vector could
 access invalid memory in the reinit() method if the MPI communicator was
 deleted before termination of the program. This usually happened when using
 vectors from GrowingVectorMemory where a pool keeps vector alive. This has
 been fixed.
 <br>
 (Martin Kronbichler, 2016/04/23)
 </li>

 <li> Fixed: The methods TrilinosWrappers::SparseMatrix::(T)mmult previously
 produced invalid matrix sizes if the final matrix was non-square. This has
 been fixed.
 <br>
 (Martin Kronbichler, Daniel Jodlbauer, 2016/04/21)
 </li>

 <li> New: Added an optional string parameter to the ParameterHandler::read_input ()
 and ParameterHandler::read_input_from_string() functions.
 When a line which equals this string is encountered, the parsing of parameters
 is terminated.
 <br>
 (Denis Davydov, 2016/04/20)
 </li>

 <li> New: Added move operations to IndexSet.
 <br>
 (Daniel Shapero, 2016/04/19)
 </li>

 <li> Improved: MeshWorker treats periodic faces as interior faces.
 <br>
 (Daniel Arndt, 2016/04/18)
 </li>

 <li> Improved: The parallel loops in the deal.II Vector class for
 vector-vector operations have been revised for performance. This includes
 adjusting the minimum parallel grain size to 4096 vector entries and using an
 affinity partitioner provided by Threading Building Blocks for better data
 locality, especially on multi-socket systems.
 <br>
 (Martin Kronbichler, 2016/04/14)
 </li>

 <li> New: added ReinitHelper for PETSc. This is required by LinearOperator
 class to reinit vectors.
 <br>
 (Mauro Bardelloni, 2016/04/13)
 </li>

 <li> New: Added TensorProductManifold to create new manifolds from two
 ChartManifold objects. This can be used, for example, to combine a
 2d manifold with a flat manifold for an extruded mesh.
 <br>
 (Timo Heister, 2016/04/12)
 </li>

 <li> Fixed and improved: Fix algorithm for incomplete assignment of level
 subdomain ids for parallel geometric multigrid. Also optimize algorithms
 used for assignment and DoF communication.
 <br>
 (Timo Heister, Martin Kronbichler, 2016/04/12)
 </li>

 <li> New: Added New option in the read_ucd function of the GridIn class.
      A flag can now be assigned to the function, to decide wether the
      indicators specified in a UCD file should be interpreted as
      boundary_ids or as manifold_ids. This is particularly useful
      when the indicators refer to internal faces, for which
      boundary_ids cannot be used.
 <br>
 (Andrea Mola, 2016/04/11)
 </li>

 <li> Improved: DoFRenumbering::compute_Cuthill_McKee when used with
 distributed triangulations contained parts that scaled as the global problem
 size, rather than the processor-local size. This prevented its use with more
 than a few hundred cores when hanging node constraints were activated. This
 has been fixed.
 <br>
 (Martin Kronbichler, 2016/04/11)
 </li>

 <li> New: added hessenberg_signal and krylov_space_signal to SolverGMRES.
 These signals allow to retrieve the Hessenberg matrix and the basis vectors
 generated by the Arnoldi algorithm.
 <br>
 (Giuseppe Pitton, Luca Heltai, 2016/04/11)
 </li>

 <li> New: Added CompositionManifold to create new manifolds from two
 ChartManifold objects. This can be used, for example, to rotate a
 cylindrical Manifold, or to make a cylinders with parabolic sides.
 <br>
 (Luca Heltai, 2016/04/09)
 </li>

 <li> Fixed: Meshworker::Assembler::ResidualSimple now also works for
 multiple blocks if no constraints are given.
 <br>
 (Daniel Arndt, 2016/04/08)
 </li>

 <li> New: A move constructor has been added to Quadrature.
 <br>
 (Daniel Shapero, 2016/04/08)
 </li>

 <li> Fixed: The multigrid transfer performed invalid data accesses on
 multigrid hierarchies that define the coarse level as a level larger than
 0. This has been fixed.
 <br>
 (Martin Kronbichler, 2016/04/03)
 </li>

 <li> New: Add GridTools::remove_hanging_nodes() and
 GridTools::remove_anisotropy() in GridTools. GridTools::remove_hanging_nodes()
 detects cells with hanging nodes and refines the neighbours in the direction
 that removes hanging nodes or in every directions.
 GridTools::remove_anisotropy() refines a mesh until the resulting mesh is
 composed by cells with ratio between the extension in each coordinate
 direction lower than a fixed value.
 <br>
 (Mauro Bardelloni, 2016/03/28)
 </li>

 <li> New: When using C++11, a move constructor and assignment operator has
 been added to SparseMatrix, so that these objects can be returned from
 functions and packed into pairs and tuples.
 <br>
 (Daniel Shapero, 2016/03/27)
 </li>

 <li> New: The product of a rank-1 tensor (a vector) and a rank-2
 symmetric tensor (a symmetric matrix) is now defined and yields
 a rank-1 tensor (a vector). The opposite product was previously
 already defined.
 <br>
 (Wolfgang Bangerth, 2016/03/25)
 </li>

 <li> New: Triangulation::add_periodicity allows for accessing neighbors across
 periodic boundaries via new functions in TriaAccessor.
 <br>
 (Daniel Arndt, Ali Samii, 2016/03/23)
 </li>

 <li> New: Added GridGenerator::torus() to generate the volume mesh of a
 torus in three dimensions and a manifold description TorusManifold to
 go with it.
 <br>
 (Timo Heister, 2016/03/21)
 </li>

 <li> Fixed: DoFHandler::locally_owned_dofs() could create a segmentation
 fault in cases where some processors do not own any cells. This was caused
 by an incorrect computation in DoFTools::locally_owned_dofs_per_subdomain().
 <br>
 (Wolfgang Bangerth, 2016/03/20)
 </li>

 <li> Improved: The distribution of degrees of freedom on multigrid levels,
 DoFHandler::distribute_mg_dofs(), contained a few steps that scaled
 quadratically in the number of local cells for certain configurations. These
 steps have been replaced by linear complexity calls.
 <br>
 (Martin Kronbichler, 2016/03/18)
 </li>

 <li> New: Added GridTools::rotate() in three space dimensions.
 <br>
 (Timo Heister, 2016/03/18)
 </li>

 <li> New: Added custom target "relocate" to Mac OS X builds, that runs
 a script to make all paths absolute in the shared libraries included
 in the deal.II package (only enabled when building a package, and when
 including external libraries to the package)
 <br>
 (Luca Heltai, 2016/03/14)
 </li>

 <li> New: Added unit tests for complex-valued PETSc and SLEPc.
 <br>
 (Toby D. Young, Denis Davydov, 2016/03/11)
 </li>

 <li> New: Add NURBSPatchManifold. This class is a child of ChartManifold and
 implements a manifold descriptor for the face of a CAD imported usign
 OpenCASCADE.
 <br>
 (Mauro Bardelloni, 2016/03/09)
 </li>

 <li> New: When using C++11, there is now a function Threads::new_task()
 that can take as an argument either a lambda function, or the result
 of a std::bind expression, or anything else that can be called as in a
 function call. There is also a similar function Threads::new_thread()
 that takes the same kind of argument.
 <br>
 (Wolfgang Bangerth, 2016/03/07)
 </li>

 <li> New: Added another scaling factor to Kelly error estimator, namely h_K.
 <br>
 (Denis Davydov, 2016/03/05)
 </li>

 <li> New: When using C++11, the function filter_iterators() allows to filter a
 range of iterators using predicates (like the ones defined in IteratorFilter).
 <br>
 (Bruno Turcksin, 2016/03/04)
 </li>

 <li> Fixed: The OpenCASCADE::push_forward_and_differential_forms()
 function is now able to change the direction of the normal vector
 according to Orientation() method.
 <br>
 (Mauro Bardelloni, 2016/03/02)
 </li>

 <li> Fixed: The function IndexSet::make_trilinos_map() now works if some
 processors have a contiguous range of indices and others do not.
 <br>
 (Bruno Turcksin, 2016/02/17)
 </li>

 <li> Updated: step-44 has been been expressed in a more dimension independent
 manner, and can be now run in both 2-d and 3-d.
 <br>
 (Jean-Paul Pelteret, 2016/02/17)
 </li>

 <li> Fixed: FE_Nedelec elements up to polynomial order 12 can now be
 constructed.
 <br>
 (Jean-Paul Pelteret, 2016/02/12)
 </li>

 <li> Fixed: The GridTools::build_triangulation_from_patches() function now
 also copies the locations of vertices from the cells of the source
 triangulation to the triangulation that is built from the list of patch cells.
 <br>
 (Spencer Patty, 2016/02/11)
 </li>
</ol>

*/
