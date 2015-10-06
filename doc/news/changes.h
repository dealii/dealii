// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2015 by the deal.II authors
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
@page changes_after_8_3 Changes after Version 8.3.0

<p>
This is the list of changes made after the release of deal.II version
8.3.0. All entries are signed with the names of the authors.
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
  <li> Removed: The <code>UpdateFlags</code> flags
  <code>update_support_points</code>, <code>update_support_jacobians</code>,
  and <code>update_support_inverse_jacobians</code> have been removed.
  <code>update_support_points</code> was deprecated in 2013 and has not done
  anything in a long time (see the commit message for more information). The
  other two appeared in 2007 and were never implemented.
  <br>
  (David Wells, 2015/09/16)
  </li>

  <li> Cleanup: The two argument variant of cross_product() that returned
  the result by reference as first argument has been removed. Use the
  function cross_product_2d(), or cross_product_3d(), that directly returns
  the result instead. Further, the exception
  Tensor<rank,dim,Number>::ExcInvalidTensorContractionIndex has been
  removed
  <br>
  (Matthias Maier, 2015/09/14 - 2015/09/17)
  </li>

  <li> Cleanup: The following functions in tensor.h have been deprecated:
  <br>
  - double_contract(). Use the generic contraction function
    contract() instead.
  - The four and five argument variants of contract() that return the
    result by reference as first argument and take the contraction indices as
    arguments. use the generic contraction function contract() instead.
  - The three argument variants of contract() that return the result by
    reference as first argument. Use <code>operator*</code> instead.
  - The three argument variant of cross_product() that returns the result
    by reference as first argument. Use the cross_product() function that
    directly returns the result instead.
  - The three argument variants of <code>outer_product</code> that return
    the result by reference as first argument. Use the function that
    directly returns the result instead.
  - determinant(dealii::Tensor<rank,1,Number>)
  <br>
  (Matthias Maier, 2015/09/14 - 2015/09/17)
  </li>

  <li> Removed: Tensor<rank,dim,Number> as well as Point<dim,Number> no
  longer have a constructor taking a boolean argument. Those were replaced
  by a default constructor will always initialize underlying values with
  zero.
  <br>
  (Matthias Maier, 2015/09/07)
  </li>

  <li> Removed: The testsuite no longer supports compiler constraints of
  the form "<code>.compiler=[NAME]...</code>".
  <br>
  (Matthias Maier, 2015/08/21)
  </li>

  <li> Changed: The parameter @p first_vector_components has been removed
  from GridTools::collect_periodic_faces(). Instead,
  DoFTools::make_periodicity_constraints() now accepts a parameter
  @p first_vector_components in all (supported) variants.
  <br>
  (Matthias Maier, 2015/08/21)
  </li>

  <li> Changed: FEValues::normal_vector() for historical reasons returned a
  value of type Point, though a normal vector is more adequately described
  as a Tensor@<1,dim@>. Many similar cases were already clarified in deal.II
  8.3. The current case has now also been changed: FEValues::normal_vector()
  now returns a Tensor, rather than a Point.
  <br>
  In a similar spirit, the FEValues::get_normal_vectors() function that
  still returns a vector of Points has been deprecated and a new function,
  FEValues::get_all_normal_vectors(), that returns a vector of tensors,
  has been added. This was necessary since there is no way to change the
  return type of the existing function in a backward compatible way. The
  old function will be removed in the next version, and the new function
  will then be renamed to the old name.
  <br>
  (Wolfgang Bangerth, 2015/08/20)
  </li>

  <li> Changed: The mesh_converter program has been removed from the
  contrib folder. The equivalent functionality can now be found in
  the GridIn class.
  <br>
  (Jean-Paul Pelteret, 2015/08/12)
  </li>

  <li> Changed: The signature of the FiniteElement::fill_fe_values(),
  FiniteElement::fill_fe_face_values(), and FiniteElement::fill_fe_subface_values()
  functions has been changed, in an effort to clarify which of these contain
  input information and which contain output information for these functions.
  The same has been done for the corresponding functions in the Mapping
  class hierarchy. As part of a general overhaul, the FEValuesData class
  has also been removed.
  <br>
  (Wolfgang Bangerth, 2015/07/20-2015/08/13)
  </li>

  <li> Changed: The functions update_once() and update_each() in the
  Mapping classes computed information that was, in essence, only of use
  internally. No external piece of code actually needed to know which
  pieces of information a mapping could compute once and which they needed
  to compute on every cell. Consequently, these two functions have been
  removed and have been replaced by Mapping::requires_update_flags().
  <br>
  (Wolfgang Bangerth, 2015/07/20-2015/08/13)
  </li>

  <li> Changed: The function DoFRenumbering::random() now produces different
  numberings than it did before, but in return has now acquired the property
  that its results are predictable and repeatable.
  <br>
  (Wolfgang Bangerth, 2015/07/21)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>
<ol>
  <li> New: Triangulation::ghost_owners() returns the set of MPI ranks of the
  ghost cells. Similarly ::level_ghost_owners() for level ghosts.
  <br>
  (Timo Heister, 2015/09/30)
  </li>

  <li> New: FunctionParser now supports pow(a,b).
  <br>
  (Timo Heister, 2015/09/30)
  </li>

  <li> New: MGTransferPrebuilt can now be used with parallel::distributed::Vector
  and TrilinosWrappers::SparseMatrix as a transfer matrix.
  <br>
  (Martin Kronbichler, 2015/09/22)
  </li>

  <li> Fixed: parallel::distributed::Vector is now fully functional for
  indices larger than 4 billion.
  <br>
  (Martin Kronbichler, 2015/09/22)
  </li>

  <li> New: PArpackSolver eigensolver interface class.
  <br>
  (Denis Davydov, 2015/09/17)
  </li>

  <li> Changed: All doxygen-generated pages now contain a link to the
  tutorial in their top-level menus.
  <br>
  (Wolfgang Bangerth, 2015/09/13)
  </li>

  <li>New: A new namespace TensorAccessors is introduced that contains
  generic algorithms for tensorial objects, i.e., objects that allow
  repeated access via the index operator <code>operator[](unsigned int)</code>.
  The methods in TensorAccessors is primarily meant to replace old internal
  code in <code>tensor.h</code>, but it might also proof useful otherwise.
  <br>
  (Matthias Maier, 2015/09/11)
  </li>

  <li> New: A python script (including instructions) for enabling pretty
  printing with GDB is now available in
  <tt>/contrib/utilities/dotgdbinit.py</tt>.
  <br>
  (Wolfgang Bangerth, David Wells, 2015/09/11)
  </li>

  <li> Improved: When available, deal.II now uses the "gold" linker, a
  reimplementation of the traditional Unix "ld" linker that is substantially
  faster. This reduces build and, in particular, test turnaround times.
  <br>
  (Wolfgang Bangerth, Matthias Maier, 2015/09/06)
  </li>

  <li> Improved: Allow continuation lines in ParameterHandler.
  <br>
  (Alberto Sartori, 2015/09/04)
  </li>

  <li> Cleanup: The interface of Tensor<rank,dim,Number> has been cleaned
  up (a lot of unnecessary partial template specializations have been
  removed). The specialization Tensor<1,dim,Number> has been removed.
  <br>
  (Matthias Maier, 2015/09/02)
  </li>

  <li> Fixed: VectorTools::integrate_difference for VectorTools::Hdiv_seminorm
  was computed incorrectly.
  <br>
  (Timo Heister, 2015/08/31)
  </li>

  <li> Improved: The testsuite now supports multiple comparison files.
  Apart from the main comparison file that ends in
  <code>[...].output</code> all files of the form
  <code>[...].output.[string]</code> are considered for comparison.
  <br>
  (Matthias Maier, 2015/08/29)
  </li>

  <li> New: A class BlockLinearOperator has been introduced that extends
  the LinearOperator concept to block structures. A BlockLinearOperator can
  be sliced back to a LinearOperator.
  <br>
  (Matthias Maier, 2015/08/27)
  </li>

  <li> Improved: Support for complex number types throughout the library.
  Several parts of the library have been reorganized to support complex
  number types.
  <br>
  <em>Classes that are now instantiated for complex number types:</em>
  - FunctionTime
  - Function
  - TensorFunction

  <br>
  <em>Classes with fixed interface that now fully support complex number
  types (pure template classes without explicit instantiations in the
  library):</em>
  - LinearOperator
  - PackagedOperation
  - Tensor
  <br>
  (Matthias Maier, 2015/08/25 - XXX)
  </li>

  <li> Fixed: The testsuite now properly supports version constraints for
  features. Those are annotated by
  <code>.with_FEATURE(&lt;=|&gt;=|=|&lt;|&gt;)VERSION.</code>.
  <br>
  (Matthias Maier, 2015/08/25)
  </li>

  <li> Fixed: The GridIn class was not instantiated for the
  <code>dim==1,spacedim==3</code> case. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2015/08/25)
  </li>

  <li> Fixed: In 1d, GridIn::read_msh() ignored boundary indicators
  associated with vertices. This is now fixed.
  <br>
  (Jan Stebel, Wolfgang Bangerth, 2015/08/25)
  </li>

  <li> Improved: The interface and documentation for periodic boundary
  conditions have been restructured. A
  @ref GlossPeriodicConstraints "glossary entry" has been written.
  <br>
  (Daniel Arndt, Matthias Maier, 2015/08/01-2015/08/21)
  </li>

  <li> New: There is a new documentation module on
  @ref FE_vs_Mapping_vs_FEValues "How Mapping, FiniteElement, and FEValues work together".
  <br>
  (Wolfgang Bangerth, 2015/08/20)
  </li>

  <li> New: parallel::shared::Triangulation class which extends
  Triangulation class to automatically partition triangulation when run
  with MPI. Identical functionality between parallel::shared::Triangulation and
  parallel::distributed::Triangulation is grouped in the parent class
  parallel::Triangulation.
  <br>
  (Denis Davydov, 2015/08/14)
  </li>

  <li> New: The online documentation of all functions now includes
  links to the file and line where that function is implemented. Both
  are clickable to provide immediate access to the source code of a
  function.
  <br>
  (Jason Sheldon, Wolfgang Bangerth, 2015/08/13)
  </li>

  <li> New: implemented the gradient method for
  InterpolatedTensorProductGridData
  <br>
  (Daniel Shapero, 2015/08/12)
  </li>


  <li> New: FE_RannacherTurek describes a discontinuous FiniteElement
  with vanishing mean values of jumps across faces.
  <br>
  (Patrick Esser, 2015/08/17)
  </li>

  <li> New: FE_Q_Bubbles describes a FiniteElement based on FE_Q
  enriched by bubble functions.
  <br>
  (Daniel Arndt, 2015/08/12)
  </li>

  <li> New: The testsuite now runs in a mode in which we abort programs for
  floating point exceptions due to divisions by zero or other invalid arithmetic.
  <br>
  (Wolfgang Bangerth, 2015/07/29)
  </li>

  <li> New: MultithreadInfo::set_thread_limit() can now be called more than
  once and the environment variable DEAL_II_NUM_THREADS will be respected
  even if user code never calls it.
  <br>
  (Timo Heister, 2015/07/26)
  </li>

  <li> New: IndexSet now implements iterators.
  <br>
  (Timo Heister, 2015/07/12)
  </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>



<ol>

  <li> New: GridIn::read_unv() can now read more element codes that
  are used in typical meshes.
  <br>
  (Aslan Kosakian, 2015/10/06)
  </li>

  <li> New: DoFTools::locally_relevant_dofs_per_subdomain() can be used
  to extract an IndexSet of locally relevant DoFs for a Triangulation
  partitioned using METIS or with a parallel::shared::Triangulation .
  <br>
  (Jean-Paul Pelteret, 2015/09/24)
  </li>

  <li> Fixed: hp::SolutionTransfer could get confused when dealing with
  FE_Nothing elements. This is now fixed.
  <br>
  (Claire Bruna-Rosso, Wolfgang Bangerth, 2015/09/23)
  </li>

  <li> Improved: The construction of the non-local graph for quick data
  exchange of TrilinosWrappers::SparseMatrix became very slow for a few
  thousand processors. This has been fixed.
  <br>
  (Martin Kronbichler, 2015/09/22)
  </li>

  <li> Improved: Initializing a TrilinosWrappers::SparseMatrix from a
  DynamicSparsityPattern included some O(global_size) operations. These have
  been replaced by operations only on the local range.
  <br>
  (Martin Kronbichler, 2015/09/22)
  </li>

  <li> Changed: All doxygen-generated pages now contain a link to the
  tutorial in their top-level menus.
  <br>
  (Wolfgang Bangerth, 2015/09/13)
  </li>

  <li>Cleanup: Constructors of AdditionalData in various linear solvers are now marked
  explicit. This avoid bugs with implicit conversions like the one fixed in step-40.
  <br>
  (Timo Heister, Lei Qiao, 2015/09/09)
  </li>

  <li>New: Introduced third-order derivatives of the shape functions, which
  can now be accessed through FEValues and FEValuesViews using similar interfaces
  as shape_values, shape_derivatives and shape_hessians.
  (Maien Hamed, 2015/09/08)
  </li>

  <li>Cleanup: TableIndices<N> can now be used (constructed and accessed)
  with N > 7.
  <br>
  (Matthias Maier, 2015/09/08)
  </li>

  <li>New: std::begin and std::end are now available within the std_cxx11
  namespace through <base/std_cxx11/iterator.h>
  <br>
  (Matthias Maier, 2015/09/08)
  </li>

  <li>New: MappingQ1Eulerian was not instantiated for the various
  Trilinos vector types. It is now instantiated for the same
  vector types as MappingQEulerian is.
  <br>
  (Wolfgang Bangerth, 2015/09/08)
  </li>

  <li>New: Introduced Hessian-related functions to the Function class.
  <br>
  (Denis Davydov, 2015/09/08)
  </li>

  <li>New: Memory consumption during compilation has been reduced by splitting
  instantiation files. For this make_instantiations now supports additional
  logic to split the the instantiations in .inst files into groups. This is
  used in fe_values.cc, error_estimator.cc, and others.
  <br>
  (Timo Heister, 2015/09/05)
  </li>

  <li> New: There is now a function SparsityPattern::print_svg() which prints the sparsity of the matrix
  in a .svg file which can be opened in a web browser.
  <br>
  (Conrad Clevenger, 2015/09/03)
  </li>

  <li> Openmp SIMD support is now enabled for Clang version 3.6, or newer
  (or the equivalent XCode version). Further, openmp support is not any
  more falsely activated for very old clang versions.
  <br>
  (Matthias Maier, 2015/09/03)
  </li>

  <li> Improved: DoFTools::make_hanging_node_constraints() now supports
  hp-refinement cases when neither_element_dominates. To that end we look for
  a least face dominating FE inside FECollection.
  <br>
  (Denis Davydov, 2015/09/02)
  </li>

  <li> Changed: FEValues::transform() has been deprecated. The functionality
  of this function is a (small) subset of what the Mapping classes
  already provide.
  <br>
  (Wolfgang Bangerth, 2015/09/02)
  </li>

  <li> New: introduced hp::FECollection::find_least_face_dominating_fe(const std::set<unsigned int> &fes)
  which aims to find the least dominating finite element w.r.t. those provided
  as fe_indices in @p fes.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2015/08/31)
  </li>

  <li> New: step-6 now has an additional subsection in the
  "Possibilities for extensions" section that discusses how
  to create a better mesh.
  <br>
  (Konstantin Ladutenko, Wolfgang Bangerth, 2015/08/31)
  </li>

  <li> New: Introduce an option for FE_Nothing to dominate any other FE.
  Therefore at interfaces where, for example, a Q1 meets an FE_Nothing,
  we will force the traces of the two functions to be the same. Because the
  FE_Nothing encodes a space that is zero everywhere, this means that the Q1
  field will be forced to become zero at this interface.
  <br>
  (Denis Davydov, 2015/08/31)
  </li>

  <li> New: Jacobian second and third derivatives are now computed by the mapping classes and can be
  accessed through FEValues in much the same way as the Jacobian and Jacobian gradient.
  <br>
  (Maien Hamed, 2015/08/28-2015/08/31)
  </li>

  <li> New: There are now a collection of functions named GridTools::compute_active_cell_halo_layer()
  that determine which cells form a layer around a specified subdomain. There is also a function
  GridTools::compute_ghost_cell_halo_layer() that returns the smallest layer of ghost cells around
  all locally relevant cells.
  <br>
  (Jean-Paul Pelteret, Denis Davydov, Wolfgang Bangerth, 2015/08/21)
  </li>

  <li> Documentation: How to set up a testsuite in a user project is now
  properly documented.
  <br>
  (Matthias Maier, 2015/08/01 - 2015/08/20)
  </li>

  <li> Fixed: The computation of gradients in FE_PolyTensor and its derived classes (in particular in
  FE_RaviartThomas and FE_Nedelec) forgot to account for terms that appear on non-affine
  cells. Consequently, the computed gradients did not match the actual derivatives of the values
  these elements report. This is now fixed.
  <br>
  (Maien Hamed, 2015/08/18-2015/08/20)
  </li>

  <li> Improved: Generalized conversion between Tensor<order+1,dim> and
  DerivativeForm<order,dim,dim> to general order using converting constructor
  and assignment operator.
  <br>
  (Maien Hamed, 2015/08/01-2015/08/09)
  </li>

  <li> Changed: The function Vector::add() that adds a scalar number to all
  elements of a vector has been deprecated. The same is true for the
  Vector::ratio() function, and for the corresponding functions in other
  vector classes.
  <br>
  (Wolfgang Bangerth, Bruno Turcksin, 2015/08/13)
  </li>

  <li> New: Direct support for Abaqus mesh files has been added to the GridIn
  class through the function GridIn::read_abaqus().
  <br>
  (Jean-Paul Pelteret, Timo Heister,  Krzysztof Bzowski, 2015/08/12)
  </li>

  <li> Improved: Finite elements now compute hessians analytically rather than
  by finite differencing.
  <br>
  (Maien Hamed, 2015/08/01-2015/08/09)
  </li>

  <li> New: There is now a function Mapping::project_real_point_to_unit_point_on_face()
  that calls Mapping::transform_real_to_unit_cell() and then projects the
  result to a provided face.
  <br>
  (Jason Sheldon, 2015/08/11)
  </li>

  <li> New: FEFaceValues and FESubfaceValues can now also compute
  gradients of the Jacobian of the transformation from unit to real cell,
  controlled by update_jacobian_grads.
  <br>
  (Martin Kronbichler, 2015/08/08)
  </li>

  <li> New: There is now a function MemoryConsumption::memory_consumption()
  for std_cxx11::unique_ptr arguments.
  <br>
  (Wolfgang Bangerth, 2015/08/07)
  </li>

  <li> Improved: CMake configuration: The DEAL_II_ADD_TEST now also
  supports unit tests writing to stdout and stderr. Further, a second test
  type consisting of an internal executable target, a configuration and a
  comparison file is now supported.
  <br>
  (Matthias Maier, 2015/08/03)
  </li>

  <li> New: VtkFlags now stores a parameter describing the compression level
  zlib uses when writing compressed output. For small problems, the flag
  ZlibCompressionLevel::best_speed can make the call to write_vtu many times
  faster.
  <br>
  (David Wells, 2015/08/03)
  </li>

  <li> Improved: The conversion Epetra_Map -> IndexSet is now an O(1)
  operation for contiguous index ranges, improving over the old O(N) behavior.
  <br>
  (Martin Kronbichler, 2015/07/30)
  </li>

  <li> Changed: The initialization methods of TrilinosWrappers::SparseMatrix,
  TrilinosWrappers::BlockSparseMatrix, TrilinosWrappers::SparsityPattern, and
  TrilinosWrappers::BlockSparsityPattern with Epetra_Map arguments have been
  marked as deprecated. Use the functions with IndexSet argument instead.
  <br>
  (Martin Kronbichler, Luca Heltai, 2015/07/30)
  </li>

  <li> New: FESystem now does some work in parallel if your system
  has multiple processors.
  <br>
  (Wolfgang Bangerth, 2015/07/19)
  </li>

  <li> Fixed: When using FESystem with base elements that require
  information other than the determinant of the Jacobian (e.g.,
  elements that require the Jacobian itself), then this information
  was not passed down to FiniteElement::fill_fe_values of the
  base element. This is now fixed.
  <br>
  (Wolfgang Bangerth, Zhen Tao, 2015/07/17)
  </li>

  <li> New: The parallel::distributed::Triangulation can now be told to
  partition the cells so that the sum of certain weights associated with each
  cell, rather than the number of cells, is roughly constant between processors.
  This is done by passing a vector of weights to the function that repartitions
  the triangulation, parallel::distributed::Triangulation::repartition().
  <br>
  (Wolfgang Bangerth, 2015/07/14)
  </li>

  <li> New: DataOutBase::TecplotFlags now takes a third argument for solution
  time which is useful to visualize transient data. If a user sets a non-negative
  time, it will be saved into the tecplot file.
  <br>
  (Praveen Chandrashekar, 2015/08/30)
  </li>

</ol>

*/
