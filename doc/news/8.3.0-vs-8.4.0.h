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
@page changes_between_8_3_and_8_4 Changes between Version 8.3.0 and 8.4.0

<p>
This is the list of changes made between the release of deal.II version
8.3.0 and that of 8.4.0. All entries are signed with the names of the
author.
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
  <li> Changed: GridGenerator::hyper_rectangle and
  GridGenerator::subdivided_hyper_rectangle now take points with @p dim
  components to correctly handle meshes embedded in higher dimensional spaces.
  <br>
  (Timo Heister, 2016/02/04)
  </li>

  <li> Changed: The constructor of FiniteElementData had a last argument
  <code>n_blocks</code> that was not actually used by the class to
  initialize anything. It has been removed. In addition, the default
  constructor of FiniteElementData has also been removed given that it
  only creates a dysfunctional element.
  <br>
  (Wolfgang Bangerth, 2016/01/23).
  </li>

  <li> Rework: SLEPcWrappers were reworked to allow usage of PETSc solvers
  and preconditioners inside SLEPc's eigensolvers. To that end extra methods
  were introduced to PETSc wrappers. Moreover, initialization of the
  underlying SLEPc objects is now done inside constructors of the wrapper
  classes. As a result, one has to provide an MPI communicator to the constructors of
  spectral transformation classes.
  <br>
  (Denis Davydov, 2015/12/29).
  </li>

  <li> Removed: The deprecated Operator class in the Algorithms namespace has been
  removed.
  <br>
  (Timo Heister, 2015/12/21)
  </li>

  <li> Changed: deallog console depth is now 0 by default, causing no
  output to the screen from solvers and other places in the library.
  <br>
  (Timo Heister, 2015/12/06)
  </li>

  <li> Changed: The function Utilities::trim() now removes general
  white space characters, such as '<tt>\\r</tt>' and '<tt>\\n</tt>', as well as
  space characters.
  <br>
  (David Wells, 2015/12/05)
  </li>

  <li> Removed: The previously deprecated global instance
  <code>multithread_info</code> of
  MultithreadInfo has been removed (all members of this class are static
  so there is no reason to use/create an instance). The deprecated
  MultithreadInfo::n_cpus member also got removed in favor of
  MultithreadInfo::n_cores().
  <br>
  (Timo Heister, 2015/11/19)
  </li>

  <li> Removed: The <code>UpdateFlags</code> flags
  <code>update_support_points</code>, <code>update_support_jacobians</code>,
  and <code>update_support_inverse_jacobians</code> have been removed.
  <code>update_support_points</code> was deprecated in 2013 and has not done
  anything in a long time. The
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

  <li> Removed: The Tensor and Point classes no
  longer have a constructor taking a boolean argument. Those were replaced
  by a default constructor that will always initialize underlying values with
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
  still returns a vector of Points has been deprecated, and a new function,
  FEValues::get_all_normal_vectors() that returns a vector of tensors,
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
  class hierarchy.
  <br>
  Likewise the signature of FiniteElement::get_data() has been changed.
  <br>
  As part of a general overhaul, the FEValuesData class
  has also been removed.
  <br>
  (Wolfgang Bangerth, 2015/07/20-2015/08/13)
  </li>

  <li> Changed: The functions update_once() and update_each() in the
  Mapping classes computed information that was, in essence, only of use
  internally. No external code actually needed to know which
  pieces of information a mapping could compute once and which they needed
  to compute on every cell. Consequently, these two functions have been
  removed and have been replaced by Mapping::requires_update_flags().
  <br>
  A similar change has been applied to the FiniteElement class.
  <br>
  (Wolfgang Bangerth, 2015/07/20-2015/12/01)
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
  <li> New: A variant for GridGenerator::subdivided_parallelepiped() was added
  that supports meshes embedded in higher dimesional spaces.
  <br>
  (Timo Heister, 2016/02/04)
  </li>

  <li> Fixed: Partitioning using METIS now works correctly with more
  domains than cells.
  <br>
  (Timo Heister, 2016/01/26)
  </li>

  <li> New: The documentation of step-17 has been completely rewritten,
  and many aspects of how one has to think when writing parallel programs
  have been much better documented now.
  <br>
  (Wolfgang Bangerth, 2016/01/07)
  </li>

  <li> New: deal.II now provides a string <code>DEAL_II_ALWAYS_INLINE</code>
  that, when supported by the compiler, can be used to annotate functions
  to ensure that the compiler always inlines them.
  <br>
  (Matthias Maier, Wolfgang Bangerth, 2016/01/07)
  </li>

  <li> New: There is a new documentation module, @ref Concepts, which describes the meaning
  behind template parameter type names.
  <br>
  (David Wells, 2015/12/09)
  </li>

  <li> Changed: The template type name arguments of some classes no longer
  shadow class names. Additionally, template type names are now much more
  consistent across deal.II.
  <br>
  (David Wells, 2015/10/18 - 2016/01/23)
  </li>

  <li> New: The WorkStream class's design and implementation are now much
  better documented in the form of a @ref workstream_paper "preprint".
  <br>
  (Wolfgang Bangerth, 2015/11/29)
  </li>

  <li> New: There is now much more documentation for the FiniteElement class,
  in particular detailing what one needs to implement when writing finite
  element descriptions in derived classes.
  <br>
  (Wolfgang Bangerth, 2015/11/29)
  </li>

  <li> New: We now experimentally support Microsoft Visual C++ compiler under
  Windows.
  <br>
  (Timo Heister, 2015/11/26)
  </li>

  <li> New: There is now a function template numbers::signaling_nan() that
  is used to create invalid floating point objects. These objects can either
  be scalars, or of type Tensor, SymmetricTensor, or DerivativeForm. The
  content of these objects is a "signaling NaN" ("NaN" stands for "not a
  number", and "signaling" implies that at least on platforms where this
  is supported, any arithmetic operation using them terminates the program).
  The purpose of this is to use them as markers for uninitialized objects
  and arrays that are required to be filled in other places, and to trigger
  an error when this later initialization does not happen before the first
  use.
  <br>
  (Wolfgang Bangerth, Timo Heister, 2015/11/24)
  </li>

  <li> Changed: The function FE_DGPNonparametric::shape_value() and similar
  functions in the same class returned values and derivatives of shape
  functions on the reference cell. However, this element is not defined
  through mapping of shape functions from the reference cell, and consequently
  it makes no sense to ask for this information. These functions have therefore
  been changed to throw an exception instead, as documented in
  FiniteElement::shape_value().
  <br>
  (Wolfgang Bangerth, 2015/11/20)
  </li>

  <li> Changed: The functionality to distribute cells across processes
  according to a vector of cell weights that was passed in a call to
  parallel::distributed::Triangulation::repartition()
  was replaced by a cell-wise signal. This signal is called during
  parallel::distributed::Triangulation::execute_coarsening_and_refinement() and
  parallel::distributed::Triangulation::repartition()
  if any function is connected to it. It allows to connect a function that
  takes the current cell iterator and a status argument that indicates whether
  this cell will be refined, coarsened or remains unchanged and returns a
  cell weight, which will be used to distribute cells across processes in a
  way that keeps the sum of weights across each individual process
  approximately equal.
  <br>
  (Rene Gassmoeller, 2015/11/02)
  </li>

  <li> New: Preliminary support for parallel, adaptive, geometric multigrid is
  now in place with changes to MGConstrainedDoFs (many new functions), MGTransfer,
  MGTools::extract_inner_interface_dofs, MGTransferPrebuilt,
  DoFTools::extract_locally_relevant_level_dofs.
  <br>
  (Timo Heister, Guido Kanschat, 2015/10/26)
  </li>

  <li> New: Two cell level signals are added to class Triangulation, namely
  pre_coarsening_on_cell and post_refinement_on_cell.
  <br>
  (Lei Qiao, 2015/10/22)
  </li>

  <li> New: parallel::distributed::Triangulation::ghost_owners()
  returns the set of MPI ranks of the ghost cells. Similarly
  parallel::distributed::Triangulation::level_ghost_owners() for level
  ghosts.
  <br>
  (Timo Heister, 2015/09/30)
  </li>

  <li> Improved: The interfaces to all deal.II type solvers and
  preconditioners have been updated such that they function as expected
  with the LinearOperator class and its associated functions (i.e.,
  linear_operator(), transpose_operator() and inverse_operator()).
  These preconditioners can now be wrapped as a LinearOperator,
  facilitating the construction of approximate matrix inverses such as in
  the development of a block matrix preconditioner. An example of this
  functionality can be found in
  <code>tests/lac/linear_operator_08.cc</code>.
  <br>
  (Jean-Paul Pelteret, 2015/09/24 - 2015/10/19)
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

  <li> Cleanup: The interface of Tensor<rank,dim,Number> has been cleaned
  up (a lot of unnecessary partial template specializations have been
  removed). The specialization Tensor<1,dim,Number> has been removed.
  <br>
  (Matthias Maier, 2015/09/02)
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
  (Matthias Maier, 2015/08/25)
  </li>

  <li> Fixed: The testsuite now properly supports version constraints for
  features. Those are annotated by
  <code>.with_FEATURE(&lt;=|&gt;=|=|&lt;|&gt;)VERSION.</code>.
  <br>
  (Matthias Maier, 2015/08/25)
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
  with MPI. Common functionality between parallel::shared::Triangulation and
  parallel::distributed::Triangulation is implemented in the parent class
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
  <li> New: FunctionParser now supports <code>rand()</code> and
  <code>rand_seed(number)</code>, which return a random value in the
  range [0,1].
  <br>
  (Wolfgang Bangerth, Luca Heltai, Alberto Sartori, 2016/02/09)
  </li>

  <li> Fixed: FullMatrix::TmTmult for matrix multiplication used to compute
  wrong results for larger matrix sizes where external BLAS is called. This
  has been fixed.
  <br>
  (Martin Kronbichler, 2016/02/02)
  </li>

  <li> Fixed: parallel::distributed::Triangulation with periodic boundary
  conditions did not respect 2:1 balance over vertices on periodic
  boundaries. This lead to incomplete ghost layers on multigrid levels. This
  has been fixed.
  <br>
  (Martin Kronbichler, Timo Heister, 2016/01/27)
  </li>

  <li> Fixed: SparseVanka now really uses second-order couplings for the
  right-hand side of the local problems.
  <br>
  (Florian Sonner, 2016/01/27)
  </li>

  <li> Fixed: A bug in the Neumann boundary handling of KellyErrorEstimator
  in 1d has been fixed and KellyErrorEstimator now correctly handles
  codimension one problems by using the correct normals from the manifold
  inside the gradient jump computation.
  <br>
  (Andrea Bonito, Timo Heister, 2016/01/21)
  </li>

  <li> New: The new class MGTransferMatrixFree implements multigrid level
  transfer using local polynomial embedding and restriction with tensor
  product evaluation techniques. This is a faster and less memory-demanding
  alternative to MGTransferPrebuilt.
  <br>
  (Martin Kronbichler, 2016/01/20)
  </li>

  <li> New: hp::FECollection now has constructors which take
  multiple finite elements as arguments.
  <br>
  (Angel Rodriguez, 2016/01/18)
  </li>

  <li> New: The glossary now contains a long entry describing what
  the term "scalability" means in the context of finite element codes.
  See @ref GlossParallelScaling.
  <br>
  (Wolfgang Bangerth, 2016/01/11)
  </li>

  <li> Fixed: Tensor::operator[] that takes TableIndices as a parameter no
  longer returns by value, but rather by reference. Tensor::operator<< for
  dim==0 now accesses values by reference instead of making a copy. This is
  useful when non-trivial number types are stored.
  <br>
  (Jean-Paul Pelteret, 2016/01/08)
  </li>

  <li> New: constrained_linear_operator() and constrained_right_hand_side()
  provide a generic mechanism of applying constraints to a LinearOperator.
  A detailed explanation with example code is given in the @ref constraints
  module.
  <br>
  (Mauro Bardelloni, Matthias Maier, 2015/10/25 - 2015/12/27)
  </li>

  <li> New: OpenCASCADE::read_IGES() and OpenCASCADE::read_STEP() have
  been unified in behaviour, and now they allow to extract *all* elements of
  the IGES and STEP files instead of only the faces. This allows the
  use of iges files describing edges only to be used as input for some of
  the OpenCASCADE Manifold wrappers.
  <br>
  (Luca Heltai, 2015/12/13)
  </li>

  <li> New: A new linear operator representing the Schur complement,
  namely schur_complement(), has been implemented. Some auxiliary functions
  that are often used in conjunction with the Schur complement
  (condense_schur_rhs() and postprocess_schur_solution()) are also provided
  as a PackagedOperation.
  An example of this functionality can be found in
  <code>tests/lac/schur_complement_01.cc</code>.
  The solution of a multi-component problem (namely step-22) using the
  schur_complement can be found in
  <code>tests/lac/schur_complement_03.cc</code> .
  <br>
  (Jean-Paul Pelteret, Matthias Maier, Martin Kronbichler, 2015/12/07)
  </li>

  <li> New: There is now a function Utilities::to_string that works like
  int_to_string, but is more safe for long integers, negative integers, and
  also handles floating point numbers. The implementation of int_to_string
  was changed to simply call to_string. int_to_string is kept for
  compatibility, but should only be used for unsigned integers.
  <br>
  (Rene Gassmoeller, 2015/12/09)
  </li>

  <li> Fixed: GridOut::write_msh() and GridOut::write_ucd() used the same
  geometric element numbers for lines and faces. This caused visualization
  programs to ignore parts with repeated geometric element numbers. This is now
  fixed.
  <br>
  (David Wells, 2016/01/16)
  </li>

  <li> New: DoFTools::extract_dofs() are now instantiated also for
  codimension different from zero.
  <br>
  (Alberto Sartori, 2016/01/13)
  </li>

  <li> Fixed: The DataOutFaces class should now also work with triangulations
  of type parallel::distributed::Triangulation.
  <br>
  (Heikki Virtanen, Wolfgang Bangerth, 2016/01/11)
  </li>

  <li> Fixed: AlignedVector<T>::fill() (and thus, Table<N,T>::reinit) did not
  correctly call the destructor of T() and could leak memory for complicated
  class types that depend on their constructor to free memory.
  <br>
  (Martin Kronbichler, 2016/01/08)
  </li>

  <li> Fixed: inverse_operator() now populates <code>Tvmult</code> and
  <code>Tvmult_add</code> correctly.
  <br>
  (Jean-Paul Pelteret, David Wells, Matthias Maier, 2015/12/30)
  </li>

  <li> New: MGTransferPrebuilt with parallel adaptive refinement has been
  finalized for parallel::distributed::Vector.
  <br>
  (Martin Kronbichler, 2015/12/23)
  </li>

  <li> Fixed: Now all members in the class SparseMatrixEZ are initialized
  correctly in the constructor. This was causing random crashes before.
  <br>
  (Timo Heister, 2015/12/21)
  </li>

  <li> New: There is now a new class ArrayView that presents a chunk of
  memory as if it was an array of fixed size. This is eventually going
  to replace the VectorSlice class which suffers from the defect that
  its template argument does not encode the type of objects it points
  to, but instead the type of the underlying container; consequently,
  where the VectorSlice class is used as a function argument, it
  automatically ties the type of object the function can be called
  with (i.e., the underlying container) even if the called function
  has no actual use for this kind of information.
  <br>
  (Wolfgang Bangerth, 2015/12/20)
  </li>

  <li> Fixed: Handling of constraints in step-26 was incorrect (hanging nodes
  were condensed twice) leading to garbage solutions. This is now fixed.
  <br>
  (Timo Heister, 2015/12/20)
  </li>

  <li> Fixed: The implementation of ShiftedMatrixGeneralized contained several
  errors that prevented it from being compiled. These have now been fixed.
  <br>
  (David Wells, 2015/12/18)
  </li>

  <li> New: There is now a function Triangulation::get_triangulation() that
  allows writing code to get at the underlying triangulation for
  everything that looks like a container, i.e., both Triangulation
  or DoFHandler objects.
  <br>
  (Wolfgang Bangerth, 2015/12/10)
  </li>

  <li> Deprecated: The functions DoFHandler::get_tria() and
  hp::DoFHandler::get_tria() were deprecated. UseDoFHandler::get_triangulation() and
  hp::DoFHandler::get_triangulation() instead.
  <br>
  (Wolfgang Bangerth, 2015/12/10)
  </li>

  <li> New: parallel::distributed::Vector has now a method to return a shared
  pointer to the underlying partitioner object.
  <br>
  (Martin Kronbichler, 2015/12/07)
  </li>

  <li> Improved: Many more functions in namespace GridTools and class
  InterGridMap are now consistely instantiated also for types
  parallel::distributed::Triangulation and parallel::shared::Triangulation.
  <br>
  (Gennadiy Rishin, Wolfgang Bangerth, 2015/12/07)
  </li>

  <li> Improved: Both versions of SparsityTools::distribute_sparsity_pattern()
  are now plain, not
  template, functions. This is not a breaking change because each function was
  instantiated for exactly one template argument.
  <br>
  (David Wells, 2015/12/06)
  </li>

  <li> Improved: The method
  parallel::distributed::Triangulation::fill_vertices_with_ghost_neighbors()
  that is used for distributing DoFs on parallel triangulations previously
  exhibited quadratic complexity in the number of coarse grid cells. This has
  been changed into linear complexity calls (apart from a few issues
  inside p4est).
  <br>
  (Martin Kronbichler, 2015/12/05)
  </li>

  <li> New: There are now new functions
  GridTools::build_triangulation_from_patch() and
  GridTools::get_cells_at_coarsest_common_level() that help build patches
  around individual cells.
  <br>
  (Arezou Ghesmati, 2015/12/02)
  </li>

  <li> Fixed: The GridTools::copy_boundary_to_manifold_id() function
  only copied boundary indicators from faces, but in 3d forgot about
  edges. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2015/11/30)
  </li>

  <li> Fixed: The constructor of SymmetricTensor that takes an array
  of initializing elements led to a compiler error. This is now
  fixed.
  <br>
  (Wolfgang Bangerth, 2015/11/28)
  </li>

  <li> Fixed: parallel::distributed::Vector now detects if the size of MPI
  messages exceeds 2GB or if the local range exceeds the size of 32-bit
  integers and throws an exception informing about the unsupported sizes.
  <br>
  (Martin Kronbichler, 2015/11/26)
  </li>

  <li> New: In 3d, GridGenerator::extract_boundary_mesh() now copies the
  manifold ids of edges of the volume mesh to the manifold ids of the edges
  of the extracted surface mesh.
  <br>
  (Wolfgang Bangerth, 2015/11/25)
  </li>

  <li> New: Triangulation::create_triangulation() now accepts subcell-data
  objects that may include information about interior edges and faces, to
  facilitate setting manifold indicators on interior edges and faces.
  <br>
  (Wolfgang Bangerth, 2015/11/25)
  </li>

  <li> Fixed: GridGenerator::extract_boundary_mesh() in 3d could generate
  surface cells that did not uniformly had a right- or left-handed coordinate
  system associated with them when viewed from one side of the surface. This
  has been fixed: they now all have a right-handed coordinate system when seen
  from one side of the surface, and a left-handed one when viewed from the
  other side.
  <br>
  (Daniel Weygand, Wolfgang Bangerth, 2015/11/22)
  </li>

  <li> Fixed: Trilinos ML preconditioner is now deterministic when using
  version 12.4 or newer.
  <br>
  (Timo Heister, 2015/11/16)
  </li>

  <li> New: Extra parameters to GD and Lanczos SLEPc solvers. Also added unit tests.
  <br>
  (Denis Davydov, 2015/11/09)
  </li>

  <li> Fixed: FETools::project_dg was adding the vector projection to
  the output vector. Now is the output vector initialized to zero.
  <br>
  (Adam Kosik, 2015/11/09)
  </li>

  <li> Fixed: A compilation issue with DEAL_II_INCLUDE_DIRS not used for
  compiling bundled boost.
  <br>
  (Lukas Korous, 2015/11/01)
  </li>

  <li> New: 2nd derivatives are implemented for PolynomialsBDM in 3D.
  <br>
  (Alistair Bentley, 2015/10/27)
  </li>

  <li> Fixed: PolynomialsBDM::degree() now returns the correct value.
  <br>
  (Alistair Bentley, 2015/10/24)
  </li>

  <li> New: Triangulation::set_all_manifold_ids_on_boundary(boundary_id, manifold_id)
  which sets the manifold_id for all parts of the boundary with a given boundary_id.
  <br>
  (Alberto Sartori, 2015/10/22)
  </li>

  <li> Fixed: The range vectors in the construction of an
  inverse_operator() is now reinitialised before solve calls. This ensures
  a consistent starting point for the solver.
  <br>
  (Jean-Paul Pelteret, 2015/10/19)
  </li>

  <li> New: Ghost cells for the multigrid levels in
  parallel::distributed::Triangulation are now correctly created also for
  periodic boundary conditions.
  <br>
  (Martin Kronbichler, 2015/10/18)
  </li>

  <li> Fixed: GridGenerator::subdivided_parallelepiped() produced
  invalid, unconnected meshes and wrong boundary indicators.
  <br>
  (Timo Heister, 2015/10/13)
  </li>

  <li> Improved: DoFTools::compute_intergrid_transfer_representation
  can now be used with a fine grid given by a parallel::Triangulation.
  <br>
  (Alexander Grayver, 2015/10/09)
  </li>

  <li> New: GridIn::read_unv() can now read more element codes that
  are used in typical meshes.
  <br>
  (Aslan Kosakian, 2015/10/06)
  </li>

  <li> New: FunctionParser now supports <code>pow(a,b)</code>.
  <br>
  (Timo Heister, 2015/09/30)
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

  <li> Improved: Allow continuation lines in ParameterHandler. Any line in a
  parameter file ending with a <tt>\\</tt> will now be combined with the next
  line; see ParameterHandler's documentation for more information.
  <br>
  (Alberto Sartori, 2015/09/04, David Wells, 2016/01/18-2016/01/28)
  </li>

  <li> New: There is now a function SparsityPattern::print_svg() which prints
  the sparsity of the matrix in a .svg file which can be opened in a web
  browser.
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

  <li> New: introduced hp::FECollection::find_least_face_dominating_fe(const
  std::set<unsigned int> &fes) which aims to find the least dominating finite
  element w.r.t. those provided as fe_indices in @p fes.
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

  <li> Fixed: VectorTools::integrate_difference() for VectorTools::Hdiv_seminorm
  was computed incorrectly.
  <br>
  (Timo Heister, 2015/08/31)
  </li>

  <li> New: Jacobian second and third derivatives are now computed by the mapping classes and can be
  accessed through FEValues in much the same way as the Jacobian and Jacobian gradient.
  <br>
  (Maien Hamed, 2015/08/28-2015/08/31)
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

  <li> Changed: The function Vector::ratio() and the corresponding
  functions in other vector classes have been deprecated.
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

  <li> New: The InterpolatedTensorProductGridData::gradient() function
  is now implemented.
  <br>
  (Daniel Shapero, 2015/08/12)
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
