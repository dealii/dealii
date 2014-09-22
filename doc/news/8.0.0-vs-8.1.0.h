// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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
 * @page changes_between_8_0_and_8_1 Changes between Version 8.0 and 8.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
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
  <li>Changed: During the implementation of the 64-bit features for deal.II
  8.0, many linear algebra classes obtained a local typedef
  <code>size_type</code> indicating the integer type that is used to index
  into them. In some instances, this was accidentally set to
  <code>types::global_dof_index</code> (which may be a 64-bit data type)
  even in cases where this is clearly not going to work, for example for
  FullMatrix::size_type, since we will not be able to store full matrix
  objects of sizes for which a 32-bit index type is not sufficient. In
  these cases, the typedef was reverted to just <code>unsigned int</code>.
  <br>
  (Wolfgang Bangerth, 2013/12/04)
  </li>

  <li> Removed: With the switch of the testsuite to CMake, the old report_features
  and build test facilities are removed.
  <br>
  (Matthias Maier, 2013/12/03)
  </li>

  <li>
  Changed: The kinds of template arguments for the VectorTools::interpolate
  function taking a Mapping as first argument has changed. This was done to
  work around a bug in the Intel ICC compiler which led to linker errors. Since
  the actual function argument list remains unchanged, the only way you will
  notice this change is if you <i>explicitly</i> specify template arguments.
  The only place one would typically do that is if you take the address of
  a template function. Since this is not a common operation, the impact of this
  change is probably limited.
  <br>
  (Wolfgang Bangerth, 2013/11/27)
  </li>

  <li>
  Changed: The ghost handling of the parallel::distributed::Vector class has
  been reworked: The vector now carries a global state that stores whether
  ghost elements have been updated or not. If a vector has ghost elements, it
  does not allow calls to compress() any more. Instead, a compress operation
  can now only be done when the ghost entries have been cleared before by
  calling zero_out_ghosts() or operator=0. The state can be queried by the new
  method has_ghost_elements(). This change avoids spurious entries to be
  inserted with compress(), but requires some change in user codes. The
  behavior of a ghosted vector is now very similar to ghosted PETSc and
  Trilinos vectors. The only difference is that the <i>same</i> vector can
  also be used as a non-ghosted vector which is designed for use in assembly
  routines.
  <br>
  (Martin Kronbichler, 2013/10/18)
  </li>

  <li>
  Removed: GridTools::collect_periodic_face_pairs. This function is superseded
  by GridTools::collect_periodic_faces which exports an
  std::vector<PeriodicFacepair<...>> instead.
  <br>
  (Matthias Maier, 2013/09/30)
  </li>

  <li>
  Removed: The member function face_to_equivalent_cell_index() in
  FiniteElementData has been removed. It had been deprecated a while
  back already. Please use FiniteElement::face_to_cell_index() instead.
  <br>
  (Wolfgang Bangerth, 2013/08/09)
  </li>

  <li>
  Changed: The typedefs DataOut::cell_iterator and
  DataOut::active_cell_iterator were previously defined as
  DoFHandler::(active)_cell_iterator, while they are now
  Triangulation::(active)_cell_iterator. This is necessary to support DataOut
  on multiple DoFHandler objects. This affects possible overloading of
  DataOut::next_cell(cell_iterator). Use the typedef
  DataOut::(active)_cell_iterator as argument type instead.
  <br>
  (Martin Kronbichler, 2013/07/24)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
  <li> New: step-26 fills a long-standing gap: There was no tutorial
          program solving the heat equation. There was also no tutorial
          showing in relatively easy terms how to do adaptive meshes
          in time dependent problems. This program fills both of these
          needs.
  <br>
  (Wolfgang Bangerth, 2013/12/18)
  </li>

  <li> Improved: The build system now supports usage of the library
  out of the build directory without prior installation. This is done by
  exporting an additional project configuration just for the build directory.
  Furthermore, a bunch of convenience targets get now defined that just build
  individual components (such as just the documentation or the libraries), and
  if <tt>CMAKE_INSTALL_PREFIX</tt> is set, also install that specific component.
  <br>
  (Matthias Maier, Luca Heltai, 2013/12/03)
  </li>

  <li> Fixed: Missing instantiations of SparseDirectMUMPS have been added.
  <br>
  (Timo Heister, 2013/11/25)
  </li>

  <li> New: introduced "make test" that runs a minimal set of tests. We
  encourage every user to run this, especially if they run in to problems.
  The tests are automatically picked depending on the configuration and
  will be shipped with every release.
  <br>
  (Timo Heister, Matthias Maier, 2013/11/08)
  </li>

  <li> Changed: It is now possible to restore a parallel Triangulation
  (and solutions) with a different number of processors it was saved with
  using Triangulation::save() and Triangulation::load().
  <br>
  (Timo Heister, 2013/11/02)
  </li>

  <li> Added support for Windows: It is now possible again to use gcc on Windows
  in order to compile the library. We support gcc-4.8.1 on Cygwin64 and MinGW-w64.
  <br>
  (Matthias Maier, 2013/11/01)
  </li>

  <li> Changed: step-9, step-13 and step-14 have been converted to use the
  more modern WorkStream concept for assembling linear systems and computing
  error indicators in parallel.
  <br>
  (Bruno Turcksin, Wolfgang Bangerth, 2013/10/26)
  </li>

  <li> New: The testsuite is now ported to <a href="http://www.cmake.org/">
  CMake</a> and uses CTest as test driver.
  <br>
  (Wolfgang Bangerth, Timo Heister, Matthias Maier, Bruno Turcksin, 2013/10/20)
  </li>

  <li>
  Changed: multithreadinfo::n_default_threads is now deprecated. Use the
  new n_threads() function instead, which works correctly with TBB.
  <br>
  (Timo Heister, 2013/10/02)
  </li>

  <li>
  Changed: if configured with TBB but the number of threads is set to 1,
  do not bother to use TBB in workstream.
  <br>
  (Timo Heister, 2013/10/02)
  </li>

  <li>
  New: step-51 demonstrates the use of hybridized discontinuous Galerkin
  methods in deal.II, using the face elements FE_FaceQ. The programs solves a
  scalar convection-diffusion equation.
  <br>
  (Martin Kronbichler and Scott Miller, 2013/10/01)
  </li>

  <li>
  New: There is now an element FE_FaceP that can be combined with FE_DGP in
  hybridized DG methods.
  <br>
  (Martin Kronbichler, 2013/09/17)
  </li>

  <li>
  Fixed: The DataOutBase::XDMFEntry class now has a proper serialization
  function to allow for checkpointing.
  <br>
  (Eric Heien, 2013/09/27)
  </li>

  <li>
  New: DataOutBase::DataOutFilter provides a way to remove duplicate vertices
  and values from a solution vector when generating output. Currently it only
  supports HDF5/XDMF output.
  <br>
  (Eric Heien, 2013/09/27)
  </li>

  <li>
  Removed: DataOutBase::HDF5MemStream was removed and the functionality replaced
  by DataOutBase::DataOutFilter. The user only manipulates these through
  DataOutBase::write_hdf5_parallel so this change should be transparent.
  <br>
  (Eric Heien, 2013/09/27)
  </li>

  <li>
  New: Like the usual DoFHandler class, the hp::DoFHandler class now also
  has a cache that makes operations such as <code>cell-@>get_dof_indices(...)</code>
  faster. This should accelerate many parts of the library that deal with
  hp finite elements.
  <br>
  (Wolfgang Bangerth, 2013/09/10)
  </li>

  <li>
  New: parallel::distributed::Triangulation now supports periodic boundaries.
  DoFTools::make_periodicity_constraints and similar functions are now working
  on parallel::distributed::Triangulation objects.
  <br>
  (Tobin Isaac, Craig Michoski, Daniel Arndt, 2013/09/06)
  </li>

  <li>
  New: It is now possible to compile and link deal.II against LLVM's libcxx. For
  this, a few issues with C++ standard violations are resolved.
  <br>
  (Matthias Maier, 2013/08/09)
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
  <li> Fixed: The DerivativeApproximation class did not work for
  parallel programs. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/12/18)
  </li>

  <li> Fixed: Move the implementation of Subscriptor::(un)subscribe() to
  the .cc file so that it is possible to link against the debug library
  without specifying <code>-DDEBUG</code>
  <br>
  (Wolfgang Bangerth, 2013/12/13)
  </li>

  <li> Fixed: Since the introduction of ThreadLocalStorage in version 8.0, the
  way in which FEValues objects visit cells in a parallel assembly loop is no
  longer deterministic. Therefore, the detection of CellSimilarity that can
  speed up computations of certain geometric quantities (shape gradients) on
  cells that are translations is disabled when the number of threads is
  greater than one. This produces somewhat slower code (usually not more than
  a few percent) but ensures exact reproducibility of results.
  <br>
  (Martin Kronbichler, Wolfgang Bangerth, 2013/12/09)
  </li>

  <li> Fixed: Several functions in namespace GridTools were not instantiated
  for parallel::distributed::Triangulation objects. This is now fixed.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2013/12/01)
  </li>

  <li> Improved: The methods ConstraintMatrix::distribute_local_to_global
  now use scratch data that is private to each thread instead of allocating
  it for every cell anew. This gives better performance, in particular in
  parallel, of these operations, while maintaining thread-safety (when
  accessing non-overlapping rows, no race condition can exist).
  <br>
  (Martin Kronbichler, 2013/12/03)
  </li>

  <li> Improved: When attempting operations such as FEValues::get_function_values()
  or FEValues::shape_value(), the FEValues object needs to know that what these
  functions return has been computed previously. What is computed is specified
  by the update flags that are passed to the constructor of all FEValues, FEFaceValues
  and FESubfaceValues objects. If a user attempts an operation for which the
  corresponding flag was not specified, an exception is generated. This exception
  did say previously what the cause was, but it was not overly explicit.
  The exception now generates a message that says exactly what is going wrong.
  <br>
  (Wolfgang Bangerth, 2013/12/01)
  </li>

  <li> Fixed: GridGenerator::truncated_cone() failed if half_length < 0.5*radius in 3d.
  <br>
  (Timo Heister, 2013/11/25)
  </li>

  <li> Fixed: make_hanging_node_constraints failed with an exception in a
  parallel::distributed computation if the element is
  RaviartThomas (and probably others).
  <br>
  (Timo Heister, 2013/11/23)
  </li>

  <li> Improved: CMake: Added a configuration check for incompatible ninja
  + icc setup, fixed several setup and performance issues with the
  testsuite.
  <br>
  (Matthias Maier, 2013/11/20)
  </li>

  <li> Changed: when a dealii::Exception is thrown, defer the symbol lookup of the
  stack trace to when it is needed. This improves performance if what() is never
  called.
  <br>
  (Timo Heister, 2013/11/17)
  </li>

  <li> Fixed: GridGenerator::parallelogram was not instantiated properly
  when using intel compilers.
  <br>
  (Timo Heister, 2013/11/17)
  </li>

  <li>
  Fixed: MappingQ1::transform_real_to_unit_cell() could fail in some
  cases with very elongated and twisted cells. This should now be fixed
  with an algorithm that uses a better method of computing the Newton
  convergence.
  <br>
  (Wolfgang Bangerth, 2013/11/17)
  </li>

  <li>
  Fixed: VectorTools::compute_no_normal_flux_constraints had a bug that
  only appeared in rare cases at vertices of the domain if one adjacent
  cell had two boundary indicators selected for no normal flux and another
  had only one. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/11/17)
  </li>

  <li> Fixed: dealii::FETools::interpolation_difference was
  not working for TrilinosWrappers::MPI::Vectors with ghost
  entries. The TrilinosWrappers::VectorBase class has now a
  get_mpi_communicator method similar to the PETSc vector
  classes.
  <br>
  (Martin Steigemann, Martin Kronbichler, 2013/11/17)
  </li>

  <li> Fixed: Bundled fparser is now compiled with FP_USE_THREAD_SAFE_EVAL in
  case of enabled threading support so that it is thread safe.
  <br>
  (Matthias Maier, reported by Francesco Cattoglio 2013/11/16)
  </li>

  <li> Fixed: The CellData class now has a default constructor that
  sets the material and boundary indicators to zero. This fixes certain
  internal errors with the Intel ICC compiler.
  <br>
  (Wolfgang Bangerth, 2013/11/13)
  </li>

  <li> Cleanup: Removed obsolete files and files with unknown licensing
  status from the source tree. Along the way, parameter_gui now uses
  default icons from the desktop environment instead of bundled ones.
  <br>
  (Matthias Maier, 2013/11/11)
  </li>

  <li> New: There is now a framework for coloring graphs, with functions
  in namespace GraphColoring.
  <br>
  (Bruno Turcksin, Martin Kronbichler, 2013/11/06)
  </li>

  <li>
  Fixed: the DerivativeApproximation class was not working correctly when
  used with parallel vectors.
  (Timo Heister, 2013/10/28)
  </li>

  <li>
  ~Subscriptor and ~GrowingVectorMemory no longer throw an exception (the
  former if disable_abort_on_exception was called) to be compatible with the
  C++11 standard which otherwise requires the program to immediately call
  std::terminate. This was done with a new macro "AssertNothrow".
  <br>
  (Wolfgang Bangerth, Matthias Maier, Bruno Turcksin 2013/10/22)
  </li>

  <li>
  dealii::SolverControl::NoConvergence now inherits dealii::ExceptionBase and
  is thrown via AssertThrow(false, ... ).
  <br>
  (Matthias Maier, 2013/10/20)
  </li>

  <li>
  New: parallel::distributed::BlockVector has now methods update_ghost_values,
  compress, zero_out_ghosts, and has_ghost_elements that do the respective
  operation on each block of parallel::distributed::Vector.
  <br>
  (Martin Kronbichler, 2013/10/18)
  </li>

  <li>
  Fixed: When deriving from DataOut to filter the cells where output is generated, there were two different bugs that result in segmentation faults or wrong cells written (example, step-18).
  <br>
  (Timo Heister, 2013/10/16)
  </li>

  <li>
  New: GridIn::read_vtk() reads 2d and 3d meshes in VTK format.
  <br>
  (Mayank Sabharwal, Andreas Putz, 2013/10/07)
  </li>

  <li>
  Fixed: ConstraintMatrix would not compress() the IndexSet in the constructor
  leading to crashes that only happen in release mode. This is now fixed.
  <br>
  (Timo Heister, 2013/09/27)
  </li>

  <li>
  Fixed: PetscWrappers::MatrixBase::row_length() no longer worked after recent changes
  to PETSc (around PETSc release 3.4). This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/09/24)
  </li>

  <li>
  New: Added write_visit_record that allows writing .visit files with multiple blocks
  and multiple time steps.
  <br>
  (Fahad Alrashed, 2013/09/21)
  </li>

  <li>
  Changed: GridTools::have_same_coarse_mesh was only instantiated for
  MGDoFHandler arguments in debug mode. This is now fixed.
  <br>
  (Timo Heister, 2013/09/20)
  </li>

  <li>
  Changed: GridTools::find_active_cell_around_point now throws the exception
  GridTools::ExcPointNotFound
  if the point is outside the Triangulation. This exception can be caught.
  <br>
  (Timo Heister, 2013/09/18)
  </li>

  <li>
  Changed: we now call MPI_Init_thread instead of MPI_Init.
  <br>
  (Timo Heister, 2013/09/17)
  </li>

  <li>
  Enhancement: It is now possible to use the build directory directly without
  the need to install first. For this, a second copy of all necessary project
  configuration files (deal.IIConfig.cmake, ..., Make.global_options) are
  generated and deployed in the build directory. (This is fully compatible with
  the old possibility to install into the build dir.)
  <br>
  (Matthias Maier, 2013/09/15)
  </li>

  <li>
  Fixed: DoFTools::extract_locally_*_dofs now instantiated for hp::DofHandler.
  <br>
  (Jean-Paul Pelteret, 2013/09/11)
  </li>

  <li>
  Changed: distributed::parallel:BlockVector::operator= now allows importing
  of ghost values like all other vector types. Also added some new constructors
  for BlockVector and Vector using IndexSets to mirror the other linear algebra
  classes.
  <br>
  (Timo Heister, 2013/09/04)
  </li>

  <li>
  Fixed: VectorTools::compute_no_normal_flux_constraints had a bug that
  only manifested on complex meshes. This is now fixed.
  <br>
  (Chih-Che Chueh, Wolfgang Bangerth, 2013/09/04)
  </li>

  <li>
  New: All vector classes now have functions <code>extract_subvector_to()</code>
  that allow extracting not just a single value but a whole set.
  <br>
  (Fahad Alrashed, 2013/09/02)
  </li>

  <li>
  Fixed: <code>common/Make.global_options</code> now exports enable-threads
  correctly, furthermore, <code>lib-suffix</code>, <code>shared-lib-suffix</code>
  and <code>static-lib-suffix</code> are now exported as well for better legacy
  support.
  <br>
  (Matthias Maier, 2013/08/30)
  </li>

  <li>
  New: The ParameterHandler class can now deal with including one parameter
  file from another.
  <br>
  (Wolfgang Bangerth, 2013/08/25)
  </li>

  <li>
  New: The method VectorTools::compute_normal_flux_constraints can be used to
  force a vector finite element function to be normal to the boundary.
  <br>
  (Martin Kronbichler, 2013/08/23)
  </li>

  <li>
  Improved: MappingQ now uses the points of the Gauss-Lobatto quadrature
  formula as support points instead of equispaced ones. This allows its use
  for high polynomial orders and also gives better interpolation of circular
  boundaries. Beware that mappings of order three and higher will behave
  slightly differently now (usually better).
  <br>
  (Martin Kronbichler, 2013/08/23)
  </li>

  <li>
  Improved: Several .cc files in the deal.II directory have been split in
  order to better utilize multiple processors when compiling in parallel and
  reduce memory requirements of the compilation stage.
  <br>
  (Martin Kronbichler, 2013/08/22)
  </li>

  <li>
  Fixed: The ParameterHandler::declare_entry() did not check that the
  default value of a parameter indeed satisfies the pattern given for this
  parameter (despite a statement in the documentation that this checking
  would happen). This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/08/21)
  </li>

  <li>
  New: Patterns::List and Patterns::Map now accept a string
  different than the default comma that denotes the separator
  between entries of the list or map.
  <br>
  (Wolfgang Bangerth, 2013/08/21)
  </li>

  <li>
  Fixed: Some operations in the MappingQ class are now done in higher
  precision arithmetic to mitigate the ill-conditioning that appears
  when using mappings of high order (say, order 6 or 8 or 10).
  <br>
  (Juan Carlos Araujo Cabarcas, 2013/08/20)
  </li>

  <li>
  Fixed: The SLEPcWrappers classes could not be compiled for 64-bit
  indices. This is now fixed.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2013/08/20)
  </li>

  <li>
  Fixed: SolutionTransfer used to crash whenever one transfered in the hp
  context between cells that use FE_Nothing and FE_Q. This is now fixed.
  <br>
  (Krzyszof Bzowski, Wolfgang Bangerth, 2013/08/18)
  </li>

  <li>
  Fixed: Under some circumstances (see http://code.google.com/p/dealii/issues/detail?id=82)
  the DoFTools::make_periodicity_constraints() function could create cycles in
  the ConstraintMatrix object. This is now fixed.
  <br>
  (David Emerson, Wolfgang Bangerth, 2013/08/16)
  </li>

  <li>
  New: There is now a function ConstraintMatrix::are_identity_constrained().
  <br>
  (Wolfgang Bangerth, 2013/08/16)
  </li>

  <li>
  New: TableHandler::write_text() now also supports output in
  org-mode (http://orgmode.org/) format via a new entry in the
  TableHandler::TextOutputFormat enumeration.
  <br>
  (Oleh Krehel, 2013/08/15)
  </li>

  <li>
  New: There are now global functions <code>scalar_product</code>
  that compute the scalar product (double contraction) between
  tensors of rank 2.
  <br>
  (Scott Miller, 2013/08/14)
  </li>

  <li>
  Fixed: Creating objects of type MappingQ was previously only possible
  for low order polynomials. For orders higher than around 6, one ran
  into assertions that tested for internal consistency. These assertions
  have now been appropriately relaxes for the growth of round-off errors
  with growing polynomial degrees.
  <br>
  (Juan Carlos Araujo Cabarcas, Wolfgang Bangerth, 2013/08/14)
  </li>

  <li>
  New: MappingQEulerian is now also instantiated for vector elements
  of type TrilinosWrappers::Vector as well as the MPI and block
  variants.
  <br>
  (Armin Ghajar Jazi, 2013/08/14)
  </li>

  <li>
  Fixed: The FiniteElement::face_to_cell_index() function had a bug
  that made it work incorrectly for elements that have more than one
  degree of freedom per line (in 2d) or per quad (in 3d). This is now
  fixed for the most common cases, namely the FE_Q elements as well
  as elements composed of FESystem elements. For all other cases, an
  exception is generated reporting that this case is not implemented.
  If you run into this, let us know.
  <br>
  (Wolfgang Bangerth, 2013/08/10)
  </li>

  <li>
  New: DataOutBase::VtkFlags now has a flag
  DataOutBase::VtkFlags::print_date_and_time that can be used to suppress output
  of date and time in output files. This is useful in test suites where a newer
  run at a different time produces differences against previously stored files,
  even though the actual data is exactly the same.
  <br>
  (Oleh Krehel, 2013/08/06)
  </li>

  <li>
  Fixed: The various block matrix classes are all derived from BlockMatrixBase
  which had race conditions when the set() or add() functions were called from
  different threads. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/08/05)
  </li>

  <li>
  Fixed: various fixes with assignment and reinit of PETScWrappers::MPI::Vector.
  <br>
  (Timo Heister, 2013/08/05)
  </li>

  <li>Fixed: An assertion wrongly triggered in
  DoFTools::make_hanging_node_constraints when used with a particular
  combination of FESystem elements containing FE_Nothing. This is now fixed.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2013/08/01)
  </li>

  <li>
  New: Add has_ghost_elements() for PETScWrappers::MPI::BlockVector and
  TrilinosWrappers::MPI::BlockVector.
  <br>
  (Timo Heister, 2013/08/01)
  </li>

  <li>
  SparsityTools::distribute_sparsity_pattern did not work correctly for
  block systems, this has been fixed (function has a different signature).
  <br>
  (Timo Heister, 2013/07/31)
  </li>

  <li>Fixed: When typing <code>make run</code> in the step-32 directory,
  the program was executed with <code>mpirun -np 2 ./step-32</code>. This
  assumes that a program <code>mpirun</code> exists, but also does that
  deal.II was in fact compiled with MPI support on. Neither was intended.
  This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>

  <li>New: The DataOut, DataOutFaces, and DataOutRotation classes now allow
  the output of data vectors using different DoFHandler objects (based on the
  same triangulation), by new functions add_data_vector. This is used in the
  step-31 tutorial program which avoids creating a joint DoFHandler just for
  output.
  <br>
  (Martin Kronbichler, 2013/07/24)
  </li>

  <li>Changed: GridGenerator used to be a class with only static members
  but is now a namespace, like all other similar constructs in deal.II.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>

  <li>Changed: In GridGenerator, several functions had erroneously been changed
  to take an argument of type <code>size_type</code> rather than <code>unsigned
  int</code>. <code>GridGenerator::size_type</code> was a typedef to
  types::global_dof_index, which for most users was <code>unsigned int</code>
  anyway, but could also be set to be a 64-bit integer type. In any case, the
  change has been reverted and these functions take just a regular
  <code>unsigned int</code> again.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>
</ol>


*/
