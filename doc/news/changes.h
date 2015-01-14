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
@page changes_after_8_2_1 Changes after Version 8.2.1

<p>
This is the list of changes made after the release of deal.II version
8.2.1. All entries are signed with the names of the authors.
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
  <li> Removed: TrilinosWrappers::SparseMatrix copy constructor got removed
  to be in line with PETSc and dealii::SparseMatrix. You can use reinit()
  and copy_from().
  <br>
  (Timo Heister, 2015/01/12)
  </li>

  <li> Removed: The following compatibility definitions were removed from
  <code>include/deal.II/base/config.h.in</code> (replacement in brackets):
  - DEAL_II_CAN_USE_CXX11 (new: DEAL_II_WITH_CXX11)
  - DEAL_II_CAN_USE_CXX1X (new: DEAL_II_WITH_CXX11)
  - DEAL_II_COMPILER_SUPPORTS_MPI (new: DEAL_II_WITH_MPI)
  - DEAL_II_MAJOR (new: DEAL_II_VERSION_MAJOR)
  - DEAL_II_MINOR (new: DEAL_II_VERSION_MINOR)
  - DEAL_II_USE_ARPACK (new: DEAL_II_WITH_ARPACK)
  - DEAL_II_USE_CXX11 (new: DEAL_II_WITH_CXX11)
  - DEAL_II_USE_METIS (new: DEAL_II_WITH_METIS)
  - DEAL_II_USE_MT (new: DEAL_II_WITH_THREADS)
  - DEAL_II_USE_P4EST (new: DEAL_II_WITH_P4EST)
  - DEAL_II_USE_PETSC (new: DEAL_II_WITH_PETSC)
  - DEAL_II_USE_SLEPC (new: DEAL_II_WITH_SLEPC)
  - DEAL_II_USE_TRILINOS (new: DEAL_II_WITH_TRILINOS)
  <br>
  (Matthias Maier, 2015/01/12)
  </li>

  <li> Removed: The direct Mumps interface through
  <code>SparseDirectMUMPS</code> has been removed. The MUMPS solver is
  still available through the Trilinos or PETSc interfaces. Alternatively,
  there is <code>SparseDirectUMFPACK</code>, which has a similar interface.
  <br>
  (Matthias Maier, 2015/01/11)
  </li>

  <li> Removed: This release removes a number of functions that have long
  been deprecated and that were previously already marked as
  deprecated (i.e., they would have yielded warnings by the compiler whenever
  you tried to use them). In almost all cases, there is a function with same
  name but different argument list that should be used instead.
  Specifically, the removed functions and classes are:
  - TimeDependent::end_sweep (with an argument).
  - PointValueHistory::mark_locations.
  - The DataPostprocessor::compute_derived_quantities_scalar and
    DataPostprocessor::compute_derived_quantities_vector functions without
    evaluation points. If you have
    data postprocessor classes implemented in your program that overload these
    functions, you will have to change it in a way that they overload the
    functions of same name but with the evaluation point argument instead.
  - The constructors of classes MGSmoother, MGSmootherRelaxation and
    MGSmootherPrecondition that take a VectorMemory object.
  - Deprecated variants of MeshWorker::loop and MeshWorker::integration_loop.
  - ThreadManagement::spawn.
  - Threads::ThreadCondition and Threads::ThreadMutex.
  - GridGenerator::laplace_transformation.
  - The version of GridGenerator::parallelogram where the corners are given
    as a rank-2 tensor rather than as an array of points.
  - DataOutBase::create_xdmf_entry with 3 arguments.
  - DataOutBase::write_hdf5_parallel with 2 arguments.
  - Algorithms::ThetaTimestepping::operator().
  - Algorithms::ThetaTimestepping::initialize.
  - Algorithms::Newton::initialize.
  - MGLevelObject::get_minlevel and MGLevelObject::get_maxlevel.
  - The versions of FunctionParser::initialize that took a
    <code>use_degrees</code> or <code>constants</code> argument.
    The implementation as it is now no longer supports either of
    these two concepts (since we switched from the FunctionParser
    library to the muparser library after the deal.II 8.1 release).
  - DoFRenumbering::downstream_dg.
  - DoFTools::count_dofs_per_component.
  - DoFTools::make_sparsity_pattern with a vector-of-vector mask.
  - GridOutFlags::XFig::level_color.
  - class BlockList.
  - MGConstrainedDoFs::non_refinement_edge_index
  - MGConstrainedDoFs::at_refinement_edge_boundary
  - The refinement listener concept of the Triangulation class. This
    approach to getting notified about what happens to triangulations
    has been superseded by the signals defined by the triangulation
    class.
  - The deprecated constructor of MPI_InitFinalize
  - The MPI support functions in namespace Utilities and Utilities::System.
  - Deprecated members of namespace types.
  - Namespace deal_II_numbers.
  - Triangulation::distort_random.
  - Triangulation::clear_user_pointers.
  - The deprecated constructor of SparseILU.
  - SparseILU::apply_decomposition.
  - The deprecated constructor of SparseMIC.
  - The compress() functions without argument in the various vector
    classes. You should use the versions with a VectorOperation
    argument instead.
  - In FEValues and related classes, the functions that contain the
    term <code>2nd_derivatives</code> were removed in favor of those
    with names containing <code>hessian</code>. Similarly, functions
    with names including <code>function_grads</code> were removed in
    favor of those called <code>function_gradients</code>. Finally,
    the <code>cell_normal_vector</code> functions were replaced by
    <code>normal_vector</code> ones. In all cases, the new functions
    have been around for a while.
  - Vector::scale.
  - TrilinosWrappers::*Vector*::compress with an Epetra_CombineMode
    argument
  <br>
  This release also removes the deprecated class MGDoFHandler. The
  functionality of this class had previously been incorporated into
  the DoFHandler class. Unlike the changes above, if you were still
  using this class, you will need to do the following changes to
  your code:
  - Where you called <code>mg_dof_handler.distribute_dofs()</code>
    you now also need to explicitly call
    <code>mg_dof_handler.distribute_mg_dofs()</code>.
  - If you called <code>mg_dof_handler.begin(level)</code>, you
    will now have to write this as
    <code>mg_dof_handler.begin_mg(level)</code> to make clear that
    you are not just interested in an iterator to a cell on a given
    level, but in fact to a cell that can access the degrees of
    freedom on a particular level of a multigrid hierarchy.
  - The type previously referred to as
    <code>MGDoFHandler::cell_iterator</code> now corresponds to
    <code>MGDoFHandler::level_cell_iterator</code>.
  - Where you previously called DoFRenumbering::component_wise
    for the entire MGDoFHandler object, you now need to call
    this function for the DoFHandler object, and then call the
    same function with the <code>level</code> argument for each
    of the levels of the triangulation individually.
  <br>
  (Wolfgang Bangerth, 2014/12/29-2015/01/11)
  </li>

  <li> Removed: The config.h file no longer exports HAVE_* definitions.
  Those are either entirely removed (for the blas/lapack symbols) or
  renamed to DEAL_II_HAVE_*. This change is done in order to avoid clashes
  with external projects also exporting HAVE_* definitions in their header
  files.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
  <li> New: The build system now queries for git branch name and
  revision sha1 (and automatically reconfigures if necessary). This
  information is used to annotate summary.log and detailed.log with the
  current revision sha1. Further, a header file <deal.II/base/revision.h>
  is now available that exports the macros: DEAL_II_GIT_BRANCH,
  DEAL_II_GIT_REVISION, DEAL_II_GIT_REVISION_SHORT.
  <br>
  (Matthias Maier, 2015/01/02)
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
  <li> New: dealii::multithread_info.n_cpus returns the correct number of CPU 
  on FreeBSD.
  <br>
  (Bruno Turcksin, 2015/01/14)
  </li>

  <li> Improved: MPI collective operations such as MPI::sum, MPI::max now
  check for job_supports_mpi() internally, which allows running them also
  without a call to MPI_Init.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> Changed: The method job_supports_mpi() now resides in the namespace
  Utilities::MPI instead of Utilities::System for consistency with other MPI
  methods. The old method has been marked deprecated and will be removed in
  a future version.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> Fixed: The update of ghost values in parallel::distributed::Vector when
  calling the assignment operator is now active when one of the two vector had
  its ghost values updated before or when the layout of the right hand side
  vector is one-to-one, more consistent with parallel PETSc and Trilinos
  vectors.
  <br>
  (Martin Kronbichler, 2015/01/13)
  </li>

  <li> New: PETScWrappers::MPI::SparseMatrix::reinit(other) copies
  the layout of another matrix. TrilinosWrappers::SparseMatrix
  operator= and copy constructor are now disabled. This brings
  functionality between PETSc and Trilinos in line.
  <br>
  (Timo Heister, 2015/01/12)
  </li>

  <li> New: Triangulation::set_all_manifold_ids() and
  Triangulation::set_all_manifold_ids_on_boundaries()
  set all manifold ids on every object or on every
  boundary object respectively.
  <br>
  (Luca Heltai, 2015/01/12)
  </li>

  <li> New: GridTools::copy_boundary_to_manifold_id() and
  GridTools::copy_material_to_manifold_id() copy
  boundary_ids and material_ids to manifold_ids for
  faces on the boundary and for cells respectively.
  <br>
  (Luca Heltai, 2015/01/09)
  </li>

  <li> Fixed: Utilities::int_to_string() produced wrong results if
  the number of digits specified was ten or greater.
  <br>
  (David Wells, 2015/01/08)
  </li>

  <li> Fixed: VectorTools::interpolate_to_different_mesh() was accidentally
  only instantiated for dealii::Vector arguments, rather than all vector
  classes. This is now fixed.
  <br>
  (Benjamin Brands, Wolfgang Bangerth, 2014/12/29)
  </li>

  <li> Fixed: Use CASROOT environment variable as additional hint for
  opencasacade.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>

  <li> Fixed: Update the run_testsuite.cmake script to also pick up
  muparser and opencascade configuration.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>

  <li> Fixed: Update several places in the documentation that were not
  updated from functionparser to muparser. Add several forgotten
  DEAL_II_WITH_* variables to certain places in the documentation.
  <br>
  (Matthias Maier, 2014/12/29)
  </li>
</ol>

*/
