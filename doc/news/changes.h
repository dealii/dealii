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
@page changes_after_8_2 Changes after Version 8.2

<p>
This is the list of changes made after the release of deal.II version
8.2.0. All entries are signed with the names of the authors.
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
</ol>

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
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
  <li> Fixed: CMake now also handles and exports the subminor version
  number correctly ("pre" and "rc?" are replaced by "0").
  <br>
  (Matthias Maier, 2015/01/02)
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
