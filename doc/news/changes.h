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

  <li> Changed: The signature of the FiniteElement::fill_fe_values(),
  FiniteElement::fill_fe_face_values(), and FiniteElement::fill_fe_subface_values()
  functions has been changed, in an effort to clarify which of these contain
  input information and which contain output information for these functions.
  The same has been done for the corresponding functions in the Mapping
  class hierarchy. As part of a general overhaul, the FEValuesData class
  has also been removed.
  <br>
  (Wolfgang Bangerth, 2015/07/20-2015/08/06)
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

</ol>

*/
