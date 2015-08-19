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

  <li> Improved: Some finite elements compute hessians analytically rather than
  by finite differencing. Namely, these are finite elements that are subclasses
  of FEPoly as well as FESystem with those as base elements.
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

</ol>

*/
