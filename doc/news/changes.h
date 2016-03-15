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
@page changes_after_8_4 Changes after Version 8.4.0

<p>
This is the list of changes made after the release of deal.II version
8.4.0. All entries are signed with the names of the authors.
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
 <li> New: Added unit tests for complex-valued PETSc and SLEPc.
 <br>
 (Toby D. Young, Denis Davydov, 2016/03/11)
 </li>

 <li> New: Added another scaling factor to Kelly error estimator, namely h_K.
 <br>
 (Denis Davydov, 2016/03/05)
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
