// ---------------------------------------------------------------------
// $Id$
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
// * @page changes_after_8_1 Changes after Version 8.1

<p>
This is the list of changes made after the release of
deal.II version 8.1.0.
All entries are signed with the names of the authors.
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
  <li> Fixed: The DerivativeApproximation class did not work for
  parallel programs. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/12/18)
  </li>

  <li> Fixed: VectorTools::project_boundary_values could not deal with
  function values close to (but not exactly equal to) zero. This is now fixed.
  <br>
  (Martin Kronbichler, 2013/12/16)
  </li>

  <li> New: It is now possible to select between different smoothers and coarse
  solvers in the Trilinos AMG preconditioners by a string to the smoother's name.
  <br>
  (Andrew Baker, 2013/12/14)
  </li>

</ol>


*/
