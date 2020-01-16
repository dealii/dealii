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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/**
@page changes_between_8_4_1_and_8_4_2 Changes between Version 8.4.1 and 8.4.2

<p>
This is the list of changes made between the release of deal.II version
8.4.1 and that of 8.4.2. All entries are signed with the names of the
author.
</p>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
 <li> Fixed: Fix MPI_InitFinalize by correctly initializing and destroying
   all p4est/libsc related objects by calls to sc_init(), p4est_init(), and
   sc_finalize(); compatibility with p4est versions >1.1.
 <br>
 (Jonathan Perry-Houts, 2016/08/31)
 </li>

 <li> Fixed: The build system now uses -fPIC instead of -fpic
 <br>
 (Matthias Maier, 2016/08/31)
 </li>

 <li> Fixed: Linker error under VS2013
 <br>
 (Alex Kokomov, 2016/08/09)
 </li>

 <li> New: The library is now compatible with PETSc 3.7.0. Part of this change
 included adding a new header, <tt>petsc_compatibility.h</tt>, which provides
 some version-independent functions for using common PETSc functions.
 <br>
 (David Wells, 2016/07/07)
 </li>

 <li> Fixed: Fix various warnings introduced by GCC 6
 <br>
 (David Wells, 2016/05/07)
 </li>

</ol>

*/
