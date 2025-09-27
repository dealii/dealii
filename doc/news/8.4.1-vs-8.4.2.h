// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
@page changes_between_8_4_1_and_8_4_2 Changes between Version 8.4.1 and 8.4.2

<p>
This is the list of changes made between the release of deal.II version
8.4.1 and that of 8.4.2. All entries are signed with the names of the
author.
</p>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="841-842-specific"></a>
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
