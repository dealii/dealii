// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2023 by the deal.II authors
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
@page changes_between_9_5_0_and_9_5_2 Changes between Version 9.5.0 and 9.5.2

<p>
This is the list of changes made between the release of deal.II version
9.5.0 and that of 9.5.2. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: various configuration, compilation and instantiation issues have
  been addressed:
  <ol>
    <li>fix a compilation issue with bundled boost and modern clang compiler</li>
    <li>instantiate various MPI helper functions for signed long long int</li>
    <li>add missing codimension-one instantiation for some DoFTools functions</li>
    <li>fix compiling with PETSc with complex scalar type</li>
    <li>allow compilation with PETSc but without MPI</li>
    <li>fix compilation of bundled tbb with gcc-13</li>
    <li>detect Trilinos NOX support and guard Trilinos NOX usage</li>
  </ol>
  <br />
  (Daniel Arndt, Wolfgang Bangerth, Matthias Maier, 2024/01/31)
 </li>

</ol>

*/
