// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
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
@page changes_between_9_5_0_and_9_5_2 Changes between Version 9.5.0 and 9.5.2

<p>
This is the list of changes made between the release of deal.II version
9.5.0 and that of 9.5.2. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

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
