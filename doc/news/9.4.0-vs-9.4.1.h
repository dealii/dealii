// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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
@page changes_between_9_4_0_and_9_4_1 Changes between Version 9.4.0 and 9.4.1

<p>
This is the list of changes made between the release of deal.II version
9.4.0 and that of 9.4.1. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: Various compatiblity issues and minor bugfixes have been resolved:
  <ol>
  <li>cmake: always export compile_commands.json in deal.II and user projects</li>
  <li>doxygen: fix various errors in formulas</li>
  <li>doxygen: fix SymmetricTensor friends</li>
  <li>cmake: fix PETSc version detection</li>
  <li>base: fix some VectorizedArrayTypes for non-default vectorization</li>
  <li>gitignore: ignore clangd files and directories</li>
  <li>change ConsensusAlgorithm deprecations to early deprecated</li>
  <li>step-81: Mention example step in the tutorial lists</li>
  <li>bugfix: use correct tolerance in MappingCartesian check</li>
  </ol>
  <br>
  (David Wells, Matthias Maier, Rene Gassmoeller, Sebastian Proell, Timo Heister, 2022/12/01)
 </li>

</ol>

*/
