// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
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
@page changes_between_9_7_0_and_9_7_1 Changes between Version 9.7.0 and 9.7.1

<p>
This is the list of changes made since the last release of deal.II.
All entries are signed with the names of the authors.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="970-971-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: FESystem now populates the
  adjust_line_dof_index_for_line_orientation_table in 2d, which in turn fixes
  higher-order simplex elements in 2d.
  <br>
  (Arthur Bawin, David Wells, 2025/08/20)
 </li>

 <li>
  Fixed: The setup of faces of MatrixFree for `dim = 1` across
  periodic boundaries added some faces twice. This is
  now fixed.
  <br>
  (Martin Kronbichler, Peter Munch, Magdalena Schreter-Fleischhacker,
  Andreas Koch, 2025/07/25)
 </li>

</ol>

*/
