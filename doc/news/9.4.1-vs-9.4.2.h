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
@page changes_between_9_4_1_and_9_4_2 Changes between Version 9.4.1 and 9.4.2

<p>
This is the list of changes made between the release of deal.II version
9.4.1 and that of 9.4.2. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: Various compatiblity issues and minor bugfixes have been resolved:
  <ol>
  <li>A compilation issue with step-70 has been resolved</li>
  <li>CMake: prefer <code>-pthread</code> for posix thread support</li>
  <li>Fix a type mismatch for suitesparse that leads to compilation failures on certain platforms</li>
  <li>A number of Microsoft Visual Code compatibility fixes concerning extern declarations</li>
  </ol>
  <br>
  (Wolfgang Bangerth, David Wells, Matthias Maier, 2023/02/06)
 </li>

</ol>

*/
