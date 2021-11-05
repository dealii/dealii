// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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
@page changes_between_9_3_1_and_9_3_2 Changes between Version 9.3.1 and 9.3.2

<p>
This is the list of changes made between the release of deal.II version
9.3.1 and that of 9.3.2. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: Various compatiblity issues and minor bugfixes have been resolved:
  <ol>
  <li>fixes a Sundials/Kinsol issue</li>
  <li>fixes a bug for plain_copy in MGTransferMatrixFree</li>
  <li>fixes a bug for Checkpointing/Restore of large triangulations greater than 4GB</li>
  <li>fixes compatibility issues with Intel 18/19 compiler</li>
  <li>fixes a compilation issue with step-79 and Intel 18 compiler</li>
  </ol>
  <br>
  (Daniel Arndt, Marc Fehling, Timo Heister, Matthias Maier, Peter Munch, David Wells, 2021/11/04)
 </li>

</ol>

*/
