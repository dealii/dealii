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
@page changes_between_9_3_2_and_9_3_3 Changes between Version 9.3.2 and 9.3.3

<p>
This is the list of changes made between the release of deal.II version
9.3.2 and that of 9.3.3. All entries are signed with the names of the
author.
</p>

<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: use correct update flags in VectorTools::point_values().
  <br />
  (Maximilian Bergbauer, 2022/01/24)
 </li>

 <li>
  Fixed: use portable "cmake -E remove" in quick_tests.
  <br />
  (Matthias Maier, 2021/12/26)
 </li>

 <li>
  Fixed: Compatibility issues with Boost 1.78 have been resolved.
  <br />
  (David wells, 2022/01/04)
 </li>

 <li>
  Fixed: deal.II has been updated to support the new OneAPI api interface
  introduced by the Intel Threading Building Blocks Library. Note that the
  threading support of the matrix-free backend will be disabled in this case:
  For now, the MatrixFree::AdditionalData::tasks_parallel_scheme control has
  no effect and simply defaults to the serial loop.
  <br />
  (Wolfgang Bangerth, Matthias Maier, 2022/01/06)
 </li>

 <li>
  Fixed: A compilation issue with sundials has been resolved.
  <br />
  (Timo Heister, 2021/12/25)
 </li>

 <li>
  Fixed: The shell script `p4est-setup.sh` for installing p4est has been
  updated for its latest releases (p4est>=2.3).
  <br />
  (Marc Fehling, 2021/06/10)
 </li>

</ol>

*/
