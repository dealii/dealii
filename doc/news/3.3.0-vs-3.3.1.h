// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
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
 * @page changes_between_3_3_0_and_3_3_1 Changes between Version 3.3.0 and 3.3.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


<ol>
  <li> <p>
       Fixed: In 3d, the function <code
       class="member">DoFTools::make_hanging_node_constraints</code>
       contained an assertion that failed erroneously for finite
       elements that do not have degrees of freedom on vertices. This
       is now fixed.
       <br>
       (WB 2002/02/21)
       </p>

  <li> <p>
       Fixed: <code>TriaAccessor<3,3>::measure</code>
       sometimes computed a negative value. This is now fixed.
       <br>
       (WB 2002/02/21)
       </p>
</ol>

*/
