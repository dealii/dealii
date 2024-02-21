// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
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
 * @page changes_between_3_1_1_and_3_1_2 Changes between Version 3.1.1 and 3.1.2

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>

<ol>
  <li> <p>
       The only is change is that the configuration detects the
       presence of gcc3.0 and works around a bug in it by compiling
       the file <code>tria_accessor.cc(3d)</code> without optimization
       flags even for optimized libraries. The gcc bug is documented
       at the <a
       href="https://gcc.gnu.org/cgi-bin/gnatsweb.pl?database=gcc&user=guest&password=guest&cmd=login"
       target="_top">gcc bug tracking system page</a>, as bugs reports c++/615 and
       optimization/2938.
       <br>
       (WB 2001/06/27)
       </p>
</ol>


*/
