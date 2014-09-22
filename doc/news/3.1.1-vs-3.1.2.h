// ---------------------------------------------------------------------
//
// Copyright (C) 2013, 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

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
       href="http://gcc.gnu.org/cgi-bin/gnatsweb.pl?database=gcc&user=guest&password=guest&cmd=login"
       target="_top">gcc bug tracking system page</a>, as bugs reports c++/615 and
       optimization/2938.
       <br>
       (WB 2001/06/27)
       </p>
</ol>


*/
