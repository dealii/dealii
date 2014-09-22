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
 * @page changes_between_3_0_0_and_3_0_1 Changes between Version 3.0.0 and 3.0.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>

<ol>
  <li> <p>
       Fix: in three space dimensions, the triangulation class
       over-estimated the necessary amount of memory needed upon
       refinement and allocated too much.
       <br>
       (WB 2000/04/19)
       </p>

  <li> <p>
       Fix: getting the support points from system finite elements
       (class <code>FESystem</code>) resulted in an
       exception. 
       <br>
       (WB 2000/05/10)
       </p>

  <li> <p>
       Fix: <code>FullMatrix::Tmmult</code> and
       <code>FullMatrix::Tvmult</code> were broken.
       <br>
       (WB 2000/05/08 and 2000/05/26)
       </p>

  <li> <p>
       Fix: slight bug in 
       <code>DataOut::build_patches</code>
       in multithreaded mode fixed.
       <br>
       (Ralf Hartmann, 2000/05/29)
       </p>

  <li> <p>
       Fix: 
       <code>SparsityPattern::print_gnuplot</code>
       wrote rows and columns exchanged. Since most matrices have
       symmetric sparsity patterns, this has gone unnoticed by now.
       <br>
       (WB 2000/05/30)
       </p>

  <li> <p>
       Fixed: the 
       <code class="program">common/scripts/make_dependencies.pl</code> 
       script that sets up the dependencies for the make files had a
       problem when the path to the library included special characters
       such as `+'. This is now fixed.
       <br>
       (WB 2000/06/15)
       </p>
</ol>

*/
