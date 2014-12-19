// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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
 * @page changes_between_6_3_0_and_6_3_1 Changes between Version 6.3.0 and 6.3.1

<p>
This is the list of changes made between the deal.II releases listed above.
</p>

<p>
All entries are signed with the names of the author. Regular
contributor's names are abbreviated by WB (Wolfgang Bangerth), GK
(Guido Kanschat), RH (Ralf Hartmann).
</p>


<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>

<ol>
  <li>None.
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li>
  <p>Fixed: The step-33 tutorial program can not be built with GCC versions
  4.5.x. There are in fact two problems, one that pertains to uses of
  <code>std::make_pair</code> that don't work any more with the upcoming
  C++ 1x standard that GCC 4.5.x already follows, and some in which the
  Trilinos package Sacado is incompatible with GCC 4.5.x, at least up to
  Trilinos version 10.4.0. While the latter problem can only be fixed in
  future Trilinos versions, at least the former problem is solved in step-33.
  <br>
  (WB 2010/07/18)
  </p>

  <li>
  <p>Fixed: GCC version 3.4.0 failed to compile the file
  <code>deal.II/source/numerics/matrices.cc</code> with
  an internal compiler error. This has
  now been worked around.
  <br>
  (WB 2010/07/15)
  </p>

  <li>
  <p>Fixed: A problem in the Makefiles caused error messages when
  building under CygWin.
  <br>
  (GK 2010/07/12)
  </p>

  <li>
  <p>Fixed: GCC version 3.3.x failed to compile the files
  <code>lac/include/lac/precondition_block.h</code>,
  <code>deal.II/source/multigrid/mg_dof_handler.cc</code> and
  <code>examples/step-34/step-34.cc</code>. These problems have
  now been worked around.
  <br>
  (WB 2010/07/12)
  </p>

  <li>
  <p>Fixed: Some older 3.x versions of GCC crashed compiling the functions in
  namespace DoFRenumbering::boost. There is now a configuration time test
  that checks that the compiler accepts the constructs in question. If the
  compiler does not, then these functions are disabled.
  <br>
  (WB 2010/07/01)
  </p>

  <li>
  <p>Fixed: Linking with more than one of the deal.II 1d, 2d, or 3d libraries
  when using static libraries did not work. This is now fixed. However, due to
  GCC bug <a href="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=10591"
  target="_top">10591</a>, GCC versions prior to and including 4.1.x will
  still not work. Working with shared libraries was not and is not affected
  by this problem.
  <br>
  (WB 2010/07/01)
  </p>

  <li>
  <p>Fixed: GCC version 4.0.1 had a bug that prevented it from compiling
  release 6.3.0 because it apparently had an infinite loop allocating
  memory when compiling <code>fe_values.cc</code> in optimized mode. This
  problem had been fixed in GCC 4.0.2, but some versions of Mac OS X still use
  this GCC version in their Xcode environment. In any case, the code in
  deal.II has been changed to avoid this problem.
  <br>
  (WB 2010/06/30)
  </p>

  <li>
  <p>Fixed: Configuring with an external BOOST version did not work when
  using shared libraries since the test ran in the wrong order with respect
  to another configure test. This is now fixed.
  <br>
  (Bradley Froehle 2010/06/29)
  </p>

  <li>
  <p>
  Fixed: deal.II release 6.3.0 did not compile with Trilinos versions 9.x and
  10.0. This is now fixed.
  <br>
  (Martin Kronbichler, WB 2010/06/28)
  </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li>None so far.
</ol>


<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li>None so far.
</ol>


<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li>
  <p>
  Fixed: The FEValues::get_cell() function was unusable from user code
  since its implementation used a class that was only forward declared
  and not visible at the point of instantiations in user code. This is now
  fixed.
  <br>
  (WB 2010/07/16)
  </p>

  <li>
  <p>
  Fixed: On some systems and compilers, the library could not be compiled
  because of a duplicate symbol in <code>MeshWorker::LocalResults</code>.
  This is now fixed.
  <br>
  (WB 2010/06/28)
  </p>
</ol>


*/
