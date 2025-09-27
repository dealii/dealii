// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2014 by the deal.II authors
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
 * @page changes_between_6_2_0_and_6_2_1 Changes between Version 6.2.0 and 6.2.1

<p>
This is the list of changes made between the deal.II releases listed above.
</p>


<ol>
  <li>
  <p>
  A trivial mistake made deal.II unable to compile against any PETSc
  version prior to 3.0.0. This is now fixed.
  </p>
  </li>

  <li>
  <p>
  When running in parallel, the step-18 tutorial program
  produced an error indicating that resetting user pointers was not
  possible. This is now fixed.
  </p>
  </li>

  <li>
  <p>
  The documentation tar-ball we provide for those who do not want to re-build
  their own documentation locally using doxygen, did not include any typeset
  formulas (an oversight: we used a machine without a latex installation to
  build this package). The 6.2.1 package gets this right.
  </p>
  </li>

  <li>
  <p>
  Some versions of gcc 3.3.x had a bug that showed with code recently introduced
  into our sources in that it erroneously warned about perfectly legitimate
  constructs. This is now fixed.
  </p>
  </li>
</ol>


*/
