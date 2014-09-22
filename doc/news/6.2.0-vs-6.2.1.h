// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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
