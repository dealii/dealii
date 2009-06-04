/**
 * @page changes_between_6.2.0_and_6.2.1 Changes between Version 6.2.0 and 6.2.1

<p>
This is the list of changes made between the release of 
deal.II version 6.2.0 and version 6.2.1:
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
  When running in parallel, the @ref step_18 "step-18" tutorial program
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
