/**
 * @page changes_after_7_3 Changes after Version 7.3

<p>
This is the list of changes made after the release of
deal.II version 7.3.0.
All entries are signed with the names of the authors.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

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

<li> Changed: The TimeDependent::end_sweep function with an argument indicating
the number of threads has been removed. Use the corresponding function without
an argument. Since the argument had a default value, few users will have used
this function.
<br>
(Wolfgang Bangerth, 2013/03/17)

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
  <li> New: deal.II now uses <a href="http://www.cmake.org/">CMake</a>
  as its configuration and build tool. Please read through the
  readme and other installation files for information about how the
  installation process has changed.
  <br>
  Because this touches the configuration of every external package we
  interact with, there are a number of other changes as a result:
  <ul>
    <li>The minimum supported version for Trilinos is now 10.8.x.
    <li>We no longer link with different versions of the p4est library
        in debug and optimized mode. Rather, we now link with the same
	library in both modes. The p4est installation instructions have
	been updated. <blink>TODO TODO TODO</blink>
  </ul>
  <br>
  (Matthias Maier, 2013/03/07)

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>

<li> Added TimerOutput::reset to remove the collected information so far and
added a new frequency TimerOutput::never to only output information if
triggered by print_summary().
<br>
(Timo Heister, 2013/03/20)

<li> Changed: FEValuesExtractors::Scalar, FEValuesExtractors::Vector and
FEValuesExtractors::SymmetricTensor could not be default constructed, and
consequently one could not easily put them into arrays (where they would
be default constructed when changing the size, and later assigned useful
values). These classes can now be default constructed to invalid
values, but can of course not be used in any useful way.
<br>
(Wolfgang Bangerth, 2013/03/15)

<li> Fixed: FETools::interpolation_difference did not work with PETSc.
This is now fixed.
<br>
(Timo Heister, 2013/03/01)

</ol>


*/
