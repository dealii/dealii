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
  <li> New: The number of threads used by deal.II/TBB can now be limited at
  run time. Using MPI based code using PETSc/Trilinos no longer requires you
  to compile the library without threads. See MPI_InitFinalize and
  MultithreadInfo::set_thread_limit for details.
  </li>
  (Timo Heister, 2013/03/26)

  <li> New: The results section of step-36 now explains how to use ARPACK
  as an alternative to SLEPc as eigenvalue solver.
  <br>
  (Juan Carlos Araujo Cabarcas, 2013/03/25)

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

<li> New: MPI_InitFinalize can also initialize PETSc/Slepc when
not compiling with MPI. This is now the preferred way to initialize
MPI/PETSc/Slepc in all cases.
<br>
(Timo Heister, 2013/3/26)
</li>

<li> Added/fixed: IterativeInverse::vmult() can now handle vectors
using a different number type than the matrix type. As usual, the
number types must be compatible. Addtitionally, the initial guess is
always set to zero, since starting with the incoming vector makes no
sense.
<br>
(Guido Kanschat, 2013/03/21)
</li>

<li> Added GridOut::write_svg to allow for the output of two-dimensional
triangulations in two space dimensions in the SVG format (Scalable Vector
Graphics, an XML-based vector image format recommended by the World
Wide Web Consortium W3C). This function also provides cell coloring
and cell labeling for the visualization of basic cell properties.
<br>
(Christian Wülker, 2013/03/21)

<li> Added TimerOutput::reset to remove the collected information so far and
added a new frequency TimerOutput::never to only output information if
triggered by print_summary().
<br>
(Timo Heister, 2013/03/20)
</li>

<li> Changed: FEValuesExtractors::Scalar, FEValuesExtractors::Vector and
FEValuesExtractors::SymmetricTensor could not be default constructed, and
consequently one could not easily put them into arrays (where they would
be default constructed when changing the size, and later assigned useful
values). These classes can now be default constructed to invalid
values, but can of course not be used in any useful way.
<br>
(Wolfgang Bangerth, 2013/03/15)
</li>

<li> Fixed: FETools::interpolation_difference did not work with PETSc.
This is now fixed.
<br>
(Timo Heister, 2013/03/01)
</li>

</ol>


*/
