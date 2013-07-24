/**
 * @page changes_after_8_0 Changes after Version 8.0

<p>
This is the list of changes made after the release of
deal.II version 8.0.0.
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

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
  <li>Fixed: When typing <code>make run</code> in the step-32 directory,
  the program was executed with <code>mpirun -np 2 ./step-32</code>. This
  assumes that a program <code>mpirun</code> exists, but also does that
  deal.II was in fact compiled with MPI support on. Neither was intended.
  This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>

  <li>New: The DataOut, DataOutFaces, and DataOutRotation classes now allow
  the output of data vectors using different DoFHandler objects (based on the
  same triangulation), by new functions add_data_vector. This is used in the
  step-31 tutorial program which avoids creating a joint DoFHandler just for
  output.
  <br>
  (Martin Kronbichler, 2013/07/24)
  </li>

  <li>Changed: GridGenerator used to be a class with only static members
  but is now a namespace, like all other similar constructs in deal.II.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>

  <li>Changed: In GridGenerator, several functions had erroneously been changed
  to take an argument of type <code>size_type</code> rather than <code>unsigned
  int</code>. <code>GridGenerator::size_type</code> was a typedef to
  types::global_dof_index, which for most users was <code>unsigned int</code>
  anyway, but could also be set to be a 64-bit integer type. In any case, the
  change has been reverted and these functions take just a regular
  <code>unsigned int</code> again.
  <br>
  (Wolfgang Bangerth, 2013/07/24)
  </li>
</ol>


*/
