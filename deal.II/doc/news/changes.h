/**
// * @page changes_after_8_0 Changes after Version 8.0

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
  <li>
  Removed: The member function face_to_equivalent_cell_index() in
  FiniteElementData has been removed. It had been deprecated a while
  back already. Please use FiniteElement::face_to_cell_index() instead.
  <br>
  (Wolfgang Bangerth, 2013/08/09)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
  <li>
  New: It is now possible to compile and link deal.II against LLVM's libcxx. For
  this, a few issues with C++ standard violations are resolved.
  <br>
  (Matthias Maier, 2013/08/09)
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
  <li>
  New: MappingQEulerian is now also instantiated for vector elements
  of type TrilinosWrappers::Vector as well as the MPI and block
  variants.
  <br>
  (Armin Ghajar Jazi, 2013/08/14)
  </li>

  <li>
  Fixed: The FiniteElement::face_to_cell_index() function had a bug
  that made it work incorrectly for elements that have more than one
  degree of freedom per line (in 2d) or per quad (in 3d). This is now
  fixed for the most common cases, namely the FE_Q elements as well
  as elements composed of FESystem elements. For all other cases, an
  exception is generated reporting that this case is not implemented.
  If you run into this, let us know.
  <br>
  (Wolfgang Bangerth, 2013/08/10)
  </li>

  <li>
  New: DataOutBase::VtkFlags now has a flag
  DataOutBase::VtkFlags::print_date_and_time that can be used to suppress output
  of date and time in output files. This is useful in test suites where a newer
  run at a different time produces differences against previously stored files,
  even though the actual data is exactly the same.
  <br>
  (Oleh Krehel, 2013/08/06)
  </li>

  <li>
  Fixed: The various block matrix classes are all derived from BlockMatrixBase
  which had race conditions when the set() or add() functions were called from
  different threads. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2013/08/05)
  </li>

  <li>
  Fixed: various fixes with assignment and reinit of PETScWrappers::MPI::Vector.
  <br>
  (Timo Heister, 2013/08/05)
  </li>

  <li>Fixed: An assertion wrongly triggered in
  DoFTools::make_hanging_node_constraints when used with a particular
  combination of FESystem elements containing FE_Nothing. This is now fixed.
  <br>
  (Denis Davydov, Wolfgang Bangerth, 2013/08/01)
  </li>

  <li>
  New: Add has_ghost_elements() for PETScWrappers::MPI::BlockVector and
  TrilinosWrappers::MPI::BlockVector.
  <br>
  (Timo Heister, 2013/08/01)
  </li>

  <li>
  SparsityTools::distribute_sparsity_pattern did not work correctly for
  block systems, this has been fixed (function has a different signature).
  <br>
  (Timo Heister, 2013/07/31)
  </li>

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
