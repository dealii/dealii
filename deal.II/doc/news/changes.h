/**
 * @page changes_after_7_1 Changes after Version 7.1

<p>
This is the list of changes made after the release of
deal.II version 7.1.0.
All entries are signed with the names of the author.
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
<li> None so far.
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>
<li> Fixed: On some systems, <code>mpiCC</code> turns out to alias the
C compiler, not the C++ compiler as expected. Consequently, try to use
<code>mpic++</code> or <code>mpicxx</code> before <code>mpiCC</code> as
these should really be unambiguous.
<br>
(Wolfgang Bangerth, 2011/11/05)

<li> Fixed: Intel's ICC compiler identifies itself as <code>icpc version
12.1.0 (gcc version 4.2.1 compatibility)</code> which we mistook as being
GCC version 4.2. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/05)

<li> Fixed: deal.II could not be compiled with gcc 4.6.1 when MPI is
enabled due to a missing include file in file
<code>source/base/utilities.cc</code>. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/25)
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> New: Implementation of an alternative handling of
inhomogeneous constraints in ConstraintMatrix. This is controlled with
a new parameter use_inhomogeneities_for_rhs in
distribute_local_to_global() and determines whether the correct or
zero values (this was the case before and still is the default) are
kept in the linear system during the solution process.
<br>
(Jörg Frohne, 2011/11/01)

<ol>
<li> Fixed: SparseMatrix::mmult and SpareMatrix::Tmmult had a number of
issues that are now fixed: (i) rebuilding the sparsity pattern was allowed
even if several of the matrices involved in these operations shared a
sparsity pattern; (ii) the functions had a vector argument that had a default
value but the default value could not be used because it wasn't used in a
template context deducible by the compiler.
<br>
(Wolfgang Bangerth, 2011/10/30)

<li> New:
parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning
setting which is necessary for save()/load() to be deterministic. Otherwise
the matrix assembly is done in a different order depending on the order of old
refinements.
<br>
(Timo Heister, 2011/10/26)

<li> New: TriaAccessor<>::minimum_vertex_distance().
<br>
(Timo Heister, 2011/10/25)

<li> New: TableHandler::print_text now supports not only printing column
keys above their own column, but also in a separate header, to make it simpler
for external plotting programs to skip this line.
<br>
(Wolfgang Bangerth, 2011/10/22)

<li> Fixed: Trying to write a TableHandler object that is empty resulted
in a segmentation fault. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/21)

<li> New: The TableHandler class can now pad columns that have only been
partially filled. See the documentation of the class for a description.
<br>
(Wolfgang Bangerth, 2011/10/18)

<li> Fixed: In TableHandler::print_text, it can happen that the function
wants to print an empty string as the element of the table to be printed.
This can confuse machine readers of this table, for example for visualization,
since they then do not see this column in that row. To prevent this, we now
print <code>""</code> in such places.
<br>
(Wolfgang Bangerth, 2011/10/18)

<li> Fixed: Using Trilinos versions 10.4 and later on Debian failed to
configure due to a different naming scheme of Trilinos libraries on
Debian. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/17)

<li> Changed: The TableHandler class has been changed significantly
internally.  It could previously store arbitrary values (though in practice,
only int, unsigned int, double and std::string were implemented). The class is
now restricted to this particular set of types. On the other hand, the
TableHandler class can now be serialized.
<br>
(Wolfgang Bangerth, 2011/10/17)

<li> Fixed: searching in the doxygen documentation.
<br>
(Timo Heister, 2011/10/13)

<li> New: parallel::distributed::Triangulation::save()/load() to store
the refinement information to disk. Also supports saving solution vectors
using the SolutionTransfer class.
<br>
(Timo Heister, 2011/10/12)

<li> Fixed: The DataOut_DoFData::merge_patches did not compile with newer compilers.
This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/11)
</ol>


*/
