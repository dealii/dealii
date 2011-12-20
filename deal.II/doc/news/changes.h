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
<li> Fixed: The Intel compiler complains that it can't copy Trilinos vector
reference objects, preventing the compiling of step-32. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/09)

<li> Fixed: Intel ICC 12.1 gets into trouble with BOOST because BOOST
believes that the compiler supports C++0x but one then still has to
specify the corresponding flag on the command line to avoid compiler
errors. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/06)

<li> Fixed: On some systems, <code>mpiCC</code> turns out to alias the
C compiler, not the C++ compiler as expected. Consequently, try to use
<code>mpic++</code> or <code>mpicxx</code> before <code>mpiCC</code> as
these should really be unambiguous.
<br>
(Wolfgang Bangerth, 2011/11/05)

<li> Fixed: Intel's ICC compiler identifies itself as <code>icpc version
12.1.0 (gcc version 4.2.1 compatibility)</code> which we mistook as being
GCC version 4.2. The same is true for the Intel C compiler. This is now fixed.
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
<li> New: The GridIn class now also reads the GMSH format 2.2 as written by
GMSH 2.5.
<br>
(Vijay Mahadevan, Wolfgang Bangerth, 2011/12/19)

<li> Fixed: When <code>./configure</code> does not detect the presence
of <code>zlib</code>, writing output in VTU format failed to produce
a valid output file.
<br>
(Timo Heister, 2011/12/16)

<li> Improved: The GridRefinement::refine_and_coarsen_optimize function
assumed that the expected convergence order was 2. It has now gotten an
argument by which the user can prescribe a different value. A bug has also
been fixed in which the function incorrectly assumed in its algorithm that
the mesh was two-dimensional.
<br>
(Christian Goll, 2011/12/16)

<li> Fixed: Restriction and prolongation didn't work for elements of
kind FE_Nothing. Consequently, many other parts of the library also
didn't work, such as the SolutionTransfer class. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/12/15)

<li> Fixed: The DerivativeApproximation class now works for distributed computations.
<br>
(Timo Heister, 2011/12/15)

<li> Changed: The ExcMessage exception class took an argument of
type <code>char*</code> that was displayed when the exception
was raised. However, character pointers are awkward to work with
because (i) they can not easily be composed to contain things like
file names that are only known at run-time, and (ii) the string
pointed to by the pointer had to live longer than the local expression
in which the exception is generated when using the AssertThrow macro
(because, when we create an exception, the exception object is passed
up the call stack until we find a catch-clause; at that time, however,
the scope in which the exception object was created has long been left).
This restriction made it impossible to construct the message using std::string
and then just do something like <code>(std::string("file: ")+filename).c_str()</code>.
<br>
To remedy this flaw, the type of the argument to ExcMessage has been
changed to std::string since objects of this type are readily copyable
and therefore live long enough.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> New: Setting up a class derived from DataPostprocessor required some
pretty mechanical steps in which one has to overload four member functions.
For common cases where a postprocessor only computes a single scalar or
a single vector, this is tedious and unnecessary. For these cases, the
new classes DataPostprocessorScalar and DataPostprocessorVector provide
short cuts that make life simpler.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> Changed: The DataPostprocessor class previously required users of this
class to overload DataPostprocessor::get_names(),
DataPostprocessor::get_data_component_interpretation()
and DataPostprocessor::n_output_variables(). The latter function is redundant
since its output must equal the length of the arrays returned by the
first two of these functions. It has therefore been removed.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> Improved: Objects of the type LogStream::Prefix can now be used
as a safe implementation of the push and pop mechanism for log
prefices.
<br>
(Guido Kanschat, 2011/12/09)

<li> New: IndexSet::add_indices(IndexSet).
<br>
(Timo Heister, 2011/12/09)

<li> Fixed: Finite element Hessians get computed for codimension one,
at least for FE_Poly derived classes.
<br>
(Guido Kanschat, 2011/12/07)

<li> Changed: Output of ParameterHandler::print_parameters with argument
ParameterHandler::LaTeX was not particularly readable. The output has
therefore been rewritten to be more structured and readable.
<br>
(Wolfgang Bangerth, 2011/11/28)

<li> Fixed: The TimerOutput class set the alignment of output to right-aligned
under some circumstances, but didn't reset this to the previous value at the
end of output. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/28)

<li> New: The copy constructor of the SparseMatrixEZ function now works
(rather than throwing an exception) if the copied matrix has size zero.
This is in accordance to the other matrix classes.
<br>
(Wolfgang Bangerth, 2011/11/15)

<li> New: The class ScalarFunctionFromFunctionObject provides a quick way to
convert an arbitrary function into a function object that can be passed
to part of the library that require the Function@<dim@> interface.
The VectorFunctionFromScalarFunctionObject class does a similar thing.
<br>
(Wolfgang Bangerth, 2011/11/15)

<li> New: The DoFTools::count_dofs_per_block function now also works
for objects of type hp::DoFHandler.
<br>
(Jason Sheldon, 2011/11/13)

<li> New: FETools::get_fe_from_name() can now return objects of type FE_Nothing.
<br>
(Scott Miller, Jonathan Pitt, 2011/11/10)

<li> New: Implementation of an alternative handling of
inhomogeneous constraints in ConstraintMatrix. This is controlled with
a new parameter use_inhomogeneities_for_rhs in
distribute_local_to_global() and determines whether the correct or
zero values (this was the case before and still is the default) are
kept in the linear system during the solution process.
<br>
(J&ouml;rg Frohne, 2011/11/01)

<li> Fixed: SparseMatrix::mmult and SpareMatrix::Tmmult had a number of
issues that are now fixed: (i) rebuilding the sparsity pattern was allowed
even if several of the matrices involved in these operations shared a
sparsity pattern; (ii) the functions had a vector argument that had a default
value but the default value could not be used because it wasn't used in a
template context deducible by the compiler.
<br>
(Wolfgang Bangerth, 2011/10/30)

<li> New:
parallel::distributed::Triangulation::mesh_reconstruction_after_repartitioning
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
