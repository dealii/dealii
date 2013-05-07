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

<li> Changed: Internal structures of ExceptionBase are now thread safe. The
Assert macro does not print an exception to deallog any more prior to
throwing if deal_II_exceptions::abort_on_exception==false. Removed: A
number of obsolete Exceptions that are not used in the library any more.
<br>
(Matthias Maier, 2013/04/16)

<li> Removed: A number of header files that have been deprecated a long time
ago have been removed. All of them had previously only included the header file
that had superseded them. To upgrade, simply include the currently used
header file. This is also backward compatible with deal.II 7.3.
<br>
(Guido Kanschat, 2013/04/12)

<li> Removed: The interfaces to the obsolete direct solvers MA27 and MA47 from
the Harwell Subroutine Library. Support for the HSL routines were not ported to
the new build system. However, the sparse direct solver UMFPACK remains to be
supported and is provided as part of the standard deal.II distribution, unlike
the HSL functions.
<br>
(Matthias Maier, 2013/04/01)

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
  <li> Fixed: setting values in TrilinosWrappers::SparseMatrix
  in parallel was adding the values instead.
  <br>
  (Bruno Turcksin, Timo Heister, 2013/05/03)
  </li>

  <li> New: The step-49 tutorial program now also has a discussion on
  what to do once you have a coarse mesh and want to refine it.
  <br>
  (Wolfgang Bangerth, 2013/04/03)
  </li>

  <li> New: The number of threads used by deal.II/TBB can now be limited at
  run time. Using MPI based code using PETSc/Trilinos no longer requires you
  to compile the library without threads. See MPI_InitFinalize and
  MultithreadInfo::set_thread_limit for details.
  <br>
  (Timo Heister, 2013/03/26)
  </li>

  <li> New: The results section of step-36 now explains how to use ARPACK
  as an alternative to SLEPc as eigenvalue solver.
  <br>
  (Juan Carlos Araujo Cabarcas, 2013/03/25)
  </li>

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
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>

<li> New: The matrix-vector product ChunkSparseMatrix::vmult now runs in
parallel in shared memory.
<br>
(Martin Kronbichler, 2013/05/07)
</li>

<li> New: The class ChunkSparseMatrix and the associated
ChunkSparsityPattern now offer iterator classes to iterate over rows or the
whole matrix in an STL-like way.
<br>
(Martin Kronbichler, 2013/05/07)
</li>

<li> Fixed: The stopping criterion for early exit in SolverBicgstab did not
work properly for systems with large values, leading to premature exit. This
is now fixed.
<br>
(Martin Kronbichler, 2013/05/07)
</li>

<li> Changed: The SolverGMRES implementation previously applied two
iterations of the modified Gram&ndash;Schmidt algorithm for
orthogonalization. In many situations one iteration is enough. The algorithm
can now detect loss of orthogonality and enables re-orthogonalization only if
necessary. The second iteration can be forced by the flag
SolverGMRES::AdditionalData::force_re_orthogonalization, though.
<br>
(Martin Kronbichler, 2013/05/06)
</li>

<li> Changed: FETools::interpolate is now instantiated for all
vector types, not just dealii::Vector and dealii::BlockVector.
<br>
(Wolfgang Bangerth, 2013/05/06)
</li>

<li> Fixed: Generate an error if the users tries to refine a cell
that is already on the maximum level in a distributed triangulation.
<br>
(Timo Heister, 2013/05/01)
</li>

<li> Fixed: The version of ParameterHandler::set that takes a boolean
as second argument was broken and did not work. This is now fixed.
<br>
(Ashkan Dorostkar, Wolfgang Bangerth, 2013/04/30)
</li>

<li> Fixed: PETScWrappers::VectorBase::print now saves and restores
the precision
and width associated with the stream it prints to around setting
the values passed as arguments.
<br>
(Fahad Alrashed, 2013/04/22)
</li>

<li> Fixed: FullMatrix::print now saves and restores the precision
and width associated with the stream it prints to around setting
the values passed as arguments.
<br>
(Fahad Alrashed, 2013/04/22)
</li>

<li> New: LogStream now has member functions LogStream::width,
LogStream::precision and LogStream::flags that make it look more
like normal objects of type <code>std::ostream</code>.
<br>
(Fahad Alrashed, 2013/04/22)
</li>

<li> New: SparseDirectUMFPACK has long had the ability to work with
BlockSparseMatrix objects, but couldn't deal with BlockVector objects.
This is now fixed.
<br>
(Wolfgang Bangerth, 2013/04/21)
</li>

<li> New: Class TimerOutput::Scope does automatic scope based enter/
exit_section of a TimerOutput object.
<br>
(Timo Heister, 2013/04/18)
</li>

<li> Fixed: TimerOutput constructed with an MPI_COMM in wall_time
mode now constructs synchronized Timer objects. This gives reliable
parallel benchmark timings.
<br>
(Timo Heister, 2013/04/18)
</li>

<li> Improved and Fixed: LogStream (and deallog) now respect std::flush in
addition to std::endl to write out content to the console/file.
Furthermore, LogStream::push(...) and LogStream::pop() now work in a thread
safe manner.
<br>
(Matthias Maier, 2013/04/18)
</li>

<li> Fixed: The HalfHyperShellBoundary class got refining
the edges that sit at the perimeter of the circular face of the domain
wrong. This is now fixed.
<br>
(Wolfgang Bangerth, J&ouml;rg Frohne, 2013/04/17)
</li>

<li> New: Functions::FEFieldFunction can now deal with
parallel::distributed::Triangulation objects.
<br>
(Wolfgang Bangerth, 2013/04/15)
</li>

<li> New: There is now a version of SparseMatrix::copy_from that can copy
from TrilinosWrappers::SparseMatrix.
<br>
(Wolfgang Bangerth, J&ouml;rg Frohne, 2013/04/15)
</li>

<li> Improved: The SolverCG implementation now uses only three auxiliary
vectors, down from previously four. Also, there are some shortcuts in case
PreconditionIdentity is used that improve the solver's performance.
<br>
(Martin Kronbichler, 2013/04/11)
</li>

<li> Fixed: The results section of step-23 did not show the movie in release 7.3
due to a poor HTML markup. This is now fixed.
<br>
(Wolfgang Bangerth, 2013/04/10)
</li>

<li> Fixed: It is now possible to use the MeshWorker framework in 1d as well.
<br>
(Wolfgang Bangerth, Scott Miller, 2013/04/09)
</li>

<li> Fixed: It was not possible to create a default-constructed object of
type Triangulation<1>::face_iterator. This is now fixed.
<br>
(Wolfgang Bangerth, Scott Miller, 2013/04/09)
</li>

<li> New: VectorTools::subtract_mean_value can now be called without a
boolean mask. The vector type is templatified and instantiated for all
non distributed vectors.
<br>
(Matthias Maier, 2013/04/08)
</li>

<li> Fixed: It is now possible to call ConvergenceTable::evaluate_convergence_rates
multiple times.
<br>
(Matthias Maier, 2013/04/08)
</li>

<li> Fixed: GridTools::distort_random (previously called Triangulation::distort_random)
had a bug where points were only ever moved in <i>positive</i> coordinate
directions rather than with uniform probability in either direction. The 1d
implementation also had the problem that it did not move vertices if the
<i>cell</i> they were on was at the boundary, even if the <i>vertex</i>
itself was not. All of these problems are now fixed.
<br>
(Wolfgang Bangerth, 2013/04/05)
</li>

<li> New: There is a class VectorFunctionFromTensorFunction that converts
between objects of type TensorFunction and Function.
<br>
(Spencer Patty, 2013/4/2)
</li>

<li> Fixed: The ParameterHandler class could not deal with parameters named
<code>"value"</code> (and a few other names). This is now fixed.
<br>
(Denis Davydov, Matthias Maier, Wolfgang Bangerth, 2013/3/31)
</li>

<li> Changed: TimerOutput no longer assumes that sections are not nested
when outputting percentage and total run time.
<br>
(Timo Heister, 2013/3/28)
</li>

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
