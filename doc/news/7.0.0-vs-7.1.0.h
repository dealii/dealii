// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
 * @page changes_between_7_0_and_7_1 Changes between Version 7.0 and 7.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


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
<li> Changed: GridTools, DoFTools, MGTools, VectorTools, MatrixCreator
and MatrixTools are now namespaces. They have long
been classes that had only public, static member functions, making
the end result semantically exactly equivalent to a namespace, which is
also how it was used. This is now also reflected in the actual code.
<br>
(Wolfgang Bangerth, 2011/04/27, 2011/09/14)

<li> Changed: The PETScWrappers::VectorBase and PETScWrappers::MatrixBase
classes tried to keep track of
whether the last operation done on a vector was to add to an element or to
write into one. If the previous such operation was of a different kind
than the current one, we would flush buffers (see the description in
@ref GlossCompress). However, trying to do this automatically turned
out to be an endless source of hard-to-find bugs in %parallel programs.
The scheme has therefore now been changed to the following: the classes
keep track of the previous operation and if it differs from the
current one, reports an error stating that the user needs to call
PETScWrapper::VectorBase::compress() or
PETScWrapper::MatrixBase::compress() instead.
<br>
(Wolfgang Bangerth, 2011/08/03, 2011/08/30)

<li> Changed: The classes Tensor, SymmetricTensor and Point now have an
additional template argument for the number type. While a default template
value of <code>double</code> ensures that all old code is still valid, this
change invalidates forward declarations of the form <code>template
@<int dim@> class Point</code> that might be present in user-defined header
files. Now forward declarations need to specify the type as well, i.e.,
<code>template @<int dim, typename Number@> class Point</code>. However,
nothing changes if the full declarations in <code>deal.II/base/tensor.h,
deal.II/base/symmetric_tensor.h</code> and <code>deal.II/base/point.h</code>
are included.
<br>
(Martin Kronbichler, 2011/08/02)

<li> Removed: deal.II no longer supports Trilinos versions prior to 10.0.
<br>
(Wolfgang Bangerth, 2011/06/29)

<li> Changed: deal.II has a namespace std_cxx1x that was used to
import classes from BOOST that are part of the upcoming C++ 1x standard. On
the other hand, if your compiler supported a sufficiently large subset
of C++ 1x, we had code that simply did
@code
  namespace std_cxx1x = std;
@endcode
allowing you to refer to everything that was part of the compiler's namespace
<code>std</code> under the alternative name. This turned out to be untenable
in connection to the changed outlined below for _1, _2, etc. Consequently,
if the compiler used supports C++ 1x, we now selectively import elements of the
compiler's namespace std into namespace std_cxx1x as well. This may lead to
incompatibilities if you are already using elements of the C++ 1x
standard by referring to them through the std_cxx1x namespace and these elements
are not on the list of selectively imported ones.
<br>
(Wolfgang Bangerth, 2011/05/29)

<li> Changed: Previously, placeholder arguments like _1, _2, etc that are used
in conjunction with the std_cxx1x::bind function could be referenced as if
they are part of the global namespace. This was achieved by importing the
corresponding elements of namespace std::placeholders into the global namespace
if your compiler supported this part of the C++ 1x standard, or otherwise using
the BOOST counterparts which are already in the global namespace. However,
this leads to a conflict if one has a C++ 1x enabled compiler (e.g. GCC 4.6)
<i>and</i> includes certain BOOST headers, since the importation of symbols
into the global namespace now leads to ambiguous names. The only solution to
the problem is to not import names into the global namespace, but rather
import the names from either BOOST or namespace std into the deal.II namespace
std_cxx1x. The downside is that all code that uses _1, _2, etc needs to be
changed to use std_cxx1x::_1, std_cxx1x::_2, etc from now on.
<br>
(Wolfgang Bangerth, 2011/05/29)
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>
<li> New: Long in the making, the new tutorial program step-32 is a massively
%parallel simulator for the Boussinesq equations that describe thermal convection.
<br>
(Martin Kronbichler, Timo Heister, Wolfgang Bangerth, 2011/10/06)

<li> New: There is now a namespace Utilities::MPI that holds some of the MPI-related
functions that were previously part of Utilities::System. Specifically, the following
functions were moved and in part renamed: Utilities::System::get_n_mpi_processes
is now Utilities::MPI::n_mpi_processes; Utilities::System::get_this_mpi_process
is now Utilities::MPI::this_mpi_process;
Utilities::System::compute_point_to_point_communication_pattern
is now Utilities::MPI::compute_point_to_point_communication_pattern;
Utilities::System::duplicate_communicator
is now Utilities::MPI::duplicate_communicator;
Utilities::System::calculate_collective_mpi_min_max_avg
is now Utilities::MPI::min_max_avg;
Utilities::System::MPI_InitFinalize
is now Utilities::MPI::MPI_InitFinalize.
In addition, some of the arguments of these functions or classes
have changed.
<br>
The previous functions should still be available, though their use
is now deprecated.
<br>
(Wolfgang Bangerth, 2011/09/26)

<li> Removed: Utilities::System::program_uses_mpi does exactly the same thing
as Utilities::System::job_supports_mpi. The former has therefore been
deprecated.
<br>
(Wolfgang Bangerth, 2011/09/26)

<li> New: When using a new enough version of GCC, debug sections in
object files are now compressed using the <code>-Wa,--compress-debug-sections</code>
flag, resulting in savings in disk space on the order of 230 MB.
<br>
(Wolfgang Bangerth, 2011/09/22)

<li> New: deal.II can now be configured and built with the
<a href="http://clang.llvm.org">Clang C++ frontend of the LLVM compiler</a>.
<br>
(Timo Heister, Wolfgang Bangerth, 2011/08/20)

<li> Changed: Several of the tutorial programs used the same class
names, often <code>LaplaceProblem</code> in the earlier programs.
There is nothing inherently wrong with this, since these are
entirely separate programs, but this sometimes confused
integrated development environments because searching for
individual symbols would turn up several occurrences in different
files. To make this a bit simpler, the main classes for step-3
through step-6 were renamed to <code>StepX</code> (steps 1 and 2
do not have a main class). Starting with step-7, everything specific
to a tutorial program has been moved into a namespace called
<code>StepX</code> to make the fully qualified names unique across
different tutorial programs.
<br>
(Wolfgang Bangerth, Guido Kanschat 2011/08/16)

<li> Extended: Many operations on objects of
type Point@<0@>, Quadrature@<0@>, etc
(including creation) were previously forbidden since such objects do not make
much sense. However, this prevented a lot of code that could otherwise work
in a dimension independent way, from working in 1d, e.g. integration on
faces. Many of these places have now been cleaned up and work.
<br>
(Wolfgang Bangerth, 2011/08/12)

<li> Extended: The classes Tensor, SymmetricTensor and Point now have an
additional template argument for the number type. It is now possible to base
these classes on any abstract data type that implements basic arithmetic
operations, like <code>Tensor<1,dim,std::complex<double> ></code>. deal.II
uses a default template argument <code>double</code> that ensures that all
code using e.g. <code>Tensor<1,dim></code> remains valid.
<br>
(Martin Kronbichler, 2011/08/02)

<li> Fixed: deal.II can link with Trilinos but previously it required a
very specific set of Trilinos sub-libraries; if Trilinos had been compiled
with a larger set of sub-libraries, linking would sometimes fail. This
has now been made more generic and deal.II obtains the proper set of
libraries from Trilinos.
<br>
(Wolfgang Bangerth, 2011/06/29)

<li> Fixed: On Mac OS X, linking with some external libraries such as Trilinos
sometimes failed due to a misconfiguration of linker flags. This should now be
fixed.
<br>
(Praveen C, Martin Kronbichler, Wolfgang Bangerth, 2011/06/23)

<li> Changed: Doing <code>make clean</code> was supposed to only remove object
files but not libraries; however, it also removed the TBB libraries and a
few executables. This has now been changed: <code>make clean</code> now only
removes stuff that isn't needed to run executables, i.e. it leaves the TBB
and other libraries alone. As before, the target <code>make distclean</code>
is responsible for removing everything.
<br>
(Max Jensen, Wolfgang Bangerth, 2011/06/14)

<li> New: The Triangulation and DoFHandler classes, together with many
smaller classes can now be serialized, i.e. their data can be written
to an output stream and later retrieved to restore the state of the program.
<br>
(Wolfgang Bangerth, 2011/06/13)

<li> New/deprecated: The Triangulation class offers ways to get informed
whenever the triangulation changes. Previously, the mechanism doing this
was through the Triangulation::RefinementListener class. This has been
deprecated and has been superceded by a BOOST signals based mechanism
that is generally more powerful and does not rely on overloading
particular virtual functions inherited from a base class.

While the old mechanism should continue to work, you should consider
upgrading. For more information on the signals mechanism, see the
documentation of the Triangulation class.

In addition to the change above, the new implementation now offers two
more signals one can subscribe to: Triangulation::Signals::clead for
when the triangulation is cleared, and Triangulation::Signals::any_change
that can be used for any operation that changes the mesh. Furthermore,
in a change from previous behavior, the Triangulations::Signal::create
signal is now also triggered when another triangulation is copied to
the one that owns the signal.
<br>
(Wolfgang Bangerth, 2011/06/01)

<li> Removed: The <code>./configure</code> script allowed configuring
for the GNU Scientific Library (GSL) in version 7.0 but didn't actually
use any of the GSL functions. The corresponding code has therefore been
removed again.
<br>
(Wolfgang Bangerth, 2011/05/22)

<li> Changed: Traditionally, include directories were set through the
<code>-I</code> flag in make files in such a way that one would do
@code
  #include <base/quadrature.h>
@endcode
In preparation for future changes that will make possible installing
header files in a directory under <code>/usr/include</code> it seemed
useful to install everything under <code>/usr/include/deal.II</code>
and include them as
@code
  #include <deal.II/base/quadrature.h>
@endcode
This change has been made throughout the library and tutorial programs.
However, the old way of using include directories will continue to work
for at least one release for backward compatibility.
<br>
(Wolfgang Bangerth, 2011/05/16)

<li> Changed: The version of BOOST we ship with deal.II has been upgraded
to 1.46.1. BOOST now also resides in the directory <code>contrib/boost-1.46.1</code>
instead of an unversioned directory.
<br>
(Wolfgang Bangerth, 2011/05/16)

<li> New: The SparseDirectUMFPACK class can now also deal with matrices
provided in SparseMatrixEZ format.
<br>
(Martin Genet, 2011/05/04)

<li> New: The new tutorial program step-46 shows how to couple different
models defined on subsets of
the domain, in this case Stokes flow around an elastic solid. The
trick here is that variables (here the flow velocity and pressure,
and the solid displacement) do not live on the entire domain, but
only on a part. The point of the program is how to represent this in
source code.
<br>
(Wolfgang Bangerth, 2011/04/30)

<li> Fixed: On Debian, the Trilinos packages use a different layout
of include files and library names. The <code>./configure</code>
script can now deal with this.
<br>
(Walter Landry, 2011/02/22)

<li> Improved: Linking the deal.II libraries on file systems that
are mounted remotely from a file server took painfully long. This
is now fixed by linking everything on the local file system
and only subsequently moving the file into its final location.
<br>
(Wolfgang Bangerth, 2011/01/28)

<li> Changed: Most classes in deal.II have a member function
<code>memory_consumption</code> that used to return an unsigned int.
However, on most 64-bit systems, unsigned int is still only 32-bit
wide, and consequently the return type does not provide enough
precision to return the size of very large objects. The return types
of all of these functions has been changed to std::size_t, which is
defined to be a type that can hold the sizes of all objects possible
on any system.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: When using the <code>--enable-mpi</code> to
<code>./configure</code>, the script only tried <code>mpiCC</code>
as the MPI C++ compiler. However, on some systems, it is called
<code>mpicxx</code>. The script now tries that as well.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: When using Trilinos and using the Intel C++ compiler,
we accidentally used invalid compiler flags that led to a warning
every time we compiled a file.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: At the bottom of the page of tutorial programs we show a "plain"
version of the tutorial program. However, the script that generates this plain
version was broken and sometimes truncated the file. This
should be fixed now.
<br>
(Wolfgang Bangerth, 2011/01/18)

<li> Extended: Several missing instantiations of functions for triangulations
and DoF handlers embedded in higher dimensional space have been added.
<br>
(Wolfgang Bangerth, 2011/01/15)
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> New: There is a new predicate IteratorFilters::LocallyOwnedCell
for filtered iterators.
<br>
(Wolfgang Bangerth, 2011/10/03)

<li> Improved: The class BlockList used in RelaxationBlock has been replaced
by SparsityPattern, since it only reproduced its functionality.
<br>
(Guido Kanschat, 2011/09/26)

<li> New: SparsityPattern::row_position() finds a column index in a row and returns
its "local" index or numbers::invalid_unsigned_int.
<br>
(Guido Kanschat, 2011/09/26)

<li> New: The functions Utilities::MPI::sum and Utilities::MPI::max
function can be used to compute the sum or maximum of values over a number of
MPI processes without having to deal with the underlying MPI functions.
<br>
(Wolfgang Bangerth, 2011/09/25)

<li> New: The CellAccessor::is_locally_owned function is a shortcut
for the <code>!cell-@>is_ghost() && !cell-@>is_artificial()</code> pattern
found in many places when dealing with distributed meshes.
<br>
(Wolfgang Bangerth, 2011/09/10)

<li> New: The GridGenerator::torus function and TorusBoundary class can
create and describe the surface of a torus.
<br>
(Daniel Castanon Quiroz, 2011/09/08)

<li> Fixed: FEFieldFunction class was not thread-safe because it keeps a
cache on the side that was invalidated when different
threads kept pouncing on it. This is now fixed.
<br>
(Patrick Sodr&eacute;, 2011/09/07)

<li> New: There is now a function GridTools::volume() computing the volume
of a triangulation.
<br>
(Wolfgang Bangerth, 2011/09/02)

<li> New: Code like <code>cell-@>face(1)-@>set_boundary_indicator(42);</code>
now also works in 1d.
<br>
(Wolfgang Bangerth, 2011/08/30)

<li> Fixed: The TimerOutput::print_summary() function changed the
precision of output on the stream it prints to, but didn't restore
the previous value. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/08/30)

<li> New: There are now two new functions Utilities::string_to_double.
<br>
(Wolfgang Bangerth, 2011/08/25)

<li> Changed: The function VectorTools::compute_no_normal_flux_constraints
used to compute its constraints by evaluating the normal vector to the
surface as described by the mapping, rather than using the normal to
the surface described by the Boundary object associated with this face.
(Note that the Mapping computes its approximation by polynomial interpolation
of the surface described by the Boundary object.) This has now been changed:
the normal vector is now obtained from the Boundary object directly, at
points computed by the Mapping.
<br>
(Wolfgang Bangerth, 2011/08/25)

<li> New: The Boundary base class now has a function Boundary::normal_vector
that returns the normal vector to the surface at a given location. Derived
classes need to implement it, of course, if they want it to be used.
<br>
(Wolfgang Bangerth, 2011/08/25)

<li> Fixed: The function VectorTools::compute_no_normal_flux_constraints had
a problem that led to extremely difficult to pin down bugs when running
with sufficiently many processors. Basically, the constraints computed
by different processors did not agree which should be the independent
degrees of freedom and which should be the constrained ones. The result
were constraints that did not lead to a consistent linear system.
<br>
(Martin Kronbichler, 2011/08/24)

<li> New: Added GridRefinement::hierarchical() to reorder the degrees of freedom
by going through the cells in hierarchical order. This ensures consistent
DoF numbering in parallel computations.
<br>
(Timo Heister, 2011/08/24)

<li> Changed: Triangulation<dim>::get_boundary_indicators() returned
wrong data for dim=1.
<br>
(Sebastian Pauletti 2011/08/17)

<li> Improved: The function LogStream::timestamp() outputs all results of the Posix
function times, namely wall time, user time and system time.
<br>
(Guido Kanschat, 2011/08/18)

<li> Extended: GridGenerator::half_hyper_shell() got the option
<code>colorize</code>, which assigns different boundary indicators to the
different parts of the boundary. Added GridGenerator::quarter_hyper_shell()
with the same options.
<br>
(Timo Heister, 2011/08/15)

<li> Fixed: The functions VectorTools::create_boundary_right_hand_side()
called with an empty set of boundary_indicators (the default), did not apply
any boundary conditions. The empty set now applies it to all boundaries.
<br>
(Timo Heister, Sebastian Pauletti, 2011/08/15)

<li> Fixed: The function VectorTools::compute_no_normal_flux_constraints had
a bug that led to an exception whenever we were computing constraints for
vector fields located on edges shared between two faces of a 3d cell if those
faces were not perpendicular. This is now fixed.
<br>
(Wolfgang Bangerth, Thomas Geenen, Timo Heister, 2011/08/10)

<li> New: The function FullMatrix::triple_product() adds triple products
like Schur complements to existing matrices.
<br>
(Guido Kanschat, 2011/08/05)

<li> Improved: The PETScWrapper::VectorBase class was rather generous in
calling the PETSc <code>VecAssembleBegin/End</code> functions that incur
communication in the %parallel case and are therefore causes of potential
slowdowns. This has been improved.
<br>
(Wolfgang Bangerth, 2011/08/03)

<li> Fixed: The function VectorTools::create_right_hand_side now also works
for objects of type hp::DoFHandler with different finite elements.
<br>
(Daniel Gerecht, 2011/07/20)

<li> Improved: Evaluation of Lagrangian basis functions has been made stable
by exchanging polynomial evaluation from the standard form
$a_n x^n+\ldots+a_1 x + a_0$ to a product of linear factors,
$c (x - x_0) (x-x_1)\ldots (x-x_n)$. This ensures accurate evaluation up to
very high order and avoids inaccuracies when using high order finite elements.
<br>
(Martin Kronbichler 2011/07/26)

<li> Improved: The internal functions in the constructor of the FE_Q element
have been improved for high order elements. Especially when the element is
constructed for a 1D quadrature formula, the initialization is now much faster.
E.g. the initialization up to order 12 in three dimension completes in less
than a second, whereas it took hundreds of seconds before.
<br>
(Martin Kronbichler 2011/07/26)

<li> New: There is now a class Threads::ThreadLocalStorage that allows threads
to have their own copy of an object without having to fear interference from
other threads in accessing this object.
<br>
(Wolfgang Bangerth 2011/07/07)

<li> Fixed: The 2d grid reordering algorithm that is used by all grid readers had
a component that was quadratic in its complexity, sometimes leading to cases
where reading in a mesh in debug mode could take minutes for just a few tens
of thousands of cells. This has now been fixed.
<br>
(Wolfgang Bangerth 2011/07/07)

<li> New: The function DoFTools::count_dofs_per_component now also works
for objects of type hp::DoFHandler.
<br>
(Christian Goll, Wolfgang Bangerth 2011/07/06)

<li> Fixed: Under some circumstances, Threads::Thread::join() could only be
called once and would generate a system exception when called a second time.
Since it is often useful to not track whether this function had already been
called, this is now worked around in such a way that one can always call
the function multiple times.
<br>
(Wolfgang Bangerth 2011/07/03)

<li> New: The Threads::Thread::join() function can now also be called even
if no thread has been assigned to this thread object. The function then simply
does nothing.
<br>
(Wolfgang Bangerth 2011/07/03)

<li> New: There is now a new function Threads::Thread::valid that can be used
to query whether the thread object has been assigned a thread.
<br>
(Wolfgang Bangerth 2011/07/01)

<li> New: The new function GridGenerator::merge_triangulations can be used to compose
coarse meshes from simpler ones by merging their cells into a single
triangulation object.
<br>
(Wolfgang Bangerth 2011/06/17)

<li> Fixed: If an FEValues object was kept around until after the triangulation
on which it works has been refined or coarsened, and is then reinitialized
with a cell from the refined triangulation, it could compute wrong results or
crash outright. This has now been fixed.
<br>
(Wolfgang Bangerth 2011/06/02)

<li> Changed: The TrilinosWrappers::SparsityPattern::trilinos_sparsity_pattern()
function returned a reference to an object of kind Epetra_CrsMatrix. However, the
actual object pointed to is of derived class Epetra_FECrsMatrix. The function
has now been changed to return a reference to the latter type. Since derived
references can be assigned to references to base, this change should not
result in any incompatibilities.
<br>
(Wolfgang Bangerth 2011/05/27)

<li> New: The class RelaxationBlockJacobi has been added to the relaxation classes.
<br> (Guido Kanschat, 2011/05/19)

<li> New: discontinuous Galerkin versions of vector-valued elements have been
implemented: FE_DGBDM, FE_DGNedelec, and FE_DGRaviartThomas.
<br> (Guido Kanschat, 2011/05/19)

<li> New: Mapping<dim,spacedim>::transform_real_to_unit_cell  now
works also in the codimension one case, where it performs the normal
projection of the point on the codimension one surface.
<br> (Luca Heltai, 2011/05/17)

<li> New: The PersistentTriangulation class now works also in
the codimension one case.
<br>
(Luca Heltai, 2011/05/16)

<li> Fixed: The TrilinosWrappers::SparseMatrix::print() function
didn't get column indices right. This is now fixed.
<br>
(Habib Talavatifard, Wolfgang Bangerth 2011/05/10)

<li> Fixed: The TrilinosWrappers::SparseMatrix::operator() and
TrilinosWrappers::SparseMatrix::el() functions sometimes produced
wrong results for rectangular matrices. The same is true for
TrilinosWrappers::SparsityPattern::exists(). This is now fixed.
<br>
(Habib Talavatifard, Wolfgang Bangerth 2011/05/09, 2011/05/27)

<li> New: The version of DoFTools::make_flux_sparsity_pattern that takes
the coupling masks is now also available for hp::DoFHandler objects.
<br>
(Wolfgang Bangerth, 2011/04/27)

<li> Fixed: If Triangulation::create_triangulation is called after an
hp::DoFHandler object is attached to the triangulation object, setting active
FE indices leads to a crash. The problem did not happen if the mesh was
refined before setting the FE indices. This is now fixed. In the process, the
Triangulation::RefinementListener::create_notification function was
introduced.
<br>
(Wolfgang Bangerth, 2011/04/22)

<li> Fixed: The function FEValuesViews::SymmetricTensor::divergence had a bug.
This is now fixed.
<br>
(Wolfgang Bangerth, Feifei Cheng, Venkat Vallala 2011/04/21)

<li> Fixed: Under some conditions, FEFaceValues applied to an FESystem element
that contained elements of type FE_Nothing would receive an erroneous
exception. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/04/17)

<li> New: There is now an operator* for the multiplication of a <code>SymmetricTensor@<2,dim@></code>
and a <code>Tensor@<1,dim@></code>.
<br>
(Wolfgang Bangerth, 2011/04/12)

<li> Fixed: Added some instantiations to make KellyErrorEstimator and SolutionTransfer
work in  codimension one. Fixed some dim in spacedim.
<br>
(Luca Heltai, 2011/04/11)

<li> Fixed: Added some instantiations to make anisotropic refinement work
in codimension one.
<br>
(Luca Heltai, 2011/03/31)

<li> Fixed: Corrections in the creation of the face and subface
interpolation matrices in the class FE_Nedelec.
<br>
(Markus B&uuml;rg, 2011/03/17)

<li> Fixed: In step-21, the inner iteration would sometimes not converge for
very coarse meshes because of numerical roundoff. This is now fixed by allowing
more than <code>rhs.size()</code> CG iterations if the number of degrees of freedom
is very small.
<br>
(Jichao Yin, Wolfgang Bangerth, 2011/04/06)

<li> New: There is now a new function ConditionalOStream::get_stream().
<br>
(Wolfgang Bangerth, 2011/03/09)

<li> Fixed: FESystem::get_unit_face_support_points would refuse to return
anything if one of the base elements did not have support points. This
condition has been relaxed: it now only doesn't return anything if this
base element has no support points and also has degrees of freedom on
the face.
<br>
(Wolfgang Bangerth, 2011/03/07)

<li> Fixed: Objects of type FE_Nothing could be generated with multiple vector components
by passing an argument to the constructor. Yet, the FE_Nothing::get_name() function
always just returned the string <code>FE_Nothing@<dim@>()</code> independently of the
number of components. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/03/07)

<li> Fixed: PETScWrappers:MPI:SparseMatrix and apply_boundary_values() produced an error in debug mode about non-existant SparsityPattern entries. Reason: clear_rows() also eliminated the whole row in the PETSc-internal SparsityPattern, which resulted in an error in the next assembly process.
<br>
(Timo Heister, 2011/02/23)

<li> Fixed: It wasn't possible to use the FE_Nothing element inside an FESystem
object and hand the result over to an FEValues object. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/02/18)

<li> New: There is now a function DataOutBase::write_visit_record that does
the equivalent for VisIt that DataOutBase::write_pvtu_record does for ParaView:
generate a file that contains a list of all other VTK or VTU files of which the
current parallel simulation consists.
<br>
(Wolfgang Bangerth, 2011/02/16)

<li> New: There is now a function TrilinosWrappers::VectorBase::minimal_value.
<br>
(Wolfgang Bangerth, 2011/02/16)

<li> Fixed: TableBase::operator= could not be compiled if the type of the
elements of the table was <code>bool</code>. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/02/16)

<li> Improved: The GridGenerator::hyper_shell function generated meshes in 3d
that are valid but of poor quality upon refinement. There is now an additional
option to generate a coarse mesh of 96 cells that has a much better quality.
<br>
(Wolfgang Bangerth, 2011/02/12)

<li> Fixed: There are systems where the <code>libz</code> library is installed
but the <code>zlib.h</code> header file is not available. Since the latter
condition was not tested, this would result in compiler errors. This is now
fixed.
<br>
(Wolfgang Bangerth, 2011/02/09)

<li> Fixed: Prolongation and restriction matrices were not computed at all
for elements of type FE_DGQ if <code>dim@<spacedim</code>. Consequently,
consumers of this information, such as the SolutionTransfer class or
the DoFCellAccess::set_dof_values_by_interpolation function did not
work either and simply returned zero results. This is now fixed.
<br>
(M. Sebastian Pauletti, Wolfgang Bangerth, 2011/02/09)

<li> Fixed: When refining a mesh with codimension one, edges were refined using
the same manifold description as adjacent cells, but this ignored that a
boundary indicator might have been purposefully set for edges that are truly at
the boundary of the mesh. For such edges, the boundary indicator is now honored.
<br>
(M. Sebastian Pauletti, Wolfgang Bangerth, 2011/02/09)

<li> Fixed: The functions VectorTools::compute_mean_value and
VectorTools::integrate_difference now also work for distributed
triangulations of type parallel::distributed::Triangulation.
<br>
(Wolfgang Bangerth, 2011/02/07)

<li> Changed: If the <code>libz</code> library was detected during library
configuration, the function DataOutBase::write_vtu now writes data in compressed
format, saving a good fraction of disk space (80-90% for big output files).
<br>
(Wolfgang Bangerth, 2011/01/28)

<li> New: Trilinos and PETSc vectors now have a function has_ghost_elements().
<br>
(Timo Heister, 2011/01/26)

<li> Changed: The TrilinosWrappers::MPI::BlockVector::compress function now takes an
argument (with a default value) in exactly the same way as the
TrilinosWrappers::MPI::Vector::compress function already did.
<br>
(Wolfgang Bangerth, 2011/01/21)

<li> Fixed: When calling DataOut::build_patches with a mapping, requesting more
than one subdivision, and when <code>dim@<spacedim</code>, then some cells
were not properly mapped. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/01/18)

<li> New: Restructured the internals of <code>PETScWrappers::Precondition*</code>
to allow a PETSc PC object to exist without a solver. New: Use
<code>Precondition*::%vmult()</code> to apply the preconditioner once.
Preconditioners now have a default constructor and an <code>initialize()</code>
function and are no longer initialized in the solver call,
but in the constructor or <code>initialize()</code>.
<br>
(Timo Heister, 2011/01/17)

<li> Fixed: Boundary conditions in the step-23 tutorial program are now
applied correctly. Matrix columns get eliminated with the used method
and introduce some contribution to the right hand side coming from
inhomogeneous boundary values. The old implementation did not reset the
matrix columns before applying new boundary values.
<br>
(Martin Stoll, Martin Kronbichler, 2011/01/14)

<li> Extended: <code>base/tensor.h</code> has an additional collection of
contractions between three tensors (<i>ie</i>. <code>contract3</code>).
This can be useful for writing matrix/vector assembly in a more compact
form than before.
<br>
(Toby D. Young, 2011/01/12)

</ol>


*/
