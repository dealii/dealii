// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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
 * @page changes_between_7_3_and_8_0 Changes between Version 7.3 and 8.0

<p>
This is the list of changes made between the deal.II releases listed above.
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

<li> Removed: it was possible to call get_dof_indices(), get_dof_values(),
set_dof_values(), and distribute_local_to_global() for cells that were not
active, if the finite element only had DoFs on vertices (i.e. Q1). This is
no longer allowed.
<br>
(Timo Heister, 2013/06/02)

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
  <li> Improvements on the CMake build system: A Working C (or Fortran) compiler
  is now fully optional; Improved support for static linkage with a toggle
  "DEAL_II_PREFER_STATIC_LIBS" that will prefer static archives for non system
  libraries and "DEAL_II_STATIC_EXECUTABLE" that will switch the complete link
  interface to static linkage.
  <br>
  (Matthias Maier, 2013/07/16)
  </li>

  <li> New: various functions for parallel computations got introduced to
  make Trilinos and PETSc interfaces similar. Now step-40 can be used with
  PETSc or Trilinos with just a few changes. This patch also introduces
  better support for block systems in the PETSc interfaces.
  <br>
  (Timo Heister, 2013/07/15)
  </li>

  <li> New: deal.II can now be compiled to 64-bit global dof indices. To turn
  this feature on, use the cmake option -DDEAL_II_WITH_64BIT_INDICES=ON. If
  PETSc and/or Trilinos are used, they must be compiled to support 64-bit
  indices. To write a code that can use 32-bit and 64-bit indices depending on
  deal.II compilation option, use types::global_dof_index for all the global
  dof indices.
  <br>
  (Kainan Wang and Bruno Turcksin, 2013/06/05)
  </li>

  <li> New: All vector classes now have a member function
  <code>locally_owned_elements</code> that returns an index
  set indicating which elements of this vector the current
  processor owns.
  <br>
  (Wolfgang Bangerth, 2013/05/24)
  </li>


  <li> New: A new element FE_Q_iso_Q1 has been implemented that is defined by
  a subdivision of the element into smaller Q1 elements. An element of order
  @p p is similar to FE_Q of degree @p p with the same numbering of degrees of
  freedom. The element is useful e.g. for defining a sparser preconditioner
  matrix for AMG at higher order FE_Q elements or for representing a component
  of a system of PDEs where higher resolution is preferred over high order.
  <br>
  (Martin Kronbichler, 2013/05/14)
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
	been updated.
  </ul>
  <br>
  (Matthias Maier, 2013/03/07)
  </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>

<li>New: VectorTools::project and a whole host of similar functions
are now also available for objects of type hp::DoFHandler.
<br>
(Wolfgang Bangerth, 2013/07/21)
</li>

<li>Fixed: hp::DoFHandler::n_boundary_dofs() had a bug that always led
to a failed assertion. This is now fixed.
<br>
(Wolfgang Bangerth, 2013/07/21)
</li>

<li>Fixed: VectorTools::project has an option to first project onto the
boundary. However, the implementation of this option ignored the mapping
that is provided to the function. This is now fixed.
<br>
(Wolfgang Bangerth, 2013/07/16)
</li>

<li>Improved: The WorkStream class used throughout deal.II is now using
thread local variables and initializes temporary variables on the thread
that uses them, leading to better cache locality.
<br>
(Wolfgang Bangerth, 2013/07/16)
</li>

<li>Improved: The "pure" functions in MeshWorker::LocalIntegrator are now implemented and throw
an exception if not overloaded.
<br>
(Guido Kanschat, 2013/07/16)
</li>

<li> New: The function SparseDirectUMFPACK::Tvmult is now implemented.
<br>
(Matthias Maier, 2013/07/03)
</li>

<li> New: In addition to the FEValuesExtractors::Scalar,
FEValuesExtractors::Vector, and FEValuesExtractors::SymmetricTensor classes,
there are now also fully featured FEValuesExtractors::Tensor extractors
for non-symmetric tensors of rank 2.
<br>
(Denis Davydov, 2013/07/02)
</li>

<li> New: There are now functions Tensor::component_to_unrolled_index()
and Tensor::unrolled_to_component_indices() in the same way as they
already exist for the SymmetricTensor class.
<br>
(Denis Davydov, 2013/07/02)
</li>

<li> New: There is now a read-write version of TableIndices::operator[].
<br>
(Denis Davydov, 2013/07/02)
</li>

<li> New: The function parallel::distributed::Triangulation::copy_triangulation() is
now implemented.
<br>
(Martin Steigemann, 2013/07/02)
</li>

<li> New: TriaRawIterator::operator < (TriaRawIterator&) now implements a total ordering
relation for cells even on distributed::parallel::Triangulation across processors.
Additionally, TriaRawAccessor and CellAccessor now have an ordering relation.
<br>
(Guido Kanschat, 2013/06/24)
</li>

<li> New: CellAccessor::id() that returns a unique CellId that
also works in parallel computations (where level and index is not
useful).
<br>
(Timo Heister, 2013/06/24)
</li>

<li> New: added ConstantTensorFunction<rank, dim> and ZeroTensorFunction<rank, dim> to provide
a tensor valued analogue to ConstantFunction and ZeroFunction.
<br>
(Matthias Maier, 2013/06/20)
</li>

<li> Fixed: BlockSparsityPattern::column_number was returning
wrong values.
<br>
(Timo Heister, 2013/06/16)
</li>

<li> Fixed: The stabilization parameter for the artificial diffusion
in the step-31 tutorial program has been increased slightly to avoid
instabilities at later times (<i>t</i> > 60).
<br>
(Martin Kronbichler, 2013/06/04)
</li>

<li> Fixed: If an exception was generated on a task created by
Threads::new_task, the program would terminate with a segmentation
fault, leaving little trace of what had happened. This is now handled
more gracefully.
<br>
(Wolfgang Bangerth, 2013/06/02)
</li>

<li> Changed: subdomain ids can now only be queried/set on active cells.
Consequently, is_artificial(), is_ghost(), and is_locally_owned() is
now restricted to active cells.
<br>
(Timo Heister, 2013/05/31)
</li>

<li> Improved: Triangulation::begin(level) and Triangulation::end(level) now return an
empty iterator range if the level is larger than the maximal locally owned level,
but still in the global level range of a distributed Triangulation.
<br>
(Timo Heister and Guido Kanschat, 2013/05/26)
</li>

<li> New: The IndexSet::add_indices function that takes another IndexSet
object now has an additional argument <code>offset</code> that can be used
to offset the indices of first argument.
<br>
(Wolfgang Bangerth, 2013/05/25)
</li>

<li> New: ConstraintMatrix::distribute is now also implemented for
arguments of type PETScWrappers::MPI::BlockVector.
<br>
(Wolfgang Bangerth, 2013/05/25)
</li>

<li> Fixed: IndexSet::operator== returned the wrong results
in some cases.
<br>
(Wolfgang Bangerth, 2013/05/25)
</li>

<li> New: The global function <code>complete_index_set()</code>
creates and returns an index set of given size that contains
every single index with this range.
<br>
(Wolfgang Bangerth, 2013/05/24)
</li>

<li> New: All vector classes now have a static member variable
<code>supports_distributed_data</code> that indicates whether the
vector class supports data that is distributed across multiple
processors. This variable is provided as a <i>traits variable</i>
to allow generic algorithms working on general vector types to
query the capabilities of the vector class at compile time.
<br>
(Wolfgang Bangerth, 2013/05/23)
</li>

<li> Fixed: FETools::back_interpolate has been revised to work correctly
also with parallel::distributed::Vector.
<br>
(Martin Steigemann, 2013/05/23)
</li>

<li> Removed: The file <code>mesh_worker/vector_info.h</code> was unused and
untested. It has thus been removed.
<br>
(Wolfgang Bangerth, Guido Kanschat, 2013/05/21)
</li>

<li> Fixed: The method parallel::distributed::Vector::compress
(VectorOperation::insert) previously set the elements of ghost elements
unconditionally on the owning processor, even if they had not been touched.
This led to a problem in certain library functions where vector entries became
zero in a spurious way. This is now fixed by discarding the elements in ghost
entries for the VectorOperation::insert operation. This is legitimate since we
assume consistency of set elements across processors, so the owning processor
sets the element already.
<br>
(Martin Kronbichler, 2013/05/21)
</li>

<li> Improved: DoFTools::make_periodicity_constraints now also works
for meshes where the refinement level of the two sides of the domain
is not the same, i.e., one side is more refined than the other.
<br>
(Wolfgang Bangerth, 2013/05/20)
</li>

<li> Improved: Through the fields DataOutBase::VtkFlags::time and
DataOutBase::VtkFlags::cycle, it is now possible to encode the time and/or
cycle within a nonlinear or other iteration in VTK and VTU files written
via DataOutBase::write_vtk and DataOutBase::write_vtu.
<br>
(Wolfgang Bangerth, 2013/05/12)
</li>

<li> Fixed: The method ConvergenceTable::evaluate_convergence_rates with
 reference column did not take the dimension of the reference column into
 account, leading to wrong logarithmic rates for dim!=2. This can now be fixed
 by specifying the dimension as a last argument.
<br>
(Martin Kronbichler, 2013/05/10)
</li>

<li> Improved: The functions MatrixTools::create_mass_matrix and
MatrixTools::create_laplace_matrix take now an optional ConstraintMatrix
argument that allows to directly apply the constraints. This also helps
VectorTools::project. Note that not providing constraints remains the default
and recommended way to ensure consistency when several matrices are combined.
<br>
(Martin Kronbichler, 2013/05/08)
</li>

<li> New: The classes TrilinosWrappers::SparseMatrix and
TrilinosWrappers::BlockSparseMatrix now fully implement vmult and Tvmult with
deal.II's own vector classes Vector<double> and
parallel::distributed::Vector<double>.
<br>
(Martin Kronbichler, 2013/05/08)
</li>

<li> Improved: The matrix-vector product ChunkSparseMatrix::vmult now runs in
parallel in shared memory.
<br>
(Martin Kronbichler, 2013/05/07)
</li>

<li> New: The class ChunkSparseMatrix and the associated
ChunkSparsityPattern now offer iterator classes to iterate over rows of the
matrix in an STL-like way.
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
necessary. The second iteration (and, hence, old behavior) can be forced by
the flag SolverGMRES::AdditionalData::force_re_orthogonalization.
<br>
(Martin Kronbichler, 2013/05/06)
</li>

<li> Changed: FETools::interpolate is now instantiated for all
vector types, not just dealii::Vector and dealii::BlockVector.
<br>
(Wolfgang Bangerth, 2013/05/06)
</li>

<li> Fixed: setting values in TrilinosWrappers::SparseMatrix
in parallel was adding the values instead.
<br>
(Bruno Turcksin, Timo Heister, 2013/05/03)
</li>

<li> Fixed: Generate an error if the user tries to refine a cell
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
safe manner. Also allow to pop() the prefix "DEAL".
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

<li> Added GridOut::write_svg() to allow for the output of
two-dimensional triangulations in two space dimensions in the SVG
format (Scalable Vector Graphics, an generic XML-based vector image
format developed and maintained by the World Wide Web Consortium W3C).
This function also provides cell coloring and cell labeling for the
visualization of basic cell properties. Pespective view is further
possible and the cell level number may be converted into altitude,
revealing the inactive cells lying below.
<br>
(Christian Wülker, 2013/03/21)
</li>

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
