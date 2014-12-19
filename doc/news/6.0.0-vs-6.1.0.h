// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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
 * @page changes_between_6_0_and_6_1 Changes between Version 6.0 and 6.1

<p>
This is the list of changes made between the deal.II releases listed above.
made to the three sub-libraries <a href="#base">base</a>,
<a href="#lac">lac</a>, and <a href="#deal.II">deal.II</a>, as well as
changes to the <a href="#general">general infrastructure,
documentation, etc</a>.
</p>

<p>
All entries are signed with the names of the author. Regular
contributor's names are abbreviated by WB (Wolfgang Bangerth), GK
(Guido Kanschat), RH (Ralf Hartmann).
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
  <li> <p>Removed: The <code>FullMatrix::add_diag</code> function was removed. It
  offered functionality of questionable use in the first place, and its
  implementation was correct only for matrices of size $3 \times 3$, $4 \times 4$
  and $8 \times 8$ without anyone noticing for 10 years. Consequently, it can't
  have been used very frequently.
  <br>
  (WB 2007/12/09)
  </p>

  <li> <p>Changed: The namespace deal_II_numbers has been renamed dealii::numbers.
  The old name stemmed from a time when not everything was already in
  namespace <code>dealii</code>. The old name is retained via a namespace
  alias but is deprecated and will eventually be removed.
  <br>
  (WB 2007/11/02)
  </p>

  <li> <p>Changed: When writing output files in UCD format using either the
  DataOutBase or the GridOut class, we used to write a preamble at the
  beginning of the file that included the date and time the file was created.
  However, several visualization programs get confused when confronted with
  comments at the beginning of UCD files. Rather than printing a sensible
  error message, they usually simply refuse to show any output, making it
  very hard to track down the actual cause.
  <br>
  The classes mentioned above previously allowed to suppress writing a preamble
  by setting a flag in the DataOutBase::UcdFlags or GridOutFlags::Ucd
  structures. Given how complicated it is to find the actual
  source of trouble, the default for these flags has been changed to not
  produce this preamble any more. If desired, the flag can be used to still
  produce a preamble.
  <br>
  (WB 2007/10/30)
  </p>

  <li> <p>Changed: The version number of the deal.II intermediate format written by
  DataOutBase::write_deal_II_intermediate has been increased to 3 to accommodate the fact that
  we now support writing vector-valued data to output files in at least some output formats.
  (Previously, vector-valued date was written as a collection of scalar fields.) Since
  we can only read files written in intermediate format that has the same number as the
  files we generate, this means that files written with previous format numbers can now
  no longer be read.
  <br>
  (WB 2007/10/11)
  </p>

  <li> <p>Changed: FilteredMatrix now can be applied to any matrix having the standard
  set of <code>vmult</code> functions. In order to achieve this, its interface
  had to be overhauled.
  Only the <code>VECTOR</code> template argument remains. Furthermore, instead of
  PreconditionJacobi being applied to FilteredMatrix, FilteredMatrix
  can now be applied to any preconditioner.
  <br>
  (GK 2007/09/25)
  </p>

  <li> <p>Changed: The deprecated typedefs
  <code>internal::Triangulation::Line</code>,
  <code>internal::Triangulation::Quad</code>, and
  <code>internal::Triangulation::Hexahedron</code> have been removed.
  <br>
  (WB 2007/09/07)
  </p>
</ol>



<a name="general"></a>
<h3>General</h3>

<ol>
<li> <p>New: step-33, written mainly by David Neckels,
  illustrates how the Euler equations of compressible inviscid gas
  dynamics can be solved. More generally, the program also provides
  a framework for the solution of arbitrary hyperbolic, nonlinear
  conservation laws.
  <br>
  (WB, 2008/05/23)
  </p>

<li> <p>Improved: step-16 has been updated to show the
  use of Multigrid with Dirichlet boundary conditions.
  <br>
  (GK, 2007/05/13)
  </p>

<li> <p>New: step-22, written by Martin Kronbichler and Wolfgang
     Bangerth, shows how to solve, and in particular how to solve
     efficiently, the Stokes equations that describe viscous
     incompressible fluid flow.
  <br>
  (WB, 2008/03/28)
  </p>

<li> <p>New: step-29 demonstrates how problems involving
  complex %numbers can be implemented by viewing real and imaginary parts
  of a complex-valued solution as the two components of a vector-valued
  function.
  <br>
  (Moritz Allmaras, 2007/10/31)
  </p>

<li> <p>New: A significantly simpler way to code the assembly of linear
  systems for problems with more than one solution variable has been
  implemented. This is explained in detail in the report on @ref vector_valued
  and tutorial programs step-20 and step-21
  have been converted as well.
  <br>
  (WB, 2008/02/23)
  </p>

  <li> <p>Improved: On Mac OS X, the operating system provides for
  "frameworks", which are essentially collections of shared libraries.
  We now link with the "Accelerate" framework (from Mac OS X 10.4
  onwards) or the "vecLib" framework (for previous versions) instead
  of the individual BLAS and LAPACK libraries if they are needed. This
  insulates us from having to use the actual name of these libraries,
  which may be subject to change, and it may also link with optimized
  or vectorized libraries if they are available.
  <br>
  (Eh Tan, WB 2007/10/22)
  </p>

</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
<li> <p> New: The FunctionParser::default_variable_names() function returns
default names for variables in the given space dimension. For example, in
2d, this would be "x,y".
<br>
(WB 2008/05/12)
</p> </li>

<li> <p>New: ConvergenceTable::evaluate_convergence_rates has been incomplete
for a long time. An implementation for reduction_rate and
reduction_rate_log2 has been added in the case where one specifies
a reference column that allows to compute error rates also for local
refinement.
<br>
(Luca Heltai 2008/04/16)
</p></li>

<li> <p> Fixed: The Utilities::int_to_string function did not work for 5 or 6
digits, while its counterpart Utilities::needed_digits does. This is now fixed,
and the two are coherent.
<br>
(Luca Heltai 2008/04/16)
</p> </li>

<li> <p> New: The SymmetricTensor class now has a constructor that creates
an object from an array consisting its independent components.
<br>
(WB 2008/02/21)
</p> </li>

<li> <p> New: There are now output operators (i.e. <code>operator@<@<</code>)
for the SymmetricTensor class.
<br>
(WB 2008/02/12)
</p> </li>

<li> <p> Improved: LogStream is now thread safe and output lines from different
threads are separated now. Additionally, a function LogStream::log_thread_id()
has been added to log the id of the thread printing the message. The semantics
has slightly changed in that the header is now printed at the time when the line
is finished.
<br>
(GK 2008/01/22)
</p> </li>

<li> <p> New: There is a function Threads::this_thread_id() now returning
an integer id of the current thread.
<br>
(GK 2008/01/22)
</p> </li>

<li> <p> Improved: Quadrature now has an assignment operator.
<br>
(GK 2007/12/27)
</p> </li>

<li> <p>Fixed: MultithreadInfo::n_cpus now gives the correct result
  also on my Intel MacBook Pro.
  <br>
  (Luca Heltai 2007/11/12)
  </p> </li>

  <li> <p>New: There is now a template numbers::NumberTraits that provides
  the means to implement linear algebra and other algorithms for both real
  and complex data types.
  <br>
  (WB 2007/11/03)
  </p> </li>

  <li> <p>New: The DataOutBase class and derived classes can now deal with data that is
  logically vector-valued, i.e. derived classes can pass down information that some of
  the components of the output should be grouped together as vectors, instead of being
  treated as a number of separate scalar fields. This allows visualization programs to
  display these components as vector-fields right away, without the need to manually group
  them together into a vector.
  <p>
  While all output format writers receive the information necessary to do this, currently
  only the VTK reader in DataOutBase::write_vtk as well as the deal.II
  intermediate format writer make use of it.
  <p>
  The use of this ability is explained in the step-22 and
  step-33 tutorial programs.
  <br>
  (WB 2007/10/11)
  </p> </li>

  <li> <p>New: Macro #AssertDimension introduced for easier handling of
  ExcDimensionMismatch.
  <br>
  (GK 2007/09/13)
  </p> </li>

</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p>Changed: A few simple PETSc wrappers for helper functions
  PETScWrappers::MatrixBase::transpose,
  PETScWrappers::MatrixBase::is_symmetric,
  PETScWrappers::MatrixBase::is_hermitian, and
  PETScWrappers::MatrixBase::write_ascii have been added.
  <br>
  (Toby D. Young 2008/04/30)
  </p>

<li> Fixed: The SolverQMRS class could not deal with block vectors
and block matrices. This is now fixed.
<br>
(WB 2008/3/20)
</li>

<li> Fixed: The implementation of SparseILU::decompose was rather
inefficient in that it accessed random elements of the matrix in its
inner loop. It has been replaced by the algorithm given in the book
by Yves Saad: "Iterative methods for sparse linear systems", second
edition, in section 10.3.2. The code is now about 8 times faster than
before.
<br>
(WB 2008/2/29)
</li>

<li> Fixed: The implementation of SparseILU::vmult very needlessly called
SparseMatrix::vmult just to throw away the (nonsensical) result away
immediately. This is now fixed.
<br>
(WB 2008/2/27)
</li>

<li> New: The functionality print_gnuplot is now also available for sparsity
patterns derived from the class BlockSparsityPatternBase, e.g. for
BlockSparsityPattern in the same way as for general sparsity patterns.
<br>
(Martin Kronbichler 2008/2/26)
</li>

<li> Fixed: The PETScWrappers::MatrixBase::clear_row and
PETScWrappers::MatrixBase::clear_rows functions had a memory leak that made
PETSc become progressively slower over time. This is now fixed.
<br>
(Daniel Goldberg 2008/2/20)
</li>

<li> Improved: All GrowingVectorMemory objects that have the same template
argument will access the same memory pool. Therefore, it is now not a
crime any more to just create a memory pool and discard later. Furthermore,
logging of statistics has been switched off by default, such that linear
solvers remain silent.
<br>
(GK 2007/12/16)
</li>

<li> Fixed: Vector::operator/= can't work when the scaling factor is zero,
but it happened to check whether the factor was positive. That's of course
bogus, the check should have been whether it is non-zero. This has now been
fixed.
<br>
(WB 2007/11/03)
</li>

<li> New: A class ScaledMatrix was introduced which combines the vector operations of
an underlying matrix with a scaling.
<br>
(GK 2007/10/30)
</li>

<li> Improved: FilteredMatrix has an iterator now that allows users to access the
constraints individually.
<br>
(GK 2007/10/30)
</li>


</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>

  <li> <p> Improved: GridOutFlags::Gnuplot now takes an additional curved_inner_cells
  parameter. Similar to curved_inner_cells in DataOut this allows the output of interior
  cells with curved boundaries. This is useful for MappingQEulerian or similar mappings.
  <br>
  (RH 2008/05/12)
  </p></li>

  <li> <p>Fixed: Functions::FEFieldFunction<dim>.set_active_cell was declared
  inline in the declaration but not in the implementation. This is now fixed.
  <br>
  (Luca Heltai 2008/04/16)
  </p></li>

  <li> <p>New: The FE_RaviartThomasNodal now supports hp functionality.
  <br>
  (Zhu Liang, 2008/04/04)
  </p></li>

  <li> <p>Improved: The StraightBoundary object Triangulation::straight_boundary
  is now made public. This allows to use Triangulation::straight_boundary instead
  of recreating a StraightBoundary object when needed.
  <br>
  (RH, 2008/04/04)
  </p></li>

  <li> <p>Extendend: SolutionTransfer has an additional template parameters DH
  for the dof handler class to use, which defaults to DoFHandler. However, an
  instantiation is also provided for hp::DoFHandler.
  <br>
  (Tobias Leicht, 2008/04/03)
  </p></li>

  <li> <p>Extended: DataOut has an extended interface now which enables a user
  decision, which cells shall be written as curved cells using the provided
  mapping: none, only cells at the boundary, or all cells. The last choice can
  be useful when a mapping like MappingQEulerian is employed.
  <br>
  (Tobias Leicht, 2008/03/20)
  </p></li>

  <li> <p>Changed: From now on the SolutionTransfer class does not use any user
  pointers of the underlying Triangulation object. Thus several SolutionTransfer
  instances with different DoFHandler objects can be used simultaneously. It is
  no longer necessary to assemble all solution vectors into a new one with a
  combined DoFHandler.
  <br>
  (Tobias Leicht, 2008/03/18)
  </p></li>

  <li> <p>New: There is now a namespace DoFRenumbering::boost that contains the
  implementation of interfaces to three reordering strategies provided by the
  Boost Graph Library, namely DoFRenumbering::boost::Cuthill_McKee,
  DoFRenumbering::boost::king_ordering, and DoFRenumbering::boost::minimum_degree.
  <br>
  (WB 2008/03/07)
  </p></li>

  <li> <p>Changed: The DoFRenumbering class has been converted into a namespace.
  Since all its members were static, this hasn't changed anything fundamental
  and use of the functions therein should work in exactly the same way as before.
  <br>
  (WB 2008/03/07)
  </p></li>

  <li> <p>New: The new TriaObjectAccessor::set_all_boundary_indicators function
  does not only set the boundary indicator of the current object but of all that
  bound it as well. For example, in 3d, if TriaObjectAccessor::set_boundary_indicator
  is called on a face, then the boundary indicator of the 4 edges that bound the face
  remain unchanged. On the other hand, the boundary indicators of face and edges are
  all set at the same time using the new function.
  <br>
  (WB 2008/02/26)
  </p></li>

  <li> <p>Improved: The tria pointer and the GridIn::debug_output_grid function
  are changed from private to protected to be accessible in derived classes.
  <br>
  (RH 2008/02/24)
  </p></li>

  <li> <p>New: The new Boundary::get_new_point_on_face and
  get_intermediate_points_on_face functions offer an dimension independent
  interface to the Boundary class.
  <br>
  (RH 2008/02/24)
  </p></li>

  <li> <p>New: There is now a FEFaceValuesBase::get_face_index function which
  returns the index of the face selected the last time the reinit() function
  was called.
  <br>
  (RH 2008/02/15)
  </p></li>

  <li> <p>New: The MappingQEulerian class provides an arbitrary order Eulerian
  mapping where the nodes of a mesh are displaced by a previously computed
  vector field.
  <br>
  (Joshua White 2008/02/05)
  </p></li>

  <li> <p>New: The MappingQ class constructor now takes a parameter that determines
  whether a higher order mapping should also be used for interior cells. By default,
  the higher order mapping is only used for cells on the boundary; cells in the
  interior are bounded by straight edges described by a (bi-, tri-)linear mapping.
  <br>
  (WB 2008/02/05)
  </p></li>

  <li> <p>New: The function VectorTools::compute_no_normal_flux_constraints computes
  the constraints that correspond to boundary conditions of the
  form $\vec u \cdot \vec n = 0$.
  <br>
  (WB 2008/01/23)
  </p></li>

  <li> <p>Fixed: Neither ConstraintMatrix::print nor ConstraintMatrix::write_dot
  produced any output for constraints of the form $x_i=0$, i.e. where the right
  hand side is a trivial linear combination of other degrees of freedom. This
  is now fixed.
  <br>
  (WB 2008/01/23)
  </p></li>

  <li> <p>Improved: GridGenerator::subdivided_hyper_rectangle now also colorizes
  cells according to the octant they are in.
  <br>
  (GK 2007/12/20)
  </p></li>

  <li> <p>New: The DataOut, DataOutRotation, DataOutStack, and DataOutFaces
  can now deal with vector-valued data if the functions in DataOutBase
  that write in a particular graphical output format can deal with it.
  Previously, if a finite element field had more than one component,
  they were all output as logically independent scalar components;
  most visualization programs then allowed to display vector fields
  by composing them of individual scalar fields for each vector component.
  </p><p>
  With the new scheme, the DataOut_DoFData::add_data_vector() functions
  inherited by the classes listed above take an additional parameter
  that may be used to indicate that certain components of the data
  logically form a vector field. The output formats for which this
  is presently implemented then indicate this fact in the output file. The
  mechanism is shown in action in step-22 and
  step-33.
  <br>
  (WB 2007/10/14)
  </p></li>

  <li> <p>Improved: UpdateFlags::update_q_points has been renamed to
  UpdateFlags::update_quadrature_points. Additional update flags for support
  points have been added without functionality, yet.
  <br>
  (GK 2007/10/12)
  </p></li>

  <li> <p>Improved: The number of blocks of an FESystem was properly defined and the
  constructors changed accordingly. At least non of the test programs noticed the change.
  <br>
  (GK 2007/10/03)
  </p></li>

  <li> <p>Improved: In an effort to make names more consistent, second
  derivatives in FEValuesBase and UpdateFlags have been renamed to
  Hessians. Thus, the clash between the forms <tt>2nd</tt> and
  <tt>second</tt> has been removed. Old function and enum names are
  available for compatibility but have been marked deprecated.
  <br>
  (GK 2007/09/12)
  </p></li>

  <li> <p>Fixed+New: The GridOutFlags::Ucd and
  GridOutFlags::Msh structures now take a new parameter
  <code>write_lines</code> to output lines with boundary id different
  from 0 in three-dimensional meshes. This fixes an annoying bug for
  which meshes with ids different from zero where not written in a
  compatible way, and if re-read with the corresponding
  GridIn functions, would not yield the same mesh upon
  refinement.
  <br>
  (Luca Heltai 2007/09/10)
  </p>
  </li>

  <li> <p>Extended: The possibilities of graphical output via the DataOut,
  DataOutFaces and DataOutRotation classes have been extended by the
  ability to perform a postprocessing step before a given data
  vector is written to the output. This way, derived variables can be
  written instead of the original data. To this end, there is a new
  version of DataOut_DoFData::add_data_vector taking a
  data vector and a DataPostprocessor, which performs the actual
  calculation of new data based on the values and possibly derivatives
  of the original data at each point of the patch.
  <br>
  (Tobias Leicht 2007/09/10)
  </p>
  </li>

</ol>


*/
