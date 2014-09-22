// ---------------------------------------------------------------------
//
// Copyright (C) 2013, 2014 by the deal.II authors
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
 * @page changes_between_5_1_and_5_2 Changes between Version 5.1 and 5.2

<p>
This is the list of changes made between the deal.II releases listed above.
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
  <li> <p>
       Changed: The way boundary constraints were handled in the
       step-17 tutorial program was conceptually flawed. We tried to
       eliminate boundary nodes locally on the cell level, and hanging
       node constraints subsequently when transferring into the global
       matrix. However, this doesn't work in general: the elimination
       of hanging node constraints could re-populate rows and columns
       that had already been vacated during boundary node
       elimination. At the end of a long thought process, we came to
       the conclusion that it is impossible to revert the traditional
       order of operations: first eliminate all hanging node
       constraints, then eliminate all boundary nodes. This leads to a
       situation where the <code
       class="member">MatrixTools::local_apply_boundary_values</code>
       function is not very useful any more, except for the special
       case where there are no hanging nodes. The step-17 example
       program has therefore been changed to use the <code
       class="member">MatrixTools::apply_boundary_values</code>
       function again, though hanging node elimination still happens
       during transfer into the global matrix.
       <br>
       (WB 2005/05/05)
       </p>

  <li> <p>
       Changed: The class <code
       class="class">MGCoarseGridLACIteration</code> lost two template
       arguments. Since the matrix and preconditioner are now stored
       in form of <code>PointerMatrix</code> objects,
       handling of the class is much simpler, in particular when
       exchanging preconditioners.
       <br>
       (GK 2005/05/03)
       </p>

  <li> <p>
       Changed: The template argument of <code
       class="class">BlockMatrixArray</code> was changed to
       <tt>number</tt>, equal to the same argument of the used <code
       class="class">BlockVector</code>. Furthermore, its constructor
       requires an additional argument of type <code
       class="class">VectorMemory&lt;Vector&lt;number&gt; &gt;</code>
       providing space for auxiliary vectors. Since the entries are
       now of type <code>PointerMatrixBase</code>, even
       matrices with blocks of different types can be constructed.
       <br>
       (GK 2005/03/21)
       </p>

  <li> <p>
       Changed: The <code>GeometryInfo</code>::<code
       class="member">vertices_adjacent_to_line</code> function has
       been renamed to <code>GeometryInfo</code>::<code
       class="member">line_to_cell_vertices</code> to be named
       analogous to the <code
       class="member">face_to_cell_vertices</code>, <code
       class="member">face_to_cell_lines</code> and <code
       class="member">child_to_cell_coordinates</code> functions.
       <br>
       (RH 2005/03/14)
       </p>

  <li> <p>
       Changed: the file <tt>multigrid/mg_dof_tools.h</tt> was renamed
       to <tt>multigrid/mg_tools.h</tt> in order to match the name of
       the class contained.
       <br>
       (GK 2005/03/09)
       </p>

  <li> <p>
       Changed: <code>DoFTools</code>::<code
       class="member">make_flux_sparsity_pattern</code>, <code
       class="class">MGTools</code>::<code
       class="member">make_flux_sparsity_pattern</code> and similar
       functions in these classes do not receive arguments of type
       <code>FullMatrix&lt;double&gt;</code>
       anymore. Instead, they get a <code
       class="class">Table&lt;2,DoFTools::Coupling&gt;</code>, which
       contains more meaningful enums.
       <br>
       (GK 2005/03/09)
       </p>

  <li> <p>
       Changed: <code>Multigrid</code>::<code
       class="member">Multigrid</code> receives an argument
       determining the type of multigrid cycle instead of the minimum
       and maximum levels. The latter were rarely used anyway and can
       be modified by <code>set_minlevel()</code> and
       <code>set_maxlevel()</code>
       <br>
       (GK 2005/03/09)
       </p>
</ol>

<a name="general"></a>
<h3>General</h3>

<ol>
  <li> <p> New:
       The step-19 tutorial demonstrates handling of parameters from a
       parameter file, as well as some simple techniques for merging output
       from parallel computations. 
       <br>
       (WB 2005/09/09)
       </p>

  <li> <p>
       Fixed: The linking of PETSc libraries into the static
       <tt>libpetsc.g.a</tt> and <tt>libpetsc.a</tt> libraries
       (located in <tt>deal.II/lib</tt>) has now been fixed.
       <br>
       (RH 2005/08/26)
       </p>

  <li> <p> New:
       The step-18 example program shows how to solve a quasi-static elasticity
       problem with moving meshes, and this all in parallel.
       <br>
       (WB 2005/08/11)
       </p>

  <li> <p>
       Changed: The files
       <tt>deal.II/include/grid/geometry_info.h</tt> and
       <tt>deal.II/source/grid/geometry_info.cc</tt> have been moved
       to <tt>base/include/base</tt> and <tt>base/source</tt>,
       respectively. A redirection header file has been left in the
       old location for compatibility.
       <br>
       (GK, 2005/07/19)
       </p>

  <li> <p>
       Changed: Cygwin has problems linking against deal.II DLLs if you
       link against more than one dimension dependent library. The linker
       might issue error messages about multiple defined symbols. This is now
       detected, and the linker is forced to ignore these errors. However, if
       you accidentally defined symbols twice you might have a hard time
       debugging now... In this case remove the <code>-Xlinker 
       --allow-multiple-definition</code> flag from Make.global_options.
       <br>
       (Ralf B. Schulz, 2005/07/13)
       </p>

  <li> <p>
       Changed: Under Cygwin, linking against LAPACK and UMFPACK did
       not work. This is now fixed in the Makefiles. Changes only
       affect Cygwin.
       <br>
       (Ralf B. Schulz, 2005/07/02)
       </p>

  <li> <p>
       New: The latest version of UMFPACK, version 4.4, has been imported.
       <br>
       (WB, 2005/06/30)
       </p>

  <li> <p>
       Fixed: The documentation pages did not contain information
       about the <code>Threads::spawn</code> functions. This is now fixed.
       <br>
       (WB, 2005/05/12)
       </p>

  <li> <p>
       Fixed: The step-17 example program had a stupid mistake in that it
       initializes the <code>pcout</code> variable with an uninitialized value
       of <code>mpi_communicator</code>. While this seems to have worked for
       some older MPI implementations, this is no longer the case on some
       other systems. It is now fixed.
       <br>
       (WB, 2005/04/18)
       </p>

  <li> <p>
       New: If there are references to other example programs in any of the
       step-XX programs, then they will now show up as crosslinks in the HTML
       version for simpler navigation of the tutorial.
       <br>
       (WB, 2005/04/04)
       </p>

  <li> <p>
       New: The configuration script can now enable special optimizations for
       x86-like processors if called with the <code>--with-cpu=...</code> option.
       This can decrease computation time by up to 30% on certain systems.
       The configuration flag has already existed for PowerPC64 processors, but
       was not documented. Documentation has been added.
       <br>
       (Ralf B. Schulz, 2005/03/10)
       </p>

  <li> <p>
       Fixed: When compiling shared libraries on CygWin systems, warnings
       concerning the <code>-fPIC</code> option were issued by the compiler.
       This is now fixed.       
       Also, configure now issues a message at the end that you should include
       the DLL-path in your <code>.bash_profile</code> file on these systems.
       <br>
       (Ralf B. Schulz, 2005/03/10)
       </p>

  <li> <p>
       Fixed: On certain systems running CygWin, a call to socket functions
       like <code>gethostname()</code> could cause deal.II to hang. The reason
       seems to be that the call disables long double computations which
       revert to just double precision. This then leads to an endless loop in
       computing weights for quadrature points because the convergence
       criterion cannot be reached with simple double precision. This is now
       fixed by disabling the call to <code>gethostname()</code> on these
       systems. A new preprocessor variable, <code>DEAL_II_BROKEN_SOCKETS</code>
       has been added to <code>base/config.h</code> which is defined on
       affected systems.
       <br> 
       (Ralf B. Schulz, WB, 2005/03/02)
       </p>

  <li> <p>
       Fixed: The step-16 example program wasn't listed in the
       navigation bar of the  tutorial section, although it was in the
       table of contents. This is fixed now.
       <br> 
       (WB, 2005/02/09)
       </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p>
       New: Class <code>PathSearch</code> allows to find
       a file in a sequence of directories and by appending suffixes. The class
       generates the complete file name including directory prefix and suffix
       such that it can be used for subsequently opening the file.
       <br> 
       (GK, 2005/09/14)
       </p>

  <li> <p>
       Extended: The <code>ParameterHandler</code> class now does a much better
       job generating output in the <code>print_parameters()</code> function if
       parameters have documentation strings attached to them. See the step-19
       example program output for this.
       <br>
       (WB, 2005/09/12)
       </p>

  <li> <p>
       New: There is now a namespace <code>Utilities</code> that contains
       functions that are not particularly specific to finite element codes,
       but are needed in various contexts of application programs anyway.
       <br>
       (WB, 2005/08/11)
       </p>

  <li> <p>
       Changed: <code>ParameterHandler::get_bool()</code> only accepted
       "true" or "false" as boolean values, now it also considers "yes"
       and "no" as "true" and "false", respectively.
       <br>
       (Ralf B. Schulz, 2005/07/19)
       </p>
       
  <li> <p>
       Removed: The <code>write_multigrid</code> flag in <code
       class="member">DataOutBase::DXFlags</code> 
       has been removed, since it wasn't used anywhere.
       <br>
       (WB, 2005/07/14)
       </p>

  <li> <p>
       New: A function <code>
       deal_II_exceptions::disable_abort_on_exception</code> now allows
       to disable program abortion when an assertion fails. This is used
       in regression tests.
       <br>
       (Ralf B. Schulz, 2005/07/13)
       </p>
       
  <li> <p>
       Improved: The <code>QProjector</code> now has
       functions <code>project_to_face</code> and <code
       class="class">project_to_subface</code> returning a <code
       class="class">Quadrature</code> object
       <br>
       (GK, 2005/07/12)
       </p>

  <li> <p>
       Changed: The <code
       class="member">ParameterHandler::print_parameters</code>
       function now puts two newlines between regular entries and
       subsection listings for better readability.
       <br>
       (WB, 2005/07/05)
       </p>

  <li> <p>
       Changed: The <code
       class="member">ParameterHandler::print_parameters</code>
       function now omits printing subsections if they and the contained
       sub-subsections contain no entries.
       <br>
       (WB, 2005/07/05)
       </p>

  <li> <p>
       Changed: Some of the dimension independent functions in the <code
       class="class">DataOutInterface</code> class have been moved into the
       <code>DataOutBase</code> class that does not depend on any
       dimension template parameters. Since the latter is a base class of the
       former, there should be no problems for application programs.
       <br>
       (WB, 2005/07/01)
       </p>

  <li> <p>
       Fixed: Several of the functions in class <code
       class="class">ConstantFunction</code> did not check whether the
       component handed to them was valid. In essence, even a scalar
       function could be evaluated for component 42, and the function
       ignored how many components it was declared to have. This is
       now fixed.
       <br>
       (WB, 2005/06/30)
       </p>

  <li> <p>
       Changed: The functions in class <code
       class="class">ExceptionBase</code> have been renamed to conform
       to our coding conventions that function names be all lower-case
       with words separated by underscores. This is technically an
       incompatible change, since these functions were publicly
       visible, but were usually hidden behind macros and not used
       anywhere explicitly in the library; we hope that this also
       holds for external code.
       <br>
       (WB, 2005/06/30)
       </p>

  <li> <p>
       New: The GNU C library supports access to the call stack. For
       systems using this library, the <code
       class="class">ExceptionBase</code> class now uses these
       functions to outputs the callstack when the exception is
       created.
       <br>
       (Ralf B. Schulz, 2005/06/30)
       </p>

  <li> <p>
       Improved: The class <code>Quadrature</code> has a new
       constructor generating the <tt>dim</tt>-dimensional tensor product of a
       one-dimensonal formula directly.
       <br>
       (GK, 2005/06/10)
       </p>

  <li> <p>
       New: The class <code
       class="class">Polynomials::RaviartThomas</code> implements
       Raviart-Thomas polynomials of arbitrary order in 2d and 3d.
       <br>
       (GK, 2005/05/30)
       </p>

  <li> <p>
       New: The <code>Polynomials::Lagrange</code> class
       allows generation of Lagrange interpolation polynomials for
       arbitrary point sets. In order to get better interpolation for
       higher order polynomials, <code
       class="class">QGaussLobatto</code> has been added to produce
       interpolation point sets clustering at the boundary of a cell
       (weights are not computed yet).
       <br>
       (GK, 2005/05/30)
       </p>

  <li> <p>
       New: The <code>QAnisotropic</code> class allows
       generation of quadrature rules of different order in each
       space dimension.
       <br>
       (GK, 2005/05/24)
       </p>

  <li> <p>
       New: The <code>Tensor</code> classes now
       have member functions that compute the Frobenius norm and its
       square. There are also global <code>operator*</code> functions
       that compute the contraction over a pair of indices of tensors.
       <br>
       (WB, 2005/04/06)
       </p>

  <li> <p>
       New: The <code>DataOutBase</code> class now
       allows to write data in a new intermediate format that
       basically dumps the raw information available in patches.
       <br>
       (WB, 2005/04/05)
       </p>

  <li> <p>
       New: The <code>DataOutReader</code> class allows
       to read data back in from a file into which it has been written
       using intermediate format. From there, it can then be converted
       into any of the supported graphics formats.
       <br>
       (WB, 2005/04/05)
       </p>

  <li> <p>
       New: There is now a <code
       class="member">operator &gt;&gt; (std::istream &amp;in,
       Point&lt;dim&gt; &amp;p)</code>, a function that had apparently been
       missing for a long time already.
       <br>
       (WB, 2005/04/05)
       </p>

  <li> <p>
       New: The new function <code
       class="member">double_contract</code> contracts two tensors of
       rank 4 and 2 into a tensor of rank 2, by contracting over two
       indices at the same time.
       <br> 
       (WB, 2005/03/29)
       </p>

  <li> <p>
       New: The new function <code
       class="member">TableIndicesBase::sort</code> allows to sort the indices
       in ascending order. This is useful for cases where the order of indices
       is irrelevant (for example in symmetric tables).
       <br> 
       (WB, 2005/03/29)
       </p>

  <li> <p>
       New: There is a new class <code>SymmetricTensor</code>
       that provides storage and operations for symmetric tensors.
       <br> 
       (WB, 2005/03/28)
       </p>

  <li> <p>
       New: Class <code>Subscriptor</code> receives a
       text argument allowing to identify the violating pointer more
       easily. The additional feature is being incorporated into <code
       class="class">SmartPointer</code> constructors throughout the
       library.
       <br> 
       (GK, 2005/03/16)
       </p>

  <li> <p>
       New: Class <code>FunctionParser</code>. Wrapper
       class for the fparser library (see 
       <a href="http://warp.povusers.org/FunctionParser/">
       http://warp.povusers.org/FunctionParser/</a>).
       <br> 
       (Luca Heltai, 2005/03/07).
       </p>

  <li> <p>
       Fixed: The class <code
       class="class">MultipleParameterLoop::UserClass</code> had only 
       virtual abstract functions but no virtual destructor. This caused
       warnings with some compilers, and is generally bad practice
       anyway. This is now fixed. The same holds with respect to the class
       <code>DataOutInterface</code>.
       <br> 
       (WB, 2005/02/20)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p>
       New: The new function <code
       class="member">PetscWrappers::VectorBase::set</code> allows to set
       several vector elements at once.
       <br> 
       (WB, 2005/08/10)
       </p>

  <li> <p>
       New: There are now functions <code
       class="member">PetscWrappers::MatrixBase::clear_row</code> and <code
       class="member">clear_rows</code> that allow to set the elements of a row
       of a matrix to zero, without having to traverse the individual
       elements.
       <br> 
       (WB, 2005/08/10)
       </p>

  <li> <p>
       Fixed: <code>SparsityPattern</code>::<code
       class="member">block_write</code> neither wrote the number of
       columns of the pattern, nor the flag for optimizing the
       diagonal, and <code>SparsityPattern</code>::<code
       class="member">block_read</code> did not attempt to read this
       data. Both is fixed now.
       <br>
       (GK, 2005/06/24)
       </p>

  <li> <p>
       Improved: <code>SparseDirectUMFPACK</code> now
       has the complete interface of a preconditioner
       class. Therefore, it can be used in iterative solvers and in
       multigrid smoothers.
       <br> 
       (GK, 2005/05/10)
       </p>

  <li> <p>
       Fixed: The PETSc matrix iterators had trouble when some rows of a
       matrix were empty. This has now been mostly fixed.
       <br> 
       (WB, 2005/05/02)
       </p>

  <li> <p>
       Fixed: The <code>PETScWrappers::MPI::SparseMatrix</code>
       class had functions that allow to take a pre-existing sparsity pattern
       as the basis for entry allocation. These functions had an option to
       allow pre-setting these entries in the underlying data structures, but
       this option was disabled by default because of unresolved questions
       about its effectiveness. This has now been fixed: The code now properly
       initializes these elements, and makes the resulting matrix much faster
       to use.
       <br> 
       (WB, 2005/04/04)
       </p>

  <li> <p>
       New: The <code>ProductSparseMatrix</code>
       implements the product of two rectangular sparse matrices with
       the same <code>value_type</code>
       <br> 
       (GK, 2005/03/11)
       </p>

  <li> <p>
       New: The <code>PreconditionRichardson</code>
       implements a Richardson preconditioner.
       <br> 
       (GK, 2005/03/10)
       </p>

  <li> <p> Fixed: The <code>BlockSparseMatrix</code>
       class had no local typedef <code
       class="member">value_type</code> like all other classes, which
       made it a little awkward to use in some places. This has now
       been fixed.

       <br> 
       (WB, 2005/03/03)
       </p>

  <li> <p> Fixed: The <code>PETScWrappers</code>::<code
       class="member">MatrixBase</code> class documented that adding
       or setting a value that hasn't been in the sparsity pattern
       before will lead to an exception being thrown. This is of
       course wrong: PETSc allocates matrix entries dynamically, as
       needed. The documentation is now fixed.
       <br> (WB, 2005/03/03)
       </p>

  <li> <p> New: The <code>SparseMatrix</code> iterators
       had no <code>operator &gt;</code>, only an <code
       class="member">operator &lt;</code>. The missing operator is
       now implemented. The same holds for the <code
       class="class">FullMatrix</code> class.
       <br> (WB, 2005/03/01)
       </p>

  <li> <p> Fixed: The <code>SparseMatrix</code>
       iterators could not be compared using <code
       class="member">operator &lt;</code>: the compiler complained
       that a private member was accessed. This is now fixed.
       <br> (WB, 2005/03/01)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>

  <li> <p> Removed: The class <code
  class="class">FiniteElementBase</code> has been removed and all its
  functions are now in <code>FiniteElement</code>.
  <br>
  (GK, 2005/08/25)
  </p>

  <li> <p> New: class <code>DoFTools</code> now has two
  functions <code>compute_row_length_vector</code>, one
  for equations and one for systems. These give a much finer estimate
  for the size of the sparsity pattern required than <code
  class="class">DoFHandler</code>::<code
  class="member">max_couplings_between_dofs</code>. This is
  particularly true for systems with few couplings.
  <br>
  (GK, 2005/08/24)
  </p>

  <li> <p>
       Remove: Due to the reimplementation of the 3d reordering
       algorithm in <code>GridReordering</code> two
       years ago, the <code>CellData::rotate</code>
       function got obsolete. <code
       class="member">CellData::rotate</code> is now removed.
       <br>
       (RH, 2005/08/09)
       </p>

  <li> <p>
       Fixed: It was possible to compare iterators into different
       <code>DoFHandler</code> objects for equality and
       inequality as long as the underlying triangulation object was
       the same. This doesn't make much sense, however, and even
       sometimes declared to iterators equal if only they had the same
       triangulation, refinement level, and index. Comparing iterators
       into different containers is now completely disallowed.
       <br> 
       (WB, 2005/08/08)
       </p>

  <li> <p>
       New: The <code>DataOutStack::add_data_vector</code>
       function now also allows to give only a single name even if the data is
       vector-valued. It then makes up names for the individual components by
       adding numbers to the name, just like the <code
       class="class">DataOut</code> class already does for a long time.
       <br> 
       (WB, 2005/07/29)
       </p>

  <li> <p>
       New: The <code>DataOutStack</code> class can now also be
       used to stack two-dimensional time or parameter dependent data into a 3d
       output.
       <br> 
       (WB, 2005/07/25)
       </p>

  <li> <p>
       New: Function <code>FETools</code>::<code
       class="member">compute_face_embedding_matrices</code> allows for
       automatic computation of embedding matrices on faces under the sole
       assumption that the coarse grid space is embedded into the fine
       grid space. In particular, no interpolation points are
       required.<p>
       
       <p>Using this function, constraint matrices can be computed in
       a general way.
       <br> 
       (GK, 2005/07/14)
       </p>

  <li> <p>
       Improved: All <code>GridIn</code>::<code
       class="member">read_*</code> functions now call <code
       class="class">GridReordering</code>::<code
       class="member">invert_all_cells_of_negative_grid</code>. This
       way, also misoriented grids are represented in the right
       orientation within deal.II.
       <br>
       (RH 2005/07/06)
       </p>

  <li> <p>
       Improved: <code
       class="class">Triangulation&lt;3&gt;</code>::<code
       class="member">create_triangulation</code> now checks that all
       cells have positive volume. If not, an exception is thrown. In
       that case use the <code
       class="class">GridReordering</code>::<code
       class="member">invert_all_cells_of_negative_grid</code>
       function, see below, to fix this.
       <br>
       (RH 2005/07/06)
       </p>

  <li> <p>
       New: There is now a <code
       class="class">GridReordering</code>::<code
       class="member">invert_all_cells_of_negative_grid</code>
       function implemented in 3d which checks if some/all cells have
       negative volumes. If all cell have negative volume then the
       whole grid is reoriented. An assertion is thrown if only a
       subset of cells have negative volumes, as then the grid might
       be broken.
       <br>
       (RH 2005/07/06)
       </p>

  <li> <p>
       New: There is now a new <code
       class="class">GridTools</code>::<code
       class="member">cell_measure</code> function. However, it is
       mostly for internal use. Use <code>cell->measure()</code>
       instead.
       <br>
       (RH 2005/07/06)
       </p>

  <li> <p>
       Improved: <code>cell->measure()</code> used to give the
       absolute value of the cell measure. It now gives the measure
       with the correct sign. This is useful to find wrongly oriented
       cells in 3d with negative volumes.
       <br>
       (RH 2005/07/05)
       </p>

  <li> <p>
       Improved: The class <code
       class="class">FiniteElementData</code> now stores information
       on the Sobolev space a finite element space conforms with.
       <br>
       (GK 2005/06/29)
       </p>
  
  <li> <p>
       New: Added function <code
       class="class">FiniteElementBase::component_to_system_index</code>
       which was referred to in the documentation, but was missing.
       <br>
       (Ralf B. Schulz, 2005/06/30)
       </p>

  <li> <p>
       Improved: The function <code
       class="class">DoFTools</code>::<code
       class="member">count_dofs_per_component</code> and its
       counterpart in <code>MGTools</code> got an
       additional argument. This argument allows to count the degrees
       of freedom of nonprimitive vector valued elements only once,
       instead of in every component.  Although this argument defaults
       to the previous behavior, it had to be put into the argument
       list ahead of <tt>target_component</tt> in order to make use of
       default arguments more efficiently. The old order of arguments
       can still be used through a wrapper function.
       <br>
       (GK 2005/06/22)
       </p>

  <li> <p>
       Improved: The <code>GeometryInfo</code>::<code
       class="member">child_cell_on_face</code>, <code
       class="member">face_to_cell_vertices</code> and <code
       class="member">face_to_cell_lines</code> now have an additional
       <code>face_orientation</code> argument, which has no effect in
       2d and which defaults to <code>true</code> (standard
       orientation) in 3d. Now these functions return the right values
       also for the case of non-standard oriented faces in 3d. This
       avoids several awful pieces of code including questioning
       face_orientation and using child_switch_tables.
       <br> 
       (RH, 2005/06/16)
       </p>

  <li> <p>
       Changed: The method <code>FETools</code>::<code
       class="member">lexicographic_to_hierarchic_numbering</code> now
       takes a <code>FiniteElementData</code> instead of
       a <code>FE_Q</code> object. Now this function can
       also be called by the <code>FE_Q</code>
       constructor which avoids code duplication.
       <br> 
       (RH, 2005/06/13)
       </p>

  <li> <p>
       New: The method <code>create_mass_matrix</code>
       in class <code>MatrixCreator</code> can now handle 
       vector valued finite elements. A similar change was applied
       to <code>create_right_hand_side</code> and
       <code>create_boundary_right_hand_side</code>
       in <code>VectorTools</code>. These two changes
       now make the <code>project</code> function work also
       for Raviart-Thomas elements and other vector valued FEs. This
       is very useful, if some initial conditions have to be specified.
       <br> 
       (Oliver Kayser-Herold, 2005/06/03)
       </p>

  <li> <p>
       New: The <code>DataOut</code> class
       now supports Eulerian Mappings. If a solution is computed
       on a deformed mesh, the output file generated by the DataOut
       now shows the solution also on the deformed mesh. This
       is an option and requires the mapping to be specified as
       additional parameter to <code>build_patches</code>.
       <br> 
       (Oliver Kayser-Herold, 2005/05/31)
       </p>

  <li> <p>
       New: The <code>MappingQ1Eulerian</code> class can
       now cope with different vector types for the Euler vector.
       This is useful if it should be used with the PETSc wrapper
       classes.  The desired vector type can be specified as template
       parameter.
       <br>
       (Oliver Kayser-Herold, 2005/05/31)
       </p>

  <li> <p>
       New: The <code>FE_RaviartThomasNodal</code>
       implements Raviart-Thomas elements using function values in
       Gauss quadrature points on edges and in the interior for its
       node values. The implementation is restricted to Cartesian mesh
       cells right now, but works in 2D and 3D.
       <br> 
       (GK, 2005/05/30)
       </p>

  <li> <p>
       Improved: The <code>Mapping::transform_*</code>
       functions accept <code>VectorSlice</code> instead
       of <code>Vector</code>, thus allowing more flexibility.
       <br> 
       (GK, 2005/05/24)
       </p>

  <li> <p>
       New: The <code>MatrixTools::apply_boundary_values</code>
       function now also works for PETSc sequential and parallel matrices.
       <br> 
       (WB, 2005/05/05)
       </p>

  <li> <p>
       Improved: The function <code>GridIn</code>::<code
       class="member">read</code> now searches for files using the
       mechanism provided by the class <code
       class="class">PathSearch</code>. Furthermore, a library of input
       meshes has been started in <code>lib/meshes</code>.
       <br> 
       (GK, 2005/05/03)
       </p>

  <li> <p>
       Fixed: The <code>DataOut</code> class did not work
       properly if the <code>DataOut::first_cell</code> and
       <code>DataOut::next_cell</code> functions were
       overloaded and cell data was to be output; in that case, data from the
       wrong cells was written out. This is now fixed. In contrast to this,
       nodal data was never affected.
       <br> 
       (WB, 2005/04/20)
       </p>

  <li> <p>
       Improved: By employing the new <code
       class="class">GeometryInfo</code>::<code
       class="member">line_to_cell_vertices</code>, <code
       class="class">GeometryInfo</code>::<code
       class="member">face_to_cell_vertices</code> and <code
       class="class">GeometryInfo</code>::<code
       class="member">face_to_cell_lines</code> functions the <code
       class="class">Triangulation</code>::<code
       class="member">create_triangulation</code> functions are now
       implemented independent of given conventions for the numbering
       of vertices, lines and faces.
       <br>
       (RH, 2005/03/11)
       </p>

  <li> <p>
       New: The new <code>GeometryInfo</code>::<code
       class="member">line_to_cell_vertices</code> function maps line
       vertex numbers to cell vertex numbers.
       <br> 
       (RH, 2005/03/11)
       </p>

  <li> <p>
       New: The new <code>GeometryInfo</code>::<code
       class="member">face_to_cell_lines</code> function maps face
       line numbers to cell line numbers.
       <br> 
       (RH, 2005/03/11)
       </p>

  <li> <p>
       New: The new <code>GeometryInfo</code>::<code
       class="member">face_to_cell_vertices</code> function maps face
       vertex numbers to cell vertex numbers.
       <br> 
       (RH, 2005/03/11)
       </p>

  <li> <p>
       Fixed: There was a bug in the <code
       class="class">KellyErrorEstimator</code> class that resulted in
       assertions being thrown when run with multithreading
       enabled. This is now fixed.
       <br> 
       (WB, 2005/03/10)
       </p>

  <li> <p>
       Changed: The <code>Triangulation<2></code>::<code
       class="member">execute_refinement</code> function has been
       re-implemented to accommodate future developments. This change
       results in different ordering and indexing of lines and
       vertices. This leads to a change in the ordering of vertices in
       output files generated by <code>GridOut</code>
       for refined grids.
       <br> 
       (RH, 2005/03/09)
       </p>

  <li> <p>
       Fixed: The class <code>MGDoFHandler</code> had trouble
       when it encountered a triangulation that had unused vertices, something
       that happens when one coarsens an existing triangulation. In
       that case, it would throw unjustified exceptions. This is now
       fixed.
       <br> 
       (WB, 2005/03/01)
       </p>

  <li> <p>
       Fixed: The class <code>Triangulation::RefinementListener</code> had only
       virtual abstract functions but no virtual destructor. This caused
       warnings with some compilers, and is generally bad practice
       anyway. This is now fixed.
       <br> 
       (WB, 2005/02/22)
       </p>

  <li> <p>
       New: Function <code>FETools</code>::<code
       class="member">compute_embedding_matrices</code> allows for
       automatic computation of embedding matrices under the sole
       assumption that the coarse grid space is embedded into the fine
       grid space. In particular, no interpolation points are required.
       <br> 
       (GK, 2005/02/08)
       </p>

  <li> <p>
       Fixed: Several wrong assertions in the Raviart-Thomas finite element
       class have been fixed, that were triggered when integrating face terms.
       <br> 
       (Oliver Kayser-Herold, 2005/01/24; 
        WB, 2005/01/31)
       </p>
</ol>


*/
