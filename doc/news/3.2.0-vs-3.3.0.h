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
 * @page changes_between_3_2_and_3_3 Changes between Version 3.2 and 3.3

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li> <p>
       New: Output for 
       <a href="http://www.amtec.org" target="_top">Tecplot</a> has
       been added. It can be used by choosing output format «tecplot».
       <br>
       (<a href="mailto:benkirk@cfdlab.ae.utexas.edu">Benjamin Shelton Kirk</a> 2002/01/29)
       </p>

  <li> <p> New: configuration detects whether the compiler has the
       include file <code>&lt;ostream&gt;</code>. Most files in the
       library then include this file over
       <code>&lt;iostream&gt;</code> to save compile time.
       <br>
       (WB 2002/01/25)
       </p>

  <li> <p> Fixed: All example and test programs as well as a number of
       large applications have been checked against the memory checker
       "purify". Only three memory leaks were found and fixed. We
       believe that no major leaks exist in the library any more.
       <br>
       (WB 2001/12/07)
       </p>

  <li> <p>
       New: Output for 
       <a href="http://www.opendx.org" target="_top">OpenDX</a> has
       been added. It can be used by choosing output format «dx» (not
       yet for grid output). The data format is very basic now, but it
       is planned to improve this to make use of the excellent
       capabilites of OpenDX.
       </p>
       <p>
       Additionally, a directory <tt>contrib/dx</tt> has been added,
       containing visual programs for OpenDX. Programs for the
       existing output of data as a single vector are found in the
       subdirectory <tt>single</tt>.
       <br>
       (GK 2001/12/07)
       </p>

  <li> <p> Fixed: Previously, the \$(INCLUDE) variable in Makefiles
       included the values of the \$INCLUDE environment variable. This
       is not desirable, since the compiler evaluates that variable
       anyway and the Makefile variable has <code>-I</code> prefixed
       to all paths while the environment variable has not.
       <br>
       (WB 2001/11/29)
       </p>

  <li> <p> Removed: the option to generate printable documentation was
       removed. Since this comprised approximately 2000 pages and
       since we believe that the online documentation is rather good,
       this is probably no big loss.
       <br>
       (WB 2001/11/29)
       </p>

  <li> <p> New: example program step-12. Discontinuous Galerkin
       discretization.
       <br>
       (RH 2001/11/21)
       </p>

  <li> <p>
       New: deal.II now uses a file
       <code>config.h</code> for most global preprocessor defines,
       instead of an overly long list of compiler flags given on the
       command line.
       <br>
       (WB 2001/10/27)
       </p>

  <li> <p>
       Changed: If available, the library now uses the C++ standard
       classes <code>istringstream</code> and <code
       class="class">ostringstream</code> over the old classes
       <code>i/ostrstream</code>. The ./configure script
       finds out whether the new classes exist, or whether the
       backward compatibility classes are to be used.
       <br>
       (WB 2001/10/25)
       </p>

  <li> <p>
       New: the ./configure script now recognizes gcc3.1
       (i.e. presently prereleases of it) and sets compilation flags
       accordingly.
       <br>
       (WB 2001/10/25)
       </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p>
       Fixed: For Mac OS X, the file <code>base/source/log.cc</code>
       needed to include <code>sys/time.h</code> in addition to
       <code>sys/resource.h</code>.
       <br>
       (Alexis Herault, WB 2002/01/28)
       </p>

  <li> <p> Fixed: Private definitions of the copy constructor and
       assignment operator of the <code
       class="class">ParameterHandler</code> class are now included in
       order to inhibit the (unintentional) use of default copy
       constructors.
       <br>
       (RH 2002/01/22)
       </p>

  <li> <p>
       Improved: The cut-off functions <code
       class="class">Functios::CutOffFunctionLinfty</code>, <code
       class="class">Functios::CutOffFunctionW1</code>, and <code
       class="class">Functios::CutOffFunctionCinfty</code> can be
       vector-valued now and optionally only a single componente can
       be selected.
       <br>
       (GK 2002/01/10)
       </p>

  <li> <p>
       New: the <code
       class="member">deal_II_exceptions::set_additional_assert_output</code>
       function allows to set additional output to be printed upon
       triggering an <code>Assert()</code> call. This
       is helpful for parallel applications where you only see the
       text of the message but do not know from which cluster node it
       stems.
       <br>
       (WB 2002/01/10)
       </p>

  <li> <p>
       Changed: when an assertion fails in the <code
       class="member">Assert()</code> function, the program is usually
       aborted. Don't abort it any more if there is an active C++
       exception somewhere since we would lose its message if we
       aborted the program. In that case only report the error and
       write out an indication why we do not abort the program any
       more. On the other hand, also suppress output of further failed
       <code>Assert()</code> calls, since they often
       are follow-ups of the first one.
       <br>
       (WB 2002/01/09)
       </p>

  <li> <p>
       New: <code>ExcFileNotOpen</code> can be used
       after initializing an <code>fstream</code>
       object. This allows to avoid some cryptic <code
       class="class">ExcIO</code>s.
       <br>
       (GK 2001/12/18)
       </p>

  <li> <p>
       Changed: The <code>OutputStyle</code> enum used
       to indicate the output format has been moved into the
       <code>ParameterHandler</code> class.
       <br>
       (WB 2001/11/30)
       </p>

  <li> <p>
       Fixed: In the <code>ParameterHandler</code>
       class, we leaked 8 or 12 bytes of memory per declared
       parameter. This is now fixed.
       <br>
       (WB 2001/11/28)
       </p>

  <li> <p>
       New: <code>Functions::CutOffFunctionCinfty</code>,
       <code>Functions::CutOffFunctionW1</code>, and
       <code>Functions::CutOffFunctionLinfty</code>
       implement functions with support in an arbitrary ball and
       differentiability as indicated by their name
       <br>
       (GK 2001/10/24)
       </p>

  <li> <p>
       Fixed: The <code>DataOutBase::EpsFlags</code>
       class forgot to declare the reverse grey scale function as one
       possible input for the color function for the
       <code>ParameterHandler</code> class. This is now
       possible.
       <br>
       (WB 2001/10/10)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p>

       Improved: all sparsity pattern classes have a function <code
       class="member">exists</code>, allowing you to check whether a
       certain index pair has been allocated in the pattern.
       <br>
       (GK 2002/02/01)
       </p>

  <li> <p>
       Fixed: Allocation of temporary vectors in <code
       class="member">FilteredMatrix::allocate_tmp_vector</code>
       is now faster since it does no more copy the value of the
       template vector.
       <br>
       (WB 2001/11/22)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">FilteredMatrix::allocate_tmp_vector</code>
       function had a bug with block vectors.
       <br>
       (WB 2001/11/22)
       </p>

  <li> <p>
       Improved: reinit function of <code>Vector</code>
       and <code>BlockVector</code> allows use of a
       vector with different number type.
       <br>
       (GK 2001/11/21)
       </p>

  <li> <p>
       Fixed: when checking for convergence in linear solvers in
       <code>SolverControl::check</code>, we first
       checked whether the maximal iteration count was reached, and
       only then whether the target residual was achieved. In cases,
       where the target residual was only reached in the very last
       allowed iteration, this led to a failure notice of the linear
       solver, rather than to a success message. This is now fixed.
       <br>
       (WB 2001/11/19)
       </p>

  <li> <p>
       New: the <code>SparseDirectMA27/47</code> classes
       now provide access to Mutex locks for external
       synchronisation. 
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Fixed: an error in the definition of the <code
       class="member">SolverMinRes::solve</code> function prevented
       its compilation under some circumstances.
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Fixed: upon breakdown, the <code
       class="class">SolverBicgstab</code> forgot to increment the
       iteration counter for the breakdown cycle. This is now fixed.
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Improved: class <code>SolverGMRES</code> accepts
       a parameter in <code>AdditionalData</code>,
       allowing for right preconditioning.
       <br>
       (GK 2001/11/09)
       </p>

  <li> <p>
       Fixed: class <code>SparsityPattern</code> can
       handle rows of length zero now. For quadratic matrices, these
       rows are changed to length one, since a diagonal element must
       be stored.
       <br>
       (GK 2001/10/11)
       </p>

  <li> <p>
       New: The <code>BlockVector</code> now has a
       full-fledged random access iterator type, working in exactly
       the same way as the iterators of the C++ standard containers.
       <br>
       (WB 2001/09/28)
       </p>

  <li> <p> New: <code>Vector</code>::<code
       class="member">operator *</code> is now templatized, allowing
       for scalar products of vectors with different underlying types.
       <br>
       (WB 2001/09/27)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p>
       Changed: The classes <code>FEQ1</code>-<code
       class="class">FEQ4</code>, <code
       class="class">FEDG_Q1</code><code>FEDG_Q4</code>
       as well as the files with their definitions,
       <tt>fe/fe_lib.lagrange.h</tt> and <tt>fe/fe_lib.dg.h</tt>
       ceased to exist. They had been left in for backward
       compatibility in an earlier version, but their existence is
       more confusing than helpful. Please change your code to use the
       classes <code>FE_Q</code> and <code
       class="class">FE_DGQ</code>, respectively.
       <br>
       (GK 2002/02/01)
       </p>

  <li> <p>
       New: The <code>FilteredIterator</code> class
       provides a view on ranges of iterators by iterating over only
       those objects that satisfy a certain predicate.
       <br>
       (WB 2002/01/07)
       </p>

  <li> <p>
       Improved: It is now possible to read in unconnected domains
       through the <code>GridIn</code> class, since
       renumbering the cells now also works for these domains.
       <br>
       (<a href="mailto:ms@biomech.tu-graz.ac.at">Michael Stadler</a> 2001/12/14)
       </p>

  <li> <p>
       Improved: Both functions <code
       class="class">VectorTools</code>::<code
       class="member">compute_mean_value</code> take ingoing and
       outgoing vector types as template arguments. This allows
       applying them to <code>BlockVector</code>.
       <br>
       (GK 2001/12/07)
       </p>

  <li> <p>
       New: <code>GridGenerator</code> has a function
       <code>cylinder</code> for cylinders in three
       space dimensions. Accoridngly, a class <code
       class="class">CylinderBoundary</code> has been created.
       <br>
       (GK 2001/12/07)
       </p>

  <li> <p>
       New: <code>FiniteElement</code>::<code
       class="member">has_support_on_face</code> allows to check
       whether a shape function has non-zero values on a certain face
       of a cell.
       <br>
       (GK 2001/12/04)
       </p>

  <li> <p>
       Changed: The <code>IteratorState</code> enum used
       to indicate the state in which an iterator can be is now
       enclosed in a namespace of the same name, to take its members
       out of the global namespace. When using one of these members,
       you now have to prefix it by <code
       class="class">IteratorState::</code>.
       <br>
       (WB 2001/11/30)
       </p>

  <li> <p>
       Changed: The <code>NormType</code> enum used to
       indicate the norm the <code
       class="member">VectorTools::integrate_difference</code>
       function shall integrate is moved from the global namespace
       into the <code>VectorTools</code> class. You
       therefore have to prefix the members of this enum by the
       respective class name.
       <br>
       (WB 2001/11/29)
       </p>

  <li> <p>
       Fixed: The functions <code
       class="member">Mapping::transform_unit_to_real_cell</code>
       leaked some memory. This is now fixed.
       <br>
       (RH, WB 2001/11/28)
       </p>

  <li> <p>
       Fixed: The <code>DoFHandler</code> class had a
       memory leak. This is now fixed. Likewise for the <code
       class="class">MGDoFHandler</code> class.
       <br>
       (WB 2001/11/28)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">GridRefinement::refine</code>
       function failed when the threshold was zero. This is now fixed.
       <br>
       (RH 2001/11/26)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">MappingQ<dim>::transform_real_to_unit_cell</code>
       function failed on a very unusual cell. This is now fixed.
       <br>
       (RH 2001/11/26)
       </p>

  <li> <p>
       New: The new <code
       class="member">CellAccessor<dim>::neighbor_of_coarser_neighbor</code>
       function returns where to find the present cell from a coarser
       neighbor.
       <br>
       (RH 2001/11/21)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">GridRefinement::refine_fixed_fraction</code>
       function sometimes had problems when indicators vary over
       several orders of magnitude, due to roundoff. This is now
       fixed. 
       <br>
       (WB 2001/11/05)
       </p>

  <li> <p> 
       New: <code
       class="member">DoFTools::extract_subdomain_dofs</code> selects
       those degrees of freedom which are located on cells with a
       specified subdomain id.
       <br>
       (WB 2001/10/27)
       </p>

  <li> <p>
       New: Cells now have an additional property
       <em>subdomain_id</em> which can be used in parallel
       computations to identify which cells are handled on which
       processor. These flags are read and set using the functions
       <code>cell->subdomain_id()</code> and <code
       class="member">cell->set_subdomain_id(new_id)</code>. The
       subdomain ids are unsigned integers, so should be sufficiently
       large also for larger numbers of subdomains.
       <br>
       (WB 2001/10/27)
       </p>

  <li> <p>
       Fixed: the <code
       class="member">GridGenerator::hyper_rectangle</code> function
       was broken in 3d.
       <br>
       (WB 2001/10/26)
       </p>

  <li> <p>
       Improved: Both functions <code
       class="class">DataOut_DoFData</code>::<code
       class="member">add_data_vector</code> accepts <code
       class="class">BlockVector</code> as argument.
       <br>
       (GK 2001/10/12)
       </p>

  <li> <p>
       Improved: Both functions <code
       class="class">VectorTools</code>::<code
       class="member">integrate_difference</code> take ingoing and
       outgoing vector types as template arguments. This allows
       applying them to <code>BlockVector</code> and of
       outputting a vector of doubles suitable for <code
       class="class">DataOut</code>.
       <br>
       (GK 2001/10/12)
       </p>

  <li> <p>
       Fixed: Functions creating sparsity patterns for DG elements in
       <code>DoFTools</code> get the pattern type as
       template argument, too..
       <br>
       (GK 2001/10/01)
       </p>

  <li> <p>
       Fixed: the iterator category template base class of grid
       iterators was incorrectly set.
       <br>
       (WB 2001/09/28)
       </p>
</ol>

*/
