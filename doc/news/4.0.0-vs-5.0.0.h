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
 * @page changes_between_4_0_and_5_0 Changes between Version 4.0 and 5.0
 
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
  <li> <p> Removed: All the matrix classes have functions <code>reinit</code> that are used to resize the
       matrix. However, the sparse matrix classes had an equally named
       function without arguments that basically left the size of the matrix
       unchanged but set all elements of the matrix to zero. It could also be
       abused to notify the matrix that the sparsity pattern on which it is
       built has changed, an inherently risky procedure. The no-argument code
       <code>reinit</code> function has therefore been removed to
       avoid confusion with the <code>reinit</code> functions
       that take arguments. Instead, you should now use <code>matrix=0</code> to simply set all elements of the
       matrix to zero. If you want to notify a sparse matrix that its sparsity
       pattern has changed, use the <code>reinit(SparsityPattern)</code> function.
       <br> 
       (WB 2004/05/10)
       </p>

  <li> <p> Removed: All the vector and block vector classes as well as
       the <code>FullMatrix</code> class (the latter
       through its <code>Table</code> base class) had a
       member function <code>clear</code> which simply
       resets all values of the vector or matrix to zero. It did not
       change the size of the object, though. This was confusing,
       since the standard C++ container classes implement the
       semantics that the <code>clear</code> functions
       delete all entries of the containers, i.e. resize it to zero,
       and we implemented similar semantics also for the <code>SparseMatrix</code>, <code>DoFHandler</code>, <code>ConstraintMatrix</code> and various other
       classes.
       <br>
       To avoid this confusion in the future, the <code>clear</code> functions have been dropped from
       the interface of the vector and full matrix classes, and the
       remaining instances where deal.II classes have a function of
       that name really mean that the object is reset to its virginal
       state. To set all
       elements of a matrix or vector to zero without changing its size, the
       recommended way is to use the more obvious notation <code>vector=0</code> and <code>matrix=0</code>. To
       reset the elements of a table over arbitrary objects, use
       <code>Table&lt;T&gt;::reset_values()</code>.
       <br> 
       (WB 2004/05/10)
       </p>

  <li> <p> Removed: The <code>SparseLUDecomposition::reinit</code> and <code>SparseMIC::reinit</code> functions without
       argument had been deprecated before, and have now been removed.
       <br> 
       (WB 2004/05/10)
       </p>

  <li> <p> Changed: the template parameter of <code>MGTransferPrebuilt</code> is now the complete vector
       type, not just the value type of the vector. This allows to operate
       on <code>Vector</code> as well as on <code
       class="class">BlockVector</code>. Unfortunately, the untested class
       <code>MGTransferBlock</code> underwent some more
       changes without testing, such that it should be used with high
       caution.
       <br>
       (GK 2004/04/01)
       </p>

  <li> <p> Changed: The <code>FiniteElement</code> classes had a
       function <code>restrict</code> that returns the
       restriction matrix from children to the mother cell. Unfortunately,
       <code>restrict</code> has become a keyword in recent standards of the C
       language, and some C++ compilers have picked this up. The function has
       therefore been renamed <code
       class="member">get_restriction_matrix</code>, which also better
       reflects what it is actually doing. To keep consistent, we have also
       rename the corresponding function <code>prolongate</code> to <code>get_prolongation_matrix</code>. 
       <br>
       (WB 2004/02/29)
       </p>

  <li> <p>
       Fixed and changed: The <code>SolutionTransfer</code><code>::(refine_)interpolate(const Vector &in, Vector
       &out)</code> functions now require the <code>in</code> and <code>out</code>
       vectors being already of right sizes,
       i.e. <code>in.size()=n_dofs_old</code> and
       <code>out.size()=n_dofs_refined</code>. Furthermore, the <code>SolutionTransfer</code><code>::(refine_)interpolate(const
       vector&lt;Vector&gt; &all_in, vector&lt;Vector&gt;
       &all_out)</code> now check that the number of in and output
       vectors are the same, i.e.
       <code>all_in.size()=all_out.size()</code>.
       <br>
       (RH 2003/10/24)
       </p>

  <li> <p>
       Changed: The <code>QProjector</code> has functions that
       project a lower-dimensional quadrature formula onto all faces or
       subfaces of a cell. In 3d, it now does this but also adds projections of
       these quadrature formula onto the faces from the other side. We need
       this due to the fact that we now support faces in 3d that have a normal
       vector opposite to the standard direction.
       <br> 
       (WB 2003/10/17)
       </p>

  <li> <p> Moved and changed: The header file
       <tt>include/numerics/dof_renumbering.h</tt> has been moved to the 
       directory <tt>include/dofs</tt>, where it logically
       belongs. Furthermore, the sorting parameter of the function <code>DoFRenumbering</code><code>::component_wise</code> has changed its meaning. See
       the reference documentation for details.
       <br>
       (GK 2003/07/03)
       </p>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>

  <li> <p> New: After the documentation tool for deal.II has been
  changed to <a href="http://www.doxygen.org">Doxygen</a>, it is delivered in two
  tar-files. Additional to the traditional source tarball, the preprocessed
  documentation is available ready for reading with a web browser.
  <br>
  (GK 2004/05/04)
  </p>

  <li> <p> New:
       The step-15 example is about solving a nonlinear 1d problem, and
       dealing with transporting a solution across mesh refinement. Step-16 is
       still not finished.
       <br>
       (WB 2004/04/17)
       </p>

  <li> <p> New:
       The step-17 example program shows how to use the new PETSc wrapper
       classes, and how to do large-scale computations with up to many
       millions of unknowns on a cluster of computers. This program shows
       that deal.II is well-suited for truly massive parallel
       computations. The step-15 and step-16 programs have not yet been
       finished (despite having been started long before step-17), which
       explains the holes in the numbering of example programs.
       <br>
       (WB 2004/04/12)
       </p>

  <li> <p> New: deal.II is now able to interface to the 
       <a href="http://www-users.cs.umn.edu/~karypis/metis/index.html"
       target="_top">METIS</a> library to generate domain partitionings. This
       is enabled if a METIS installation is detected, which either happens
       automatically in <code>./configure</code>, or
       by passing a value to the switch <code>--with-metis</code> to configure
       the path of a METIS installation. For more information see the README
       file. 
       <br>
       (WB 2004/03/08)
       </p>

  <li> <p>
       New: We now support MIPSpro compiler version 7.41 from SGI. deal.II
       now runs on IRIX64 machines in 64-bit mode.
       <br>
       Please note, that we cannot support earlier MIPSpro compilers
       because the compiler version 7.3 was not C++ standard
       conforming. Version 7.4 is standard conforming but still too
       buggy.
       <br>
       (RH 2004/03/08)
       </p>

  <li> <p> New: deal.II now comes with a complete set of
       wrappers classes for <a href="http://www.mcs.anl.gov/petsc/"
       target="_top">PETSc</a> vectors, matrices, linear solvers and 
       preconditioners. Many of the algorithms in deal.II have also been
       updated to make use of these wrappers. All of this is only enabled if a
       PETSc installation is detected. This either happens automatically in
       <code>./configure</code>, or
       by passing values to the switches <code>--with-petsc</code> and
       <code>--with-petsc-arch</code> to configure path and architecture of an
       existing PETSc installation. If these switches are not used, then
       environment variables are searched for this information. For more
       information see the README file.
       <br>
       (WB 2004/03/01)
       </p>

  <li> <p>
       Changed: The part of the boost library in the <tt>contrib</tt>
       directory is now updated to boost version 1-30.2. We include
       only a minimal part (about 3% in size) of boost which is needed
       to compile deal.II. For the case the compilation of deal.II on
       some compiler/platforms requires boost files in addition to
       those included in the <tt>contrib</tt> directory please report
       the missing boost files to the deal.II mailing list.
       <br> 
       (RH 2004/02/16)
       </p>

  <li> <p>
       Changed: We don't support compilation by Intel'c icc compiler version 5
       anymore. Since the present version of this compiler is 8, this is
       probably not a real problem.
       <br> 
       (WB 2003/12/20)
       </p>

  <li> <p>
       Fixed: <code>step-9</code> had the computation of the value of the
       solution in the mid-point of the neighbor cell wrong. This is now
       fixed. Consequently, the resulting mesh also looks much nicer now (much
       smoother). 
       <br> 
       (<a href="mailto:werner.scholz@tuwien.ac.at">Werner Scholz</a>
        2003/12/11)
       </p>

  <li> <p>
       New: The <code>config.h</code> file now declares a variable <code
       class="member">deal_II_numbers::invalid_unsigned_int</code>.
       It is a representation of the largest number that can be put into an
       unsigned integer. This value is widely used throughout the library as a
       marker for an invalid unsigned integer value, such as an invalid array
       index, an invalid array size, and the like.
       <br> 
       (WB 2003/10/24)
       </p>

  <li> <p>
       Augmented: The <code
       class="member">GeometryInfo::child_cell_on_face</code> 
       results in a result that might not be what you expect in 3d in some
       cases. A warning has been added to the documentation, as well as a
       reference to the function to call instead.
       <br> 
       (WB 2003/10/17)
       </p>

  <li> <p> Fixed: The step-14 program had a bug in the rare case that
       there are more CPUs in a machine than there are cells. This is
       now fixed.
       <br>
       (WB 2003/09/23)
       </p>

  <li> <p> Fixed: In the step-14 example program, overly conservative
       quadrature formulas were chosen (with 2N+1 quadrature points per space
       dimension, where N is the degree of polynomials). This is unnecessary,
       and now fixed.
       <br>
       (WB 2003/09/23)
       </p>

  <li> <p> Fixed: On AIX systems, the xlf Fortran77 compiler wasn't recognized 
       due to the fact that when called with -v it generates many pages
       of output, later leading to a shell error. This is now fixed.
       <br>
       (WB 2003/09/19)
       </p>

  <li> <p> Fixed: The elastic example program, step-8, had a bug in the way
       the system matrix is assembled. In effect, we were solving the
       wrong equation. This should now be fixed.
       <br>
       (WB 2003/09/11)
       </p>

  <li> <p> Fixed: When building with both sparse direct solver MA27 and the
       TECPLOT libraries, the <code>detached_ma27</code> would not
       link properly. This is now fixed.
       <br>
       (WB 2003/09/05)
       </p>

  <li> <p> Improved: The script that builds the dependency lists for Makefiles
       has been rewritten in C++, since the previous perl script
       became just too slow after the addition of parts of
       BOOST. Using the old perl script should still work, although it
       simply forwards to the new program. In order to use the new
       one, simply replace the line
       <code><pre>
         \$(PERL) \$D/common/scripts/make_dependencies.pl ...
       </pre></code>
       by
       <code><pre>
         \$D/common/scripts/make_dependencies ...
       </pre></code>
       i.e. call the program directly without the perl interpreter and
       without the file extension for a perl program.
       <br>
       (WB 2003/08/19)
       </p>

  <li> <p> New: First steps to a migration of the documentation from
       <tt>kdoc</tt> to <a href="http://www.doxygen.org">Doxygen</a> have
       been done. It can be generated after installing <a
       href="http://www.doxygen.org">Doxygen</a> by calling <tt>make</tt>
       in <tt>doc/doxygen</tt> and using the preliminary link page <a
       href="../doxygen/index.html">index.html</a> in that directory.
       <br>
       (GK 2003/08/02)
       </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p> 
       New: There is now a new <code>PolynomialsP</code>
       class which is based on <code>Polynomials::Monomial</code> and <code>PolynomialSpace</code>. In contrast to the
       default ordering of the polynomials in <code>PolynomialSpace</code>, (e.g. for degree=2) <i>1,
       x, x<sup>2</sup>, y, xy, y<sup>2</sup></i>, the <code>PolynomialsP</code> class now gives the
       (natural?!)  ordering <i>1, x, y, xy, x<sup>2</sup>,
       y<sup>2</sup></i>.
       <br>
       (RH 2004/03/11)
       </p>

  <li> <p> 
       New: The classes <code>PolynomialSpace</code> and
       <code>TensorProductPolynomials</code> now have
       new <code>set_numbering</code> functions which
       allow to change the ordering of the polynomials. The ordering
       and the indices of the polynomials kann be checked by using the
       new <code>output_indices</code> functions.
       <br>
       (RH 2004/03/11)
       </p>

  <li> <p> New: The class <code>PolynomialsBDM</code> implements BDM polynomials in
       two dimensions on the unit square. It was implemented as is
       according to some urgent need, but should be suitable to be fit
       into a <code>FiniteElement</code> similar to
       Raviart/Thomas.
       <br>
       (GK 2004/01/05)
       </p>

  <li> <p> New: Objects of type <code>Polynomial</code>
       can now be added to and subtracted from each other through
       operators <code>+=</code> and <code>-=</code>.
       <br>
       (GK 2003/12/16)
       </p>

  <li> <p> New: There is now a class <code>QuadratureSelector</code> that allows to select a
       quadrature formula based on a string argument and possibly a
       number indicating the order of the formula.
       <br>
       (Ralf B. Schulz 2003/10/29)
       </p>

  <li> <p> Fixed: The constructor of the <code>QGauss</code> class
       computed positions and weights of quadrature points in long double accuracy.
       However, on machines where long double is the same as double, it 
       never reached the requested accuracy, in effect leading to an infinite loop.
       This is now fixed.
       <br>
       (WB 2003/09/19)
       </p>

  <li> <p> New: The <code>Function</code> class now
       exports the value of its template argument through the static
       member variable <code>dimension</code>.
       <br>
       (WB 2003/09/15)
       </p>

  <li> <p> Changed: The <code>ParameterHandler::declare_entry</code> function
       now allows to redeclare an entry that has already been
       declared. This can be used to override a default value
       previously set.
       <br>
       (WB 2003/09/03)
       </p>

  <li> <p> Improved: The <code>ParameterHandler::declare_entry</code> function now takes an
       additional parameter (defaulting to the empty string) that can be used
       to document the intent of a parameter. This string, together with the
       default value of an entry, is written to the output by the <code>ParameterHandler::print_parameters</code> function that
       can be used to generate a virginial parameter file, or one that contains
       the settings of all parameters used in a computation.
       <br>
       (WB 2003/08/13)
       </p>

  <li> <p> Changed: The <code>ParameterHandler::declare_entry</code> previously
       returned a value indicating whether the just-declared entry didn't
       already existed and that the default value matches the given
       pattern. However, this value could only always be true since these two
       conditions were already guarded by assertions in the implementation at
       least in debug mode, so the return value was meaningless. The function
       has now no return type any more.
       <br>
       (WB 2003/08/13)
       </p>

  <li> <p> Improved: <code>Logstream</code>::<code>depth_console</code>, <code>Logstream</code>::<code>depth_file</code>, <code>Logstream</code>::<code>log_execution_time</code> and <code>Logstream</code>::<code>log_time_differences</code> return the previous value.
       <br>
       (GK 2003/06/22)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p> Improved: The matrix-vector operations of <code>SparseMatrix</code> accept arguments of type <code>BlockVector</code>.
  <br>
  (GK/2004/03/31)
  </p>

  <li> <p> Fixed: The <code>SparseMatrix</code> iterator classes
       had various interesting bugs when rows of the matrix were completely
       empty. These should now be fixed.
       <br>
       (WB 2004/03/30)
       </p>

  <li> <p> New: The <code>SparsityPattern</code> class now also
       has an iterator class that allows to walk over all nonzero entries of a
       sparsity pattern.
       <br>
       (WB 2004/03/29)
       </p>

  <li> <p> New: The iterator classes into <code>SparseMatrix</code> have been rearranged and extended, so
       that it is now also possible to write to individual matrix entries
       through these iterators.
       <br>
       (WB 2004/03/29)
       </p>

  <li> <p> New: The <code>Vector</code> and <code>BlockVector</code> classes now have member functions
       <code>is_non_negative</code> that check whether a vector
       has no negative entries.
       <br>
       (WB 2004/02/29)
       </p>

  <li> <p> Fixed: The <code>SolverMinRes</code> class had a nasty bug where we were
       inadvertently copying vectors; this could also have led to a memory
       corruption bug. This is now fixed.
       <br>
       (WB 2004/02/26)
       </p>

  <li> <p> New: There is now a function <code>FullMatrix::add_scaled</code>. It replaces the old
       function <code>FullMatrix::add</code> which did the same,
       but had a name that was incompatible with respective functions in the
       other matrix classes.
       <br>
       (WB 2004/02/23)
       </p>

  <li> <p> New: <code>FullMatrix</code> has new functions <code>add</code> and ,<code>Tadd</code>
       allowing to add to a selected block of the matrix.
       <br>
       (GK 2004/02/12)
       </p>

  <li> <p> New: The <code>Vector</code> class now has operators to compare for
       equality and inequality.
       <br>
       (WB 2004/02/09)
       </p>

  <li> <p> New: The <code>SparseMatrix::operator()</code> generated an assertion
       failure if the requested entry in the matrix isn't there. This has been
       changed so that it actually throws an exception instead, also in
       optimized mode.
       <br>
       (WB 2004/02/06)
       </p>

  <li> <p> New: There is now a function <code>SparseMatrix::frobenius_norm</code> that computes the
       Frobenius norm of a sparse matrix.
       <br>
       (WB 2004/02/06)
       </p>

  <li> <p> Changed: In matrix-vector operations of the <code>Full/SparseMatrix</code> classes, source and destination
       cannot be the same. We now also check that this is indeed the case.
       <br>
       (WB 2004/01/26)
       </p>
  
  <li> <p> Improved: Initialization routines of class <code>SparseMatrix</code> have an additional parameter
       controlling the storage of diagonal entries.
       <br>
       (GK 2003/11/18)
       </p>

  <li> <p> New: 
       <code>SolverFGMRES</code> implements the flexible
       GMRES method with varying preconditioner from the right. It is
       also accessible in <code>SolverSelector</code> by choosing <tt>fgmres</tt>.
       <br>
       (GK 2003/10/07)
       </p>

  <li> <p> Changed: The <code>SparseDirectMA27</code>
       class used to store a pointer to the sparsity pattern of the
       matrix. It now releases this as soon as it doesn't need it any
       more.
       <br>
       (WB 2003/09/09)
       </p>

  <li> <p> New: Some of the member matrix-vector functions of the
       <code>BlockSparseMatrix</code> class that
       previously could only be used with arguments of type <code>BlockVector</code> can now also be used with the
       usual <code>Vector</code> class provided the
       block matrix has only one block row or column.
       <br>
       (WB 2003/09/09)
       </p>

  <li> <p> Fixed: <code>FullMatrix</code>::<code>copy_from</code> didn't compile when copying
       from a sparse matrix. This is now fixed.
       <br>
       (Ralf B. Schulz 2003/09/04)
       </p>

  <li> <p> New: The classes <code>FullMatrix</code> and
       <code>PreconditionBlockJacobi</code> have a <code>const_iterator</code>.
       <br>
       (GK 2003/07/18)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p> Improved: The <code>DoFTools::compute_Cuthill_McKee</code> function
       needs to build a sparsity pattern for its operations, and uses
       the <code>DoFHandler::max_couplings_per_dof</code>
       function for this. However, the estimates returned by the
       latter function are rather bad in 3d, leading to excessive
       memory allocation in the Cuthill McKee algorithm in 3d. This is
       fixed by using an intermediate compressed sparsity pattern
       instead if we are in 3d.
       <br>
       (WB 2004/05/18)
  </p>

  <li> <p> Improved: <code>Triangulation</code> has
       functions <code>n_faces</code> and <code>n_active_faces</code>, globally as well as by level,
       similar to <code>n_cells</code>.
       <br>
       (GK 2004/05/18)
  </p>

  <li> <p>
       New: Added support for <a href="http://www.geuz.org/gmsh/"
       target="_top">gmsh</a> mesh format in <code>GridIn::read_msh</code>.
       <br>
       (Luca Heltai 2004/04/21)
       </p>
              
  <li> <p>
       New: The function <code>GridGenerator::cylinder_shell</code> generates a domain
       of the type which the name suggests.
       <br>
       (WB 2004/04/19)
       </p>
              
  <li> <p>
       Changed: The <code>KellyErrorEstimator::estimate</code> function takes an
       additional parameter that lets it only compute error indicators for a
       certain subdomain. This is meant to allow for a better parallelization
       of efforts in parallel programs.
       <br>
       (GK 2004/04/01)
       </p>

  <li> <p>
       Changed: <code>MGTransferSelect</code> uses target components
       correctly. Unfortunately, the untested class <code>MGTransferBlock</code> does not work anymore. Since its
       usefulness was not clear anyway, this state may continue for a while.
       <br>
       (GK 2004/04/01)
       </p>

  <li> <p>
       New: There is now a new <code>FE_Poly</code>
       class which is templatized for polynomial spaces like <code>TensorProductPolynomials</code>, <code>PolynomialSpace</code> or <code>PolynomialsP</code>. Many finite element classes
       are now derived from this class and the implementation of all
       common functionality is now moved from these finite element
       classes to <code>FE_Poly</code>.
       <br>
       (RH 2004/03/18)
       </p>

  <li> <p> New: The new function <code>MatrixTools::local_apply_boundary_values</code> folds
       boundary value information into local matrices and vectors before they
       are written into the global matrix and vector. This way, the final call
       to  <code>MatrixTools::apply_boundary_values</code> can
       be avoided.
       <br>
       (WB 2004/03/16)
       </p>

  <li> <p> New: There are now functions <code>ConstraintMatrix::distribute_local_to_global</code> that
       take a local matrix or vector and distribute it to a global one, but
       taking care of constrained degrees of freedom; in that case, the
       respective entries are already distributed to the final place in the
       global matrix or vector. That way, the final call to the <code>ConstraintMatrix::condense</code> function can be
       avoided.
       <br>
       (WB 2004/03/15)
       </p>

  <li> <p> New: The new functions <code>SparsityPattern::partition</code>, <code>GridTools::partition_triangulation</code>, <code>DoFTools::get_subdomain_association</code>, <code>DoFTools::count_dofs_with_subdomain_association</code>,
       <code>GridTools::get_subdomain_association</code>, <code>GridTools::count_cells_with_subdomain_association</code>, 
       and <code>DoFRenumbering::subdomain_wise</code> can now
       be used to generate partitions of a triangulation and its associated
       degrees of freedom suitable for parallel computations with PETSc.
       <br>
       (WB 2004/03/08)
       </p>

  <li> <p> Improved: When eliminating nodes from a matrix using the <code>ConstraintMatrix::condense</code> functions, the
       diagonal entry was set to one. It is now set to an entry that more
       resembles the size of the other diagonal entries, so that we don't run
       into scaling problems for applications that have very large or small
       domains.
       <br>
       (WB 2004/03/02)
       </p>

  <li> <p> Changed: The classes <code>DataOut*</code> and <code>KellyErrorEstimator</code> have been generalized to take
       more and different vector types as input parameters. In particular,
       they now take normal and block vectors over doubles and floats, as well
       as PETSc vectors if PETSc support is detected during configuration of
       the library.
       <br>
       (WB 2004/03/01)
       </p>

  <li> <p> Changed: The template parameter of the functions in the <code>GridRefinement</code> class have been changed. Where they
       previously denoted the type over which the <code>Vector</code> class is to be templated, they now mean the
       actual vector class. Thus, they can be any other template class as long
       as there are suitable operations defined over them. In addition,
       the documentation stated that they must be vectors of floats; this
       hasn't been true any more for quite a while already, and is duly
       removed from the documentation.
       <br>
       (WB 2004/02/28)
       </p>

  <li> <p>
       New: The function
       <code>FETools::project_dg</code>
       performs <i>L<sup>2</sup></i>-projections between finite element spaces
       of different degrees on the same mesh.
       <br>
       (GK 2003/11/28)
       </p>

  <li> <p>
       Improved: <code>FiniteElementData</code> has a function
       <code>tensor_degree()</code>, returning the degree of the
       polynomial space suitable for choosing a tensor product quadrature
       formula. 
       <br>
       (GK 2003/11/28)
       </p>

  <li> <p>
       New: Long requested but never implemented before in the
       library: there is now a function <code>GridTool::find_active_cell_around_point</code>
       that, given a point, finds the active cell in which this point
       lies.
       <br>
       (WB 2003/10/30)
       </p>

  <li> <p>
       New: <code>MGCoarseGridHouseholder</code>
       implements a coarse grid solver using QR-factorization.
       <br>
       (GK 2003/10/17)
       </p>

  <li> <p>
       Fixed: The <code>FEFaceValuesBase::boundary_form</code>
       function was declared but not implemented. This is now fixed.
       <br>
       (J&ouml;rg R. Weimar 2003/10/22)
       </p>

  <li> <p>
       Improved: The <code>MatrixCreator::create_mass_matrix</code>
       functions are now templatized also on the template argument of
       the <code>SparseMatrix</code> class. This allows
       invoking this function for <code>SparseMatrix&lt;double&gt;</code> and <code>SparseMatrix&lt;float&gt;</code> objects.
       <br>
       (RH 2003/10/22)
       </p>

  <li> <p>
       New: There is now also a function <code>MGDoFCellAccessor::neighbor_child_on_subface</code>
       that returns the result of the <code>CellAccessor::neighbor_child_on_subface</code>
       function but converts it so that one can also access MGDoF
       data.
       <br>
       (RH 2003/10/22)
       </p>

  <li> <p>
       New: There are now functions <code>CellAccessor::neighbor_child_on_subface</code> and <code>DoFCellAccessor::neighbor_child_on_subface</code>
       that should be called instead of using <code>GeometryInfo::child_cell_on_face</code> in most cases.
       <br> 
       (WB 2003/10/17)
       </p>

  <li> <p>
       New: <code>GridGenerator</code> has a new
       function <code>subdivided_hyper_rectangle</code> 
       which generates a rectangle with given corner points and possibly 
       different numbers of subdivisions in different directions.
       Use it, e.g., to generate a domain of 1*4 length units
       with square cells.
       <br> 
       (Joerg Weimar 2003/09/09)
       </p>

  <li> <p>
       Improved: The 3d grid reordering code in the <code>GridReordering</code> class now uses an algorithm
       that is linear in the number of elements. The old code was
       exponential, so this is a vast improvement.
       <br> 
       (Michael Anderson 2003/09/23)
       </p>

  <li> <p>
       Improved: <code>GridOut</code> has a an improved
       functionality for <code>write_eps</code> 
       to color the grid according to the refinement level.
       A corresponding option is included in 
       <code>GridOutFlags::Eps<2></code>.
       <br> 
       (Joerg Weimar 2003/09/09)
       </p>

  <li> <p> New: The <code>TriaAccessor</code>::<code>point_inside</code> function is now also
       implemented in 3d.
       <br>
       (Joerg Weimar, WB 2003/09/04)
       </p>

  <li> <p> New: The <code>TriaAccessor</code>::<code>recursively_set_material_id</code> function sets
       the material id of the present cell and of all its children,
       grandchildren, etc to the given value.
       <br>
       (WB 2003/09/04)
       </p>

  <li> <p> New: The new <code>FETools</code>::<code>get_fe_from_name</code> function can do the
       reverse of the <code>FiniteElement</code>::<code>get_name</code> function: it takes a string and
       parses it to regenerate a finite element from it. Useful for
       parsing finite element names from input files.
       <br>
       (WB 2003/07/08)
       </p>

  <li> <p> New: The <code>DataOut_DoFData</code>::<code>merge_patches</code> now takes a second
       parameter that indicates a shift for each vertex of the given
       patches to be merged. This is sometimes nice if one wants to
       generate "exploded" views of a collection of subdomains. It is
       also templatized on the first argument, so can merge in some
       other <code>DataOut_DoFData</code> that create
       the same type of patches but are otherwise different.
       <br>
       (WB 2003/06/30)
       </p>

  <li> <p> Fixed: The <code>FETools</code>::<code>extrapolate</code> function operates on patches
       of cells, but didn't check whether the grid is at least refined
       once everywhere. If this was not the case, it would generate
       wrong results. It now checks for this, and if the grid has
       unrefined coarse grid cells, an exception is generated.
       <br>
       (WB 2003/06/25)
       </p>

  <li> <p>
       Improved: <code>FEFaceValuesBase</code> has a new
       function <code>orientation</code> accessing a unique
       and consistent orientation for each face.
       <br> 
       (GK 2003/06/23)
       </p>

  <li> <p> 
       Changed: Embedding and restriction matrices for intergrid transfer are
       now computed in the constructor of most elements, rather than taken from
       precomputed and tabulated values. This removes restrictions on which
       elements are available since the old tables were only precomputed for
       certain polynomial degrees and are now available for all.
       <br>
       (WB 2003/06/09)
       </p>

  <li> <p> 
       New: Finite elements got a function <code>get_interpolation_matrix</code>, with which they can
       compute interpolations between different finite elements. Most will use
       this to compute interpolations between finite elements of the same kind
       but different polynomial degrees. The <code>FETools::get_interpolation_matrix</code> makes use of
       this function if it implements the desired interpolation, and falls back
       to the old algorithm if this is not the case.
       <br>
       (WB 2003/06/09)
       </p>

  <li> <p> 
       New: Finite elements got a function <code>get_name</code>, which can be used to identify a finite
       element by its name.
       <br>
       (WB 2003/06/09)
       </p>

  <li> <p> 
       New: Raviart-Thomas elements are now implemented in the <code>FE_RaviartThomas</code> class.
       <br>
       (WB 2003/06/09)
       </p>
</ol>


*/
