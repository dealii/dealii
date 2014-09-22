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
 * @page changes_between_3_1_0_and_3_2_0 Changes between Version 3.1.0 and 3.2.0

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li> <p>
       New: example program step-11.
       <br>
       (WB 2001/09/17)
       </p>

  <li> <p>
       New: example program step-10.
       <br>
       (WB, RH 2001/09/12)
       </p>

  <li> <p>
       New: the <code>./configure</code> script now also recognizes
       Compaq's cxx compiler. The library now also compiles cleanly
       with this compiler.
       <br>
       (WB 2001/07/06)
       </p>

  <li> <p>
       Changed: The libraries are no more linked using the C++
       compilation flags, but rather using LDFLAGS. Some compilers
       object to compilation flags on the linker line.
       <br>
       (WB 2001/05/23)
       </p>

  <li> <p>
       Changed: If in multithreaded mode, the ACE library is now
       automatically added to the <code>\$(LIBS)</code> Makefile
       variable. There is no need anymore for a special clause in your
       Makefile.
       <br>
       (WB 2001/05/23)
       </p>

  <li> <p>
       New: the <code>./configure</code> script now also recognizes
       the Intel ICC compiler that was released for Linux lately. The
       library now also compiles cleanly with this compiler.
       <br>
       (WB 2001/05/22)
       </p>

  <li> <p>
       Fixed: the program that generated HTML from the example
       programs was broken on some platforms for reasons beyond our
       knowledge. This is now fixed.
       <br>
       (Roger Young, WB 2001/03/22)
       </p>

  <li> <p>
       Improved: libraries are now all in a subdirectory lib in the
       deal.II home directory. It is now easy to
       make this directory a link to put the libraries on a faster
       disk.
       <br>
       (GK 2001/03/14)
       </p>

  <li> <a name="new_fe_mapping_design"></a>
       <i> New Design of <code>FiniteElements</code>
       and <code>Mappings</code></i>
              
       <p>
       New: The finite element
       classes have been redesigned from scratch to allow also
       elements that depend on the actual cell shape, like
       Raviart-Thomas, BDM or edge elements. We have implemented
       continuous Lagrange elements up to degree four and discontinous
       Lagrange elements up to degree 10. They can be easily extended
       to higher degree.
       </p>

       <p> 
       Furthermore we have totally reimplemented the mapping
       between unit and real cell.  In future we won't be restricted
       to Q<sub>1</sub> mapping any more, but we will have
       Q<sub>p</sub> mappings of arbitrary degree p. Presently
       implemented are mappings up to degree 10 in 2d, and degree
       three in 3d. This allows to approximate curved boundaries not
       only by straight lines (linear approximation) but also by
       quadratic, cubic, etc approximation. The finite elements will
       be able to be combined with arbitrary mappings. So in future we
       can use subparametric, iso- as well as superparametric
       elements.
       </p>

       <p>
       The new implementation uses a totally new structure of
       communication between <code>FEValues</code>, the
       <code>FiniteElements</code> and the <code
       class="class">Mappings</code>. Despite of this, the new
       structure will almost not be 'visible' to the user of
       deal.II as we tried to retain the interface
       to the user (mainly that of <code
       class="class">FEValues</code>) as far as possible.
       </p>

       <p>
       Together with this new design comes a reduction of 25000(!)
       lines of deal.II code. This elimination
       mainly affected code that previously was machine generated
       (e.g. by maple). In view of almost unchanged execution times of
       applications, the faster compilation and additional
       functionality constitutes a significant improvement of
       deal.II. Results produced by the new code
       match those of the old design up to machine precision.
       <br>
       (RH & GK 2001/03/13)
       </p>

  <li> <p>
       New: There is now some support to include and use routines from the 
       <a href="http://www.cse.clrc.ac.uk/Activity/HSL" 
       target="_top">Harwell Subroutine Library</a>.
       <br>
       (WB 2001/01/30)
       </p>

  <li> <p>
       New: The <code>./configure</code> script now checks for the
       existence of a Fortran 77 compiler, and sets its path, some
       compiler flags and the libraries to be linked in when mixing
       C++ and F77 in some variables in the file
       <code>common/Make.global_options</code>.
       <br>
       (WB 2000/12/30)
       </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p> 
       New: Color function <code
       class="member">DataOutBase::EpsFlags::reverse_grey_scale_color_function</code>.
       <br>
       (WB 2001/08/24)
       </p>

  <li> <p> 
       New: Function <code>QProjector::project_to_child</code>
       generates quadrature formulae which act on the area which a
       child would occupy.
       <br>
       (WB 2001/08/23)
       </p>

  <li> <p>
       Changed: The examples classes in the base directory are now
       moved into a namespace <code>Functions</code> of
       their own. This improves encapsulation, but also keeps the
       documentation of these functions together, as they were
       previously scrambled all over the screen in the documentation
       page of the base library.
       <br>
       (WB 2001/07/18)
       </p>

  <li> <p> 
       New: classes <code>FourierSineFunction</code>
       and <code>FourierCosineFunction</code>,
       resembling one mode of a Fourier decomposition. Classes <code
       class="class">FourierSineSum</code> and <code
       class="class">FourierCosineSum</code>, resembling sums of such
       modes of a Fourier decomposition.
       <br>
       (WB 2001/07/18)
       </p>

  <li> <p>
       New: class <code>vector2d</code> was introduced
       in analogy to STL class <code>vector</code>. The
       bew class provides a two-dimensional array and replaces the use
       of <code>FullMatrix</code> in the base library.
       <br>
       (GK 2001/05/21)
       </p>

  <li> <p>
       Improved: <code>JobIdentifier</code>::<code
       class="member">operator()</code> includes host name if
       supported by the operating system
       <br>
       (GK 2001/05/17)
       </p>

  <li> <p>
       New: There is now a new <code
       class="class">AutoDerivativeFunction</code> class that
       automatically computes the gradient of a function by employing
       numerical difference quotients. This only, if the user function
       does not provide the gradient function himself. This class can
       replace the <code>Function</code> class as base
       class for user defined <code>Function</code>
       classes.
       <br>
       This new class can be used if the user only implements the
       <code>value</code> function but wants to call
       also the <code>gradient</code> functions.
       <br>
       (RH 2001/05/15)
       </p>

  <li> <p>
       New: The <code>Quadrature</code> class now has a
       constructor that only accepts quadrature points and sets the
       weights to invalid values.
       <br>
       (GK 2001/04/26)
       </p>

  <li> <p>
       New: The function <code>Logstream</code>::<code
       class="member">get_prefix</code> allows access to the prefix
       string used for log-files.
       <br>
       (GK 2001/04/26)
       </p>

  <li> <p>
       New: There is now a global function <code
       class="function">trace</code> that computes the trace of a tensor
       of rank 2.
       <br>
       (WB 2001/04/12)
       </p>

  <li> <p>
       New: The <code>Threads</code> now has a barrier
       class in order to synchronise multiple threads. In
       multithreaded mode, this is an alias to the <code
       class="class">ACE_Barrier</code> class of the ACE library, if
       deal.II is not configured for multithreading, then the class is
       a dummy class that does nothing.
       <br>
       (WB 2001/03/19)
       </p>

  <li> <p>
       New: We now support the output format of the <a
       href="http://www.kitware.com/vtk.html"
       target="_top">Visualization Toolkit (Vtk)</a> from the <code
       class="class">DataOutBase</code> class and all derived classes.
       <br>
       (WB 2001/03/19)
       </p>

  <li> <p>
       New: The class <code
       class="class">TensorProductPolynomials&lt;dim&gt;</code>
       performs the tensor product of 1d polynomials, computes its
       values, its first and second derivatives. If <code>n</code>
       polynomials are given to the constructor of this class, it
       constructs <code>n<sup>dim</sup></code> tensor product
       polynomials.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: <code>LagrangeEquidistant</code> is a new
       class derived from <code>Polynomial</code>
       providing the 1d Lagrange interpolation of degree
       <code>n</code> of <code>n+1</code> equidistant support points
       in <code>[0,1]</code>. It is implemented up to degree 10.
       This class is used for the implementation of the continuous and
       discontinuous Lagrange elements.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: The new <code>Polynomial</code> class can
       be used as base class for all 1d polynomials. It stores the
       coefficients of the polynomial and uses the Horner scheme to
       evaluate values and all derivates.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: function <code>contract</code> for two arguments of
       <code>Tensor&lt;1,dim&gt;</code>
       <br>
       (GK 2001/03/05)
       </p>

  <li> <p> New: <code>Logstream</code>::<code
       class="member">log_time_differences (bool)</code> makes <code
       class="class">Logstream</code> print the time since the last
       log line instead of accumulated time since program start.
       <br>
       (GK 2001/03/05)
       </p>

  <li> <p> Fix: <code>Logstream</code>::<code
       class="member">pop()</code> does not perform anything on empty
       stacks.
       <br>
       (GK 2001/03/05)
       </p>

  <li> <p>
       Changed: Sort the quadrature points of each <code
       class="class">Quadrature&lt;1&gt;</code> in ascending order. This
       actually changed the order of the points of only <code
       class="class">QGauss4</code> and <code
       class="class">QGauss5</code>.
       <br>
       (Ralf Hartmann 2001/01/22)
       </p>

  <li> <p>
       New: function <code>contract</code> for two arguments of
       <code>Tensor&lt;1,dim&gt;</code>
       <br>
       (GK 2001/01/15)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p>
       New: Function <code>FullMatrix::symmetrize()</code>.
       <br>
       (WB 2001/09/20)
       </p>

  <li> <p>
       Improved: the stopping criterion of <code
       class="class">SolverBicgstab</code> without computing the exact
       residual is implemented
       <br>
       (GK 2001/09/11)
       </p>

  <li> <p>
       New: The <code>FullMatrix</code> class now has a
       function <code>operator*=</code>, which simply
       scales the matrix.
       <br>
       (WB 2001/09/02)
       </p>

  <li> <p>
       New: Function <code>BlockSparseMatrix::el()</code>,
       just like <code>SparseMatrix::el()</code>.
       <br>
       (WB 2001/08/21)
       </p>

  <li> <p>
       New: There is now a class <code>MatrixOut</code>
       which can be used to convert a matrix into a graphical output,
       for example in order to see visually that a matrix is
       diagonally dominant or not.
       <br>
       (WB 2001/08/21)
       </p>

  <li> <p>
       Changed: Base class <code>Solver</code> and all
       <code>Preconditioner</code> classes are now
       derived from <code>Subscriptor</code>.
       Class <code>PreconditionLACSolver</code> now uses
       smartpointers to the given solver and preconditioner
       objects. You will, therefore, have to derive your own
       preconditioners and solvers from <code
       class="class">Subscriptor</code> if you want to use it with
       <code>PreconditionLACSolver</code>.
       <br>
       (WB 2001/08/16)
       </p>

  <li> <p>
       New: Classes <code>Vector</code> and
       <code>BlockVector</code> now have member functions
       <code>operator*=</code> which scale the vectors
       by a constant factor. These functions are aliases for the 
       <code>scale</code> member functions except that
       they return a reference to themselves.
       <br>
       (WB 2001/08/14)
       </p>

  <li> <p>
       New: There is now a function 
       <code>FullMatrix::precondition_Jacobi</code>. The 
       <code>PreconditionJacobi</code> class is
       therefore now also applicable with the full matrix class.
       <br>
       (WB 2001/08/11)
       </p>

  <li> <p>
       New: The <code>Vector</code> and
       <code>BlockVector</code> classes can now be
       initialized using a new constructor that takes two iterators
       that denote a range of elements which are to be copied.
       <br>
       (WB 2001/08/08)
       </p>

  <li> <p>
       Changed: The <code>SolverCG</code> class now
       saves the initial matrix vector product if the initial value of
       the solution is vector is zero, as it is common practice to use
       this is starting vector.
       <br>
       (WB 2001/06/25)
       </p>

  <li> <p>
       New: <code>SparsityPattern::reinit</code> no
       more throws an error if the given maximal row lengths are all
       zero.
       <br>
       (WB 2001/06/25)
       </p>

  <li> <p>
       New: Class <code>CompressedBlockSparsityPattern</code>
       may be used as an intermediate form of the <code
       class="class">BlockSparsityPattern</code>.
       <br>
       (WB 2001/06/23)
       </p>

  <li> <p>
       New: Class <code>CompressedSparsityPattern</code>
       may be used as an intermediate form of the <code
       class="class">SparsityPattern</code> class if memory
       requirements are tight during construction of the sparsity
       pattern.
       <br>
       (WB 2001/06/22)
       </p>

  <li> <p>
       New: There are now functions 
       <code>SparsityPattern::copy_from</code> and
       <code>SparseMatrix::copy_from</code>
       that can be used to construct sparsity patterns and matrix at
       once, i.e. without setting up a sparsity pattern element by
       element, possibly after using a way too large value for the
       maximal row length, then calling 
       <code>SparsityPattern::compress</code> etc.
       <br>
       (WB 2001/05/07)
       </p>

  <li> <p>
       New: <code>BlockIndices::block_size</code>
       returns the size of a specified block.
       <br>
       (WB 2001/05/07)
       </p>

  <li> <p>
       New: There is now a (private) function <code
       class="member">SparsityPattern::optimized_lower_bound</code>
       that is used as an optimized replacement for <code
       class="member">std::lower_bound</code> for searching in the
       column number arrays. It unrolls small loops and it also seems
       that the compiler is able to optimized it better due to
       eliminated template parameters, making it about twice as fast
       as the standard implementation. In effect, it also speeds up
       the SSOR preconditioner that spent about one third of its time
       in that function by approximately 15 per cent.
       <br>
       (WB 2001/04/25)
       </p>

  <li> <p>
       New: The <code>FilteredMatrix</code> class is a
       replacement for the <code
       class="class">MatrixTools::apply_boundary_values</code>
       function for cases where you would like to solve several times
       with the same matrix, either for different right hand sides, or
       for different boundary values.
       <br>
       (WB 2001/04/27)
       </p>

  <li> <p>
       New: There is now a function <code
       class="class">Vector</code>::<code
       class="member">scale(Vector)</code>
       that scales each element of the vector by the corresponding
       element of the argument.
       <br>
       (WB 2001/04/23)
       </p>

  <li> <p> Changed: Solver functions <code>solve</code>
       return void now. If the solver has not converged within the
       maximum number of iterations or encounters a breakdown, it
       throws an exception of type <code
       class="class">SolverControl</code>::<code
       class="member">NoConvergence</code> instead of
       returning a special value.
       <br>
       (GK 2001/03/29)
       </p>

  <li> <p> 
       New: The functions <code>FullMatrix</code>::<code
       class="member">mmult</code> and <code
       class="member">Tmmult</code> now have an additional
       <code>adding</code> argument. If this flag is
       <code>true</code>, the matrix-matrix product is added to the
       resulting matrix; if <code>false</code>, the resulting matrix
       is set to (overwritten by) the matrix-matrix product. The
       latter is the behaviour these functions had before. Hence the
       default value is set to be <code>false</code> to ensure
       backward compatibility.
       <br>
       (RH 2001/03/29)
       </p>

  <li> <p>
       New: class <code>SchurMatrix</code> implements
       a Schur complement for block matrices. It provides matrix-vector
       multiplication suitable for iterative methods as well as pre- and
       post-processing of right hand side and slution, respectively. 
       <br>
       (GK 2001/03/22)
       </p>

  <li> <p>
       Removed: The explicite instantiations of <code
       class="class">SparseMatrix&lt;long double&gt;</code> are
       removed as a prerelease of gcc3.0 fails to compile it. A user
       of <code>SparseMatrix&lt;long double&gt;</code>
       needs now to include
       <code>lac/include/lac/sparse_matrix.templates.h</code> into his
       source file and to use an appropriate compiler, e.g. gcc2.95.2 or
       a future version of gcc3.0 (that will then hopefully be fixed).
       <br>
       (RH 2001/03/14)
       </p>
  
  <li> <p> 
       New: class <code
       class="class">BlockMatrixArray&lt;MATRIX&gt;</code> implements
       a block matrix based on an array of matrix pointers. Since this
       array may be incomplete, no direct access to entries is
       allowed. Matrix blocks may be scaled and transposed.
       <br>
       (GK 2001/02/13)
       </p>

  <li> <p>
       New: There is now some support to include and use routines from the 
       <a href="http://www.cse.clrc.ac.uk/Activity/HSL" 
       target="_top">Harwell Subroutine Library</a>, and support
       classes 
       <code>SparseDirectMA27</code> and
       <code>SparseDirectMA47</code>
       for the sparse direct solvers MA27 and MA47.
       <br>
       (WB 2001/01/30)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p> 
       New: Class <code>MappingQ1Eulerian</code>
       implementing an Eulerian mapping.
       <br>
       (<a href="mailto:ms@biomech.tu-graz.ac.at">Michael Stadler</a> 2001/09/24)
       </p>

  <li> <p> 
       New: <code>VectorTools::create_boundary_right_hand_side</code>
       integrates boundary forces for inhomogeneous Neumann boundary values.
       <br>
       (WB 2001/09/13)
       </p>

  <li> <p> 
       New: <code>DoFTools::make_flux_sparsity_pattern</code>
       now exists also for 1d.
       <br>
       (WB 2001/09/03)
       </p>

  <li> <p> 
       New: There are now two functions
       <code
       class="member">FETools::hierarchic_to_lexicographic_numbering</code>
       and <code
       class="member">FETools::lexicographic_to_hierarchic_numbering</code>
       which map the hierarchical numbering used in continuous finite
       element classes to a lexicographical numbering and back.
       <br>
       (WB 2001/08/31)
       </p>

  <li> <p> 
       New: <code>ConstraintMatrix::close</code>
       now simply returns instead of throwing an exception, if the
       matrix was already closed.
       <br>
       (WB 2001/08/30)
       </p>

  <li> <p>
       New: Member function
       <code>ConstraintMatrix::is_identity_constrained</code>.
       <br>
       (WB 2001/08/29)
       </p>

  <li> <p>
       Fixed: in a rather rare case, some work was done twice in the
       <code>KellyErrorEstimator</code> class when in
       multithread mode. This is now fixed.
       <br>
       (WB 2001/08/24)
       </p>

  <li> <p>
       New: There is now a class <code>GridTools</code>
       which provides algorithms working on triangulations. At
       present, it offers a function computing the diameter of a
       triangulation.
       <br>
       (WB 2001/08/16)
       </p>

  <li> <p>
       Changed: The <code>MatrixCreator</code> and <code
       class="class">MatrixTools</code> class have lost their template
       arguments. Rather, the individual functions now are templated
       on their own, so the compiler can pick out the correct space
       dimension on its own.
       <br>
       (WB 2001/08/15)
       </p>

  <li> <p>
       Extended: <code>ConstraintMatrix::merge</code>
       can now handle arguments which further constrain the present object.
       <br>
       (WB 2001/07/30)
       </p>

  <li> <p>
       New: Implement
       <code>DoFTools::make_sparsity_pattern</code>,
       <code>DoFTools::make_boundary_sparsity_pattern</code>,
       and 
       <code>ConstraintMatrix::condense</code> to work on
       the <code>CompressedSparsityPattern</code> and
       <code>CompressedBlockSparsityPattern</code> and
       classes.
       <br>
       (WB 2001/06/22)
       </p>

  <li> <p>
       New: <code>FE*Values::get_quadrature</code>
       returns a reference to the quadrature formula used by a
       FEValues object.
       <br>
       (WB 2001/06/21)
       </p>

  <li> <p>
       New:  <code
       class="member">DoFTools::compute_intergrid_transfer_representation</code>
       is a function that generates a matrix representation that can
       be used to transfer data from one grid to one component of the
       fields on another, finer, grid.
       <br>
       (WB 2001/05/24)
       </p>

  <li> <p>
       Changed: the  <code>GeometryInfo</code> class
       has been reverted from a general template that calculates the
       values of its member variables by recursion/induction to a set
       of explicitely specialized classes. This seemed necessary since
       most compilers objected to the old implementation as the
       declaration of the values of the class required knowledge of
       the elements of the class with one lower dimension which was,
       however, only being declared at present except for 1d. Also,
       allocating space for variables was difficult as that would mean
       declaring specializations after their first use. The new (and
       very old) implementation is entirely compatible to the previous
       one.
       <br>
       (WB 2001/05/23)
       </p>

  <li> <p>
       Changed: the classes that denote flags to the <code
       class="class">TimeStepBase_Tria</code> class have been move
       from local classes to a namespace of their own, names
       <code>TimeStepBase_Tria_Flags</code>.
       <br>
       (WB 2001/05/23)
       </p>

  <li> <p>
       Fixed: due to a bug in gcc, the compiler did not check that we
       did not use the set of given boundary indicators to the
       <code>DoFTools::extract_boundary_dofs</code>
       function in 1d. That set was therefore ignored. This is now
       fixed.
       <br>
       (WB 2001/05/18)
       </p>

  <li> <p>
       Changed: the flags which are given to the <code
       class="class">GridOut</code> class to modify the appearance of
       the output have been moved from local classes of the <code
       class="class">GridOut</code> class to a namespace names
       <code>GridOutFlags</code> and have lost the
       trailing <code>Flags</code> part in their name.
       This change was necessary as C++ does not allow explicit
       specialization of member classes; the previous use in the
       library was only accepted by GCC as an extension.
       <br>
       (WB 2001/05/17)
       </p>

  <li> <p>
       New: The functions <code
       class="member">DoFTools::map_dof_to_boundary_indices</code>,
       <code>DoFTools::make_boundary_sparsity_pattern</code>,
       <code>DoFHandler::n_boundary_dofs</code>,
       <code
       class="member">DoFHandler::max_couplings_between_boundary_dofs</code>,
       and <code>VectorTools::project_boundary_values</code>,
       are now also implemented in 1d.
       <br>
       (WB 2001/04/29)
       </p>

  <li> <p>
       New: There are now functions <code
       class="class">DoFTools</code>::<code
       class="member">map_dofs_to_support_points</code> and 
       <code
       class="class">DoFTools</code>::<code
       class="member">map_support_points_to_dofs</code> that
       generate a map between DoF indices and the support points of
       these DoFs, and the other way round.
       <br>
       (WB 2001/04/24)
       </p>

  <li> <p>
       New: Implement the <code
       class="member">VectorTools::interpolate_boundary_values</code>
       a second time, this time taking a <code
       class="class">FunctionMap</code> object as do all the other
       functions of this type.
       <br>
       (WB 2001/04/23)
       </p>

  <li> <p>
       Fixed: The vertex numbers written by the <code
       class="class">GridOut</code>::<code
       class="member">write_ucd_faces</code> function are now also
       1-based. This, to be consistent with the vertex numbers given
       by the <code>GridOut</code>::<code
       class="member">write_ucd</code> function.
       <br>
       (RH 2001/04/20)
       </p>

  <li> <p>
       Extended: the <code
       class="class">DerivativeApproximation</code> class now also
       works in 3d, as the spectral norm of a symmetric 3x3 matrix can
       now be computed.
       <br>
       (<a href="mailto:roger@kea.grace.cri.nz">Roger Young</a> 2001/04/12)
       </p>
       
  <li> <p>
       Fixed: the <code
       class="class">DoFRenumbering</code>::<code
       class="member">Cuthill_McKee</code> function
       did not work correctly when giving the <code>reversed_numbering</code>
       flag (off-by-one indexing). This is now fixed.
       <br>
       (<a href="mailto:or@winfos.com">Oliver Rheinbach</a> 2001/04/12)
       </p>
       
  <li> <p>
       Changed: A <code>typedef FunctionMap</code> was
       declared in at least four places previously (in classes <code
       class="class">DoFHandler</code>, <code
       class="class">MatrixTools</code>, <code
       class="class">VectorTools</code>, and <code
       class="class">KellyErrorEstimator</code>). It is now unified in
       one place and is called <code
       class="class">FunctionMap&lt;dim&gt;</code>::<code
       class="member">type</code> (it is a local
       typedef in a class called <code
       class="class">FunctionMap</code>). This type is defined in the
       file <tt>dofs/function_map.h</tt>.
       <br>
       (WB 2001/04/09)
       </p>
       
  <li> <p>
       Fixed: When using Neuman boundary functions in the 
       <code>KellyErrorEstimator</code> class, it was
       assumed that the function object had <code
       class="class">Function</code>::<code
       class="member">vector_value</code> overloaded, even
       in the scalar case. We now use <code
       class="member">Function</code>::<code
       class="member">value</code> instead.
       <br>
       (WB 2001/04/09)
       </p>
       
  <li> <p>
       New: The various <code
       class="class">MatrixCreator</code>::<code
       class="member">create_*_matrix</code> functions
       are now using multiple threads to do their work, if the library
       is configured to use multithreading.
       <br>
       (WB 2001/04/08)
       </p>
       
  <li> <p>
       New: The following functions are extended to work
       on arbitrary mappings:
       <br>
       The static <code>create_mass_matrix</code>,
       <code>create_boundary_mass_matrix</code> and
       <code>create_laplace_matrix</code> member
       functions of the class <code>MatrixCreator</code>,
       <br>
       the static <code>project</code> and <code
       class="member">project_boundary_values</code> member functions
       of the class <code>VectorTools</code>,
       <br>
       the two versions of the static <code
       class="member">estimate</code> member functions of the class
       <code>KellyErrorEstimator</code>,
       <br>
       and the static <code>approximate_gradient</code>
       and <code>approximate_second_derivative</code>
       member functions of the class <code
       class="class">DerivativeApproximation</code>.
       </p>

       <p>
       All these functions now take a <code
       class="class">Mapping</code> object as additional argument.
       <br>      
       For backward compatibility there exist second versions of
       each of these functions that still don't have a <code
       class="class">Mapping</code> argument. These functions
       implicitely use a mapping of the class <code
       class="class">MappingQ1</code>.
       <br>
       (RH 2001/04/05)
       </p>

  <li> <p>
       Removed: Along with the redesign of the <code
       class="class">FiniteElement</code> and the <code
       class="class">Mapping</code> classes the <code
       class="class">FiniteElement</code>::<code
       class="member">get_local_mass_matrix</code> function is now
       removed. This was necessary as the local mass matrix depends not only on the
       considered <code>FiniteElement</code> but also on
       the mapping between the unit to the real cell. Hence this function
       cannot be a member function of a <code
       class="class">FiniteElement</code> any more.
       </p>

       <p>
       As a consequence also the <code
       class="class">MatrixCreator</code>::<code
       class="member">create_mass_matrix</code> function with two
       arguments is removed, as it relied on the <code
       class="member">get_local_mass_matrix</code> function. If in
       future the user wants to create a mass matrix he needs to use
       one of the remaining <code
       class="class">MatrixCreator</code>::<code
       class="member">create_mass_matrix</code> functions that require
       an appropriate <code>Quadrature</code> as
       argument.
       <br>
       (RH 2001/04/04)
       </p>
       
  <li> <p>
       New/Fixed: Now there exists a new <code
       class="class">Triangulation</code>::<code
       class="member">ExcMultiplySetLineInfoOfLine</code> exception,
       that is thrown if the <code>SubCellData</code>
       that is given to <code
       class="class">Triangulation</code>::<code
       class="member">create_triangulation</code>, multiply includes
       the line info of a specific line. Before the fix the wrong
       <code>ExcInteriorLineCantBeBoundary</code>
       exception was thrown.
       <br>
       (RH 2001/04/03)
       </p>
       
  <li> <p>
       Fixed: Missing <code>ucd</code> is now added to the list of
       supported output formats returned by <code
       class="class">GridOut</code>::<code
       class="member">get_output_format_names</code>.
       <br>
       (RH 2001/04/03)
       </p>

  <li> <p>
       New/fixed: In some obscure corner cases, the detection logic in <code
       class="member">DataOut_DoFData::add_data_vector</code> would
       not have been able to detect whether something is a DoF data
       vector or a vector of cell data, and in some equally rare cases
       this would also have made a difference. This is now fixed by
       adding another argument to the function telling it either to
       automatically detect the vector type (default) or to assume
       that it has a certain type (for these corner cases).
       <br>
       (WB 2001/03/30)
       </p>

  <li> <p>
       Removed: the <code>ProblemBase</code> class,
       which has been deprecated since before the release of 
       deal.II 3.0, has finally been removed. The
       same applies for the classes
       <code>Assembler</code>,
       <code>Equation</code>,
       <code>MassMatrix</code>, and
       <code>LaplaceMatrix</code>.
       <br>
       (WB 2001/03/27)
       </p>

  <li> <p>
       New: There is now a class <code>MappingC1</code>
       that implements a continuously differentiable C<sup>1</sup>
       mapping of the boundary by using a cubic mapping with
       continuous derivatives at cell boundaries. The class presently
       only implements something for 2d and 1d (where it does nothing).
       <br>
       (WB 2001/03/27)
       </p>

  <li> <p>
       New: The static <code>interpolate</code>, <code
       class="member">create_right_hand_side</code>, <code
       class="member">interpolate_boundary_values</code>, <code
       class="member">integrate_difference</code> and <code
       class="member">compute_mean_value</code> member functions of
       the class <code>VectorTools</code> are extended
       to work on arbitrary mappings. All these functions now take a
       <code>Mapping</code> object as additional
       argument. 
       <br>
       For backward compatibility there exist second versions of
       each of these functions that still don't have a <code
       class="class">Mapping</code> argument. These functions
       implicitely use a mapping of the class <code
       class="class">MappingQ1</code>.
       <br>
       (RH 2001/03/27)
       </p>

  <li> <p>
       New: <code>Boundary</code> and derived classes
       now have a function <code
       class="member">get_normals_at_vertices</code> that returns a
       multiple of the normal vector to the boundary at the
       vertices of a given face. This is used in the construction of
       C<sup>1</sup> mappings of the boundary.
       <br>
       (WB 2001/03/23)
       </p>

  <li> <p>
       New: The new function <code
       class="class">PersistentTriangulation&lt;dim&gt;</code>::<code
       class="member">restore(unsigned int)</code> allows to restore
       the triangulation step by step.
       <br>
       New: Now there exists a function <code
       class="class">PersistentTriangulation&lt;dim&gt;</code>::<code
       class="member">clear_flags</code> to allow to re-<code
       class="member">read_flags</code> without the need of clearing
       the triangulation in advance.
       <br>
       New: The new function <code
       class="class">PersistentTriangulation&lt;dim&gt;</code>::<code
       class="member">n_refinement_steps</code> returns the number of
       refinement steps stored in <code
       class="member">refine_flags</code>.
       <br>
       (RH 2001/03/20)
       </p>

  <li> <p>
       New: <code>GridGenerator</code>::<code
       class="member">hyper_rectangle</code> creates
       coordinate-parallel rectangles in arbitrary dimension.
       <br>
       (GK 2001/03/16)
       </p>

  <li> <p>
       Changed: The syntax of the <code
       class="class">FiniteElement&lt;dim&gt;</code>::<code
       class="member">get_unit_support_points</code> function is
       changed; it returns a reference to the vector of points in lieu
       of taking this vector as argument. The unit support points are
       now computed by the constructor of the <code
       class="class">FiniteElement</code> and not on each <code
       class="class">FiniteElement&lt;dim&gt;</code>::<code
       class="member">get_unit_support_points</code> function call as
       before.
       <br>
       (WB 2001/03/14)
       </p>

  <li> <p> 
       Removed: The function <code
       class="class">FiniteElement&lt;dim&gt;</code>::<code
       class="member">get_support_points</code>
       is removed as the <code>FiniteElement</code>
       cannot know the position of the support points by itself. This is
       because the support points depend on the unit support points
       (known to the <code>FiniteElement</code> class)
       but they also depend on the mapping.

       In future the support points can be computed using a <code
       class="class">FEValues</code> object initialized with the <code
       class="class">Mapping</code> and a <code
       class="class">Quadrature</code> that uses the unit support
       points as quadrature points.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p> 
       New: The class <code>Boundary</code> has two new
       functions <code
       class="member">get_intermediate_points_on_line</code> and <code
       class="member">get_intermediate_points_on_quad</code> that
       needs to be implemented in derived classes only if <code
       class="class">GridOut</code> is used with <code
       class="class">MappingQ</code> of degree <code>p>2</code>.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>

       New: The functions <code
       class="member">GridOut::write_gnuplot</code> and <code
       class="member">GridOut::write_eps</code> now take a pointer to
       a <code>Mapping</code> object as additional
       argument. This allows to write grids in 2d whereby cells with a
       curved boundary are transformed by the given <code
       class="class">Mapping</code> class. The default mapping is
       <code>MappingQ1</code>. Note, that the grids do
       not show the `original' domain with boundaries described by the
       <code>Boundary</code> classes but the
       discretized domain whose boundary cells are transformed using the
       given mapping.
       </p>
       There are also a new <code
       class="member">GnuplotFlags::n_boundary_face_points</code> and
       <code>EpsFlags::n_boundary_face_points</code>
       variables to set the number of additional points written
       to represent the curved boundary faces.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: The constructor of the <code>FEValues</code>
       and <code>FE(Sub)FaceValues</code> classes now
       take a reference of a <code>Mapping</code> object
       as additional argument. This is the new possibility to combine
       a <code>FiniteElement</code> with an arbitrary
       <code>Mapping</code>, see also <a
       href="#new_fe_mapping_design">New FE and Mapping Design</a>.
       </p>

       <p>
       For backward compatibility there still exists a constructor of
       <code>FEValues</code> without a <code
       class="class">Mapping</code> argument that uses a <code
       class="class">MappingQ1</code> by default.
       <br>
       (RH 2001/03/14)
       </p>
       
  <li> <p>
       Changed: We now unconditionally include
       <code>deal.II/include/grid/tria_accessor.templates.h</code>
       (which containts some inline functions for triangulation
       accessor classes) into 
       <code>deal.II/include/grid/tria_accessor.h</code> to work
       around a problem with gcc3.0 which does not place a copy of
       these functions into the library. Previously we only included
       the file in optimized mode.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: The class <code>GridReordering::Cell</code> has
       now a copy constructor to work around a bug in a gcc3.0
       snapshot.
       <br>
       (RH, WB 2001/03/14)
       </p>

  <li> <p>
       Changed: when refining a hexahedron in 3d, the rules by which
       the new vertices at the centers of the faces and at the center
       of the cell are placed, are changed. They are changed in a way
       as to minimize the distortion of the resulting cells from the
       optimal shape and to make them look as much alike as possible
       to generate a smoother grid.
       <br>
       (WB 2001/03/02)
       </p>

  <li> <p>
       Fix: the <code
       class="class">DoFTools</code>::<code
       class="member">compute_intergrid_constraints</code>
       function took memory quadratic in the number of degrees of
       freedom. This is now reduced to linear behaviour, with a
       constant that depends on the number of levels by which the two
       grids differ.
       <br>
       (WB 2001/02/26)
       </p>

  <li> <p>
       Fix: in the triangulation, the <code
       class="member">straight_boundary</code> variable, which is a
       reference, was assigned the address of a temporary object. It
       is unclear how this could have worked for three years, but it
       apparently did...
       <br>
       (WB 2001/02/26)
       </p>

  <li> <p> 
       New: The <code>DoFTools</code> class now has
       a function <code>count_dofs_per_component</code>
       that counts the number of degrees of freedom in each of the
       components of the finite element, i.e. how many degrees of
       freedom there are on the global mesh for each variable (field).
       <br>
       (WB 2001/02/22)
       </p>

  <li> <p>
       New: The <code>CellAccessor</code> class now has a function
       <code>has_boundary_lines</code> that mostly has
       the same semantics that <code>at_boundary</code>
       has, but also covers the case that a hexahedron in 3d may be at
       the boundary only by one line, rather than by a whole face. In
       that case, <code>at_boundary</code> reports
       false, while <code>has_boundary_lines</code>
       reports true.
       <br>
       (WB 2001/02/21)
       </p>

  <li> <p>
       New: There is now a function
       <code>ConstraintMatrix</code>::<code
       class="member">merge</code> that merges
       the constraints represented by two constraint matrices into one
       object.
       <br>
       (WB 2001/01/31)
       </p>
</ol>

*/
