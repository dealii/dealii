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
 * @page changes_between_3_3_and_3_4 Changes between Version 3.3 and 3.4

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>

<a name="general"></a>
<h3>General</h3>

<ol>
  <li> <p>
       <strong>
       Deprecated: The functions <code
       class="member">FEValuesBase::get_shape_values</code>, <code
       class="member">FEValuesBase::get_shape_grads</code>, and <code
       class="member">FEValuesBase::get_shape_2nd_derivatives</code> are
       now deprecated as they expose too much of the internal data
       structure of their class, and interfere with plans for the
       extension of this and related classes. The functions still
       exist in this release of the library, but will be removed in
       the next version. Use the <code
       class="member">FEValuesBase::shape_value</code> and alike
       functions as a replacement.
       <br>
       For more information, read
       <a href="http://www.dealii.org/mail/msg00638.html" target="body">this mail</a>.
       </strong>
       <br>
       (WB 2002/06/03)
       </p>

  <li> <p>
       Added: The <code>configure</code> script now recognizes Intel's ECC
       compiler when run on Itanium systems with this compiler present.
       The ECC compiler is similar to the ICC compiler but it acccepts some
       different options.
       <br>
       (BK 2002/05/22)
       </p>

  <li> <p> 
       New: The step-14 example program demonstrates duality based
       error estimators and some more software design concepts.
       <br>
       (WB 2002/05/05)
       </p>

  <li> <p> 
       New: In all previous versions, deal.II used
       the <a href="http://www.cs.wustl.edu/~schmidt/ACE.html"
       target="_top">ACE (Adaptive Communications Environment)</a>
       library to support cross-platform threading
       facilities. While this is still supported, the default way
       is now to use the POSIX threading functions that are
       available on many systems. The relieves you from the need of
       installing a huge library of which the most part is not used
       anyway. However, if you use ACE for other reasons, then it is
       still supported. For installation instructions, see the 
       <a href="../../readme.html" target="body">ReadMe</a> file.
       <br>
       (WB 2002/04/30)
       </p>

  <li> <p> 
       Changed: The Makefiles for the library are now truly
       parallel. To this end, the automatic generation of the files
       <tt>forward_declarations.h</tt> in the various directories had
       to be changed. They will now be generated automatically at the
       end of the targets <tt>all</tt>, <tt>debug</tt> and
       <tt>optimized</tt>. They will not be generated while building a
       single library. In this case, <tt>make forward</tt> can be used
       to build them manually.
       <br>
       I introduced separate targets for the generation of the
       optimized versions only.
       <br>
       (GK 2002/04/17)
       </p>

  <li> <p> 
       New: The step-13 example program tells you something about
       software design things for finite element programs.
       <br>
       (WB 2002/04/16)
       </p>

  <li> <p> 
       Changed: Due to problems with undeclared functions and general
       compatibility concerns, <code>-ansi</code> is now no more part
       of the compile flags for gcc.
       <br>
       (WB 2002/04/16)
       </p>

  <li> <p> 
       Fixed: Explicit specializations of member templates are now
       conforming to the C++ standard. While most compilers accepted
       the previous form, Sun's Forte compiler wants a strictly
       conforming one.
       <br>
       (WB 2002/03/25)
       </p>

  <li> <p> 
       Fixed: For gcc versions that used <code
       class="class">ostrstream</code> instead of <code
       class="class">ostringstream</code>, it was necessary to append
       a final <code>std::ends</code> when piping text
       into the string stream. This was not previously
       conditionalized, but done for old and new classes.
       <br>
       (WB 2002/03/13)
       </p>

  <li> <p> 
       Changed: The configure machinery has been revamped
       significantly.
       <br>
       (WB 2002/03/08)
       </p>

  <li> <p> 
       Added: The top-level Makefile now supports "optimized" as a
       target that builds only optimized versions of the <code>base</code>,
       <code>lac</code>, <code>1d</code>, <code>2d</code>, and <code>3d</code>
       libraries. 
       <br>
       (BK 2002/02/19)
       </p>

  <li> <p> 
       Changed: The build system was entirely revised. Object
       files in debug mode now have the suffix <code>.g.o</code>
       instead of <code>.go</code>. All object files from the
       subdirectories are now placed into the <code>/lib</code>
       top-level directory, rather than in library directories in the
       individual subdirs.
       <br>
       (WB 2002/02/11)
       </p>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p> 
       New: The <code>vector2d</code> row accessor
       classes now have member functions <code
       class="member">begin</code> and <code>end</code>
       which allow iterating over the elements of a row of such an
       object. 
       <br>
       (WB 2002/05/30)
       </p>

  <li> <p> 
       New: The <code>Legendre</code> and
       <code>LagrangeEquidistant</code> classes now have
       static member functions <code
       class="member">generate_complete_basis</code> which returns an
       array of polynomial objects spanning the complete space up to a
       specified order in 1d. This may be used to generate the
       respective polynomial spaces in higher space dimensions.
       <br>
       (WB 2002/05/27)
       </p>

  <li> <p> 
       Changed: The <code>Polynomial</code> and
       <code>LagrangeEquidistant</code> classes have lost
       their default constructor, as that did not make much sense
       anyway.
       <br>
       (WB 2002/05/27)
       </p>

  <li> <p> 
       Fixed: When forward declaring the <code
       class="class">Tensor</code> class, we now also forward declare
       its partial specialization for a rank one tensor. Not doing so
       confused Sun's Forte compiler.
       <br>
       (WB 2002/03/22)
       </p>

  <li> <p> 
       Fixed: The class <code>TensorFunction</code>
       now uses local types <code>value_type</code> and
       <code>gradient_type</code> as return values of
       its member functions. This works around a bug in Sun's Forte
       C++ compilers.
       <br>
       (WB 2002/03/20)
       </p>

  <li> <p> 
       Improved: The <code>AssertThrow</code> macro now
       uses <code>__builtin_expect</code> if the
       compiler supports this. This indicates to the compiler that we
       expect the condition to be true and that throwing an exception
       is a rare case. By this information, the compiler can help the
       branch prediction unit of modern processors to better predict
       which direction a branch will take.
       <br>
       (WB 2002/03/13)
       </p>

  <li> <p>
       New: The <code>vector2d</code> class now not only
       allows access to elements through the <code
       class="member">operator()(unsingned int,unsigned int)</code>
       (i.e. matrix or Fortran style access), but also through nested
       brackets via an <code>operator[]</code>
       (i.e. like to a two-dimensional C-style array).
       <br>
       (WB 2002/03/08)
       </p> 

  <li> <p>
       Changed: The function <code>MultithreadInfo</code>::
       <code>get_n_cpus</code> now reports the proper number
       of CPUs when running on Silicon Graphics.
       <br>
       (BK 2002/02/19)
       </p> 

  <li> <p> 
       Changed: The quite logorrhoeic function name <code
       class="class">TensorProductPolynomials</code>::<code
       class="member">n_tensor_product_polynomials</code> was changed to
       <code>n</code> to be compliant wth the new class <code
       class="class">PolynomialSpace</code>.
       <br>
       (GK 2002/02/11)
       </p>

  <li> <p> 
       New: The class <code>PolynomialSpace</code>
       implements the space of polynomials at most a certain degree in
       arbitrary space dimensions.
       <br>
       (GK 2002/02/11)
       </p>

  <li> <p> 
       New: The function <code>DataOutBase</code>::
       <code>write_tecplot_binary</code> has been
       added.  This function will write Tecplot binary files if the
       Tecplot API is detected by ./configure.  To use this feature be
       sure that the environment variable TECHOME points to a valid
       Tecplot installation and that the files
       $TECHOME/include/TECIO.h and $TECHOME/lib/tecio.a exist.  The
       name of the file to be written is specified through the <code
       class="class">DataOutBase</code> ::<code
       class="member">TecplotFlags</code>.  <code
       class="member">tecplot_binary_file_name</code> variable. If the
       API is not available this code simply calls the existing ASCII
       output function.
       <br>
       (BK 2002/02/11)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p> Improved: <code>SolverGMRES</code> allocates
       basis vectors only, when they are needed. Therefore, it is safe
       now to ask for a basis larger than the expected number of
       iteration steps. On the other hand, memory allocation failures
       may occur during the iteration now.
       <br>
       (GK 2002/05/24)
       </p>

  <li> <p> 
       New: Function <code
       class="member">SparsityPattern::matrix_position</code> is the
       inverse function for <code
       class="member">SparsityPattern::operator()</code>.
       <br>
       (WB 2002/02/13)
       </p>

  <li> <p> 
       New: Functions <code
       class="member">SparsityPattern::copy_from</code> and <code
       class="member">SparseMatrix::copy_from</code> allow to copy a full
       matrix into a sparse matrix.
       <br>
       (WB 2002/02/06)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p> 
       New: The <code>GeometryInfo</code> class now
       provides two methods,
       <code>unit_cell_vertex</code> and <code
       class="member">vertices_adjacent_to_line</code>, that reveal
       something about the placement and numbering of vertices on the
       uni cell.
       <br>
       (GK 2002/05/29)
       </p>

  <li> <p> 
       New: The <code>GridOut::</code>
       <code>write_dx</code> function is now implemented.
       It allows to write the mesh (cells and faces) with some additional
       information that may be useful once in a while.
       <br>
       (GK 2002/05/02)
       </p>

  <li> <p> 
       Fixed: The <code>IteratorState::IteratorState</code>
       enum is now called <code
       class="class">IteratorState::IteratorStates</code>. This works
       around a bug in Sun's Forte C++ compilers which can't handle
       members of namespaces with the same name as the enclosing
       namespace.
       <br>
       (WB 2002/03/20)
       </p>

  <li> <p> 
       Fixed: In 3d, the function <code
       class="member">DoFTools::make_hanging_node_constraints</code> 
       contained an assertion that failed erroneously for finite
       elements that do not have degrees of freedom on vertices. This
       is now fixed.
       <br> 
       (WB 2002/02/21)
       </p>

  <li> <p> 
       Fixed: <code>TriaAccessor<3,3>::measure</code>
       sometimes computed a negative value. This is now fixed.
       <br> 
       (WB 2002/02/21)
       </p>

  <li> <p> 
       New: Finite element family with complete polynomial spaces
       for discontinuous Galerkin: <code>FE_DGP</code>
       <br> 
       (GK 2002/02/11)
       </p>
</ol>

*/
