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
 * @page changes_between_3_4_and_4_0 Changes between Version 3.4 and 4.0

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author. Regular
contributor's names are abbreviated by WB (Wolfgang Bangerth), GK
(Guido Kanschat), RH (Ralf Hartmann).
</p>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li> <p> 
       New: deal.II now uses a new threading
       scheme. The new scheme is simpler to use, and in particular
       more flexible in some cases where only one thread is started,
       or where a thread is completely detached, since we got rid of
       the <code>ThreadManager</code> class and now only store handles
       to individual threads (which can be discarded, or added to a
       <code>ThreadGroup</code> variable that is able to wait for a
       whole set of threads at once.
       <br>
       The new scheme also implements a much needed feature: calling
       functions on a new thread that return values. Previously, such
       functions needed to be written in a way that they return their
       return value through an additional reference parameter. This
       was inflexible if one wanted to call functions that already
       exist. This restriction is now lifted: such functions can be
       called, and the return value can be accessed once the thread
       has finished.
       <br>
       (WB 2003/02/06)
       </p>

  <li> <p> 
       New: deal.II now makes use of some parts of
       the <a href="http://www.boost.org/">boost</a> library, which is
       supposed to be a testground for the next generation C++ standard
       library. The parts which we use are now in
       <code>contrib/boost/include/boost_local/</code> and can be
       referenced from within your programs. The directory contains
       the string <code>_local</code> since you may still want to use
       another version or installation of boost in your own programs.
       </p>
       
       <p>
       Also note that boost is large -- much larger than the subset we
       have imported --, so we only took what we needed.
       <br>
       (WB 2003/02/06)
       </p>

  <li> <p> 
       Fixed: A longstanding bug in the documentation system has been fixed: if
       a namespace was ended with <code>}</code> instead of <code>};</code>
       (note the semicolon), then the documentation tool assumed that the
       following entities were still part of the namespace just left. This was
       since the closing brace alone was not accepted as ending a namespace
       (after all, structures, classes, and enums -- the other entities that
       can enclose other declarations -- need the semicolon). This led to some
       classes not showing up in the class index of the sublibraries. This is
       now fixed.
       <br>
       (WB 2003/02/02)
       </p>

  <li> <p> 
       Changed: Classes and structures that are declared inside namespaces named
       <code>internal</code> are now no longer shown in the class
       index view of the documentation of each of the sublibraries. Since they
       are supposed to be used internally only, this is no drawback for the
       general public. However, they are documented as members of these
       namespaces.
       <br>
       (WB 2003/02/02)
       </p>

  <li> <p> 
       Fixed: Some of the formulas in the step-14 tutorial were obviously 
       scrambled a little. This is now fixed.
       <br>
       (Roy Stogner, WB 2003/01/30)
       </p>

  <li> <p> 
       Changed: The main Makefile has been changed to sequentialize building 
       the base, lac, and deal.II sublibraries. We changed this, since on some
       systems (notably AIX), the latter libraries need to be linked against
       the former ones, when creating shared libraries.
       <br>
       (WB 2003/01/24)
       </p>

  <li> <p> 
       New: Changes have been made to support compiling and using 
       deal.II on AIX 5 systems.
       <br>
       (WB 2003/01/24)
       </p>

  <li> <p> 
       Removed: Thread support now relies solely on the use of POSIX
       functions. The use of the
       <a href="http://www.cs.wustl.edu/~schmidt/ACE.html" target="_top">ACE
       (Adaptive Communications Environment)</a> library for this is now no
       longer supported. However, application programs can of course still use
       ACE, but they will need to generate paths to this library in their
       makefiles themselves.
       <br>
       (WB 2003/01/11)
       </p>

  <li> <p> 
       New: Some changes have been made to support Mac OS X 10.2. Shared
       libraries are not supported on this architecture, but everything else
       should work.
       <br>
       (WB 2002/12/18)
       </p>

  <li> <p> 
       New: deal.II can be compiled with version 7.0 of
       Intel's icc compiler, which was recently released. Since this compiler
       finally supports the very restrictive flags <code>-Xc -ansi</code> that
       check for close conformance with the C++ standard, we use them (previous
       versions of icc would crash when these two flags are given). This
       requires that we distinguish between these compiler versions, and the
       corresponding Makefile variable <code>GXX-VERSION</code> now no longer
       holds the non-versioned string <code>intel_icc</code> when icc is
       detected, but rather either <code>intel_icc5</code>,
       <code>intel_icc6</code>, or <code>intel_icc7</code>, depending on what
       version of the compiler was detected.
       <br>
       (WB 2002/12/05)
       </p>

  <li> <p> 
       Changed: Previously, we just set the preprocessor variable
       <code>DEAL_II_USE_MT</code>, when <code>--with-multithreading</code> was
       given as argument to <code>./configure</code>. Tests in the code
       therefore looked like <code>#ifdef DEAL_II_USE_MT</code>. This has been
       changed so that the variable is always defined, but its value is now
       equal to <code>1</code> 
       when multithreading was requested, and zero otherwise. The reason for
       this is that you can now write <code>if (DEAL_II_USE_MT && ...)</code>
       conditions, and need not interleave if-else clauses from regular code
       and the preprocessor, if conditions involve both the state of this
       preprocessor variable and the run-time state of your program. However,
       this change requires that all appearances of <code>#ifdef
       DEAL_II_USE_MT</code> be changed to <code>#if DEAL_II_USE_MT ==
       1</code>, since the variable is now defined unconditionally.
       <br>
       (WB 2002/11/14)
       </p>

  <li> <p> 
       New: Object files are now named according to the local defaults
       on the system we are running on. On Unix systems, this is
       usually a <code>.o</code> suffix, while on Windows it is
       <code>.obj</code>. Likewise for executables, which have no
       suffix on Unix, but <code>.exe</code> on Windows.
       <br>
       (WB 2002/11/11)
       </p>

  <li> <p> 
       New: deal.II can now also be compiled with Red Hat's
       version of the gcc compiler, gcc 2.96. However, some problems remain;
       for more information.
       <br>
       (WB 2002/10/14)
       </p>

  <li> <p> 
       Fixed: On CygWin, one header files has a <code>#define quad
       quad_t</code>. This is annoying, since we have local variables and
       member functions with the name <code>quad</code>, and in fact it breaks
       compilation on those versions of CygWin that has this. Fortunately, the
       define is only active if a preprocessor variable
       <code>_POSIX_SOURCE</code> is not set. Thus, we now check for the define
       when configuring the library, and if necessary set the preprocessor
       variable. However, while this allows to compile the library on these
       systems, it may otherwise affect your code, if you use functions or
       other features of the system that are not available when the flag is
       set.
       <br>
       (Stephen Kolaroff, WB
       2002/09/28) 
       </p>

  <li> <p> 
       New: Since <code>math.h</code> only defines the values of PI or E (as
       <code>M_PI</code> and <code>M_E</code>) when certain defines are set (on
       Linux, these are <code>__USE_BSD</code> or <code>__USE_XOPEN</code>),
       portable programs usually defined these constants themselves. In
       deal.II, this happened at 6 different places. To
       avoid this in the future, <code>base/config.h</code> now exports a
       namespace <code>deal_II_numbers</code> that defines these two, and a
       number of other numerical constants.
       <br>
       (WB 2002/09/12)
       </p>

  <li> <p> 
       New: <code>base/config.h</code> now exports the
       deal.II base directory through the
       <code>DEAL_II_PATH</code> preprocessor variable.
       <br>
       (WB 2002/09/09)
       </p>

  <li> <p> 
       Removed: The forward declarations files have gone. We have never
       propagated their use in the example programs, but these files have been
       there in the base, lac, and grid include directories, and forward
       declared all classes that were present in the respective parts of the
       library. This, the idea was, enables you to use just this include file
       in your own header files, rather than including the full declarations of
       these classes. However, maintaining these forward declaration files has
       been a constant thorn in our side, be it that the timing of their
       generation was difficult when using parallel builds, or that they were
       difficult to generate at first. The latter is now the reason for their
       abolition: we had a script for their generation, but it did not take
       into account namespaces, so we got clashes when we found that we had
       used the same class name in two different namespaces, since the script
       put the forward declaration incorrectly into the global namespace where
       they conflicted. Since we do not plan to extend the script by a
       parser that can properly handle opening and closing braces of
       namespaces, we simply drop these files.
       <br>
       What you should do if you have used these forward declaration files: you
       have two possibilities - either include the respective header file in
       which the class is fully declared, or write the forward declaration into
       your headers yourself.
       <br>
       (WB 2002/09/05)
       </p>

  <li> <p> 
       New: There is now
       a new report
       on assembling matrices available from the 
       <a href="../documentation.html" target="body">documentation
       page</a>. The main focus is assembling of matrices for
       vector-valued problems, where shape functions are
       vector-valued, and may have only one or more non-zero vector
       components. 
       <br>
       (WB 2002/06/17)
       </p>

  <li> <p>
       Removed: The functions <code
       class="member">FEValuesBase::get_shape_values</code>, <code
       class="member">FEValuesBase::get_shape_grads</code>, and <code
       class="member">FEValuesBase::get_shape_2nd_derivatives</code> are
       now removed as they expose too much of the internal data
       structure of their class, and interfere with plans for the
       extension of this and related classes. These functions, which
       had been deprecated in the previous version, are thus not
       part of the release after version 3.4 any more. Use the <code
       class="member">FEValuesBase::shape_value</code> and alike
       functions as a replacement.
       <br>
       For more information, read
       <a href="http://www.dealii.org/mail/msg00638.html" target="body">this mail</a>.
       <br>
       (WB 2002/06/10)
       </p>

  <li> <p> 
       New: deal.II now also supports vector-valued
       finite elements with shape functions for which more than just
       one vector component is non-zero. Such elements are, for
       example, the Nedelec and Raviart-Thomas families. Previously,
       vector-valued elements were only supported insofar as they
       could be composed of scalar ones; in that case, each
       (vector-valued) shape function had only one non-zero vector
       component.
       <br>
       (WB 2002/06/10)
       </p>

  <li> <p> 
       New: The top-level makefile now how a target <tt>distclean</tt>.
       <tt>clean</tt> leaves the libraries now, removing
       everything that is not needed to use
       deal.II. <tt>distclean</tt> removes even the
       libraries, leaving the directory more or less in the state like
       after <tt>configure</tt>.
       <br>
       (GK 2002/06/07)
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li> <p>
       Fixed: A bug in the <code
       class="class">Patterns::MultipleSelection</code> class if more
       than two elements in a comma-separated list were given.
       <br>
       (Brian Carnes 2003/05/14)
       </p>

  <li> <p>
       Changed: The <code>Polynomials::Legendre</code>
       class lost its template argument and is now just a regular
       class. There was no real good reason for the template argument,
       it had just crept in.
       <br>
       (WB 2003/05/12)
       </p>

  <li> <p>
       New: There is now a class <code
       class="class">AnisotropicPolynomials</code> that constructs a higher
       dimensional polynomial space from a set of 1-d polynomials in each of
       the coordinate directions.
       <br>
       (WB 2003/04/20)
       </p>

  <li> <p>
       Changed: The <code>Table</code> accessor classes
       have been moved to a namespace <code
       class="class">internal</code>. Since these classes are not (or should
       not be) used directly in applications, this should not change
       compatibility. However, they will now no longer show up in the class
       overview of the documentation, which they were cluttering up.
       <br>
       (WB 2003/02/13)
       </p>

  <li> <p> 
       New: The <code>Function</code> class now has an assignment
       operator. This way, you can put function objects into
       containers. However, the assignment operator makes sure that only
       objects can be assigned that have the same number of vector components.
       <br>
       (WB 2003/02/06)
       </p>

  <li> <p> 
       New: The <code>ThreadMutex</code> classes now have a
       member class <code>ScopedLock</code> that implements the
       scoped thread-safe locking pattern of Doug Schmidt. It is also used in
       various places of the code now.
       <br>
       (WB 2003/01/28)
       </p>

  <li> <p> 
       Fixed: The <code>PosixThreadManager</code> called
       its <code>wait</code> function in the
       destructor. If this had been called before already, then the
       same threads would have been waited for twice, which invokes
       undefined behavior. This is fixed by making sure that <code
       class="member">wait</code> removes the id's of the threads it
       has already waited for, and so calling it more than once will
       not wait for threads which have already been waited for.
       <br>
       (Michael Anderson, WB 2003/01/16)
       </p>

  <li> <p> 
       Fixed: The <code>Subscriptor</code> uses a counter to
       count how many <code>SmartPointer</code> objects subscribe
       to the pointed-to object. This counter needs to be a volatile variable
       in multithreaded mode, to avoid false compiler optimizations based on
       the assumption that the variable cannot change between two subsequent
       reads.
       <br>
       (WB 2003/01/11)
       </p>

  <li> <p> 
       Fixed: In multithreaded mode, when a new thread is started, the
       arguments to the function being called need to be copied from the stack
       of the starting thread to that of the new thread. In order to
       synchronise this, mutexes were used that were acquired from one thread
       and released from another thread. On Linux this does not lead to
       problems, but POSIX functions do not guarantee that this actually works,
       and it also leads to problems when running programs under valgrind. This
       is now fixed with the help of condition variables.
       <br>
       (Michael Anderson, WB 2003/01/11)
       </p>

  <li> <p> 
       New: There are now classes <code>ThreadCondition</code>
       that implement thread condition variable operations through POSIX
       functions (in multithreaded mode) or doing nothing (in singlethreaded
       mode).
       <br>
       (WB 2003/01/11)
       </p>

  <li> <p> 
       New: Newer versions of gcc have a very nice feature: you can set
       a verbose terminate handler, that not only aborts a program
       when an exception is thrown and not caught somewhere, but
       before aborting it prints that an exception has been thrown,
       and possibly what the std::exception::what() function has to
       say. Since many people run into the trap of not having a
       catch clause in main(), they wonder where that abort may be
       coming from. The terminate handler then at least says what is
       missing in their program.
       <br>
       (WB 2002/12/19)
       </p>

  <li> <p> 
       New: There is now a <code>Patterns::List</code> pattern
       that matches a list of elements each of which has to satisfy a pattern
       that was given to the constructor.
       <br>
       (WB 2002/11/29)
       </p>

  <li> <p> 
       Changed: In POSIX mode, when the <code
       class="member">ThreadManager</code> class created a new thread through
       <code>pthread_create</code>, it only checked for the
       error code and aborted if it was non-zero. Now, it checks whether the
       error code is <code>EAGAIN</code> and simply retries the
       call if this is the case. This may, in rare cases, lead to a deadlock or
       an infinite loop, but will usually just wait until the respective
       resources for thread creation are available from the operating system
       and will then succeed.
       <br>
       (WB 2002/11/13)
       </p>

  <li> <p> 
       Fixed: The <code>write_text</code> and <code
       class="member">write_tex</code> functions of the <code
       class="class">TableHandler</code> class now check whether their
       <code>ofstream</code> arguments are in a proper state before
       using them for output.
       <br>
       (RH 2002/11/11)
       </p>

  <li> <p> 
       New: Added Hierarchical Polynomial (similar to Legendre class). Will
       eventually be used in a hierarchical FiniteElement class similar to
       FE_Q class. Included in Polynomials namespace.
       <br>
       (Brian Carnes 2002/10/15)
       </p>

  <li> <p> Changed: Because they became too many, the classes describing 1d
       polynomials are now in a <code>namespace
       Polynomials</code>.
       <br>
       (WB 2002/10/14)
       </p>

  <li> <p> Changed: When an exception is thrown but not caught in a sub-thread,
       this exception is not passed to the main thread by the operating
       system. Rather, if the exception is not caught from the function that
       was invoked by the spawning system function, the entire program is
       terminated without an additional message. The wrapper functions which
       are used to spawn new threads in the <code>Threads</code>
       namespace therefore now catch these exceptions and at least print the
       message they carry, before aborting the program. This way, at least the
       message gets displayed.
       <br>
       (WB 2002/10/13)
       </p>

  <li> <p> Changed: The class <code
       class="member">Table&lt;2&gt;::fill</code> function, which is also
       inherited from the <code>FullMatrix</code> class, used to
       work also when the size of the matrix was zero, by simply not copying
       something. This led to difficult to detect errors. It is therefore no
       more allowed to call this function when the matrix is empty. For all
       other cases, the status of copying without checking the size of the
       array copied from remains unchanged.
       <br>
       (WB 2002/09/28)
       </p>

  <li> <p> New: The classes <code
       class="class">TableIndices&lt;N&gt;</code> and <code
       class="class">Table&lt;N,T&gt;</code> are now implemented also
       for <code>N=4,5</code> and <code>6</code>. The <code
       class="class">Table&lt;N,T&gt;</code> class represents an
       <code>N</code>-dimensional array and might replace the
       N-times-nested-use of the <code
       class="class">std::vector</code> class.
       <br>
       (RH 2002/09/24)
       </p>

  <li> <p> 
       New: The <code>Threads::n_existing_threads</code>
       function returns the present number of existing threads, allowing an
       assessment whether it is useful to spawn new threads or rather perform
       the operation in question sequentially.
       <br>
       (WB 2002/09/08)
       </p>

  <li> <p> New: Global exception class <code
       class="class">ExcIteratorPastEnd</code>, which should be used if an
       iterator is incremented or decremented beyond its end value.
       <br>
       (GK 2002/09/08)
       </p>

  <li> <p> 
       Extended: Previously, the <code
       class="class">Threads::PosixThreadBarrier</code>
       class could not be used at all (they threw exceptions), if your system
       did not have the respective POSIX functions. This restriction is lifted
       for the special case that you give one as the number of parties that
       will be waiting for the barrier: in this case, a barrier is a
       no-operation, so we do not need assistence from the operating
       system. This change makes it slightly simpler to write programs in a way
       such that they run in both single- and multithreaded environments.
       <br>
       (WB 2002/09/05)
       </p>

  <li> <p> 
       New: The old class <code>vector2d</code>, implementing a
       two-dimensional array of objects is now gone. Instead, there is the new
       <code>Table</code> class that implements tables of
       arbitrary dimension. Transition is painless: where there was
       <code>vector2d&lt;type&gt;</code> before, use
       <code>Table&lt;2,type&gt;</code> now (and don't forget to update the
       name of the header file, of course). If you have a three-dimensional
       array, you can use <code>Table&lt;3,type&gt;</code> now.
       <br>
       (WB 2002/09/03)
       </p>

  <li> <p> 
       New: There are now functions returning the transpose of <code
       class="class">Tensor</code> objects of rank 2.
       <br>
       (WB 2002/08/31)
       </p>

  <li> <p> 
       New: Row accessors for the <code>vector2d</code>
       class now have a member function <code
       class="member">size</code> that returns the size of the row,
       i.e. the number of columns of the table.
       <br>
       (WB 2002/07/24)
       </p>

  <li> <p> 
       Fixed: In EPS output, colors were set to invalid values if the
       values of the field that is used for coloring are all
       equal. This happens, for example, in the very first time step
       of time dependent problems, if initial values are zero. The
       color value now used is arbitrary, but valid.
       <br>
       (WB 2002/07/24)
       </p>

  <li> <p> 
       Changed: To save disk space, color values in EPS output are
       written as grey scale with only one value instead of three RGB
       values if the color so represented is actually a grey scale.
       <br>
       (WB 2002/07/24)
       </p>

  <li> <p> 
       New: There are now operators to allow multiplication and
       division of <code>Tensor</code> objects by scalar
       factors.
       <br>
       (WB 2002/06/07)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p> New: Function 
       <code>BlockSparseMatrix::print_formatted</code>
       added, similar to the 
       <code>BlockVector::print_formatted</code> function.
       <br>
       (Brian Carnes 2003/06/02)
       </p> 

  <li> <p> New: Functions <code>SparseMatrix::operator *=</code>
       and <code>SparseMatrix::operator /=</code> as well as 
       <code>BlockSparseMatrix::operator *=</code>
       and <code>BlockSparseMatrix::operator /=</code> 
       are added.
       <br>
       (Brian Carnes 2003/06/02)
       </p> 

  <li> <p> Deprecated: The functions <code>Vector::scale</code>
       and <code>BlockVector::scale</code> are now deprecated
       and will be removed in a future version. Use <code>operator*=</code> and
       <code>operator/=</code> instead.
       <br>
       (WB 2003/05/31)
       </p> 

  <li> <p> New: <code>Vector</code>, <code
       class="class">BlockVector</code> and
       <code>FullMatrix</code> now have <code>operator/=</code>
       for scaling by a scalar.
       <br>
       (WB 2003/05/31)
       </p> 

  <li> <p> New: <code>PointerMatrix</code> now
       has <code>empty()</code> function, which returns true if the pointer is
       null; otherwise we call the pointer's <code>empty()</code> function.
       This requires the class MATRIX to have an <code>empty()</code>
       function.
       <br>
       (Brian Carnes 2003/05/22)
       </p>

  <li> <p> New: <code>SparseLUDecomposition</code> now
       has <code>empty()</code> function, which calls the inherited <code
       class="class">SparseMatrix</code> <code>empty()</code> function.
       <br>
       (Brian Carnes 2003/05/22)
       </p>

  <li> <p> New: <code>PreconditionPSOR</code> implements
       interface to permuted SOR preconditioner function in <code
       class="class">SparseMatrix</code>.
       <br>
       (GK 2003/05/12)
       </p>

  <li> <p> Improved: <code>FullMatrix</code>::<code
       class="member">fill</code> now copies the largest possible block,
       whether the destination or source matrix is bigger. Additionally, an
       offset inside the source matrix may be specified.
       <br>
       (GK 2003/03/31)
       </p>

  <li> <p>
       New/Changed: The <code>SparseILU</code>, <code
       class="class">SparseMIC</code> and <code
       class="class">SparseLUDecomposition</code> now use the same
       interface (<code>initialize</code>, <code
       class="class">vmult</code>, <code>clear</code>)
       as all <code>PreconditionBlock</code> classes. In
       virtue of an unified preconditioner interface it is now
       recommended to use the new methods. The old methods (<code
       class="class">reinit</code>, <code
       class="class">decompose</code>, <code
       class="class">apply_decomposition</code>) are now deprecated,
       and will be removed in a later version.
       <br>
       (RH 2003/02/19)
       </p>

  <li> <p>
       Changed: The <code>BlockVector</code> accessor classes
       have been moved to a namespace <code
       class="class">internal</code>. Since these classes are not (or should
       not be) used directly in applications, this should not change
       compatibility. However, they will now no longer show up in the class
       overview of the documentation, which they were cluttering up.
       <br>
       (WB 2003/02/13)
       </p>

  <li> <p> 
       New: The <code>SolverGMRES</code> now accepts an
       <code>AdditionalData</code> parameter
       <code>use_default_residual</code> whose default is
       <code>true</code>. By setting this flag to <code>false</code>,
       the stopping criterion of the left-preconditioned GMRes solver
       is not the default preconditioned residual but the normal
       (unpreconditioned) residual and visa versa for the
       right-preconditioned GMRes solver. Due to a performance loss of
       the solver this flag should be set to <code>false</code> only
       for debugging/testing purposes.
       <br>
       (RH 2003/01/31)
       </p>

  <li> <p> 
       New: <code>FullMatrix</code> has a function <code
       class="member">copy_from</code>, copying from sparse matrices.
       It uses iterators of the sparse matrix classes.
       <br>
       (GK 2003/01/08)
       </p>


  <li><p>
       Changed: In an attempt to unify the use of preconditioners a
       little, the function <code>initialize</code> of
       classes <code>PreconditionRelaxation</code> and
       <code>PreconditionBlock</code> take an argument
       of type <code>AdditionalData</code>, defined in
       the same class.  Standard behavior of <code
       class="class">PreconditionBlock</code> has been changed on this
       occasion to invert diagonal blocks during initialization.
       <br>
       (GK 2003/01/06)
       </p>

  <li> <p>
       New: The interface for sparse decompositions has been abstracted, and
       there is now an Modified Incomplete Cholesky (MIC) decomposition in
       addition to the Incomplete LU (ILU) decomposition.
       <br>
       (Stephen Kolaroff 2002/11/27)
       </p>
 
  <li> <p>
       Changed: In multithread mode, the <code
       class="class">SparseMatrix</code> would spawn
       <code>multithread_info.n_default_threads</code> threads to
       perform matrix-vector multiplications and similar
       operations. It would even do so if
       <code>multithread_info.n_default_threads</code> was equal to
       one. In that case, we now do the operation on the thread we are
       presently on, eliminating the overhead of spawning a single
       thread, and later waiting and terminating it.
       <br>
       (WB 2002/11/13)
       </p>
 
  <li> <p>
       Fixed: In the <code>SparseDirectMA27</code> class, wrapping 
       the MA27 solver written in Fortran77 into some structure amenable to C++,
       we wrote beyond the end of an array on rare inputs. This is now fixed. The
       same probably holds for the respective code for MA47.
       <br>
       (WB 2002/09/30)
       </p>

  <li> <p>
       New: Since the MA27 sparse direct solver uses Fortran common blocks, it
       was previously impossible to run several instances of this solver in
       parallel, in a multihreaded environment. To solve this problem, the
       <code>SparseDirectMA27</code> class now has a detached
       mode, in which it forks off a separate program that will do the
       computations using this solver. The actual operations are therefore
       distributed to distint programs that have separate address spaces. This
       allows to have as many concurrent instances of this solver in parallel
       as you want. For more information, read the documentation of the 
       <code>SparseDirectMA27</code> class.
       <br>
       (WB 2002/09/25)
       </p>

  <li> <p> Changed: The classes <code
       class="class">PreconditionBlock</code>, <code
       class="class">PreconditionBlockJacobi</code>, <code
       class="class">PreconditionBlockSOR</code>, and <code
       class="class">PreconditionBlockSSOR</code> have changed their
       template signature. The first template argument is now the matrix
       type, not just a number type.
       <br>
       (GK 2002/09/18)
       </p>

  <li> <p> New: Class <code>BlockVector</code> has a
       function <code>collect_sizes()</code>, very much as
       <code>BlockSparsityPattern</code>. This allows
       updating internal structures after blocks have been resized.
       <br>
       (GK 2002/09/17)
       </p>

  <li> <p> New: Class <code>SparseMatrix</code> has an
       STL-conforming <code>const_iterator</code> and
       functions <code>begin()</code> and <code
       class="member">end()</code> for looping through all existing
       entries. Furthermore, <code>begin(row)</code> and
       <code>end(row)</code> allow looping through all
       entries of a single line.
       <br>
       (GK 2002/09/11)
       </p>

  <li> <p>
       New: Classes <code>SparsityPattern</code> and <code
       class="class">SparseMatrix</code> now have functions <code 
       class="member">block_write/block_read</code>, allowing to dump the data
       of these objects into a file in binary format, and later to re-read it
       without much need for parsing.
       <br>
       (WB 2002/09/09)
       </p>

  <li> <p>
       New: <code>Vector</code> has a function <code
       class="member">lp_norm</code>, computing the <i>l<sub>p</sub></i>-norm
       of a vector for arbitrary <i>p</i>.
       <br>
       (GK 2002/08/13)
       </p>

  <li> <p> 
       New: a way of using abstract base classes for matrices has
       been implemented with <code>PointerMatrixBase</code>
       and <code>PointerMatrix</code>. Storing a matrix in
       <code>PointerMatrix</code> allows to use the base
       class in functions only templated for the vector class.
       <br>
       (GK 2002/07/18)
  </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p>
       Fixed: The restriction matrices for the Q1 element in 1d had a
       trivial bug in that one element was not set. Due to the fact
       that contributions from all child cells are taken into account,
       this did no harm, though, since all computations were done
       correctly anyway.
       <br>
       (WB 2003/05/06)
       </p>

  <li> <p>
       New: The <code>GeometryInfo</code> classes now
       have new static member functions <code
       class="member">child_cell_from_point</code> that, given a point
       in the unit cell, returns which child cell it is on; <code
       class="member">cell_to_child_coordinates</code> that transforms
       coordinates between the unit coordinate systems from the mother
       to the child cell; <code
       class="member">child_to_cell_coordinates</code> that does
       exactly the opposite; and <code
       class="member">is_inside_unit_cell</code> that tells whether a 
       given point is inside the unit cell.
       <br>
       (WB 2003/05/01)
       </p>

  <li> <p>
       New: There are now functions <code
       class="member">recursively_set/clear_user_pointer</code> that
       do much the same as  <code
       class="member">recursively_set/clear_user_flag</code> on a
       line/quad/hex and all of its descendants, but on user pointers
       rather than flags.
       <br>
       (WB 2003/04/11)
       </p>

  <li> <p>
       New: The functions in the <code
       class="class">DerivativeApproximation</code> class can now also
       work on <code>BlockVector</code>.
       <br>
       (WB 2003/04/11)
       </p>

  <li> <p>
       Changed: The cell argument to <code
       class="member">Mapping::transform_unit_to_real_cell</code> (and its
       reverse function) is now passed by reference, rather than by value, for
       efficiency reasons.
       <br>
       (WB 2003/04/06)
       </p>

  <li> <p>
       New: The <code
       class="member">GridGenerator::subdivided_hyper_cube</code> generated a
       hypercube as coarse grid that is subdivided a given number of times.
       <br>
       (WB 2003/04/05)
       </p>

  <li> <p>
       New: The <code>DataOut_DoFData::merge_patches</code>
       function allows to merge the patches of two objects, so as to create one
       output file from several DoF handlers. This is useful if one uses a
       domain decomposition algorithm where each block of the domain is
       represented by one DoF handler.
       <br>
       (WB 2003/04/05)
       </p>

  <li> <p>
       Fixed: The <code>DataOutStack</code> class had a problem
       when there were as many degrees of freedom as there were cells (i.e. if
       we were using DG0 elements). This should now be fixed.
       <br>
       (WB 2003/03/02)
       </p>

  <li> <p>
       Fixed: The <code>DataOutFaces</code> class was
       broken for cell data. It should now be correct, although the
       algorithm used is not optimal, being approximately quadratic in
       runtime.
       <br>
       (WB 2003/02/25)
       </p>

  <li> <p>
       New: The <code>ConstraintMatrix::shift</code>
       function shifts and translates the elements of the constraint
       matrix by a certain number of indices.
       <br>
       (Roy Stogner 2003/02/24)
       </p>

  <li> <p>
       New: The <code>GridReordering</code> class now uses a
       vastly better algorithm in 2d than previously. The new algorithm is
       linear in time, where it could be exponential before.
       <br>
       (Michael Anderson 2003/02/21)
       </p>

  <li> <p>
       New: There is now a function <code
       class="member">GridIn::read_xda</code> that allows reading
       grids from a file in XDA format.
       <br>
       (WB 2003/02/21)
       </p>

  <li> <p> 
       Changed: Some implementation details of the <code
       class="class">GridReordering</code> class have been moved to a
       namespace <code>internal</code>.
       <br>
       (WB 2003/02/19)
       </p>

  <li> <p>
       New: There are now functions <code
       class="member">FiniteElement::prolongation_is_implemented</code> and
       <code>FiniteElement::constraints_are_implemented</code>
       that inform the caller about whether the finite element in question in
       fact implements its prolongation matrices and hanging node constraints.
       <br>
       (WB 2003/02/15)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">ConstraintMatrix::is_identity_constrained</code> function
       would previously generate a segmentation fault if called on a constraint
       matrix object that did not contain any constraints at all. This is now
       fixed.
       <br>
       (WB 2003/02/15)
       </p>

  <li> <p>
       Fixed: Objects of type <code>FESystem</code> 
       could not be constructed if one of the elements it is to be composed of
       did not implement interface constraints for hanging nodes. This is now
       fixed: you can construct such a composed element, but it does not
       implement hanging node constraints either.
       <br>
       (WB 2003/02/15)
       </p>

  <li> <p>
       New: For each of the renumbering functions in the <code
       class="class">DoFRenumbering</code> class there is now an
       additional <code>compute_*</code>
       function. These new functions compute and return the
       corresponding renumbering vectors but do not perform the actual
       renumbering on the <code>DoFHandler</code>
       object. The behaviour of the old functions is not changed.
       <br>
       (RH 2003/02/03)
       </p>

  <li> <p> 
       Fixed: The <code>GridReordering</code> tried to be
       thread-safe in the initialization of some data, but was not due to a
       typo. This is now fixed.
       <br>
       (WB 2003/01/28)
       </p>

  <li> <p> 
       Changed: The <code>FEValues::get_cell</code> and
       <code>FEValues::get_face</code> functions have
       been removed, since they limited our ability to use this
       class for other types of DoFHandlers, for example future
       extensions for hp elements.
       <br>
       (WB 2003/01/17)
       </p>

  <li> <p> 
       New: The DoF accessor classes now have a function <code
       class="member">get_fe()</code> that returns a reference to the finite
       element being used on this cell. The result is of course identical to
       what a call to <code
       class="member">iterator-&gt;get_dof_handler().get_fe()</code> would have
       yielded.
       <br>
       (WB 2003/01/10)
       </p>

  <li> <p> 
       New: Checked in new <code>GridGenerator</code> 
       member function <code>half_hyper_ball</code>,
       derived from member <code>hyper_ball</code>.     
       The initial mesh contains four elements. This mesh will work with
       the boundary class <code>HalfHyperBallBoundary</code>.
       <br>
       (Brian Carnes 2002/12/16)
       </p>

  <li> <p> 
       New: Checked in new class <code>FE_Q_Hierarchical</code>
       derived from class <code>FiniteElement</code>.  
       This element is analogous to <code>FE_Q</code>, but 
       makes use of hierarchical shape functions, based on the 
       <code>Polynomials::Hierarchical</code> class. 
       For <code>degree&gt;1</code>, the non-nodal basis functions are "bubble"
       functions, which are not Lagrange polynomials. Therefore, the usual
       interpolation based on using unit support points will not work for
       <code>degree&gt;1</code>. It is planned to implement a different
       interpolation-projection operator, based on an hp-type interpolant.
       <br>
       The files for this element are
       <code>deal.II/include/fe/fe_q_hierarchical.h</code> and 
       <code>deal.II/source/fe/fe_q_hierarchical.cc</code>.
       <br>
       (Brian Carnes 2002/12/13)
       </p>

  <li> <p> 
       New: For finite element classes, the functions
       <code>unit_support_point</code> and
       <code>unit_face_support_point</code> return the position
       of an individual support point. This is necessary when you want to get
       information about the support points of certain components in a composed
       finite element, where not all components provide support points, and the
       composed element thus does not fill the array the 
       <code>get_unit_support_points</code> function returns.
       <br>
       (WB 2002/11/23)
       </p>

  <li> <p> 
       Fixed: Vectors could not be given as input and output vectors to the
       <code>SolutionTransfer</code> class at the same time, but
       this was not checked. An assertion has now been added to ensure this
       requirement. 
       <br>
       (WB 2002/10/28)
       </p>

  <li> <p> 
       Fixed: The <code
       class="member">DoFRenumbering::component_wise</code> function accepts a
       parameter indicating the order in which the degrees of freedom
       corresponding to one vector component are to be sorted. However, it did
       not honor this order, but always sorted them in the order in which the
       components appear. This is now fixed.
       <br>
       (WB 2002/10/17)
       </p>

  <li> <p> 
       New: The <code
       class="member">FiniteElement::system_to_base_index</code> function now
       exports the values of the <code
       class="member">FiniteElement::system_to_base_table</code> member
       variable. Likewise for the indices on faces.
       <br>
       (WB 2002/10/17)
       </p>

  <li> <p> 
       New: The <code
       class="member">FiniteElement::element_multiplicity</code> function was
       previously only available in the <code>FESystem</code>
       class, where it actually returned a non-trivial value. However, in some
       cases one would need to access this field also for general finite
       elements, even if the returned value will be equal to one in all cases
       other than a composed element.
       <br>
       (WB 2002/10/17)
       </p>

  <li> <p> 
       Fixed: The algorithm to generate neighbor information for patches from
       cells in <code>DataOut::build_patches</code> was rather
       inefficient, tripling the time for patch generation when support for
       neighbor information was added. Furthermore, the algorithm was at least
       O(N log N), where the rest was all O(N), making this particularly
       problematic when the data set was already large. This should now be back
       to the previous level, by using a more efficient algorithm.
       <br>
       (WB 2002/10/10)
       </p>

  <li> <p> 
       Changed (internals): Previously, the finite element base class
       initialized the restriction, prolongation, and face constraints matrices
       to their correct size. Derived classes had to fill these classes, and
       should have set their size back to zero in case they chose not to
       implement them. However, we found a class that forgot to resize it to
       zero, so it is now the other way round: they remain at size zero, and a
       class that chooses to implement these matrices has to set them to the
       correct size, to avoid programs that run on data that as just been
       forgotten to add. (This information only concerns programs that
       implement some finite element class on their own.)
       <br>
       (WB 2002/09/28)
       </p>

  <li> <p>Improved: The different transfer functions in <code
       class="class">FETools</code> operate on template vector arguments.
       <br>
       (GK 2002/09/24)
       </p>

  <li> <p> New: the class <code
       class="class">FE_DGPNonparametric</code> implements finite elements
       where shape functions are polynomials of order <i>k</i> on the
       actual grid cell. This is achieved by evaluating the polynomials at
       the mapped quadrature points. No grid transfer matrices are
       available for this class.
       <br>
       (GK 2002/09/19)
       </p>

  <li> <p> 
       Fixed: Some of the various instances of the <code
       class="member">VectorTools::interpolate_boundary_values</code> functions
       were not explicitly instantiated, leading to linker errors. This is now
       fixed. 
       <br>
       (WB 2002/09/19)
       </p>

  <li> <p> 
       Removed: The <code
       class="member">FiniteElement::component_to_system_index</code> function
       and its counterpart for faces is gone. This data did not make much sense
       in the case of elements that are non-zero in more than one vector
       component, such as the Nedelec element. The respective information can
       also be obtained from other sources provided by the finite element
       classes, if so necessary.
       <br>
       (WB 2002/09/16)
       </p>

  <li> <p> 
       Changed: The <code
       class="member">FiniteElement::restriction_is_additive</code> function
       used to take an argument that denoted the vector component of a finite
       element. This has become difficult with elements that are non-zero in
       more than one vector component, such as the Nedelec element. Thus, the
       semantics of the function have been changed so that the argument now
       denotes the index of a shape function, rather than a vector
       component. Since the function is probably not used in application code,
       this will most probably not lead to more serious problems.
       <br>
       (WB 2002/09/16)
       </p>

  <li> <p> 
       New: The mapping classes now also know how to transform tensors of rank
       2, not only vectors (tensors of rank 1) in a co- and contravariant way.
       <br>
       (WB 2002/09/03)
       </p>

  <li> <p> 
       Fixed: the <code>GridIn</code> class had problems
       when reading in UCD grids with comment lines that contained
       only the comment sign, but nothing else. This is now fixed.
       <br>
       (WB 2002/08/30)
       </p>

  <li> <p> 
       Improved: <code>VectorTools</code>::<code
       class="member">integrate_difference</code> can compute <i>L<sup>p</sup></i>
       and <i>W<sup>1,p</sup></i> norms for arbitrary <i>p</i>. The function
       receives an additional optional argument <i>p</i> for this. All previous
       fuctionality remains unchanged, although the code has been cleaned up a bit.
       <br>
       (GK 2002/08/01)
       </p>

  <li> <p> New: The <code>GridTools</code> class now
       offers functions to apply general transformations to all
       vertices of a triangulation, and also more specialized
       functions that shift, scale, and rotate triangulations.
       <br>
       (RH 2002/07/25)
       </p>

  <li> <p> New: The existing <code
       class="member">FETools::extrapolate</code> functions does not
       take into account hanging node constraints. Therefore, it works
       for continuous elements on globally refined grids and on
       discontinuous elements, only. Now, there is a new <code
       class="member">FETools::extrapolate</code> function with an
       additional <code>ConstraintMatrix</code> argument
       for the hanging node constraints. This new function works for
       continuous as well as for discontinuous elements. But, the old
       function is still supported.
       <br>
       (RH 2002/06/17)
       </p>

  <li> <p> New: There are now <code
       class="member">Triangulation::save/load_user_pointers(vector&lt;void
       *&gt; &amp;)</code> functions similar to the respective <code
       class="member">Triangulation::save/load_user_flags(vector&lt;bool&gt;
       &amp;)</code> functions.
       <br>
       (RH 2002/06/14)
       </p>

  <li> <p> Fixed: Bug in <code
       class="member">Triangulation<3>::load_user_flags(vector&lt;bool&gt;
       &amp;)</code> is now fixed.
       <br>
       (RH 2002/06/14)
       </p>

  <li> <p> Fixed: <code
       class="member">Triangulation::load_user_flags(vector&lt;bool&gt;
       &amp;)</code> erroneously threw an assertion for
       <code>dim==1</code>. This is now fixed.
       <br>
       (RH 2002/06/14)
       </p>

  <li> <p> Changed: <code
       class="member">Triangulation<1>::n_quads</code> now returns 0,
       instead of throwing an assertion, as before. The same holds for
       similar functions like <code
       class="member">Triangulation&lt;dim&gt;::n_hexs</code> that now
       returns 0 for <code>dim<3</code>.
       <br>
       (RH 2002/06/14)
       </p>

  <li> <p> Improved: Several functions like the <code
       class="member">DoFHandler::distribute_dofs</code> and <code
       class="member">DoFTools::make_flux_sparsity_pattern</code>
       functions altered the <code
       class="member">user_flags</code>. This was stated in the
       documentation of these functions. Nevertheless, it might have
       led to unexpected behaviour of the <code
       class="member">user_flags</code> for users who weren't aware of
       this <em>side-effect</em>. Now, these functions do not alter
       the <code>user_flags</code>, any
       more. Consequently, the users do not need to worry any more
       about the reliability of the <code
       class="member">user_flags</code> when calling any function of
       the library.
       <br>
       (RH 2002/06/14)
       </p>  

  <li> <p> Fixed: <code
       class="member">FE_Q::has_support_on_face</code> always returned
       true in 1d, partly because faces are not really an issue in
       1d. It now does so only when the support point of the
       respective shape function is actually the requested vertex. The
       same applied to <code
       class="member">FE_DGQ::has_support_on_face</code>
       <br>
       (WB 2002/06/13)
       </p>

  <li> <p> New: The existing <code
       class="member">FETools::interpolate</code>, <code
       class="member">FETools::back_interpolate</code> and <code
       class="member">FETools::interpolate_difference</code> functions
       do not take into account hanging node constraints. Therefore,
       they work for continuous elements on globally refined grids and
       on discontinuous elements, only. Now, there are new functions
       with the same names but additional <code
       class="class">ConstraintMatrix</code> arguments for the hanging
       node constraints. These new functions work for continuous as
       well as for discontinuous elements. But, the old functions are
       still supported.
       <br>
       (RH 2002/06/13)
       </p>

  <li> <p> Changed: The constructor of <code
       class="class">DoFHandler</code> now takes a reference to a
       <code>const</code> <code>Triangulation</code>.
       <br>
       (RH 2002/06/12)
       </p>

  <li> <p> Changed: The constructors of all <code
       class="class">DoFAccessor</code>, <code
       class="class">TriaAccessor</code> and <code
       class="class">TriaIterator</code> classes now take pointers to
       <tt>const</tt> <code>Triangulation</code>s.
       <br>
       (RH 2002/06/12)
       </p>

  <li> <p> Fixed: In debug mode the <code
       class="member">MappingQ1::transform_real_to_unit_cell</code>
       function erroneously threw an assertion when used in 1 or 3
       dimensions. This is now fixed.
       <br>
       (RH 2002/06/10)
       </p>

  <li> <p> 
       Fixed: The <code>get_dof_indices</code>
       functions of DoF accessor classes used to work also for
       non-active cells. However, the results were bogus except for
       the special case that we had a finite element that has all its
       degrees of freedom located in vertices. This is now fixed: the
       function throws an exception in all other cases, since there is
       no useful meaning for it then. It continues to work in the
       special case.
       <br>
       (WB 2002/06/08)
       </p>

  <li> <p> 
       New: For encapsulated postscript output of 2d grids, it is now
       possible to tell the <code>GridOut</code> class to
       write the cell numbers into each cell, as well as the numbers
       of the vertices.
       <br>
       (WB 2002/06/04)
       </p>
</ol>


*/
