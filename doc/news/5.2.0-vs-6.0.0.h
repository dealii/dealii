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
 * @page changes_between_5_2_and_6_0 Changes between Version 5.2 and 6.0

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

  <li> <p>Removed: Several functions in the linear algebra classes that have
   been deprecated for more than a year have been removed.
   <br>
   (GK 2007/08/22)
   </p>

  <li> <p>Changed: Implementing gradients for the class <code>FunctionDerivative</code>, it became evident that its enums for
  difference formulas clashed with those of <code>AutoDerivativeFunction</code>. Therefore, only the latter
  survived.
  <br>
  (GK 2007/08/02)
  </p>

  <li> <p>Changed: When new multigrid transfer classes were introduced,
  the existing class <code>MGTransferSelect</code> was
  moved to the new header file
  <code>multigrid/mg_transfer_component</code>. It also received a new
  base class <code>MGTransferComponent</code>.
  <br>
  (GK 2007/04/26)
  </p>

  <li> <p>
       Changed: Everything that is part of the deal.II library has
       been moved into a namespace <code>dealii</code>. To allow your
       code to compile, you will have to either put a <code>using
       namespace dealii;</code> declaration at the top of each file,
       or qualify each deal.II class and function explicitly with
       <code>dealii::</code> namespace access.
       <br> 
       (WB 2006/10/22)
       </p>

  <li> <p>
       Changed: Indices, such as vertex indices, are usually represented by
       unsigned integers in deal.II. One place where we didn't do this was in
       the <code>CellData</code> structure that can be used to describe cells
       when building a triangulation. This has now been rectified. 
       <br> 
       (WB 2006/09/06)
       </p>

  <li> <p>
       Changed: Lower dimensional objects have been removed from the
       hierarchical structure of levels in <code>TriaLevel</code>. Faces, i.e. lines in 2D and
       lines and quads in 3D, have no associated level
       anymore. Therefore, the level argument of some iterator
       functions have been removed. The affected functions are <code>Triangulation::begin_raw_face</code>, <code>begin_face</code>, <code>begin_active_face</code> and all <code>last_*_face</code> or <code>end_*_face</code> functions, no matter whether
       raw, used or active. Also effected are the direct iterator
       functions which are related to faces like <code>begin_line</code> in 2D and 3D or <code>begin_quad</code> in 3D. Again, the same applies
       to <code>last_*</code> and <code>end_*</code>.
       <br>
       The respective functions in <code>DoFHandler</code>, <code>hp::DoFHandler</code> and <code>MGDoFHandler</code> have been changed
       accordingly.
       <br>
       Nested loops with an outer loop over all levels and an inner
       loop over all faces on that level have to be changed to a
       single loop over all faces. In most cases, the necessary
       changes should be internal to the library.
       <br>
       (Tobias Leicht, 2006/06/13)
       </p>

  <li> <p>
       Changed: In order to facilitate the implementation of hp finite
       elements, the ordering of boundary DoFs returned by <code>DoFTools::map_dof_to_boundary_indices</code> has been
       changed. Fortunately, this is a rarely used function so that the effect
       should be limited to only a few programs.
       <br> 
       (WB 2006/04/26)
       </p>

  <li> <p>
       Changed: The <code>distribute_local_to_global</code> functions have
       been moved from the classes implementing accessors to arbitrary
       objects in the triangulation to the cell accessors, to
       facilitate the implementation of hp methods. That means
       that you can't call these functions for face or edge iterators
       any more, but only cells. Since this is what is usually
       desired, this should not be too severe a restriction.
       <br> 
       (WB 2006/04/26)
       </p>

  <li> <p>
       Changed: The template argument of the <code>InterGridMap</code> class has been changed. Code
       like e.g. <code>InterGridMap&lt;DoFHandler,2&gt;</code> must be
       replaced by <code>InterGridMap&lt;DoFHandler&lt;2&gt; &gt;</code>.
       <br> 
       (RH 2006/02/27)
       </p>

  <li> <p>
       Changed: The <code>DoFObjectAccessor</code> and derived classes
       used to take only the space dimension and the dimension of the
       object they represented as template arguments. Now the later
       argument is replaced by the type of DoFHandler they correspond
       to. This may be <code>::DoFHandler</code> or
       <code>hp::DoFHandler</code>.
       <br> 
       (RH 2006/02/10, WB 2005/12/20)
       </p>

  <li> <p>
       Removed: Support for gcc2.95 in particular, and all compilers that do
       not support the <code>std::ostringstream</code> class in general,
       has been removed.
       <br> 
       (WB 2006/02/03)
       </p>

  <li> <p>
       Changed: Several classes related to the storage of data in the
       Triangulation and DoFHandler classes have been moved into internal
       namespaces, and some have also been removed. Since these were not for
       use anywhere outside the library, it is unlikely that this poses
       problems to application programs. On the other hand, it moves common
       names like <code>Line</code> and <code>Quad</code> out of the global
       namespace.
       <br> 
       (WB 2006/01/13)
       </p>

  <li> <p>
       Changed: Previously, if <code>mg_dof_handler</code> was a
       multilevel DoF handler object, calling
       <code>DoFRenumbering::component_wise(mg_dof_handler)</code>)
       didn't quite do what was probably intended: it did an implicit
       cast of <code>mg_dof_handler</code> to its base class,
       <code>DoFHandler</code> and then renumbered the <em>global</em>
       degrees of freedom, but didn't touch the multilevel part. This
       has now been fixed: there is another function that takes a
       <code>MGDoFHandler</code> object and that renumbers the global
       part as well as all levels of the multigrid part. The change is
       incompatible since the previous call now leads to a different
       action; if you really want the old behavior back, cast the
       argument to the base class before the call, like so:
       <code>DoFRenumbering::component_wise(static_cast&lt;DoFHandler&lt;dim&gt;&&gt;(mg_dof_handler))</code>)
       <br> 
       (WB 2005/12/15)
       </p>

  <li> <p>
       Changed: The internal numbering of faces, lines and vertices
       has been reimplemented. Now the numbering scheme uses a
       lexicographic ordering (with x running fastest) wherever
       possible. For more information on the new numbering scheme, see
       the <a
       href="http://ganymed.iwr.uni-heidelberg.de/pipermail/dealii/2005/000827.html">announcement</a>
       on the mailing list.  <br> The ordering of vertices in <code>CellData</code> given to the <code>Triangulation::create_triangulation</code>
       function has been changed accordingly. For backward
       compatibility there is now a new <code>Triangulation::create_triangulation_compatibility</code>
       function which takes <code>CellData</code> in the
       old vertex ordering. However, as this function might not be
       supported forever users are advised to change their code to the
       new numbering scheme and to use the <code>Triangulation::create_triangulation</code>
       function.  
       <br> 
       (RH 2005/11/02)
       </p>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>

  <li> <p>Improved: Tests for libraries in configure can handle --without now.
  <br>
  (GK 2007/08/22)
  </p>

  <li> <p>Updated: The version of UMFPACK bundled with deal.II has been updated
       to release 5.1.0.
       <br>
       (WB 2007/08/15)
       </p>

  <li> <p>New: The step-27 tutorial program has been added. It solves a Laplace
          equation with hp finite elements 
	  and shows how to set finite element degrees, assemble matrices on
	  cells with varying polynomial degrees, and how to compute a simple
	  criterion for estimating the local Sobolev smoothness of a function.
       <br>
       (WB 2007/08/09)
       </p>

  <li> <p>
       Improved: The table of contents of tutorial programs now also shows
       subsections of introduction and results.
       <br>
       (WB 2007/08/05)
       </p>

  <li> <p>
       Extended: Up to now, in 3D only 'orientable' meshes could be used in
       deal.II, i.e. all lines are in standard orientation and the faces can be
       either in standard orientation or with a reversed normal. In order to
       enable the treatment of all purely hexahedral meshes, this interface
       (only the boolean flag <tt>face_orientation</tt> so far) has
       been extended by new flags <tt>face_flip</tt> and <tt>face_rotation</tt>,
       denoting a face which is rotated against the 'standard face' by 180 and 90
       degrees, respectively. Furthermore, there is a new flag
       <tt>line_orientation</tt> with obvious meaning. 
		 <br> These flags have to be respected during creation and refinement of a
		 triangulation, when projecting quadrature points to faces or when the dof
		 indices on a cell are extracted. Furthermore, asking for vertices and
		 lines of cells is in some cases a bit more complicated. The same applies,
       for example, 
		 to the extraction of the information, which child of a neighbor is behind
		 a given subface. However, this information is supplied by various
		 functions in <code>GeometryInfo</code>. As a rule-of-thumb:
		 if you want to use non-standard meshes, all occurrences of
		 <tt>face_orientation</tt> have to be supplemented by <tt>face_flip</tt>
		 and <tt>face_rotation</tt>.
       <br> In order to reduce the impact of possible bugs, the grid is still given to
       <code>GridReordering</code>, which now returns the original
       connectivity if the reordering was not possible.
       <br> A more detailed document concerning these topics will be available
       in the future.
		 <br>
       (Tobias Leicht 2007/02/09)
       </p>

  <li> <p>
       New: a program <tt>reconfigure</tt> has been added to the
       deal.II main directory, which reruns <tt>configure</tt> with
       the sam command line arguments as last time.
       <br>
       (GK 2006/10/24)
       </p>

  <li> <p>New: The step-25 tutorial program, demonstrating the
       solution of the nonlinear sine-Gordon equation (a variant of
       the time dependent wave equation) has been added.
       <br>
       (WB 2006/11/05)
       </p>

  <li> <p>New: The step-21 tutorial program is ready. It considers time
       dependent two-phase flow through a porous medium.
       <br>
       (WB 2006/10/30)
       </p>

  <li> <p>New: The step-23/step-24 tutorial programs have been added,
       discussing how to solve the wave equation.
       <br>
       (WB 2006/10/03)
       </p>

  <li> <p>Fixed: when building the library in parallel, using <code>make
        -jN</code> (most often with rather large values of <code>N</code>, for
	example 8), it could sometimes happen that the same file is compiled
	twice at exactly the same time. This leads to invalid object files and
	libraries that would contain illegal instructions. This is now fixed.
	<br>
	(WB 2006/08/09)
	</p> 

  <li> <p> Changed: The default number of subdivisions in the <code>build_patches</code> functions of the DataOut classes is
       now part of the parameters read by <code>parse_parameters</code>. This default value is used
       whenever patches with zero subdivisions are being built.
       <br>
       (GK 2006/06/18)
       </p>

  <li> <p>New: a configuration parameter <tt>--with-boost</tt> has been
	introduced, allowing the selection of a previously installed boost
	library instead of the contibuted.
	<br>
	(GK 2006/03/23)
	</p> 

  <li> <p>New: The step-20 tutorial program was added. It shows how to use
       Raviart-Thomas elements to solve a mixed formulation of the Laplace
       equation. It also demonstrates the use of block matrices and vectors,
       and how they can be used to define more complex solvers and
       preconditioners working on the Schur complement.
       <br>
       (WB 2006/02/12)
       </p>

  <li> <p>
       Updated: The <code>deal.II/contrib/boost</code> library has
       been updated to boost version 1.33.1.
       <br>
       (RH 2006/02/08)
       </p>

  <li> <p>
       Improved: Documentation has been greatly enhanced, both in the
       API documentation as well as in the tutorial programs. In
       particular, the step-1 through step-8 tutorial programs have
       been overhauled and in many cases the documentation has been
       rewritten in large parts. The presentation of tutorial programs
       and reference manual should also be nicer and better readable
       in many cases.
       <br>
       (WB 2006/02/07)
       </p>

  <li> <p>
       Changed: Previously, the UMFPACK code was linked into its own
       library, <code>liblac_umfpack.*</code>. Instead, it is now
       directly linked into the <code>liblac.*</code> libraries. This
       makes further linking simpler than before and also simplified
       the link process.
       <br>
       (WB 2006/01/19)
       </p>

  <li> <p>
       Improved: The order of quadrature points is now x fastest, z slowest, as
       well as the vertices of cells. Also, the nodes and subcells of patches
       in DataOutBase follow the same ordering. Output code could be simplified
       a lot, saving 900 lines of code.
       <br>
       (GK 2006/01/18)
       </p>

  <li> <p>
       Improved: Documentation of tutorial program is now run through
       Doxygen. This allows cross-referencing the tutorial programs to
       the API manual (and the other way round) and leads to a
       generally much nicer output. In particular, formulas embedded
       in the documentation of the programs are now properly displayed.
       <br> 
       (WB 2006/01/16)
       </p>

  <li> <p>
       Improved: The link to the manual now points directly to
       Doxygen's module list, and almost all classes are now grouped
       into modules that capture related functionality.
       <br> 
       (WB 2006/01/16)
       </p>

  <li> <p>
       Improved: If you do have a PETSc installation
       and have set the <code>PETSC_DIR</code> and
       <code>PETSC_ARCH</code> environment variables but do not wish
       deal.II to be configured for PETSc use, you
       should specify <code>--with-petsc=no</code> as a flag during
       configuration. 
       <br> 
       (WB 2006/01/16)
       </p>

  <li> <p>
       Improved: In all directories, including those for tutorial
       programs, we generate a file called <code>Makefile.dep</code>
       that lists the dependencies of object files on source files,
       and that the <code>make</code> uses to determine which files to
       recompile. If the generation of this file failed (for example
       if <code>LD_LIBRARY_PATH</code> hadn't been set correctly), a
       number of very strange things would happen, among others the
       attempt to link an empty object file. This was almost
       impossible to figure out if you didn't know what was going
       on. This mechanism has now been robustified and should yield
       better error messages.
       <br> 
       (WB 2005/11/23)
       </p>

  <li> <p>
       New: The dynamic libraries, which is the name of shared libs
       under Apples OSX, are now supported on these platforms and
       enabled by default. They should safe a lot of harddisk space.
       <br> 
       (Oliver Kayser-Herold 2005/10/20)
       </p>

  <li> <p>
       Fixed: The <tt>Tecplot</tt> library (<tt>tecio.a</tt>) was
       detected but not added to <tt>LIBS</tt>. This is now fixed.
       <br> 
       (RH 2005/09/23)
       </p>

  <li> <p>
       New: <tt>configure</tt> will automatically detect a
       <tt>NetCDF</tt> installation, when its path is given by the
       <tt>NETCDF_DIR</tt> environment variable. The path of the
       <tt>NetCDF</tt> installation directory can also be specified
       through the <tt>--with-netcdf=/path/to/netcdf</tt> configure
       option.
       <br> 
       (RH 2005/09/23)
       </p>

  <li> <p>
       Fixed: The PETSc libraries have been relinked each time
       <tt>make</tt> was called. Now, the PETSc libraries will only be
       relinked if necessary.
       <br> 
       (RH 2005/09/15)
       </p>

</ol>



<a name="base"></a>
<h3>base</h3>

<ol>

  <li> <p> New: There is a new class:
       <code>Functions::FEFieldFunction</code> which is a Function
       interface to a finite element solution. 
       <br> 
       (Luca Heltai 2007/08/29)

  
  <li> <p> Improved: <code>FunctionDerivative</code> is now
  derived from <code>AutoDerivativeFunction</code> and implements
  gradients as well, giving you automatic second derivatives of a function.
  <br>
  (GK 2007/08/02)
  </p>

  <li> <p> New: The function <code>Utilities::fixed_power&lt;n&gt;(q)</code>
       calculates <code>q</code> to the power of <code>n</code> where
       <code>n</code> is a constant known at compile time. It allows to
       calculate powers efficiently at compile time, most often things like a
       number to the power <code>dim</code>.
       <br> 
       (WB 2007/06/23)
       </p>

  <li> <p> New: The deal.II intermediate format has been
       changed. Since files in this format are only meant to be
       processed by exactly the same version of deal.II, this should
       be of no concern to users. However, this restriction was
       previously undocumented, even if implied. The documentation for
       this has now been updated.
       <br> 
       (WB 2007/02/08)
       </p>

  <li> <p> New: The new <code>Functions::Monomial</code> class implements
       monomials as a function object.
       <br> 
       (WB 2006/12/15)
       </p>

  <li> <p> Fixed: If no substring is found with a width smaller than the given
       threshold, the <code>Utilities::break_text_into_lines</code> function now
       returns the smallest possible substring (larger than the threshold).
       <br> 
       (Tobias Leicht 2006/12/13)
       </p>

  <li> <p> Changed: We used to set preprocessor variables
       <code>PACKAGE_NAME</code>, <code>PACKAGE_TARNAME</code>,
       <code>PACKAGE_VERSION</code>, <code>PACKAGE_STRING</code>,
       <code>PACKAGE_BUGREPORT</code> in the file
       <code>&lt;base/config.h&gt;</code>. These variables were automatically
       set by the autoconf package, with which we generate the
       <code>./configure</code> script. However, these variables would conflict
       if someone were to use a different package that also uses autoconf,
       since that package would have the same variables, just set to different
       values. The preprocessor warns about these cases.
       <br>
       We now avoid this conflict by not using these names any more. Instead,
       we use the preprocessor names as above, prefixed with
       <code>DEAL_II_</code>.
       <br>
       (WB 2006/11/10)
       </p>

  <li> <p> Extended: The <code>QGaussLobatto</code> quadrature
       rule computes Legendre-Gauss-Lobatto nodes and quadrature weights.
       <br>
       (Florian Prill 2006/11/02)
       </p>

  <li> <p> Extended: The <code>contract</code> function family contracts two
       tensors. There is now a new version of this function, that contracts a
       tensor of rank three with a second one of rank two over given indices
       <tt>index1</tt> and <tt>index2</tt> of the first and second tensor,
       respectively, resulting in a tensor of rank three. 
       <br>
       (Tobias Leicht 2006/09/27)
       </p>

  <li> <p> Fixed: The
       <code>DataOutBase::write_deal_II_intermediate</code> function
       did not check whether it can actually write to the stream
       given. It would therefore silently do nothing if the file could
       not be written to. This is now fixed: an exception is generated
       in that case, as for all the other output functions in that class.
       <br>
       (WB 2006/08/31)
       </p>

  <li> <p> Improved: Stack backtraces in multithreaded mode were often
       very long and almost unreadable because the functions in
       namespace <code>Threads</code> that are used to set up new
       threads have long and awkward signatures. The output is now
       filtered to make these backtraces easier to read.
       <br>
       (WB 2006/08/15)
       </p>

  <li> <p> New: The second argument to <code>Utilities::int_to_string</code>
       can now be omitted, leading to a string that isn't zero padded at all.
       <br>
       (WB 2006/08/02)
       </p>

  <li> <p> Changed: When there is still a <code>SmartPointer</code> object
       pointing to another object at the time it is destroyed, this would cause
       the program to be aborted. However, there are cases where this is not
       desirable, for example here:
       <pre>
       <code>
          void f() 
          {
	    Triangulation tria;
	    DoFHandler *dh = new DoFHandler(tria);
	    ...some function that throws an exception
          }
       </code>
       </pre>
       When the exception is thrown but not caught, the two local objects are
       destroyed in reverse order of their construction, i.e. first the pointer
       then the triangulation. However, only the pointer, not the
       <code>DoFHandler</code> pointed to is destroyed, triggering the abort in
       the triangulation since there is still the <code>DoFHandler</code>
       object pointing to it at the time of destruction of the
       triangulation. In such cases, one would not want to see the program
       aborted, since then one would never learn about the actual exception
       being thrown.
       <br>
       The behavior of the <code>Subscriptor</code> class as therefore been
       changed to not abort the program any more if an exception is being
       handled at the moment. Rather, only an error message is shown on
       <code>std::cerr</code>.
       <br>
       (WB 2006/08/01)
       </p>

  <li> <p> Fixed: The <code>TableHandler::write_tex</code>
       accidentally took a parameter of type
       <code>std::ofstream</code> instead of <code>std::ostream</code>
       and could therefore only write to files, not to other kinds of
       streams. This is now fixed.
       <br>
       (WB 2006/07/28)
       </p>

  <li> <p> New: <code>GeometryInfo</code> offers several new
       functions, <code>is_inside_unit_cell</code> with an epsilon parameter to
       specify allowable offsets from the actual unit cell,
       <code>distance_to_unit_cell</code> returning the infinity norm of the
       distance of a given point to the unit cell, and
       <code>project_to_unit_cell</code> returning the projection of a point
       onto the unit cell. Also, a new member <code>vertex_to_face</code> allow
       to determine to which faces of a cell a vertex belongs. 
       <br>
       (Ralf B. Schulz 2006/05/10)
       </p>

  <li> <p> Improved: <code>DataOutBase</code>::<code>OutputFormat</code> has a new value <tt>none</tt>,
       writing no output at all. This way, the writing of output files can be
       controlled more easily from parameter files.
       <br>
       (GK 2006/04/14)
       </p>

  <li> <p> Improved: <code>VectorSlice</code> has new functions
       <code>begin()</code> and <code>end()</code>,
       returning the corresponding vector iterators.
       <br>
       (GK 2006/03/31)
       </p>

  <li> <p> New: The various tensor classes can now effectively be reset to zero
       by simply writing <code>t=0;</code> as has long been allowed for
       matrices and vectors.
       <br>
       (WB 2006/03/24)
       </p>

  <li> <p> New: The new <code>deal_II_numbers::is_finite</code> function can
       be used to check whether a floating point number is finite,
       i.e. neither plus or minus infinite nor NaN (not a number);
       it is in the new include file <tt>base/numbers.h</tt>, which will
       contain non-configured inline functions of this type.
       <br>
       (WB 2006/03/15) (GK 2006/03/22)
       </p>

  <li> <p> Improved: If the C++ runtime environment allows it, we now
       demangle stacktraces generated whenever an exception is thrown.
       <br>
       (WB 2006/03/14)
       </p>

  <li> <p> Improved: The function <code>Subscriptor</code>::<code>list_subscribers</code> logs all current subscribers of an
       object to <tt>deallog</tt>.
       <br>
       (GK 2006/03/08)
       </p>

  <li> <p>Fixed: Writing a denormal "NaN" through <code>LogStream</code> objects such as 
       <code>deallog</code> erroneously printed zero, rather
       than "nan". This is now fixed.
       <br>
       (Luca Heltai, WB 2006/03/07)
       </p>

  <li> <p>Improved: The <code>TableBase</code> base class of all
       the <code>Table</code> classes had an
       <code>>operator()</code> that takes a <code>TableIndices</code> object
       that represents the entire set of indices at once. However, due to C++
       name lookup rules, this operator wasn't accessible through the
       <code>Table</code> class without explicitly specifying the base class
       name. This is now fixed by also providing a similar
       <code>>operator()</code> in the derived classes.
       <br>
       (WB 2006/02/21)
       </p>

  <li> <p>Fixed: The <code>QGauss</code> constructor
       hangs on x86 linux systems when deal.II is run from inside
       MATLAB. The reason is that while the processor offers long
       double accuracy, MATLAB switches off support for that inside
       the CPU. It therefore leaves codes that expect the additional
       accuracy high and dry. Annoyingly, these codes of course run
       just fine outside of MATLAB's environment. This behavior leads
       to an infinite loop inside the QGauss constructor, where we
       want to use the additional accuracy to compute quadrature
       points to high precision. To avoid the infinite loop, we now
       check during runtime whether the extra precision is available,
       not only at compile time.
       <br>
       (Christian Kamm, WB 2006/02/21)
       </p>

  <li> <p>Fixed: The <code>ParameterHandler::get_integer</code>
       had an erroneous check that the value given for a parameter really
       represented an integer. This check always succeeded, even in face of an
       error, and a zero value was returned every time the value in the
       parameter file happened to be not an integer.
       <br>
       (WB 2006/02/17)
       </p>

  <li> <p>Improved: The <code>ComponentSelect</code> class can
       now also handle the case that one wants to select multiple vector
       components at once, not only a single one.
       <br>
       (WB 2006/02/12)
       </p>

  <li> <p>Improved: A new function <code>TableBase</code>::<code>fill</code>, filling the whole table with the same
       element has been introduced.
       <br>
       (GK 2006/01/18)
       </p>

  <li> <p>Improved: <code>DataOutBase</code> now writes binary
       files for OpenDX.
       <br>
       (GK 2006/01/18)
       </p>

  <li> <p>
       New: There are now functions
       <code>ParameterHandler::set</code> that allow to set the value of a
       parameter to something different later on.
       <br>
       (Marc Secanell, WB 2006/1/2)
       </p>

  <li> <p>
       New: There are new functions
       <code>Utilities::match_at_string_start</code> and 
       <code>Utilities::get_integer_at_position</code>.
       <br>
       (WB 2005/12/19)
       </p>

  <li> <p>
       Fixed: The computation of <code>HMIN</code> and
       <code>HMAX</code> in <code>DataOutBase::write_povray</code> has been
       wrong. This is now fixed.
       <br>
       (RH 2005/11/02)
       </p>

  <li> <p>
       Fixed: <code>DataOutBase&lt;2,3&gt;::write_tecplot_binary</code>
       did not write the <code>z</code> coordinates. This is now
       fixed.
       <br>
       (RH 2005/09/28)
       </p>

  <li> <p>
       Fixed: The <tt>tecplot_binary</tt> <code>OutputFormat</code> has been added to <code>DataOutBase::get_output_format_names</code>. Now
       an exception will be raised if <code>write_tecplot_binary</code> is invoked without
       specifying the filename through the <code>DataOutBase::TecplotFlags</code> interface.
       <br> 
       (RH 2005/09/23)
       </p>

  <li> <p>
       New: There are now new <code>get_n_mpi_processes</code> and <code>get_this_mpi_process</code> functions in the
       <code>Utilities::System</code> namespace. In case
       of a code not running in parallel, they simply return 1 and 0,
       respectively.
       <br> 
       (RH 2005/09/20)
       </p>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li> <p>Changed: The <code>CompressedBlockSparsityPattern</code> has been renamed to
       <code>BlockCompressedSparsityPattern</code> to be
       consistent, since the "block" part builds on the "compressed sparsity
       pattern", not the other way around. The old name remains as a typedef,
       but its use is deprecated.
       <br>
       (WB 2007/06/17)
       </p>

  <li> <p>New: The <code>CompressedSetSparsityPattern</code>
       class is an alternative to the <code>CompressedSparsityPattern</code> class that appears to be
       better suited for problems that have many entries per row. There is also
       a block version, <code>BlockCompressedSetSparsityPattern</code>. The two new
       classes can used in all places where a regular compressed sparsity
       pattern can also be used.
       <br>
       (Oliver Kayser-Herold, WB 2007/06/03)
       </p>

  <li> <p>Fixed: The <code>SolverBicgstab</code> class
       did not handle hitting on the solution early very
       gracefully (it threw an exception). This is now fixed.
       <br>
       (Roger Young 2007/03/07)
       </p>

  <li> <p>Fixed: The <code>SparseDirectMA27</code> class allows
       to run the sparse direct solver as a separate program (the
       <code>detached_ma27</code> program) rather than as part of the current
       program in order to be able to run several instances of it in parallel
       during multithreaded computations. However, the destructor of the <code>SparseDirectMA27</code> class had a bug that when using
       this detached mode led to a segmentation fault. This is now fixed.
       <br>
       (WB 2007/02/09)
       </p>

  <li> <p>Improved: A simple <code>print_formatted</code> 
       function has been added to <code>SparseMatrixEZ</code>.
       <br>
       (Moritz Allmaras 2007/02/08)
       </p>

  <li> <p>Fixed: The <code>SparseDirectMA27</code> 
       class works for symmetric matrices and had a check to make sure the
       matrix is indeed symmetric. However, the check compared entries for
       equality, something one should not do for floating point numbers. It now
       checks that symmetric entries are equal up to a small relative number.
       <br>
       (WB 2006/12/21)
       </p>

  <li> <p>Fixed: The <code>FullMatrix::invert</code> 
       function would return wrong results if one tried to invert a
       matrix of dimension smaller than 5 in situ. This is now fixed.
       <br>
       (Florian Prill 2006/12/18)
       </p>

  <li> <p>Improved: The <code>SparsityPattern</code> class would produce
       segmentation faults if one tried to allocate a matrix with more
       than about 4.2 billion elements (i.e. the number that one can
       store in a 32-bit unsigned integer). This is now fixed: if you
       have enough memory, this is now possible on 64-bit
       systems. (The number of rows is still limited by the 32-bit
       limit, but this is probably going to be enough for some time to
       come.) This fix also required changing the return type of the
       <code>SparsityPattern::get_rowstart_indices</code> function
       from <code>const unsigned int *</code> to <code>const
       std::size_t *</code>; since this function should not be used
       outside of the library anyway, this is probably not a change
       that affects user code.
       <br>
       (WB 2006/12/14)
       </p>

  <li> <p>New: The class <code>PointerMatrixVector</code> implements the functions
       <code>vmult</code> and <code>Tvmult</code> for a <i>1xn</i>-matrix represented by
       a single vector.
       <br>
       (GK 2006/09/17)
       </p>

  <li> <p>Improved: The class <code>SolverCG</code> can
       now estimate the condition number of the linear system using <code>TridiagonalMatrix</code> and LAPACK.
       <br>
       (GK 2006/09/06)
       </p>  

  <li> <p>New: the class <code>TridiagonalMatrix</code>
       has been introduced together with its basic methods, including an
       eigenvalue solver for the symmetric case (LAPACK).
       <br>
       (GK 2006/09/06)
       </p>

  <li> <p> Improved: PETSc 2.3.1 (<i>build 16</i>) is now supported by the
       linear preconditioner classes. The new PETSc functions 
       <code>PCFactorSetXXX</code> are used.<br>       
       (Florian Prill, 2006/08/04)
       </p>

  <li> <p>New: The class <code>TransposeMatrix</code>
       modeled after <code>PointerMatrix</code> swaps the
       <code>vmult</code> functions such that its effect is
       the transpose of the matrix it points to.
       <br>
       (GK 2006/07/07)
       </p>  

  <li> <p>New: There is now a function <code>FullMatrix::trace</code> that does what its name suggests
       it does.
       <br>
       (WB 2006/06/16)
       </p>
  
  <li> <p>Improved: <code>PointerMatrixAux</code> now has a
       default constructor and a function for setting the <code>VectorMemory</code> object.
       <br>
       (GK 2006/06/14)
       </p>
  
  <li> <p>Fixed: <code>FullMatrix::print</code> would yield a link error
       when used with <code>std::ofstream</code> since an explicit
       instantiation was missing. The function has now been moved to
       the header file and made inline so that it is always visible,
       whatever the type of the output stream is.
       <br> 
       (WB 2006/06/07)
       </p>
  
  <li> <p> Improved: The <code>SparseDirectUMFPACK</code> solver can now also be
       used with sparse matrices with elements of type float, as well
       as with block matrices with types float and double.
       <br> 
       (WB 2006/04/25)
       </p>

  <li> <p> New: The function <code>BlockSparsityPattern::row_length</code> adds up
       the row lengths of the individual blocks of a block matrix for
       a given row.
       <br> 
       (WB 2006/04/25)
       </p>
  
  <li> <p> New: There is a new class <code>IdentityMatrix</code> that represents an
       efficient version of the matrix with ones on the diagonal and
       zeros everywhere else. The main usefulness of this matrix lies
       in the fact that most other important matrix classes have
       assignment operators and copy constructors that allow to
       conveniently fill a (sparse, full) matrix with the identity
       matrix.
       <br> 
       (WB 2006/04/24)
       </p>
  
  <li> <p> New: There are now assignment operators from <code>BlockVector</code> to <code>Vector</code> and back.
       <br>
       (WB 2006/03/30)
       </p>

  <li> <p> Improved: <code>BlockSparsityPattern</code> can be
       initialized directly using the vector generated by
       <code>DoFTools</code>::<code>compute_row_length_vector</code>.
       <br>
       (GK 2006/03/30)
       </p>

  <li> <p>
       Improved: All matrix (and some vector) classes now check whether
       entries written into them represent finite floating point
       values. This should catch some bugs earlier where one writes
       infinite or NaN values into a matrix only to realize later that
       the linear solver fails.
       <br>
       (Stephan Kramer, WB 2006/03/15)
       </p>  

  <li> <p>
       Changed: There are new <code>FullMatrix</code>::<code>equ</code> functions which assign <tt>this</tt>
       matrix to the linear combination of one, two or three
       matrices. Also there is now a new <code>Vector</code>::<code>equ</code>
       function for three vectors.
       <br>
       (RH 2006/02/21)
       </p>  

  <li> <p>
       Fixed: The <code>SolverMinRes</code> class did not work
       with the <code>BlockVector</code> class.
       This is now fixed.
       <br>
       (Luca Heltai 2006/02/19)
       </p>  

  <li> <p>
       Changed: There are now new <code>FullMatrix</code>::<code>add</code> functions which add two and three
       scaled matrices.
       <br>
       (RH 2006/02/16)
       </p>  

  <li> <p>
       Changed: The matrix classes <code>FullMatrix</code>, <code>SparseMatrix</code>, <code>SparseMatrixEZ</code> and <code>BlockMatrixBase</code> now have an <code>add</code> function analogous to the <code>add</code> of vector classes. The old <code>add_scaled</code> functions are now deprecated.
       <br>
       (RH 2006/02/16)
       </p>  

  <li> <p>       
       Improved: <code>BlockMatrixArray</code>::<code>enter_aux</code> allows using matrices without
       adding vector multiplication by using <code>PointerMatrixAux</code>.
       <br>
       (GK 2006/02/02)
       </p>  

  <li> <p>
       New: The class <code>PointerMatrixAux</code> was
       introduced for use with matrix classes lacking the adding
       vector multiplication functions. It implements these by using
       its own auxiliary vector.
       <br>
       (GK 2006/02/02)
       </p>

  <li> <p>
       Improved: The <code>FilteredMatrix</code> class
       was able to filter only the <code>SparseMatrix</code>
       matrix class.
       A new (faster) version is now in its place with the same interface 
       that is able to perform such a filter on anything that provides 
       the usual matrix vector operations.
       <br> 
       (Luca Heltai 2006/01/09)
       </p>

  <li> <p>
       New: All solver classes took an argument of type <code>VectorMemory</code> which they use to allocate
       memory for temporary vectors. A clever implementation of this
       class might allow reusing temporary vectors, and thus reduce
       the overhead of repeatedly allocating and freeing
       memory. However, almost all instances of using these classes
       used the <code>PrimitiveVectorMemory</code>
       class. The solver class have now been changed so that if one
       omits the respective argument, they fall back to using such a
       memory object. Also, all example programs that did not
       specifically use a different memory allocation class, have been
       changed to not specify anything at all, and thus fall back to
       the default.
       <br> 
       (WB 2005/11/23)
       </p>

  <li> <p>
       Fixed: The <code>SparseMatrix</code> iterators had a
       problem when one wrote code like <code>iterator->value()=0</code>
       (i.e. with a zero integer, rather than a floating point number on the
       right), in that that opened up a second conversion sequence and the
       compiler complained about an ambiguity. This is now fixed.
       <br> 
       (WB 2005/10/17)
       </p>

  <li> <p>
       New: Now the class <code>SparseILU</code> supports
       also the <code>Tvmult</code> method. This feature
       allows the straightforward use of <code>SparseILU</code>
       as a multigrid smoother within the framework.
       <br> 
       (Oliver Kayser-Herold 2005/09/15)
       </p>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li> <p> New: There is a new class:
       <code>Functions::ParsedFunction</code> which is friendly
       wrapper to the <code>FunctionParser</code> class.
       <br> 
       (Luca Heltai 2007/08/29)
  
   <li> <p>Fixed: the function
       <code>DataOut::build_patches</code> 
       had a quadratic algorithm when generatic cell-data (as opposed
       to DoF data). This algorithm became a bottleneck when
       generating output on meshes with large number of cells. This is
       now fixed.
       <br>
       (WB 2007/08/28)
       </p>

   <li> <p>New: the function
       <code>DoFTools::get_active_fe_indices</code> 
       extracts for each cell the active finite element index used on it.
       <br>
       (WB 2007/08/07)
       </p>

   <li> <p>Extended: the function
       <code>Triangulation::copy_triangulation</code> 
       copies all members of a triangulation except for the list of
       <code>RefinementListener</code>s. In most cases this is exactly the
       intended behavior. However, if a RefinementListener should be copied to
       the new triangulation, e.g. if a Persistent Triangulation is created from
       an ordinary one, it can do so now through implementing the new
       function <code>RefinementListener::copy_notification</code>.
       <br>
       (Tobias Leicht 2007/06/05)
       </p>

   <li> <p>New: the function
       <code>DoFTools::make_sparsity_pattern</code> 
       now takes an optional constraint matrix argument that can be used to
       avoid the subsequent call to condense the sparsity pattern. In effect,
       the sparsity pattern is generated condensed right away. For problems in
       3d and with many constraints, this is many times faster than separate
       creation and condensation.
       <br>
       (WB 2007/06/03)
       </p>

   <li> <p>New: the new function
       <code>ConstraintMatrix::add_entries_local_to_global</code> 
       can be used to add entries to a matrix where entries that will appear if
       the current entry corresponds to a constrained degree of freedom are
       also added.
       <br>
       (WB 2007/06/03)
       </p>

   <li> <p>Fixed: the function
       <code>GridTools::find_cells_adjacent_to_vertex</code> 
       was not detecting properly the coarse cells adjacent to
       refined cells.
       <br>
       (Luca Heltai 2007/05/15)
       </p>

   <li> <p>Fixed: the two tools
       <code>DoFTools::count_dofs_per_component</code> and
       <code>DoFTools::count_dofs_per_block</code> where changing the
       size of the destination vector. Consistently with (most of) the 
       rest of the library, now the vectors are expected to be the 
       right size before calling these functions. 
       <br>
       (Luca Heltai 2007/05/15)
       </p>

  <li> <p>New: The classes <code>MGTransferBlockSelect</code> and <code>MGTransferBlock</code> allow for transfer of multigrid
  vectors for single blocks of a system and for several blocks.
  <br>
  (GK 2007/04/26)
  </p>

   <li> <p>New: There is a new variant of
       <code>DoFTools::make_sparsity_pattern</code> that can be used to
       construct sparsity patterns corresponding to problems where one would
       like to test shape functions defined on one mesh with shape functions
       defined on a different mesh (or on the same mesh but using a different
       DoFHandler that may use a different finite element, for example).
       <br>
       (Yaqi Wang 2007/03/09)
       </p>

   <li> <p>New: The
       <code>GridTools::get_active_child_cells</code> function determines all
       active children of a cell.
       <br>
       (Yaqi Wang 2007/03/09)
       </p>

   <li> <p>New: The
       <code>ConstraintMatrix::write_dot</code> function can be used
       to produce a graphical representation of the graph of
       constraints stored in the constraint matrix. The output can be
       sent through the "dot" program to produce a number of graphical
       formats, such as postscript, png, or xfig.
       <br>
       (WB 2007/03/08)
       </p>

   <li> <p>Fixed: The
       <code>GridGenerator::cylinder</code> function in 3d properly
       set the boundary indicators of the end faces of the generated
       cylinder to a different value than the outer faces, but it
       forgot to also adjust the boundary indicators of the edges
       purely in the interior to the ends of the cylinder. This is now
       fixed.
       <br>
       (WB 2007/03/08)
       </p>

   <li> <p>Improved: The
       <code>CylinderBoundary</code> class can now describe the
       boundaries of arbitrarily oriented cylinders, not only does
       oriented parallel to the axes and going through the origin.
       <br>
       (WB 2007/03/08)
       </p>

   <li> <p>Improved: The
       <code>GridRefinement::refine_and_coarsen_fixed_number</code> and
       <code>GridRefinement::refine_and_coarsen_fixed_fraction</code> functions
       have gained an additional last argument that can be used to specify a
       maximal number of cells that we would like to use in a
       triangulation. Its default value is set to indicate that no limit is
       desired, as is the previous behavior.
       <br>
       (WB 2007/02/20)
       </p>

   <li> <p>New: Added function <code>GridGenerator</code>::<code >hyper_cube_with_cylindrical_hole</code> that produces 
	a square with a circular hole in the middle in 2d, and extrudes it
	along the z-direction between 0 and L. 
        <br>
        (Luca Heltai 2007/02/15)
        </p>

   <li> <p>Workaround: The class <code>GridOut</code>::<code >write_msh</code> produces a mesh which can be visualized
	in the Gmsh reader. A bug in Gmsh prevented the boundary indicator to
	be properly visualized. The boundary indicator was added to the
	material flag of the faces (which is ignored when reading back the
	mesh) in order for the Gmsh reader to be able to display the boundary
	indicator in a proper way. 
        <br>
        (Luca Heltai 2007/02/15)
        </p>


   <li> <p>Fixed: The function <code>DoFTools</code>::<code >distribute_cell_to_dof_vector</code> produced
       wrong results if the vector into which the result is to be
       stored contained nonzero values. This is now fixed.
       <br>
       (Rohallah Tavakoli, WB 2007/02/13)
       </p>

   <li> <p>Fixed: A local variable in
       <code>TriaAccessor&lt;3,3&gt;::measure</code> was
       erroneously marked as static could lead to
       wrong results and crashes in multithreaded mode. This is now fixed.
       <br>
       (WB 2007/02/09)
       </p>

   <li> <p>Fixed: <code>MatrixCreator</code>::<code >create_mass_matrix</code> and <code>
       create_laplace_matrix</code> computed wrong values for the right 
       hand sides. This has been fixed. 
       <br>
       (Moritz Allmaras 2007/02/08)
       </p>

   <li> <p>Extended: <code>DataOutBase::Patch</code> has been extended by
       a new boolean flag <tt>points_are_available</tt>, which
       defaults to <tt>false</tt>. It is set to <tt>true</tt> if the
       coordinates of the points defining the subdivision of a patch
       are appended to the <tt>data</tt> table contained in a
       <code>Patch</code>. This way, <code>DataOut</code>::<code>build_patches()</code> can use a <code>Mapping</code> to represent curved boundaries,
       especially for higher order elements. This change corresponds
       to an extension of the intermediate format for graphics.
       <br>
       Fixed: Using the given <code>Mapping</code> to
       obtain function values in <code>DataOut</code>::<code>build_patches()</code> also fixes a bug for
       <code>FE_RaviartThomas</code> and <code>FE_ABF</code>elements, which need to evaluate the
       function values on the real (mapped) cell.
       <br>
       (Tobias Leicht 2007/01/30)
       </p>

  <li> <p> Fixed: On faces with wrong <code>face_orientation</code> the dofs
       have to reordered in order to be combined with the correct shape
       functions. This is only relevant for continuous elements in 3D. At least for
       <code>FE_Q</code> and systems of <code>FE_Q</code> this works now, for other finite elements the
       reordering vector still has to be implemented.
       <br>
       (Tobias Leicht, 2007/01/17)
       </p>

  <li> <p>
       Fixed: The <code>Triangulation</code>::<code>execute_coarsening_and_refinement</code> function has to
       discard coarsening flags in 2d and 3d if a neighbor of a flagged cell is
       refined or will be refined, in order to avoid that we end up with
       neighboring cells that differ in refinement by two levels. The function
       was overly conservative, however, in that it didn't allow a cell to be
       coarsened if its neighbor is once refined but is also marked for
       coarsening. This is now fixed and will lead to a few more cells being
       coarsened.
       <br>
       (Yaqi Wang, WB 2006/12/27)
       </p>

  <li> <p>
       Fixed: The <code>Triangulation</code>::<code>
       MeshSmoothing</code> flag <code>patch_level_1</code> wrongly produced cells
       on level 0 during <code>coarsening_and_refinement</code>, which is now fixed.
       <br>
       (RH 2006/12/20)
       </p>

  <li> <p>
       Extended: There is now a new <code>Triangulation</code>
       ::<code>MeshSmoothing</code> flag <code>coarsest_level_1</code>
       which ensures that after <code>coarsening_and_refinement</code> there are no
       cells on level 0, i.e. the coarsest level is 1, if the triangulation had
       <code>coarsest_level_1</code> already before.
       <br>
       (RH 2006/12/20)
       </p>

  <li> <p>
       Extended: The <code>GridIn</code> class can now
       read in tecplot files in ASCII format (block and point format,
       ordered and unstructured grids, format specifiers acccording to
       Tecplot 10 and younger versions). At the moment the
       implementation is restricted to 2d grids but can easily be
       extended to 3d as well.
       <br>
       (Tobias Leicht 2006/12/14)
       </p>

  <li> <p>
       Extended: So far, the <code>GridReordering</code>::<code>invert_all_cells_of_negative_grid</code>
       function did nothing in 2d. Now, it inverts cells from
       clockwise to counterclockwise sense (in the old numbering
       scheme).
       <br>
       (Tobias Leicht 2006/12/13)
       </p>

  <li> <p>
       New: There is now a function <code>GridTools</code>::<code>delete_duplicated_vertices</code> that deletes
       duplicate vertices which can occur if structured grids are read
       into deal.II, leading to unwanted interior
       boundaries. In order to reduce the effort of the quadratic
       algorithm, a list of vertices to be considered can be supplied
       if available.
       <br>
       (Tobias Leicht 2006/12/13)
       </p>

  <li> <p>
       New: There are now two new functions <code>GridGenerator</code>::<code>subdivided_hyper_rectangle</code> that produces
       a non-uniformly subdivided rectangle, ideally suited for graded
       meshes. One of these functions is able to create meshes with holes. 
       <br>
       (Yaqi Wang 2006/11/15, 2006/12/29)
       </p>

  <li> <p> Fixed: Corrected <code>clone</code> method 
       of <code>FE_DGQ</code> class for non-equidistant 
       support points.
       <br>
       (Florian Prill 2006/10/31)
       </p>

  <li> <p> Improved: The lookup mechanism in <code>FETools</code>::<code>get_fe_from_name</code> has been changed, so
       additional custom finite elements can be added using <code>FETools</code>::<code>add_fe_name</code>. In the course of this
       change, the implementation of the lookup mechanism has been
       simplified.
       <br>
       (GK 2006/10/24)
       </p>

  <li> <p>
       New: There is a new functions <code>GridTools</code>::<code>create_union_triangulation</code>
       that generates a triangulation that contains the respectively
       finest cells of two input triangulations.
       <br> 
       (WB 2006/10/23)
       </p>

  <li> <p>
       New: The <code>ConstraintMatrix</code> class did not
       allow that one degree of freedom was constrained against another DoF
       that was itself constrained. However, this is necessary for the
       implementation of hp methods and is now allowed. Arbitrarily long chains
       of constraints are resolved at the time the
       <code>ConstraintMatrix::close()</code> function is called. The only
       thing that is not allowed are cycles in such constraints.
       <br> 
       (WB 2006/09/19)
       </p>

  <li> <p>
       New: There are new functions <code>GridTools::minimal_cell_diameter</code> and <code>GridTools::maximal_cell_diameter</code>, with obvious
       functionality. 
       <br> 
       (WB 2006/09/06)
       </p>

  <li> <p>
       Changed: The functions <code>FETools::compute_embedding_matrices</code>,
       <code>FETools::compute_face_embedding_matrices</code>, and
       <code>FETools::compute_projection_matrices</code>
       (mostly used in internal computations in setting up finite
       element objects) previously took pointers to the first element
       of an array of matrices as arguments. This isn't type-safe, and
       in particular did not allow to check for the number of matrices
       in the array. The functions now take a reference to an array of
       the correct length.
       <br> 
       (WB 2006/08/14)
       </p>

  <li> <p>
       Extended: The <code>VectorTools::project</code> functions
       are now also implemented for 1d.
       <br> 
       (WB 2006/08/08)
       </p>

  <li> <p>
       Extended: <code >DerivativeApproximation</code> now offers access to the
       full tensor of derivatives of orders one, two and three. This
       information can be requested for a single cell by means of the <code>DerivativeApproximation</code><code>::approximate_derivative_tensor</code> function. If the
       norm of this tensor is required later on, the new <code >DerivativeApproximation</code><code>::derivative_norm</code> function can be used. Note, that
       for second order derivatives, this returns the largest eigenvalue
       instead of the Frobenius norm, which is returned by the <code >Tensor&lt;rank_,dim&gt;</code><code>::norm</code> function.
       <br> 
       (Tobias Leicht 2006/08/03)
       </p>

  <li> <p>
       Fixed: <code >DerivativeApproximation</code> offers approximated
       derivatives of a discrete function. The symmetrization of the derivative
       tensor is now done at the right place, i.e. the derivative itself is
       symmetrized instead of an intermediate tensor. This should improve the
       results slightly, but cause no problems otherwise, as this is completely
       internal to the class.
       <br> 
       (Tobias Leicht 2006/08/03)
       </p>

  <li> <p>
       Fixed: The <code >DataOut::build_patches</code> and similar
       functions in the related <code>DataOut*</code> classes allowed to pass
       zero as the second argument (denoting the number of threads to use if
       multithreading is enabled). This led to no output being created at
       all. This is now fixed by throwing an exception in this case.
       <br> 
       (WB 2006/07/31)
       </p>

  <li> <p>
       New: The new function <code>FiniteElementBase</code>::<code>n_dofs_per_object</code> returns either
       <code>dofs_per_vertex</code>, <code>dofs_per_line</code>,
       <code>dofs_per_quad</code>, ..., depending on the explicitly
       specified function template argument. This is often useful for
       template trickery.
       <br> 
       (WB, 2006/07/28)
       </p>

  <li> <p>
       Fixed: <code >Triangulation&lt;dim&gt;::fix_coarsen_flags</code>
       has been modified to allow coarsening in all possible cases. Up
       to now, coarsening was forbidden, if the neighbor cell was not refined
       but had the <code>refine_flag</code> set, whereas it was allowed, if
       the neighbor cell already was refined. Now, in both cases coarsening is
       allowed (if all children are flagged for coarsening). This leeds to
       triangulations with a slightly reduced number of cells. In some cases
       older references will have to be updated.
       <br> 
       (Tobias Leicht 2006/06/22)
       </p>

  <li> <p>
       New: There are now new internal <code>TriaObjectAccessor&lt;1,dim&gt;</code>::<code>lines()</code> and <code>TriaObjectAccessor&lt;2,dim&gt;</code>::<code>quads()</code> functions. By using these
       functions, 30 function specializations could be removed,
       significantly reducing code duplication.
       <br>
       (RH 2006/06/13)
       </p>

  <li> <p>
       New: Function <code>VectorTools</code>::<code>create_point_source_vector</code> to calculate the projection
       of a Dirac pulse on the base functions. This models a point source as
       used in the calculations of Green's functions and can also be used to
       extract the value of a finite element solution in a single point.
       <br> 
       (Ralf B. Schulz, 2006/05/15)
       </p>

  <li> <p>
       Changed: Functions <code>VectorTools</code>::<code>point_value</code> and <code>VectorTools</code>::<code>point_difference</code> using the old interface
       without boundary mapping were replaced by wrapper functions
       calling the new versions.
       <br> 
       (Ralf B. Schulz, 2006/05/15)
       </p>

  <li> <p>
       Changed: The old version of <code>GridTools</code>::<code>find_active_cell_around_point</code> has been replaced
       by a wrapper function for backward compatibility. The wrapper calls the
       new version of this function, and it is highly recommended to use the new
       version as it automatically delivers also the local coordinate of the
       point (so it can save some computation time in most cases as you don't have
       to calculate that again).
       <br> 
       (Ralf B. Schulz, 2006/05/12)
       </p>

  <li> <p>
       Improved: The functions <code>VectorTools</code>::<code>point_value</code> and <code>VectorTools</code>::<code>point_difference</code> now can also use arbitrary
       mappings, using the new <code>GridTools</code>::<code>find_active_cell_around_point</code> algorithm.
       <br> 
       (Ralf B. Schulz, 2006/05/11)
       </p>

  <li> <p>
       New: In <code>GridTools</code>, several functions have been added.
       <code>GridTools::find_closest_vertex</code> searches for the vertex
       located at closest distance to a given point. <code>
       GridTools::find_cells_adjacent_to_vertex</code> allows to determine
       all cells adjacent to a given vertex. And finally, a new version of
       <code>find_active_cell_around_point</code>, which takes additionally
       a mapping as parameter, implements a new and faster algorithm to
       determine the active cell in which a given point is located. For
       points located on boundaries and edges, it is in most cases also able
       to give the finest cell.
       <br> 
       (Ralf B. Schulz 2006/05/10)
       </p>

  <li> <p>
       Changed: The <code>DoFObjectAccessor::get_dof_values</code> and
       <code>DoFObjectAccessor::set_dof_values</code> were part of the
       accessors for lines, quads, and hexes. However, they could not
       be called for these objects unless the object was actually a
       cell, i.e. one could never call this function for a line or
       face in 3d, for example. The functions have therefore been
       moved to the <code>DoFCellAccessor</code> class that provides
       access to cells (i.e. lines in 1d, quads in 2d, and hexes in
       3d) for which this operation is actually useful.
       <br> 
       (WB 2006/05/01)
       </p>
  
  <li> <p> Fixed: second derivatives where not computed correctly in <code>FEFaceValuesBase</code>, i.e. when evaluating
       second derivatives on faces of cells. This is now fixed. Using
       second derivatives evaluated at quadrature points within a cell
       was not affected by this problem.
       <br>
       (GK 2006/04/28)
       </p>

  <li> <p>
       New: The functions <code>Triangulation::clear_user_flags_line</code>,
       <code>Triangulation::clear_user_flags_quad</code>, and
       <code>Triangulation::clear_user_flags_hex</code> can be used to
       selectively clear only some of the user flags as needed.
       <br> 
       (WB 2006/04/25)
       </p>

  <li> <p>
       New: The function <code>VectorTools::project</code> functions can now
       also be used for vector arguments of type other than
       <code>Vector&lt;double&gt;</code>.
       <br> 
       (WB 2006/04/17)
       </p>

  <li> <p>
       New: The function <code>GridTools::get_finest_common_cells</code> can be
       used to find those cells that are shared by two triangulations that are
       based on the same coarse mesh. This is useful, for example, when having
       to integrate the solution on one mesh against the shape functions on
       another mesh, a situation that frequently happens in time-dependent but
       also in nonlinear problems.
       <br> 
       (WB 2006/04/14)
       </p>

  <li> <p>
       New: The function <code>GridTools::have_same_coarse_mesh</code> figures
       out whether two triangulations, two DoFHandlers, etc, are built on the
       same coarse mesh.
       <br> 
       (WB 2006/04/14)
       </p>

  <li> <p>
       New: Since calling <code>cell-&gt;get_dof_indices</code> is a
       fairly frequent operation (called about 7 times per cell in
       step-14, but much more often in more complicated programs), the
       results of this function are now cached for faster access.
       <br> 
       (WB 2006/03/29)
       </p>

  <li> <p>
       Fixed: An exception was generated when trying to interpolate
       boundary values for a primitive component of a composite finite
       element that also has non-primitive base elements (for example,
       trying to interpolate boundary values for a the last component
       of a RT x Q1 element). This is now fixed.
       <br> 
       (WB 2006/03/27)
       </p>

  <li> <p>
       New: There is a new function <code>FiniteElement</code>::<code>face_to_equivalent_cell_index</code> that can
       convert a face index of a degree of freedom into a
       corresponding cell index.
       <br> 
       (WB 2006/03/27)
       </p>

  <li> <p>
       New: There are now functions <code>VectorTools</code>::<code>point_value</code> that evaluate the value of a
       finite element function at a given point inside the domain.
       <br> 
       (WB 2006/03/06)
       </p>

  <li> <p> Improved: <code>GridOut</code> now has functions for
       declaring and parsing parameters in a <code>ParameterHandler</code>, like used in <code>DataOut_Interface</code> for a long time.
       <br>
       (GK 2006/03/06)
       </p>

  <li> <p>
       Improved: The <code>Triangulation</code>, <code>PersistentTriangulation</code>, <code>DoFHandler</code>, <code>hp::DoFHandler</code> and <code>MGDoFHandler</code> classes now all have a <code>dimension</code> variable which allows to ask
       the template argument <tt>dim</tt> by
       <code>SomeClass::dimension</code>.
       <br> 
       (RH 2006/02/27)
       </p>

  <li> <p>
       Fixed: When used, the copy constructor of <code>MappingQ</code> would lead to memory
       corruption. This is now fixed.
       <br> 
       (RH 2006/02/23)
       </p>

  <li> <p>
       New: There is now a <code>StaticMappingQ1</code>::<code>mapping</code> object of type <code>MappingQ1</code> which is used in several parts
       of the library. This way multiple creation and storage of <code>MappingQ1</code> objects can be avoided.
       <br>
       Similar to <code>StaticMappingQ1</code>::<code>mapping</code> there is now also an object <code>hp::StaticMappingQ1</code>::<code>mapping_collection</code> of type <code>MappingCollection</code>.
       <br>
       (RH 2006/02/21)
       </p>

  <li> <p>
       New: There is now a <code>GridOut</code>::<code>write_gmsh</code> function to write in the gmsh 
       format. It has the same features of the write_ucd one. 
       <br> 
       (Luca Heltai 2006/02/17)
       </p>

  <li> <p>
       Changed: The <code>hp::FECollection</code>::<code>get_fe</code>, <code>hp::QCollection</code>::<code>get_quadrature</code> and <code>hp::MappingCollection</code>::<code>get_mapping</code> functions are now renamed to
       the standard <code>operator[]</code>
       function. Similarly, the <code>FECollection</code>::<code>n_finite_elements</code>, <code>QCollection</code>::<code>n_quadratures</code> and <code>MappingCollection</code>::<code>n_mappings</code> functions are now renamed to
       the standard <code>size()</code> function.
       <br> 
       (RH 2006/02/10)
       </p>

  <li> <p>
       Improved: <code>GridGenerator</code>::<code>half_hyper_ball</code> now is implemented also 
       in three dimensions, and the boundary now is colored with 0 and one 
       respectively on the curved and plane.
       <br> 
       (Luca Heltai 2006/02/07)
       </p>

  <li> <p>
       Improved: <code>FiniteElementData</code> and
       <code>FiniteElement</code> now support the
       concept of blocks along with components. See the glossary for
       the difference.
       <br>
       (GK 2006/01/25)
       </p>

  <li> <p>
       Fixed: The <code>MGDoFHandler</code>::<code>distribute_dofs</code> and <code>renumber</code> functions could not handle
       coarsened grids (unused vertices, faces and cells). This is now
       fixed.
       <br> 
       (RH 2006/01/23)
       </p>

  <li> <p>
       New: The <code>Mapping</code> class and derived classes now have
       functions <code>clone()</code> that return a new object of the same
       (derived) type.
       <br> 
       (WB 2006/1/8)
       </p>

  <li> <p>
       Improved: The <code>KellyErrorEstimator</code> class now also allows to
       select the cells on which it is supposed to work by giving a material
       ID, instead of only a subdomain ID.
       <br> 
       (WB 2006/1/4)
       </p>

  <li> <p>
       Improved: A new <code>TriaAccessor::ExcCellHasNoChildren</code>
       exception will be raised if the <code>TriaObjectAccessor::child_index</code> function
       is invoked for cells without children.
       <br> 
       (RH 2005/12/09)
       </p>

  <li> <p>
       Fixed: We had a bug in <code>DataOut::build_patches</code>
       that, when used in multithreaded mode, caused an exception to be
       thrown. In particular, this happened when the step-14 program was run on
       dual processor machines and the library was compiled for
       multithreading. This is now fixed.
       <br> 
       (WB 2005/10/20)
       </p>

  <li> <p>
       New: There is now a new <code>Triangulation::get_boundary_indicators</code>
       function.
       <br> 
       (RH 2005/09/30)
       </p>

  <li> <p>
       New: There is now a new <code>GridTools::delete_unused_vertices</code>
       function. Previously a <tt>private</tt> function in <code>GridIn</code> it has now been moved to and made
       <tt>public</tt> in <code>GridTools</code>.
       <br> 
       (RH 2005/09/28)
       </p>

  <li> <p>
       New: The <code>GridIn&lt;dim&gt;::read_netcdf(string
       &filename)</code> function reads grids from NetCDF files. The
       only data format currently supported is the <tt>TAU grid
       format</tt>.
       <br> 
       (RH 2005/09/23)
       </p>

  <li> <p>
       Fixed: The <code>GridIn&lt;dim&gt;::read(filename, format)</code>
       function ran into an exception when called with
       <tt>Default</tt> format and a filename which does not include
       <tt>/</tt>. This is now fixed.
       <br> 
       (RH 2005/09/21)
       </p>
</ol>



*/
