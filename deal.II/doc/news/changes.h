/**
 * @page changes_after_6.2 Changes after Version 6.2

<p>
This is the list of changes made after the release of 
deal.II version 6.2. It is subdivided into changes
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
  <li>
  <p>
  Changed: Previously, the Triangulation::create_triangulation
  function silently accepted input meshes with inverted cells
  (i.e. cells with a zero or negative determinant of the Jacobian of
  the mapping from the reference cell to the real cell). This has been
  changed now: The function checks whether cells are distorted or
  inverted, and may throw an exception containing a list of cells
  for which this is the case. If you know that this is harmless, for
  example if you have cells with collapsed vertices in your mesh but
  you do not intend to integrate on them, then you can catch and
  ignore this message. In all other cases, the output of your
  computations are likely to be wrong anyway.
  <br>
  The same is true for the Triangulation::execute_coarsening_and_refinement
  function: if it creates cells that are distorted, it throws a list of cells
  whose children are distorted.
  <br>
  The whole issue is described in some detail in the entry on
  @ref GlossDistorted "distorted cells" in the glossary.
  <br>
  (WB 2009/06/29)
  </p>
  </li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
   <li>
   <p> 
   Changed: When configuring to use METIS for partitioning meshes in parallel,
   the METIS header files had to be modified by hand. In addition, with some
   MPI implementations one would get into trouble if <code>mpi.h</code>
   included <code>mpicxx.h</code>. These two problems have now been
   worked around.
   <br>
   (WB 2009/07/06)
   </p>
 
   <li>
   <p> 
   Changed: Previously, one had to give the two flags
   <code>--enable-multithreading --with-multithreading</code> to
   <code>./configure</code> to enable thread usage throughout the library.
   This has now been simplified: only the flag <code>--enable-threads</code>
   is now necessary. Furthermore, since most current machines have multiple
   cores these days, the default is now to use threads. This can be switched
   off using <code>--disable-threads</code>, however.
   <br>
   (WB 2008/09/29)
   </p>
 
   <li>
   <p>
   New: As a primary means of parallelizing programs, deal.II now uses
   a task-based, rather than thread-based approach, in which one
   uses a high-level description of <i>what</i> needs to be done,
   rather than how these jobs have to be mapped onto threads. We then
   use the <a href="http://www.threadingbuildingblocks.org">Threading
   Building Blocks (TBB) library</a> to schedule tasks onto available
   hardware resources. This new scheme of describing parallism and
   various abstractions to make programming in this framework easier
   are described in great detail in the
   @ref threads "Parallel computing with multiple processors" module.
   In addition, most of the parallelism already used within deal.II
   has been converted to use tasks, rather than threads, and so have
   some of the tutorial programs.
   <br>
   (WB 2009/01/09)
   </p>
   </li>
 
   <li>
   <p>
   Changed: The support for threading has been completely re-written. In
   particular, the Threads::spawn functions have been deprecated, and
   new functions Threads::new_threads have been introduced.
   Threading is now discussed in a lot of detail in the
   @ref threads "Parallel computing with multiple processors" module.
   <br>
   (WB 2009/01/09)
   </p>
   </li>   
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li>
  <p>
  New: There is now a specialization Tensor<0,dim> of tensors of rank 0. Since rank-0
  tensors are scalars, this class essentially acts like a scalar, but it allows for
  some neat template tricks that involve tensors of arbitrary rank.
  <br>
  (WB 2009/09/15)
  </p>
  </li>

  <li>
  <p>
  New: The GeometryInfo::jacobian_determinants_at_vertices can be used
  to investigate the degree of distortion of cells.
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  New: The GeometryInfo::d_linear_shape_function and
  GeometryInfo::d_linear_shape_function_gradient functions can be used
  to represent the $d$-linear shape functions that are frequently
  used to map the reference cell to real cells (though the
  Mapping class hierarchy also allows to use higher order mappings).
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  New: The determinant() function is now implemented for rank-2 Tensor
  arguments of all sizes. The implementation is not efficient for very large
  matrix sizes, however.
  <br>
  (WB 2009/06/28)
  </p>
  </li>

  <li>
  <p>
  Improved: The QGaussLobatto::gamma function now returns a long double
  instead of an unsigned int, otherwise we will get an overflow and thus
  meaningless weights for higher QGaussLobatto quadrature rules.  
  <br>
  (Tobias Leicht, RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: The new function Utilities::duplicate_communicator can be used
  to duplicate an MPI communicator to produce a unique duplicate.
  <br>
  (WB 2009/05/13)
  </p>
  </li>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li>
  <p>
  New: Based on work by Francisco Alvaro, the existing SLEPcWrappers now 
  have a handle on the generalized eigenvalue problem where B=I.
  <br>
  (Toby D. Young 2009/06/25)
  </p>
  </li>
</ol>


<ol>
  <li>
  <p>
  Fixed: The TrilinosWrappers::MPI::BlockVector class declares an assignment
  operator from the non-Trilinos BlockVector class but it could not be
  compiled due to an oversight. This is now fixed.
  <br>
  (WB 2009/06/29)
  </p>
  </li>

  <li>
  <p>
  Fixed: The TrilinosWrappers::BlockVector class declares an assignment
  operator from the non-Trilinos BlockVector class but it wasn't implemented.
  This is now fixed.
  <br>
  (WB 2009/06/24)
  </p>
  </li>

  <li>
  <p>
  New: The SparseMatrix class has now a function SparseMatrix::mmult that
  can multiply two sparse matrices with each other.
  <br>
  (Martin Kronbichler 2009/05/04)
  </p>
  </li>
</ol>

<ol>
  <li>
  <p>
  New: Based on work with Rickard Armiento, Francisco Alvaro, and Jose
  E. Roman, SLEPcWrappers that give a handle on some of the features
  of SLEPc (Scalable Library for Eigenvalue Problem Computations): (1)
  The SLEPcWrappers::SolverBase class can be used for specifying an
  eigenvalue problem, either in standard or generalized form, on
  serial or parallel architectures with support for a few solver
  types; and (2) The SLEPcWrappers::TransformationBase class
  encapsulates a variety of spectral transformations providing some
  functionality required for acceleration techniques based on the
  transformation of the spectrum.  
  <br>
  (Toby D. Young 2009/06/25)
  </p>
  </li>
</ol>


<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li>
  <p>
  Fixed: The CellAccessor::recursively_set_material_id function did not
  set the material id for all children, but only for the first two, which
  is obviously a bug. This should now be fixed.
  <br>
  (WB 2009/07/14)
  </p>
  </li>

  <li>
  <p>
  Fixed: The GridIn class sometimes had problems with input files that had
  whitespace at the end of lines. This should now be fixed.
  <br>
  (WB 2009/07/10)
  </p>
  </li>

  <li>
  <p>
  New: The new hp::DoFHandler::set_active_fe_indices function allows
  to distribute all active FE indices at once based on a given
  vector. This might be useful if this information is stored
  somewhere and has to be reconstructed or else if two DoFHandler
  objects with the same FE index distribution should be created.
  There is now also a corresponding
  hp::DoFHandler::get_active_fe_indices function.  
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>

  <li>
  <p>
  Fixed: The projection of quadrature points to subfaces in
  MappingQ in case of 3d anisotropic refinement did not respect
  non-standard face orientation/flip/rotation cases. This
  has now been fixed.  
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>
  
  <li>
  <p>
  New: The new Triangulation::n_raw_faces() function forwards
  to Triangulation::n_raw_lines() in 2d and
  Triangulation::n_raw_quads() in 3d.  
  <br>
  (Tobias Leicht, RH 2009/06/12)
  </p>
  </li>

  <li>
  <p>
  New: There is now a new DataOutFaces::build_patches function which
  takes a Mapping argument. For higher order mappings this allows to
  represent curved boundaries by using more subdivisions. This function
  is also useful in the context of MappingQ1Eulerian.  
  <br>
  (Tobias Leicht, RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: For empty triangulations the new Triangulation::set_mesh_smoothing
  function allows to override the MeshSmoothing given to the constructor.  
  <br>
  (RH 2009/06/05)
  </p>
  </li>

  <li>
  <p>
  New: The new function TriaAccessor::is_translation_of computes
  whether a cell, face, or edge is a translation of another.
  <br>
  (Martin Kronbichler, WB 2009/05/19)
  </p>
  </li>

  <li>
  <p>
  New: The DoFTools::make_sparsity_pattern functions have acquired a
  new paramater <code>subdomain_id</code>. If a value other than the
  default value is passed for it, the function only assembles the
  sparsity pattern on those cells that have the given subdomain id.
  This is useful, for example, in conjunction with the
  TrilinosWrappers::SparsityPattern class that can hold a sparsity
  pattern distributed across several MPI processes; in that case, it is
  not necessary that each process builds the entire sparsity pattern.
  <br>
  (WB 2009/04/29)
  </p>
  </li>
  
   <li>
   <p>
   Fixed: The DoFRenumbering::component_wise function for MGDoFHandler objects
   did a few things in parallel that weren't thread-safe. This is now fixed.
   <br>
   (WB, 2009/01/20)
   </p>
 
   <li>
   <p> 
   Changed: The two DataOut::build_patches, DataOutFaces::build_patches, and
   DataOutRotation::build_patches functions have lost the argument
   that indicated the number of threads with which they should build the
   intermediate representation. This is something that now happens
   transparently in the background and doesn't need caller input any more.
   <br>
   (WB 2008/12/16)
   </p>

   <li>
   <p>
   Changed: The KellyErrorEstimator::estimate functions had a parameter
   that indicates the number of threads to be used in the computation.
   This parameter continues to exist for compatibility, but is now ignored.
   Rather, the number of threads is determined automatically by scheduling
   the requested computations on available compute resources.
   <br>
   (WB, 2008/12/29)
   </p>
 
   <li>
   <p>
   New: The new function internal::hp::FEValuesBase::get_fe_collection function
   allows to query the finite element collection currently in used in an hp::FEValues,
   hp::FEFaceValues, or hp::FESubfaceValues object.
   <br>
   (WB 2008/09/30)
   </p> 
 
   <li>
   <p>
   New: The new function FEValuesBase::get_update_flags allows to query
   the update flags that are currently set on an FEValues, FEFaceValues, or
   FESubfaceValues object.
   <br>
   (WB 2008/09/29)
   </p> 
</ol>


*/
