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
  </p>
  </li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
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
  New: The SparseMatrix class has now a function SparseMatrix::mmult that
  can multiply two sparse matrices with each other.
  <br>
  (Martin Kronbichler 2009/05/04)
  </p>
  </li>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
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
/ol>


*/
