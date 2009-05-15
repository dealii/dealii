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
   
  <li>
   <p>
   Fixed: The DoFRenumbering::component_wise function for MGDoFHandler objects
   did a few things in parallel that weren't thread-safe. This is now fixed.
   <br>
   (WB, 2009/01/20)
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
