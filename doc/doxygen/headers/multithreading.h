// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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
 * @defgroup threads Parallel computing with multiple processors accessing shared memory
 * @ingroup Parallel
 *
 * @brief A module discussing the use of parallelism on shared memory
 * machines. See the detailed documentation and
 * @ref MTToC "Table of Contents" below the lenghty list of members
 * of this module.
 *
 * @dealiiVideoLecture{39,40}
 *
 * On machines with more than one processor (or multicore processors),
 * it is often profitable to run several parts of the computations in
 * %parallel. For example, one could have several threads running in
 * %parallel, each of which assembles the cell matrices of a subset of
 * the triangulation and then writes them into the global
 * matrix. Since assembling matrices is often an expensive operation,
 * this frequently leads to significant savings in compute time on
 * multiprocessor machines.
 *
 * deal.II supports operations running in %parallel on shared-memory
 * (SMP) machines through the functions and classes in the Threads
 * namespace. The MultithreadInfo class allows to query certain
 * properties of the system, such as the number of CPUs. These
 * facilities for %parallel computing are described in the
 * following. The step-9, step-13, step-14, step-32, step-35 and
 * step-37 tutorial programs also show their use in practice, with the
 * ones starting with step-32 using a more modern style of doing
 * things in which essentially we describe <i>what</i> can be done in
 * %parallel, while the older tutorial programs describe <i>how</i>
 * things have to be done in %parallel.
 *
 * On the other hand, programs running on distributed memory machines
 * (i.e. clusters) need a different programming model built on top of MPI and
 * PETSc or Trilinos. This is described in the step-17, step-18 and step-32
 * example programs.
 *
 * @anchor MTToC
 * <table class="tutorial" width="50%">
 * <tr><th><b>%Table of contents</b></th></tr>
 * <tr><td width="100%" valign="top">
 * <ol>
 *  <li> @ref MTTasks "Task-based parallelism"
 *  <li> @ref MTUsing "Using tasks from within deal.II"
 *  <li> @ref MTHow "How scheduling tasks works and when task-based programming is not efficient"
 *  <li> @ref MTSimpleLoops "Abstractions for tasks: Simple loops"
 *  <li> @ref MTComplexLoops "Abstractions for tasks: More complex loops"
 *  <li> @ref MTWorkStream "Abstractions for tasks: Work streams"
 *  <li> @ref MTTaskSynchronization "Tasks and synchronization"
 *  <li> @ref MTThreads "Thread-based parallelism"
 *  <li> @ref MTTaskThreads "Controlling the number of threads used for tasks"
 * </ol> </td> </tr> </table>
 *
 *
 * @anchor MTTasks
 * <h3>Task-based parallelism</h3>
 *
 * The traditional view of parallelism on shared memory machines has been to
 * decompose a program into <i>threads</i>, i.e. running different parts of
 * the program in %parallel <i>at the same time</i> (if there are more threads
 * than processor cores on your machine, the operating system will run each
 * thread round-robin for a brief amount of time before switching execution to
 * another thread, thereby simulating that threads run
 * concurrently). deal.II's facilities for threads are described below (see
 * @ref MTThreads "Thread-based parallelism"), but we would first like to
 * discuss an abstraction that is often more suitable than threads:
 * <i>tasks</i>.
 *
 * Tasks are essentially the individual parts of a program. Some of them are
 * independent, whereas others depend on previous tasks to be completed
 * first. By way of example, consider the typical layout of a part of the
 * <code>setup_dofs</code> function that most of the tutorial programs have:
 * @code
1  dof_handler.distribute_dofs (fe);
2  DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
3  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
4  hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 *
 * Here, each of the operations require a significant amount of
 * computations. But note that not all of them depend on each other: clearly
 * we can not run statements 2-4 before 1, and 4 needs to wait for the
 * completion of statements 2 and 3. But statements 2 and 3 are independent:
 * they could be run in any order, or in %parallel. In essence, we have
 * identified four <i>tasks</i>, some of which are dependent on each other,
 * whereas others are independent. In the current example, tasks are
 * identified with individual C++ statements, but often they more
 * generally coincide with entire code blocks.
 *
 * The point here is this: If we wanted to use threads to exploit the
 * independence of tasks 2 and 3, we would start two threads and run each of
 * tasks 2 and 3 on its own thread; we would then wait for the two threads to
 * finish (an operation called "joining a thread") and go on with statement
 * 4. Code to achieve this would look like this (the actual syntax is
 * explained in more detail below):
 * @code
   dof_handler.distribute_dofs (fe);

   Threads::Thread<void>
     thread_1 = Threads::new_thread (&DoFTools::make_hanging_node_constraints,
                                     dof_handler, hanging_node_constraints);
   Threads::Thread<void>
     thread_2 = Threads::new_thread (&DoFTools::make_sparsity_pattern,
                                     dof_handler, sparsity_pattern);
   thread_1.join();
   thread_2.join();
   hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 *
 * But what if
 * your computer has only one processor core, or if we have two but there is
 * already a different part of the program running in %parallel to the code
 * above? In that case, the code above would still start new threads, but the
 * program is not going to run faster since no additional compute resources
 * are available; rather, the program will run slower since threads have to be
 * created and destroyed, and the operating system has to schedule threads to
 * oversubscribed compute resources.
 *
 * A better scheme would identify independent tasks and then hand them off to
 * a scheduler that maps tasks to available compute resources. This way, the
 * program could, for example, start one thread per processor core and then
 * let threads work on tasks. Tasks would run to completion, rather than
 * concurrently, avoiding the overhead of interrupting threads to run a
 * different thread. In this model, if two processor cores are available,
 * tasks 2 and 3 above would run in %parallel; if only one is available, the
 * scheduler would first completely execute task 2 before doing task 3, or the
 * other way around. This model is able to execute much more efficiently in
 * particular if a large number of tasks is available for execution, see for
 * example the discussion below in section
 * @ref MTWorkStream "Abstractions for tasks: Work streams". In
 * essence, tasks are a high-level description of what needs to be done,
 * whereas threads are a low-level way of implementing how these tasks can be
 * completed. As in many other instances, being able to use a high-level
 * description allows to find efficient low-level implementations; in this
 * vein, it often pays off to use tasks, rather than threads, in a program.
 *
 * deal.II does not implement scheduling tasks to threads itself. For this, we
 * use the <a href="http://www.threadingbuildingblocks.org" target="_top">Threading Building
 * Blocks (TBB) library</a> for which we provide simple wrappers. TBB
 * abstracts the details of how to start or stop threads, start tasks on
 * individual threads, etc, and provides interfaces that are portable across
 * many different systems.
 *
 *
 *
 * @anchor MTUsing
 * <h3>Using tasks from within deal.II</h3>
 *
 * Ideally, the syntax to start tasks (and similarly for threads, for that
 * matter), would be something like this for the example above:
 * @code
   Threads::Task<void>
     thread
     = new_task DoFTools::make_hanging_node_constraints (dof_handler,
                                                         hanging_node_constraints);
 * @endcode
 * In other words, we would like to indicate the fact that the function call
 * should be run on a separate task by simply prefixing the call with a
 * keyword (such as <code>new_task</code> here, with a similar keyword
 * <code>new_thread</code> for threads). Prefixing a call would return a
 * handle for the task that we can use to wait for the tasks's completion and
 * that we may use to query the return value of the function called (unless it
 * is void, as it is here).
 *
 * Since C++ does not support the creation of new keywords, we have to be a
 * bit more creative. The way chosen is to introduce a function
 * <code>new_task</code> that takes as arguments the function to call as well
 * as the arguments to the call. The <code>new_task</code> function is
 * overloaded to accommodate starting tasks with functions that take no, one,
 * two, and up to 9 arguments. In deal.II, these functions live in the Threads
 * namespace. Consequently, the actual code for what we try to do above looks
 * like this:
 * @code
   Threads::Task<void>
     thread
     = Threads::new_task (&DoFTools::make_hanging_node_constraints,
                          dof_handler,
                          hanging_node_constraints);
 * @endcode
 * Note that DoFTools::make_hanging_node_constraints is a static member
 * function and so does not need an object of type DoFTools to work on.
 * (In fact, DoFTools has only static member functions and could as well be
 * a namespace instead of class; that it is a class at the time of writing
 * this is mostly a historic relic.)
 *
 * Similarly, if we want to call a member function on a different task, we can
 * do so by specifying the object on which to call the function as first
 * argument after the function pointer:
 * @code
   class C {
     public:
       double f(int);
   };

   int main () {
     C c;

     // call f(13) as usual, i.e. using the current processor:
     c.f(13);

     // call f(42) as a separate task, to be scheduled
     // whenever processor resources become available:
     Threads::Task<double>
       task = Threads::new_task (&C::f, c, 42);

     // do something else in between:
     ...;

     // having finished our other business, wait until the task
     // above has terminated and get the value returns by c.f(42):
     double result = task.return_value();
 * @endcode
 * Here, note first how we pass the object <code>c</code> (i.e. the
 * <code>this</code> pointer the function <code>C::f</code> will see) as if it
 * was the first argument to the function. Secondly, note how we can acquire
 * the value returned by the function on the separate task by calling
 * Threads::Task::return_value(). This function implies waiting for the
 * completion of the task, i.e. the last line is completely equivalent to
 * @code
     task.join ();
     double result = task.return_value();
 * @endcode
 *
 * Note also that it is entirely valid if <code>C::f</code> wants to start
 * tasks of its own:
 * @code
   class C {
     public:
       double f(int);
     private:
       double f1(int);
       double f2(int);
   };

   double C::f (int i) {
     Threads::Task<double> t1 = Threads::new_task (&C::f1, *this, i);
     Threads::Task<double> t2 = Threads::new_task (&C::f2, *this, i);
     return t1.return_value() + t2.return_value();
   }

   int main () {
     C c;

     Threads::Task<double>
       task = Threads::new_task (&C::f, c, 42);

     // do something else in between:
     ...;

     double result = task.return_value();
 * @endcode
 * Here, we let <code>C::f</code> compute its return value as
 * <code>c.f1(i)+c.f2(i)</code>. If sufficient CPU resources are available,
 * then the two parts of the addition as well as the other things in
 * <code>main()</code> will all run in %parallel. If not, then we will
 * eventually block at one of the places where the return value is needed,
 * thereby freeing up the CPU resources necessary to run all those spawned
 * tasks to completion.
 *
 *
 * In many cases, such as the introductory example of the
 * <code>setup_dofs</code> function outlined above, one can identify several
 * independent jobs that can be run as tasks, but will have to wait for all of
 * them to finish at one point. One can do so by storing the returned object
 * from all the Threads::new_task() calls, and calling Threads::Task::join()
 * on each one of them. A simpler way to do this is to put all of these task
 * objects into a Threads::TaskGroup object and waiting for all of them at
 * once. The code would then look like this:
 * @code
   dof_handler.distribute_dofs (fe);

   Threads::TaskGroup<void> task_group;
   task_group += Threads::new_task (&DoFTools::make_hanging_node_constraints,
                                    dof_handler, hanging_node_constraints);
   task_group += Threads::new_task (&DoFTools::make_sparsity_pattern,
                                    dof_handler, sparsity_pattern);
   task_group.join_all ();
   hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 *
 *
 * @anchor MTHow
 * <h3>How scheduling tasks works and when task-based programming is not efficient</h3>
 *
 * The exact details of how tasks are scheduled to run are %internal to the
 * Threading Building Blocks (TBB) library that deal.II uses for tasks. The
 * documentation of TBB gives a detailed description of how tasks are
 * scheduled to threads but is rather quiet on how many threads are actually
 * used. However, a reasonable guess is probably to assume that TBB creates as
 * many threads as there are processor cores on your system. This way, it is
 * able to fully utilize the entire system, without having too many threads
 * that the operating system will then have to interrupt regularly so that
 * other threads can run on the available processor cores.
 *
 * The point then is that the TBB scheduler takes tasks and lets threads
 * execute them. %Threads execute tasks completely, i.e. the TBB scheduler does
 * not interrupt a task half way through to make some halfway progress with
 * another task. This makes sure that caches are always hot, for example, and
 * avoids the overhead of preemptive interrupts.
 *
 * The downside is that the CPU cores are only fully utilized if the threads
 * are actually doing something, and that means that (i) there must be enough
 * tasks available, and (ii) these tasks are actually doing something. Note
 * that both conditions must be met; in particular, this means that CPU cores
 * are underutilized if we have identified a sufficient number of tasks but if
 * some of them twiddle thumbs, for example because a task is writing data to
 * disk (a process where the CPU frequently has to wait for the disk to
 * complete a transaction) or is waiting for input. Other cases are where
 * tasks block on other external events, for example synchronising with other
 * tasks or threads through a mutex. In such cases, the scheduler would let a
 * task run on a thread, but doesn't notice that that thread doesn't fully
 * utilize the CPU core.
 *
 * In cases like these, it <i>does</i> make sense to create a new thread (see
 * @ref MTThreads "Thread-based parallelism" below) that the operating system
 * can put on hold while they are waiting for something external, and let a
 * different thread (for example one running a task scheduled by TBB) use the
 * CPU at the same time.
 *
 *
 * @anchor MTSimpleLoops
 * <h3>Abstractions for tasks: Simple loops</h3>
 *
 * Some loops execute bodies on data that is completely independent
 * and that can therefore be executed in %parallel. Rather than a
 * priori split the loop into a fixed number of chunks and executing
 * them on tasks or threads, the TBB library uses the following
 * concept: the range over which the loop iterates is split into a
 * certain number of sub-ranges (for example two or three times as
 * many as there are CPU cores) and are equally distributed among
 * threads; threads then execute sub-ranges and, if they are done with
 * their work, steal entire or parts of sub-ranges from other threads
 * to keep busy. This way, work is load-balanced even if not every
 * loop iteration takes equally much work, or if some of the CPU cores fall
 * behind because the operating system interrupted them for some other
 * work.
 *
 * The TBB library primitives for this are a bit clumsy so deal.II has
 * wrapper routines for the most frequently used operations. The
 * simplest one is akin to what the std::transform does: it takes
 * one or more ranges of input operators, one output iterator, and a
 * function object. A typical implementation of std::transform would
 * look like this:
 * @code
     template <typename InputIterator1, typename InputIterator,
               typename OutputIterator, typename FunctionObject>
     void transform (const InputIterator1 &begin_in_1,
                     const InputIterator1 &end_in_1,
                     const InputIterator2 &begin_in_2,
                     const OutputIterator &begin_out,
                     FunctionObject       &function)
     {
       InputIterator1 in_1 = begin_in_1;
       InputIterator2 in_2 = begin_in_2;
       OutputIterator out  = begin_out;

       for (; in_1 != end_in_1; ++in_1, ++in_2, ++out)
         *out = function(*in_1, *in_2);
     }
 * @endcode
 *
 * In many cases, <code>function</code> has no state, and so we can
 * split this loop into several sub-ranges as explained
 * above. Consequently, deal.II has a set of functions
 * parallel::transform that look like the one above but that do their
 * work in %parallel (there are several versions with one, two, and
 * more input iterators for function objects that take one, two, or
 * more arguments). The only difference in calling these functions is
 * that they take an additional last argument that denotes the minimum
 * size of sub-ranges of <code>[begin_in_1,end_in_1)</code>; it should
 * be big enough so that we don't spend more time on scheduling
 * sub-ranges to processors but small enough that processors can be
 * efficiently load balanced. A rule of thumb appears to be that a
 * sub-range is too small if it takes less than 2000 instructions to
 * execute it.
 *
 * An example of how to use these functions are vector operations like
 * the addition in $z = x+y$ where all three objects are of type Vector:
 * @code
     parallel::transform (x.begin(), x.end(),
                          y.begin(),
                          z.begin(),
                          (boost::lambda::_1 + boost::lambda::_2),
                          1000);
 * @endcode
 *
 * In this example, we used the <a
 * href="http://www.boost.org/doc/libs/1_37_0/doc/html/lambda.html">Boost
 * Lambda</a> library to construct, on the fly, a function object that
 * takes two arguments and returns the sum of the two. This is exactly
 * what we needed when we want to add the individual elements of
 * vectors $x$ and $y$ and write the sum of the two into the elements
 * of $z$. Because of the way Boost Lambda is written, the function
 * object that we get here is completely known to the compiler and
 * when it expands the loop that results from parallel::transform will
 * be as if we had written the loop in its obvious form:
 * @code
       InputIterator1 in_1 = x.begin();
       InputIterator2 in_2 = y.begin();
       OutputIterator out  = z.begin();

       for (; in_1 != x.end(); ++in_1, ++in_2, ++out)
         *out = *in_1 + *in_2;
 * @endcode
 * The next C++ standard will contain a more elegant way to achieve the
 * same effect shown above using the Boost library, through a
 * mechanism known as <i>lambda expressions</i> and <i>closures</i>.
 *
 * Note also that we have made sure that no CPU ever gets a chunk of
 * the whole loop that is smaller than 1000 iterations (unless the
 * whole range is smaller).
 *
 *
 * @anchor MTComplexLoops
 * <h3>Abstractions for tasks: More complex loops</h3>
 *
 * The scheme shown in the previous section is effective if the
 * operation done in each iteration is such that it does not require
 * significant setup costs and can be inlined by the compiler. Lambda
 * expressions are exactly of this kind because the compiler knows
 * everything about the lambda expression and can inline it, thereby
 * eliminating the overhead of calling an external function. However,
 * there are cases where it is inefficient to call some object or
 * function within each iteration.
 *
 * An example for this case is sparse matrix-vector multiplication. If you
 * know how data is stored in compressed row format like in the SparseMatrix
 * class, then a matrix-vector product function looks like this:
 * @code
    void SparseMatrix::vmult_one_row (const Vector &src,
                                      Vector       &dst) const
    {
      const double       *val_ptr    = &values[0];
      const unsigned int *colnum_ptr = &colnums[0];
      Vector::iterator dst_ptr = dst.begin();

      for (unsigned int row=0; row<n_rows; ++row, ++dst_ptr)
        {
          double s = 0.;
          const double *const val_end_of_row = &values[rowstart[row+1]];
          while (val_ptr != val_end_of_row)
            s += *val_ptr++ * src(*colnum_ptr++);
          *dst_ptr = s;
        }
    }
 * @endcode
 * Inside the for loop, we compute the dot product of a single row of the
 * matrix with the right hand side vector <code>src</code> and write it into
 * the corresponding element of the <code>dst</code> vector. The code is made
 * more efficient by utilizing that the elements of the <i>next</i> row follow
 * the ones of the current row <i>immediately</i>, i.e. at the beginning of
 * the loop body we do not have to re-set the pointers that point to the
 * values and column %numbers of each row.
 *
 * Using the parallel::transform function above, we could in principle write
 * this code as follows:
 * @code
    void SparseMatrix::vmult (const Vector     &src,
                              Vector           &dst,
                              Vector::iterator &dst_row) const
    {
      const unsigned int  row = (dst_row - dst.begin());

      const double       *val_ptr    = &values[rowstart[row]];
      const unsigned int *colnum_ptr = &colnums[rowstart[row]];

      double s = 0.;
      const double *const val_end_of_row = &values[rowstart[row+1]];
      while (val_ptr != val_end_of_row)
        s += *val_ptr++ * src(*colnum_ptr++);
      *dst_row = s;
    }

    void SparseMatrix::vmult (const Vector &src,
                              Vector       &dst) const
    {
      parallel::transform (dst.begin(), dst.end(),
                           std_cxx11::bind (&SparseMatrix::vmult_one_row,
                                        this,
                                        std_cxx11::cref(src),
                                        std_cxx11::ref(dst),
                                        std_cxx11::_1),
                           200);
    }
 * @endcode
 * Note how we use <a
 * href="http://www.boost.org/doc/libs/1_37_0/libs/bind/bind.html">std_cxx11::bind</a>
 * to <i>bind</i> certain arguments to the <code>vmult_one_row</code>
 * function, leaving one argument open and thus allowing the
 * parallel::transform function to consider the passed function argument as
 * unary. Also note that we need to make the source and destination vectors as
 * (const) references to prevent std_cxx11::bind from passing them by value
 * (implying a copy for <code>src</code> and writing the result into a
 * temporary copy of <code>dst</code>, neither of which is what we desired).
 * Finally, notice the grainsize of a minimum of 200 rows of a matrix that
 * should be processed by an individual CPU core.
 *
 * The point is that while this is correct, it is not efficient: we have to
 * set up the <code>row, val_ptr, colnum_ptr</code> variables in each
 * iteration of the loop. Furthermore, since now the function object to be
 * called on each row is not a simple Boost Lambda expression any more, there
 * is an implied function call including argument passing in each iteration of
 * the loop.
 *
 * A more efficient way is to let TBB split the original range into
 * sub-ranges, and then call a target function not on each individual element
 * of the loop, but on the entire range. This is facilitated by the
 * parallel::apply_to_subranges function:
 * @code
    void
    SparseMatrix::vmult_on_subrange (const unsigned int  begin_row,
                                     const unsigned int  end_row,
                                     const Vector     &src,
                                     Vector           &dst)
    {
      const double       *val_ptr    = &values[rowstart[begin_row]];
      const unsigned int *colnum_ptr = &colnums[rowstart[begin_row]];
      Vector::iterator dst_ptr = dst.begin() + begin_row;

      for (unsigned int row=begin_row; row<end_row; ++row, ++dst_ptr)
        {
          double s = 0.;
          const double *const val_end_of_row = &values[rowstart[row+1]];
          while (val_ptr != val_end_of_row)
            s += *val_ptr++ * src(*colnum_ptr++);
          *dst_ptr = s;
        }
    }

    void SparseMatrix::vmult (const Vector &src,
                              Vector       &dst) const
    {
       parallel::apply_to_subranges (0, n_rows(),
                                     std_cxx11::bind (vmult_on_subrange,
                                                  this,
                                                  std_cxx11::_1, std_cxx11::_2,
                                                  std_cxx11::cref(src),
                                                  std_cxx11::ref(dst)),
                                     200);
    }
 * @endcode
 * Here, we call the <code>vmult_on_subrange</code> function on sub-ranges
 * of at least 200 elements each, so that the initial setup cost can amortize.
 *
 * A related operation is when the loops over elements each produce a
 * result that must then be accumulated (other reduction operations
 * than addition of numbers would work as well). An example is to form
 * the matrix norm $x^T M x$ (it really is only a norm if $M$ is
 * positive definite, but let's assume for a moment that it is). A
 * sequential implementation would look like this for sparse matrices:
 * @code
    double SparseMatrix::mat_norm (const Vector &x) const
    {
      const double       *val_ptr    = &values[0];
      const unsigned int *colnum_ptr = &colnums[0];

      double norm_sqr = 0;

      for (unsigned int row=0; row<n_rows; ++row, ++dst_ptr)
        {
          double s = 0.;
          const double *const val_end_of_row = &values[rowstart[row+1]];
          while (val_ptr != val_end_of_row)
            s += *val_ptr++ * x(*colnum_ptr++);
          norm_sqr += x(row) * s;
        }

      return std::sqrt (norm_sqr);
    }
 * @endcode
 *
 * It would be nice if we could split this operation over several
 * sub-ranges of rows, each of which compute their part of the square
 * of the norm, add results together from the various sub-ranges, and
 * then take the square root of the result. This is what the
 * parallel::accumulate_from_subranges function does (note that you
 * have to specify the result type as a template argument and that, as
 * usual, the minimumum number of elements of the outer loop that can
 * be scheduled on a single CPU core is given as the last argument):
 * @code
    double
    SparseMatrix::mat_norm_sqr_on_subrange (const unsigned int  begin_row,
                                            const unsigned int  end_row,
                                            const Vector     &x)
    {
      const double       *val_ptr    = &values[rowstart[begin_row]];
      const unsigned int *colnum_ptr = &colnums[rowstart[begin_row]];
      Vector::iterator dst_ptr = dst.begin() + begin_row;

      double norm_sqr = 0;

      for (unsigned int row=begin_row; row<end_row; ++row, ++dst_ptr)
        {
          double s = 0.;
          const double *const val_end_of_row = &values[rowstart[row+1]];
          while (val_ptr != val_end_of_row)
            s += *val_ptr++ * x(*colnum_ptr++);
          norm_sqr += x(row) * s;
        }

      return norm_sqr;
    }

    double SparseMatrix::mat_norm (const Vector &x) const
    {
      return
        std::sqrt
        (parallel::accumulate_from_subranges (0, n_rows(),
                                              std_cxx11::bind (mat_norm_sqr_on_subrange,
                                                           this,
                                                           std_cxx11::_1, std_cxx11::_2,
                                                           std_cxx11::cref(x)),
                                              200));
    }
 * @endcode
 *
 *
 * @anchor MTWorkStream
 * <h3>Abstractions for tasks: Work streams</h3>
 *
 * In the examples shown in the introduction we had identified a
 * number of functions that can be run as independent tasks. Ideally,
 * this number of tasks is larger than the number of CPU cores (to
 * keep them busy) but is also not exceedingly huge (so as not to
 * inundate the scheduler with millions of tasks that will then have
 * to be distributed to 2 or 4 cores, for example). There are,
 * however, cases where we have many thousands or even millions of
 * relatively independent jobs: for example, assembling local
 * contributions to the global linear system on each cell of a mesh;
 * evaluating an error estimator on each cell; or postprocessing on
 * each cell computed data for output fall into this class. These
 * cases can be treated using a software design pattern we call
 * "%WorkStream". In the following, we will walk through the rationale
 * for this pattern and its implementation; more details as well as
 * examples for the speedup that can be achieved with it are given in
 * the @ref workstream_paper .
 *
 * Code like this could then be written like this:
 * @code
   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
   { ... }

   template <int dim>
   void MyClass<dim>::assemble_system ()
   {
     Threads::TaskGroup<void> task_group;
     for (typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active();
          cell != dof_handler.end(); ++cell)
       task_group += Threads::new_task (&MyClass<dim>::assemble_on_one_cell,
                                        *this,
                                        cell);
     task_group.join_all ();
   }
 * @endcode
 * On a big mesh, with maybe a million cells, this would create a massive
 * number of tasks; while it would keep all CPU cores busy for a while, the
 * overhead of first creating so many tasks, scheduling them, and then waiting
 * for them would probably not lead to efficient code. A better strategy would
 * be if the scheduler could somehow indicate that it has available resources,
 * at which point we would feed it another newly created task, and we would do
 * so until we run out of tasks and the ones that were created have been
 * worked on.
 *
 * This is essentially what the WorkStream::run function does: You give it an iterator
 * range from which it can draw objects to work on (in the above case it is
 * the interval given by <code>dof_handler.begin_active()</code> to
 * <code>dof_handler.end()</code>), and a function that would do the work on
 * each item (the function <code>MyClass::assemble_on_one_cell</code>)
 * together with an object if it is a member function.
 *
 * In the following, let us lay out a rationale for why the functions in the
 * WorkStream namespace are implemented the way they are. More information on
 * their implementation can be found in the @ref workstream_paper .
 * To see the WorkStream class used in practice on tasks like the ones
 * outlined above, take a look at the step-9, step-13, step-14, step-32, step-35 or step-37
 * tutorial programs.
 *
 * To begin with, given the brief description above,
 * the way the <code>MyClass::assemble_system</code>
 * function could then be written is like this (note that this is not quite
 * the correct syntax, as will be described below):
 * @code
   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
   { ... }

   template <int dim>
   void MyClass<dim>::assemble_system ()
   {
     WorkStream::run (dof_handler.begin_active(),
                      dof_handler.end(),
                      *this,
                      &MyClass<dim>::assemble_on_one_cell);
   }
 * @endcode
 *
 * There are at least three problems with this, however:
 *<ul>
 *<li>First, let us take a look at how the <code>MyClass::assemble_on_one_cell</code>
 *   function likely looks:
 * @code
   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
   {
     FEValues<dim> fe_values (...);
     FullMatrix<double> cell_matrix (...);
     Vector<double>     cell_rhs (...);

     // assemble local contributions
     fe_values.reinit (cell);
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
           cell_matrix(i,j) += ...;
     ...same for cell_rhs...

     // now copy results into global system
     std::vector<unsigned int> dof_indices (...);
     cell->get_dof_indices (dof_indices);
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         system_matrix.add (dof_indices[i], dof_indices[j],
                            cell_matrix(i,j));
     ...same for rhs...
   }
 * @endcode

 *   The problem here is that several tasks, each running
 *   <code>MyClass::assemble_on_one_cell</code>, could potentially try
 *   to write into the object <code>MyClass::system_matrix</code> <i>at
 *   the same time</i>. This could be avoided by explicit synchronisation
 *   using a Threads::Mutex, for example, and would look like this:
 * @code
     // now copy results into global system
     std::vector<unsigned int> dof_indices (...);
     cell->get_dof_indices (dof_indices);

     static Threads::Mutex mutex;
     mutex.acquire ();
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         system_matrix.add (dof_indices[i], dof_indices[j],
                            cell_matrix(i,j));
     ...same for rhs...
     mutex.release ();
   }
 * @endcode

 *   By making the mutex a static variable, it exists only once globally
 *   (i.e. once for all tasks that may be running in %parallel) and only one of
 *   the tasks can enter the region protected by the acquire/release calls on
 *   the mutex. As an aside, a better way to write this code would be like
 *   this, ensuring that the mutex is released even in case an exception is
 *   thrown, and without the need to remember to write the call to
 *   Threads::Mutex::release():
 * @code
     // now copy results into global system
     static Threads::Mutex mutex;
     Threads::Mutex::ScopedLock lock (mutex);
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         system_matrix.add (dof_indices[i], dof_indices[j],
                            cell_matrix(i,j));
     ...same for rhs...
   }
 * @endcode
 *   Here, the mutex remains locked from the time the ScopedLock is created to
 *   where it is destroyed, at the end of the code block.
 *
 *   Note that although we now avoid the race condition that multiple threads
 *   could be writing to the same object, this code is not very efficient:
 *   mutexes are expensive on multicore machines, and we also block threads
 *   some of the time which is inefficient with tasks as explained above in
 *   the section on
 *   @ref MTHow "How scheduling tasks works and when task-based programming is not efficient".
 *
 *<li>A second correctness problem is that even if we do lock the global matrix
 *   and right hand side objects using a mutex, we do so in a more or less
 *   random order: while tasks are created in the order in which we traverse
 *   cells normally, there is no guarantee that by the time we get to the
 *   point where we want to copy the local into the global contributions the
 *   order is still as if we computed things sequentially. In other words, it
 *   may happen that we add the contributions of cell 1 before those of cell
 *   0. That may seem harmless because addition is commutative and
 *   associative, but in fact it
 *   is not if done in floating point arithmetic: $a+b+c \neq a+c+b$ -- take
 *   for example $a=1, b=-1, c=10^{-20}$ (because $1+10^{-20}=1$ in floating
 *   point arithmetic, using double precision).
 *
 *   As a consequence, the exact values that end up in the global matrix and
 *   right hand side will be close but may differ by amounts close to
 *   round-off depending on the order in which tasks happened to finish their
 *   job. That's not a desirable outcome, since results will not be
 *   reproducible this way.
 *
 *   As a consequence, the way the WorkStream class is designed is to use two
 *   functions: the <code>MyClass::assemble_on_one_cell</code> computes the
 *   local contributions and stores them somewhere (we'll get to that next), and
 *   a second function, say <code>MyClass::copy_local_to_global</code>, that
 *   copies the results computed on each cell into the global objects. The
 *   trick implemented in the WorkStream class is that (i) the
 *   <code>MyClass::copy_local_to_global</code> never runs more than once in
 *   %parallel, so we do not need to synchronise execution through a mutex, and
 *   (ii) it runs in exactly the same order on cells as they appear in the
 *   iterator range, i.e. we add elements into the global matrix the same way
 *   <i>every time, independently of when the computation of these element
 *   finishes</i>.
 *
 *   We now only have to discuss how the
 *   <code>MyClass::assemble_on_one_cell</code> communicates to
 *   <code>MyClass::copy_local_to_global</code> what it has computed. The way
 *   this is done is to use an object that holds all temporary data:
 * @code
   struct PerTaskData {
     FullMatrix<double>        cell_matrix;
     Vector<double>            cell_rhs;
     std::vector<unsigned int> dof_indices;
   }

   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            PerTaskData &data)
   {
     FEValues<dim> fe_values (...);

     data.cell_matrix = 0;
     data.cell_rhs    = 0;

     // assemble local contributions
     fe_values.reinit (cell);
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
           data.cell_matrix(i,j) += ...;
     ...same for cell_rhs...

     cell->get_dof_indices (data.dof_indices);
   }

   template <int dim>
   void MyClass<dim>::copy_local_to_global (const PerTaskData &data)
   {
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
       for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
         system_matrix.add (data.dof_indices[i], data.dof_indices[j],
                            data.cell_matrix(i,j));
     ...same for rhs...
   }

   template <int dim>
   void MyClass<dim>::assemble_system ()
   {
     PerTaskData per_task_data;
     ...initialize members of per_task_data to the correct sizes...

     WorkStream::run (dof_handler.begin_active(),
                      dof_handler.end(),
                      *this,
                      &MyClass<dim>::assemble_on_one_cell,
                      &MyClass<dim>::copy_local_to_global,
                      per_task_data);
   }
 * @endcode
 *
 *   The way this works is that we create a sample <code>per_task_data</code>
 *   object that the work stream object will replicate once per task that runs
 *   in %parallel. For each task, this object will be passed first to one of
 *   possibly several instances of <code>MyClass::assemble_on_one_cell</code>
 *   running in %parallel which fills it with the data obtained on a single
 *   cell, and then to a sequentially running
 *   <code>MyClass::copy_local_to_global</code> that copies data into the
 *   global object. In practice, of course, we will not generate millions of
 *   <code>per_task_data</code> objects if we have millions of cells; rather,
 *   we recycle these objects after they have been used by
 *   <code>MyClass::copy_local_to_global</code> and feed them back into
 *   another instance of <code>MyClass::assemble_on_one_cell</code>; this
 *   means that the number of such objects we actually do create is a small
 *   multiple of the number of threads the scheduler uses, which is typically
 *   about as many as there are CPU cores on a system.
 *
 * <li>The last issue that is worth addressing is that the way we wrote the
 *   <code>MyClass::assemble_on_one_cell</code> function above, we create and
 *   destroy an FEValues object every time the function is called, i.e. once
 *   for each cell in the triangulation. That's an immensely expensive
 *   operation because the FEValues class tries to do a lot of work in its
 *   constructor in an attempt to reduce the number of operations we have to
 *   do on each cell (i.e. it increases the constant in the ${\cal O}(1)$
 *   effort to initialize such an object in order to reduce the constant in
 *   the ${\cal O}(N)$ operations to call FEValues::reinit on the $N$ cells of
 *   a triangulation). Creating and destroying an FEValues object on each cell
 *   invalidates this effort.
 *
 *   The way to avoid this is to put the FEValues object into a second
 *   structure that will hold scratch data, and initialize it in the
 *   constructor:
 * @code
   struct PerTaskData {
     FullMatrix<double>        cell_matrix;
     Vector<double>            cell_rhs;
     std::vector<unsigned int> dof_indices;

     PerTaskData (const FiniteElement<dim> &fe)
                :
                cell_matrix (fe.dofs_per_cell, fe.dofs_per_cell),
                cell_rhs (fe.dofs_per_cell),
                dof_indices (fe.dofs_per_cell)
       {}
   }

   struct ScratchData {
     FEValues<dim>             fe_values;

     ScratchData (const FiniteElement<dim> &fe,
                  const Quadrature<dim>    &quadrature,
                  const UpdateFlags         update_flags)
                :
                fe_values (fe, quadrature, update_flags)
       {}

     ScratchData (const ScratchData &scratch)
                :
                fe_values (scratch.fe_values.get_fe(),
                           scratch.fe_values.get_quadrature(),
                           scratch.fe_values.get_update_flags())
       {}
   }
 * @endcode
 * and then use this FEValues object in the assemble function:
 * @code
   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            ScratchData &scratch,
                                            PerTaskData &data)
   {
     scratch.fe_values.reinit (cell);
     ...
   }
 * @endcode
 *   Just as for the <code>PerTaskData</code> structure, we will create a
 *   sample <code>ScratchData</code> object and pass it to the work stream
 *   object, which will replicate it as many times as necessary. For this
 *   to work <code>ScratchData</code> structures need to copyable. Since FEValues
 *   objects are rather complex and cannot be copied implicitly, we provided
 *   our own copy constructor for the <code>ScratchData</code> structure.
 *
 *   The same approach, putting things into the <code>ScratchData</code>
 *   data structure, should be used for everything that is expensive to
 *   construct. This holds, in particular, for everything that needs to
 *   allocate memory upon construction; for example, if the values of a
 *   function need to be evaluated at quadrature points, then this is
 *   expensive:
 * @code
   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            ScratchData &scratch,
                                            PerTaskData &data)
   {
     std::vector<double> rhs_values (fe_values.n_quadrature_points);
     rhs_function.value_list (data.fe_values.get_quadrature_points,
                              rhs_values)
     ...
   }
 * @endcode
 * whereas this is a much cheaper way:
 * @code
   struct ScratchData {
     std::vector<double>       rhs_values;
     FEValues<dim>             fe_values;

     ScratchData (const FiniteElement<dim> &fe,
                  const Quadrature<dim>    &quadrature,
                  const UpdateFlags         update_flags)
                :
                rhs_values (quadrature.size()),
                fe_values (fe, quadrature, update_flags)
       {}

      ScratchData (const ScratchData &scratch)
                :
                rhs_values (scratch.rhs_values),
                fe_values (scratch.fe_values.get_fe(),
                           scratch.fe_values.get_quadrature(),
                           scratch.fe_values.get_update_flags())
       {}
   }

   template <int dim>
   void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            ScratchData &scratch,
                                            PerTaskData &data)
   {
     rhs_function.value_list (scratch.fe_values.get_quadrature_points,
                              scratch.rhs_values)
     ...
   }
 * @endcode
 *
 * </ul>
 *
 * As a final point: What if, for some reason, my assembler and copier
 * function do not match the above signature with three and one argument,
 * respectively? That's not a problem either. The WorkStream namespace offers two
 * versions of the WorkStream::run() function: one that takes an object and
 * the addresses of two member functions, and one that simply takes two
 * function objects that can be called with three and one argument,
 * respectively. So, in other words, the following two calls are exactly
 * identical:
 * @code
     WorkStream::run (dof_handler.begin_active(),
                      dof_handler.end(),
                      *this,
                      &MyClass<dim>::assemble_on_one_cell,
                      &MyClass<dim>::copy_local_to_global,
                      per_task_data);
     // ...is the same as:
     WorkStream::run (dof_handler.begin_active(),
                      dof_handler.end(),
                      std_cxx11::bind(&MyClass<dim>::assemble_on_one_cell,
                                      *this,
                                      std_cxx11::_1,
                                      std_cxx11::_2,
                                      std_cxx11::_3),
                      std_cxx11::bind(&MyClass<dim>::copy_local_to_global,
                                      *this,
                                      std_cxx11::_1),
                      per_task_data);
 * @endcode
 * Note how <code>std_cxx11::bind</code> produces a function object that takes three
 * arguments by binding the member function to the <code>*this</code>
 * object. <code>std_cxx11::_1, std_cxx11::_2</code> and <code>std_cxx11::_3</code> are placeholders for the first,
 * second and third argument that can be specified later on. In other words, for
 * example if <code>p</code> is the result of the first call to
 * <code>std_cxx11::bind</code>, then the call <code>p(cell, scratch_data,
 * per_task_data)</code> will result in executing
 * <code>this-@>assemble_on_one_cell (cell, scratch_data, per_task_data)</code>,
 * i.e. <code>std_cxx11::bind</code> has bound the object to the function pointer
 * but left the three arguments open for later.
 *
 * Similarly, let us assume that <code>MyClass::assemble_on_one_cell</code>
 * has the following signature in the solver of a nonlinear, time-dependent problem:
 * @code
   template <int dim>
   void
   MyClass<dim>::assemble_on_one_cell (const Vector<double> &linearization_point,
                                       const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       ScratchData &scratch,
                                       PerTaskData &data,
                                       const double current_time)
   { ... }
 * @endcode
 * Because WorkStream expects to be able to call the worker function with
 * just three arguments, the first of which is the iterator and the second
 * and third the ScratchData and PerTaskData objects, we need to pass the following
 * to it:
 * @code
     WorkStream::run (dof_handler.begin_active(),
                      dof_handler.end(),
                      std_cxx11::bind(&MyClass<dim>::assemble_on_one_cell,
                                      *this,
                                      current_solution,
                                      std_cxx11::_1,
                                      std_cxx11::_2,
                                      std_cxx11::_3,
                                      previous_time+time_step),
                      std_cxx11::bind(&MyClass<dim>::copy_local_to_global,
                                      *this,
                                      std_cxx11::_1),
                      per_task_data);
 * @endcode
 * Here, we bind the object, the linearization point argument, and the
 * current time argument to the function before we hand it off to
 * WorkStream::run(). WorkStream::run() will then simply call the
 * function with the cell and scratch and per task objects which will be filled
 * in at the positions indicated by <code>std_cxx11::_1, std_cxx11::_2</code>
 * and <code>std_cxx11::_3</code>.
 *
 * There are refinements to the WorkStream::run function shown above.
 * For example, one may realize that the basic idea above can only scale
 * if the copy-local-to-global function is much quicker than the
 * local assembly function because the former has to run sequentially.
 * This limitation can only be improved upon by scheduling more work
 * in parallel. This leads to the notion of coloring the graph of
 * cells (or, more generally, iterators) we work on by recording
 * which write operations conflict with each other. Consequently, there
 * is a third version of WorkStream::run that doesn't just take a
 * range of iterators, but instead a vector of vectors consisting of
 * elements that can be worked on at the same time. This concept
 * is explained in great detail in the @ref workstream_paper , along
 * with performance evaluations for common examples.
 *
 *
 * @anchor MTTaskSynchronization
 * <h3>Tasks and synchronization</h3>
 *
 * Tasks are powerful but they do have their limitation: to make
 * things efficient, the task scheduler never interrupts tasks by
 * itself. With the exception of the situation where one calls the
 * Threads::Task::join function to wait for another task to finish,
 * the task scheduler always runs a task to completion. The downside
 * is that the scheduler does not see if a task is actually idling,
 * for example if it waits for something else to happen (file IO to
 * finish, input from the keyboard, etc). In cases like this, the task
 * scheduler could in principle run a different task, but since it
 * doesn't know what tasks are doing it doesn't. Functions that do
 * wait for external events to happen are therefore not good
 * candidates for tasks and should use threads (see below).
 *
 * However, there are cases where tasks are not only a bad abstraction
 * for a job but can actually not be used: As a matter of principle,
 * tasks can not synchronize with other tasks through the use of a
 * mutex or a condition variable (see the Threads::Mutex and
 * Threads::ConditionVariable classes). The reason is that if task A
 * needs to wait for task B to finish something, then this is only
 * going to work if there is a guarantee that task B will eventually
 * be able to run and finish the task. Now imagine that you have 2
 * processors, and tasks A1 and A2 are currently running; let's assume
 * that they have queued tasks B1 and B2, and are now waiting with a
 * mutex for these queued tasks to finish (part of) their work. Since
 * the machine has only two processors, the task scheduler will only
 * start B1 or B2 once either A1 or A2 are done -- but this isn't
 * happening since they are waiting using operating system resources
 * (a mutex) rather than task scheduler resources. The result is a
 * deadlock.
 *
 * The bottom line is that tasks can not use mutices or condition variables to
 * synchronize with other tasks. If communication between tasks is necessary,
 * you need to use threads because the operating system makes sure that all
 * threads eventually get to run, independent of the total number of threads.
 * Note however that the same is not true if you only use a Thread::Mutex on
 * each task separately to protect access to a variable that the tasks may
 * write to: this use of mutices is ok; tasks may simply not want to wait for
 * another task to do something.
 *
 *
 * @anchor MTThreads
 * <h3>Thread-based parallelism</h3>
 *
 * Even though tasks are a higher-level way to describe things, there are
 * cases that are poorly suited to a task (for a discussion of some of
 * these cases see
 * @ref MTHow "How scheduling tasks works and when task-based programming is not efficient"
 * above). Generally, jobs that are not able to fully utilize the CPU are bad
 * fits for tasks and good fits for threads.
 *
 * In a case like this, you can resort to explicitly start threads, rather
 * than tasks, using pretty much the same syntax as above. For example, if you
 * had a function in your application that generates graphical output and then
 * estimates the error to refine the mesh for the next iteration of an
 * adaptive mesh scheme, it could look like this:
 * @code
   template <int dim>
   void MyClass<dim>::output_and_estimate_error () const
   {
     DataOut<dim> data_out;
     data_out.attach_dof_handler (dof_handler);
     data_out.add_data_vector (solution, "solution");
     data_out.build_patches ();

     std::ofstream output ("solution.vtk");

     Threads::Thread<void>
       thread = Threads::new_thread (&DataOut<dim>::write_vtk, data_out, output);

     Vector<float> error_per_cell (triangulation.n_active_cells());
     KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);
     thread.join ();
 * @endcode
 *
 * Here, Threads::new_thread starts the given function that writes to the
 * output file on a new thread that can run in %parallel to everything else:
 * In %parallel to the KellyErrorEstimator::estimate() function, the
 * DataOut::write_vtk() function will run on a separate thread. This execution
 * is independent of the scheduler that takes care of tasks, but that is
 * not a problem because writing lots of data to a file is not something that
 * will keep a CPU very busy.
 *
 * Creating threads works pretty much the same way as tasks, i.e. you can wait
 * for the termination of a thread using Threads::Thread::join(), query the
 * return value of a finished thread using Threads::Thread::return_value(),
 * and you can group threads into a Threads::ThreadGroup object and wait for
 * all of them to finish.
 *
 *
 * @anchor MTTaskThreads
 * <h3>Controlling the number of threads used for tasks</h3>
 * As mentioned earlier, deal.II does not implement scheduling tasks to
 * threads or even starting threads itself. The TBB library does a good job at
 * deciding how many threads to use and they do not recommend setting the
 * number of threads explicitly. However, on large symmetric multiprocessing
 * (SMP) machines, especially ones with a resource/job manager or on systems
 * on which access to some parts of the memory is possible but very expensive
 * for processors far away (e.g. large NUMA SMP machines), it may be necessary
 * to explicitly set the number of threads to prevent the TBB from using too
 * many CPUs. In this case the following commands can be used:
 * @code
   #include <tbb/task_scheduler_init.h>

   // In the program, before any task-based parallelism is reached.
   // Early in the main method is a good place to call this:
   tbb::task_scheduler_init init(n_desired_threads + 1);
 * @endcode
 * The method of setting the number of threads relies on this call to
 * <code>task_scheduler_init</code> occurring before any other calls to the
 * TBB functions.  If this call is first, it will set the number of threads
 * used by TBB for the remainder of the program.  Notice that the call starts
 * one less thread than the argument, since the main thread is counted. If
 * this method does not seem successful (gdb can be used to see when and how
 * many threads are started), ensure that the call is early enough in the
 * program (for example right at the start of <code>main()</code>), and check
 * that there are no other libraries starting threads.
 *
 * Note that a small number of places inside deal.II also uses thread-based
 * parallelism controlled by MultithreadInfo::n_default_threads
 * generally. Under some circumstances, deal.II also calls the BLAS library
 * which may sometimes also start threads of its own. You will have to consult
 * the documentation of your BLAS installation to determine how to set the
 * number of threads for these operations.
 */
