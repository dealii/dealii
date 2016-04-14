// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2015 by the deal.II authors
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

#ifndef dealii__parallel_h
#define dealii__parallel_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/base/std_cxx11/bind.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <cstddef>

#ifdef DEAL_II_WITH_THREADS
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/partitioner.h>
#  include <tbb/blocked_range.h>
#endif


//TODO[WB]: allow calling functions to pass along a tbb::affinity_partitioner object to ensure that subsequent calls use the same cache lines

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace internal
  {
    /**
     * Helper struct to tell us if we can use SIMD instructions for the given
     * @p Number type.
     */
    template <typename Number>
    struct EnableOpenMPSimdFor
    {
      static const bool value = true;
    };

#ifdef __INTEL_COMPILER
    // Disable long double SIMD instructions on ICC. This is to work around a bug
    // that generates wrong code at least up to intel 15 (see
    // tests/lac/vector-vector, tests/lac/intel-15-bug, and the discussion at
    // https://github.com/dealii/dealii/issues/598).
    template <>
    struct EnableOpenMPSimdFor<long double>
    {
      static const bool value = false;
    };
#endif



    /**
     * Convert a function object of type F into an object that can be applied
     * to all elements of a range of synchronous iterators.
     */
    template <typename F>
    struct Body
    {
      /**
       * Constructor. Take and package the given function object.
       */
      Body (const F &f)
        :
        f (f)
      {}

      template <typename Range>
      void
      operator () (const Range &range) const
      {
        for (typename Range::const_iterator p=range.begin();
             p != range.end(); ++p)
          apply (f, p.iterators);
      }

    private:
      /**
       * The stored function object.
       */
      const F f;

      /**
       * Apply F to a set of iterators with two elements.
       */
      template <typename I1, typename I2>
      static
      void
      apply (const F &f,
             const std_cxx11::tuple<I1,I2> &p)
      {
        *std_cxx11::get<1>(p) = f (*std_cxx11::get<0>(p));
      }

      /**
       * Apply F to a set of iterators with three elements.
       */
      template <typename I1, typename I2, typename I3>
      static
      void
      apply (const F &f,
             const std_cxx11::tuple<I1,I2,I3> &p)
      {
        *std_cxx11::get<2>(p) = f (*std_cxx11::get<0>(p),
                                   *std_cxx11::get<1>(p));
      }

      /**
       * Apply F to a set of iterators with three elements.
       */
      template <typename I1, typename I2,
                typename I3, typename I4>
      static
      void
      apply (const F &f,
             const std_cxx11::tuple<I1,I2,I3,I4> &p)
      {
        *std_cxx11::get<3>(p) = f (*std_cxx11::get<0>(p),
                                   *std_cxx11::get<1>(p),
                                   *std_cxx11::get<2>(p));
      }
    };


    /**
     * Take a function object and create a Body object from it. We do this in
     * this helper function since alternatively we would have to specify the
     * actual data type of F -- which for function objects is often
     * extraordinarily complicated.
     */
    template <typename F>
    Body<F> make_body(const F &f)
    {
      return Body<F>(f);
    }
  }

  /**
   * An algorithm that performs the action <code>*out++ =
   * predicate(*in++)</code> where the <code>in</code> iterator ranges over
   * the given input range.
   *
   * This algorithm does pretty much what std::transform does. The difference
   * is that the function can run in parallel when deal.II is configured to
   * use multiple threads.
   *
   * If running in parallel, the iterator range is split into several chunks
   * that are each packaged up as a task and given to the Threading Building
   * Blocks scheduler to work on as compute resources are available. The
   * function returns once all chunks have been worked on. The last argument
   * denotes the minimum number of elements of the iterator range per task;
   * the number must be large enough to amortize the startup cost of new
   * tasks, and small enough to ensure that tasks can be reasonably load
   * balanced.
   *
   * For a discussion of the kind of problems to which this function is
   * applicable, see the
   * @ref threads "Parallel computing with multiple processors"
   * module.
   */
  template <typename InputIterator,
            typename OutputIterator,
            typename Predicate>
  void transform (const InputIterator &begin_in,
                  const InputIterator &end_in,
                  OutputIterator       out,
                  Predicate           &predicate,
                  const unsigned int   grainsize)
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in = begin_in; in != end_in;)
      *out++ = predicate (*in++);
#else
    typedef std_cxx11::tuple<InputIterator,OutputIterator> Iterators;
    typedef SynchronousIterators<Iterators> SyncIterators;
    Iterators x_begin (begin_in, out);
    Iterators x_end (end_in, OutputIterator());
    tbb::parallel_for (tbb::blocked_range<SyncIterators>(x_begin,
                                                         x_end,
                                                         grainsize),
                       internal::make_body (predicate),
                       tbb::auto_partitioner());
#endif
  }



  /**
   * An algorithm that performs the action <code>*out++ = predicate(*in1++,
   * *in2++)</code> where the <code>in1</code> iterator ranges over the given
   * input range, using the parallel for operator of tbb.
   *
   * This algorithm does pretty much what std::transform does. The difference
   * is that the function can run in parallel when deal.II is configured to
   * use multiple threads.
   *
   * If running in parallel, the iterator range is split into several chunks
   * that are each packaged up as a task and given to the Threading Building
   * Blocks scheduler to work on as compute resources are available. The
   * function returns once all chunks have been worked on. The last argument
   * denotes the minimum number of elements of the iterator range per task;
   * the number must be large enough to amortize the startup cost of new
   * tasks, and small enough to ensure that tasks can be reasonably load
   * balanced.
   *
   * For a discussion of the kind of problems to which this function is
   * applicable, see the
   * @ref threads "Parallel computing with multiple processors"
   * module.
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename OutputIterator,
            typename Predicate>
  void transform (const InputIterator1 &begin_in1,
                  const InputIterator1 &end_in1,
                  InputIterator2        in2,
                  OutputIterator        out,
                  Predicate            &predicate,
                  const unsigned int    grainsize)
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate (*in1++, *in2++);
#else
    typedef
    std_cxx11::tuple<InputIterator1,InputIterator2,OutputIterator>
    Iterators;
    typedef SynchronousIterators<Iterators> SyncIterators;
    Iterators x_begin (begin_in1, in2, out);
    Iterators x_end (end_in1, InputIterator2(), OutputIterator());
    tbb::parallel_for (tbb::blocked_range<SyncIterators>(x_begin,
                                                         x_end,
                                                         grainsize),
                       internal::make_body (predicate),
                       tbb::auto_partitioner());
#endif
  }



  /**
   * An algorithm that performs the action <code>*out++ = predicate(*in1++,
   * *in2++, *in3++)</code> where the <code>in1</code> iterator ranges over
   * the given input range.
   *
   * This algorithm does pretty much what std::transform does. The difference
   * is that the function can run in parallel when deal.II is configured to
   * use multiple threads.
   *
   * If running in parallel, the iterator range is split into several chunks
   * that are each packaged up as a task and given to the Threading Building
   * Blocks scheduler to work on as compute resources are available. The
   * function returns once all chunks have been worked on. The last argument
   * denotes the minimum number of elements of the iterator range per task;
   * the number must be large enough to amortize the startup cost of new
   * tasks, and small enough to ensure that tasks can be reasonably load
   * balanced.
   *
   * For a discussion of the kind of problems to which this function is
   * applicable, see the
   * @ref threads "Parallel computing with multiple processors"
   * module.
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename InputIterator3,
            typename OutputIterator,
            typename Predicate>
  void transform (const InputIterator1 &begin_in1,
                  const InputIterator1 &end_in1,
                  InputIterator2        in2,
                  InputIterator3        in3,
                  OutputIterator        out,
                  Predicate            &predicate,
                  const unsigned int    grainsize)
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate (*in1++, *in2++, *in3++);
#else
    typedef
    std_cxx11::tuple<InputIterator1,InputIterator2,InputIterator3,OutputIterator>
    Iterators;
    typedef SynchronousIterators<Iterators> SyncIterators;
    Iterators x_begin (begin_in1, in2, in3, out);
    Iterators x_end (end_in1, InputIterator2(),
                     InputIterator3(), OutputIterator());
    tbb::parallel_for (tbb::blocked_range<SyncIterators>(x_begin,
                                                         x_end,
                                                         grainsize),
                       internal::make_body (predicate),
                       tbb::auto_partitioner());
#endif
  }


  namespace internal
  {
#ifdef DEAL_II_WITH_THREADS
    /**
     * Take a range argument and call the given function with its begin and
     * end.
     */
    template <typename RangeType, typename Function>
    void apply_to_subranges (const tbb::blocked_range<RangeType> &range,
                             const Function  &f)
    {
      f (range.begin(), range.end());
    }
#endif
  }


  /**
   * This function applies the given function argument @p f to all elements in
   * the range <code>[begin,end)</code> and may do so in parallel.
   *
   * However, in many cases it is not efficient to call a function on each
   * element, so this function calls the given function object on sub-ranges.
   * In other words: if the given range <code>[begin,end)</code> is smaller
   * than grainsize or if multithreading is not enabled, then we call
   * <code>f(begin,end)</code>; otherwise, we may execute, possibly in
   * %parallel, a sequence of calls <code>f(b,e)</code> where
   * <code>[b,e)</code> are subintervals of <code>[begin,end)</code> and the
   * collection of calls we do to <code>f(.,.)</code> will happen on disjoint
   * subintervals that collectively cover the original interval
   * <code>[begin,end)</code>.
   *
   * Oftentimes, the called function will of course have to get additional
   * information, such as the object to work on for a given value of the
   * iterator argument. This can be achieved by <i>binding</i> certain
   * arguments. For example, here is an implementation of a matrix-vector
   * multiplication $y=Ax$ for a full matrix $A$ and vectors $x,y$:
   * @code
   *   void matrix_vector_product (const FullMatrix &A,
   *                               const Vector     &x,
   *                               Vector           &y)
   *   {
   *     parallel::apply_to_subranges
   *        (0, A.n_rows(),
   *         std_cxx11::bind (&mat_vec_on_subranges,
   *                          std_cxx11::_1, std_cxx11::_2,
   *                          std_cxx11::cref(A),
   *                          std_cxx11::cref(x),
   *                          std_cxx11::ref(y)),
   *         50);
   *   }
   *
   *   void mat_vec_on_subranges (const unsigned int begin_row,
   *                              const unsigned int end_row,
   *                              const FullMatrix &A,
   *                              const Vector     &x,
   *                              Vector           &y)
   *   {
   *     for (unsigned int row=begin_row; row!=end_row; ++row)
   *       for (unsigned int col=0; col<x.size(); ++col)
   *         y(row) += A(row,col) * x(col);
   *   }
   * @endcode
   *
   * Note how we use the <code>std_cxx11::bind</code> function to convert
   * <code>mat_vec_on_subranges</code> from a function that takes 5 arguments
   * to one taking 2 by binding the remaining arguments (the modifiers
   * <code>std_cxx11::ref</code> and <code>std_cxx11::cref</code> make sure
   * that the enclosed variables are actually passed by reference and constant
   * reference, rather than by value). The resulting function object requires
   * only two arguments, begin_row and end_row, with all other arguments
   * fixed.
   *
   * The code, if in single-thread mode, will call
   * <code>mat_vec_on_subranges</code> on the entire range
   * <code>[0,n_rows)</code> exactly once. In multi-threaded mode, however, it
   * may be called multiple times on subranges of this interval, possibly
   * allowing more than one CPU core to take care of part of the work.
   *
   * The @p grainsize argument (50 in the example above) makes sure that
   * subranges do not become too small, to avoid spending more time on
   * scheduling subranges to CPU resources than on doing actual work.
   *
   * For a discussion of the kind of problems to which this function is
   * applicable, see also the
   * @ref threads "Parallel computing with multiple processors"
   * module.
   */
  template <typename RangeType, typename Function>
  void apply_to_subranges (const RangeType                          &begin,
                           const typename identity<RangeType>::type &end,
                           const Function                           &f,
                           const unsigned int                        grainsize)
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) grainsize;

#  ifndef DEAL_II_BIND_NO_CONST_OP_PARENTHESES
    f (begin, end);
#  else
    // work around a problem with MS VC++ where there is no const
    // operator() in 'Function' if 'Function' is the result of std::bind
    Function ff = f;
    ff (begin, end);
#  endif
#else
    tbb::parallel_for (tbb::blocked_range<RangeType>
                       (begin, end, grainsize),
                       std_cxx11::bind (&internal::apply_to_subranges<RangeType,Function>,
                                        std_cxx11::_1,
                                        std_cxx11::cref(f)),
                       tbb::auto_partitioner());
#endif
  }



  /**
   * This is a class specialized to for loops with a fixed range given by
   * unsigned integers. This is an abstract base class that an actual worker
   * function is derived from. There is a public function apply that issues a
   * for loop in parallel, subdividing the work onto available processor cores
   * whenever there is enough work to be done (i.e., the number of elements is
   * larger than grain_size). Inside the function, a virtual function
   * apply_to_subrange specifying a range of two integers <tt>[lower,
   * upper)</tt> is called which needs to be defined in a derived class.
   *
   * The parallelization cases covered by this class are a subset of what is
   * possible with the function apply_to_subranges (which also covers the case
   * of more general iterators that might not be described by an integer
   * range). However, for simple integer ranges one might prefer this class,
   * like when there are many structurally similar loops, e.g., some simple
   * copy or arithmetic operations on an array of pointers. In that case,
   * apply_to_subranges will generate a lot of code (or rather, a lot of
   * symbols) because it passes the long names generated by std::bind to the
   * templated parallel for functions in TBB. This can considerably increase
   * compile times and the size of the object code. Similarly, the incorrect
   * use of std::bind often results in very cryptic error messages, which can
   * be avoided by this class (only a virtual function needs to be defined in
   * a derived class). Finally, the additional cost of a virtual function is
   * negligible in the context of parallel functions: It is much more
   * expensive to actually issue the work onto a thread, which in turn should
   * be much less than the actual work done in the for loop.
   */
  struct ParallelForInteger
  {
    /**
     * Destructor. Made virtual to ensure that derived classes also have
     * virtual destructors.
     */
    virtual ~ParallelForInteger ();

    /**
     * This function runs the for loop over the given range
     * <tt>[lower,upper)</tt>, possibly in parallel when end-begin is larger
     * than the minimum parallel grain size. This function is marked const
     * because it any operation that changes the data of a derived class will
     * inherently not be thread-safe when several threads work with the same
     * data simultaneously.
     */
    void apply_parallel (const std::size_t begin,
                         const std::size_t end,
                         const std::size_t minimum_parallel_grain_size) const;

    /**
     * Virtual function for working on subrange to be defined in a derived
     * class.  This function is marked const because it any operation that
     * changes the data of a derived class will inherently not be thread-safe
     * when several threads work with the same data simultaneously.
     */
    virtual void apply_to_subrange (const std::size_t,
                                    const std::size_t) const = 0;
  };



  namespace internal
  {
#ifdef DEAL_II_WITH_THREADS
    /**
     * A class that conforms to the Body requirements of the TBB
     * parallel_reduce function. The first template argument denotes the type
     * on which the reduction is to be done. The second denotes the type of
     * the function object that shall be called for each subrange.
     */
    template <typename ResultType,
              typename Function>
    struct ReductionOnSubranges
    {
      /**
       * A variable that will hold the result of the reduction.
       */
      ResultType result;

      /**
       * Constructor. Take the function object to call on each sub-range as
       * well as the neutral element with respect to the reduction operation.
       *
       * The second argument denotes a function object that will be used to
       * reduce the result of two computations into one number. An example if
       * we want to simply accumulate integer results would be
       * std::plus<int>().
       */
      template <typename Reductor>
      ReductionOnSubranges (const Function &f,
                            const Reductor &reductor,
                            const ResultType neutral_element = ResultType())
        :
        result (neutral_element),
        f (f),
        neutral_element (neutral_element),
        reductor (reductor)
      {}

      /**
       * Splitting constructor. See the TBB book for more details about this.
       */
      ReductionOnSubranges (const ReductionOnSubranges &r,
                            tbb::split)
        :
        result (r.neutral_element),
        f (r.f),
        neutral_element (r.neutral_element),
        reductor (r.reductor)
      {}

      /**
       * Join operation: merge the results from computations on different sub-
       * intervals.
       */
      void join (const ReductionOnSubranges &r)
      {
        result = reductor(result, r.result);
      }

      /**
       * Execute the given function on the specified range.
       */
      template <typename RangeType>
      void operator () (const tbb::blocked_range<RangeType> &range)
      {
        result = reductor(result,
                          f (range.begin(), range.end()));
      }

    private:
      /**
       * The function object to call on every sub-range.
       */
      const Function f;

      /**
       * The neutral element with respect to the reduction operation. This is
       * needed when calling the splitting constructor since we have to re-set
       * the result variable in this case.
       */
      const ResultType neutral_element;

      /**
       * The function object to be used to reduce the result of two calls into
       * one number.
       */
      const std_cxx11::function<ResultType (ResultType, ResultType)> reductor;
    };
#endif
  }


  /**
   * This function works a lot like the apply_to_subranges(), but it allows to
   * accumulate numerical results computed on each subrange into one number.
   * The type of this number is given by the ResultType template argument that
   * needs to be explicitly specified.
   *
   * An example of use of this function is to compute the value of the
   * expression $x^T A x$ for a square matrix $A$ and a vector $x$. The sum
   * over rows can be parallelized and the whole code might look like this:
   * @code
   *   void matrix_norm (const FullMatrix &A,
   *                     const Vector     &x)
   *   {
   *     return
   *      std::sqrt
   *       (parallel::accumulate_from_subranges<double>
   *        (0, A.n_rows(),
   *         std_cxx11::bind (&mat_norm_sqr_on_subranges,
   *                          std_cxx11::_1, std_cxx11::_2,
   *                          std_cxx11::cref(A),
   *                          std_cxx11::cref(x)),
   *         50);
   *   }
   *
   *   double
   *   mat_norm_sqr_on_subranges (const unsigned int begin_row,
   *                              const unsigned int end_row,
   *                              const FullMatrix &A,
   *                              const Vector     &x)
   *   {
   *     double norm_sqr = 0;
   *     for (unsigned int row=begin_row; row!=end_row; ++row)
   *       for (unsigned int col=0; col<x.size(); ++col)
   *         norm_sqr += x(row) * A(row,col) * x(col);
   *     return norm_sqr;
   *   }
   * @endcode
   *
   * Here, <code>mat_norm_sqr_on_subranges</code> is called on the entire
   * range <code>[0,A.n_rows())</code> if this range is less than the minimum
   * grainsize (above chosen as 50) or if deal.II is configured to not use
   * multithreading. Otherwise, it may be called on subsets of the given
   * range, with results from the individual subranges accumulated internally.
   *
   * @warning If ResultType is a floating point type, then accumulation is not
   * an associative operation. In other words, if the given function object is
   * called three times on three subranges, returning values $a,b,c$, then the
   * returned result of this function is $(a+b)+c$. However, depending on how
   * the three sub-tasks are distributed on available CPU resources, the
   * result may also be $(a+c)+b$ or any other permutation; because floating
   * point addition is not associative (as opposed, of course, to addition of
   * real %numbers), the result of invoking this function several times may
   * differ on the order of round-off.
   *
   * For a discussion of the kind of problems to which this function is
   * applicable, see also the
   * @ref threads "Parallel computing with multiple processors"
   * module.
   */
  template <typename ResultType, typename RangeType, typename Function>
  ResultType accumulate_from_subranges (const Function &f,
                                        const RangeType                          &begin,
                                        const typename identity<RangeType>::type &end,
                                        const unsigned int grainsize)
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) grainsize;

#  ifndef DEAL_II_BIND_NO_CONST_OP_PARENTHESES
    return f (begin, end);
#  else
    // work around a problem with MS VC++ where there is no const
    // operator() in 'Function' if 'Function' is the result of std::bind
    Function ff = f;
    return ff (begin, end);
#  endif
#else
    internal::ReductionOnSubranges<ResultType,Function>
    reductor (f, std::plus<ResultType>(), 0);
    tbb::parallel_reduce (tbb::blocked_range<RangeType>(begin, end, grainsize),
                          reductor,
                          tbb::auto_partitioner());
    return reductor.result;
#endif
  }


// --------------------- for loop affinity partitioner -----------------------

  /**
   * A class that wraps a TBB affinity partitioner in a thread-safe way. In
   * Vector, we use a shared pointer to share an affinity partitioner
   * between different vectors of the same size for improving data (and
   * NUMA) locality. However, when an outer task does multiple vector
   * operations, the shared pointer could lead to race conditions. This
   * class only allows one instance to get a partitioner. The other objects
   * cannot use that object and need to create their own copy.
   */
  namespace internal
  {
    class TBBPartitioner
    {
    public:
      /**
       * Constructor.
       */
      TBBPartitioner()
#ifdef DEAL_II_WITH_THREADS
        :
        my_partitioner(new tbb::affinity_partitioner()),
        in_use(false)
#endif
      {}

#ifdef DEAL_II_WITH_THREADS
      /**
       * Destructor. Check that the object is not in use any more, i.e., all
       * loops have been completed.
       */
      ~TBBPartitioner()
      {
        Assert(in_use == false,
               ExcInternalError("A vector partitioner goes out of scope, but "
                                "it appears to be still in use."));
      }

      /**
       * Return an affinity partitioner. In case the partitioner owned by the
       * class is free, it is returned here. In case another thread has not
       * released it yet, a new object is created. To free the partitioner
       * again, return it by the release_one_partitioner() call.
       */
      std_cxx11::shared_ptr<tbb::affinity_partitioner>
      acquire_one_partitioner()
      {
        dealii::Threads::Mutex::ScopedLock lock(mutex);
        if (in_use)
          return std_cxx11::shared_ptr<tbb::affinity_partitioner>(new tbb::affinity_partitioner());

        in_use = true;
        return my_partitioner;
      }

      /**
       * After using the partitioner in a tbb loop through
       * acquire_one_partitioner(), this call makes the partitioner available
       * again.
       */
      void release_one_partitioner(std_cxx11::shared_ptr<tbb::affinity_partitioner> &p)
      {
        if (p.get() == my_partitioner.get())
          {
            dealii::Threads::Mutex::ScopedLock lock(mutex);
            in_use = false;
          }
      }

    private:
      /**
       * The stored partitioner that can accumulate knowledge over several
       * runs of tbb::parallel_for
       */
      std_cxx11::shared_ptr<tbb::affinity_partitioner> my_partitioner;

      /**
       * A flag to indicate whether the partitioner has been acquired but not
       * released yet, i.e., it is in use somewhere else.
       */
      bool in_use;

      /**
       * A mutex to guard the access to the in_use flag.
       */
      dealii::Threads::Mutex mutex;
#endif
    };
  }
}


namespace internal
{
  namespace Vector
  {
    /**
     * If we do computations on vectors in parallel (say, we add two vectors
     * to get a third, and we do the loop over all elements in parallel), then
     * this variable determines the minimum number of elements for which it is
     * profitable to split a range of elements any further to distribute to
     * different threads.
     *
     * This variable is available as a global writable variable in order to
     * allow the testsuite to also test the parallel case. By default, it is
     * set to several thousand elements, which is a case that the testsuite
     * would not normally encounter. As a consequence, in the testsuite we set
     * it to one -- a value that's hugely unprofitable but definitely tests
     * parallel operations.
     */
    extern unsigned int minimum_parallel_grain_size;
  }


  namespace SparseMatrix
  {
    /**
     * Like internal::Vector::minimum_parallel_grain_size, but now denoting
     * the number of rows of a matrix that should be worked on as a minimum.
     */
    extern unsigned int minimum_parallel_grain_size;
  }

} // end of namespace internal


/* --------------------------- inline functions ------------------------- */

namespace parallel
{

#ifdef DEAL_II_WITH_THREADS

  namespace internal
  {
    /**
     * This is the function actually called by TBB for the ParallelForInteger
     * class.
     */
    struct ParallelForWrapper
    {
      ParallelForWrapper (const parallel::ParallelForInteger &worker)
        :
        worker_ (worker)
      {}

      void operator() (const tbb::blocked_range<std::size_t> &range) const
      {
        worker_.apply_to_subrange (range.begin(), range.end());
      }

      const parallel::ParallelForInteger &worker_;
    };
  }

#endif


  inline
  ParallelForInteger::~ParallelForInteger ()
  {}


  inline
  void
  ParallelForInteger::apply_parallel (const std::size_t begin,
                                      const std::size_t end,
                                      const std::size_t minimum_parallel_grain_size) const
  {
#ifndef DEAL_II_WITH_THREADS
    // make sure we don't get compiler
    // warnings about unused arguments
    (void) minimum_parallel_grain_size;

    apply_to_subrange (begin, end);
#else
    internal::ParallelForWrapper worker(*this);
    tbb::parallel_for (tbb::blocked_range<std::size_t>
                       (begin, end, minimum_parallel_grain_size),
                       worker,
                       tbb::auto_partitioner());
#endif
  }

} // end of namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
