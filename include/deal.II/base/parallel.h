// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_parallel_h
#define dealii_parallel_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/std_cxx20/type_traits.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/template_constraints.h>

#include <cstddef>
#include <functional>
#include <memory>
#include <tuple>

#ifdef DEAL_II_WITH_TASKFLOW
#  include <deal.II/base/multithread_info.h>

#  include <taskflow/algorithm/for_each.hpp>
#  include <taskflow/taskflow.hpp>
#endif

#ifdef DEAL_II_WITH_TBB
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/partitioner.h>
#else
#  include <boost/range/iterator_range.hpp>
#endif

#ifdef DEAL_II_HAVE_CXX20
#  include <concepts>
#endif


// TODO[WB]: allow calling functions to pass along a tbb::affinity_partitioner
// object to ensure that subsequent calls use the same cache lines

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
    // Disable long double SIMD instructions on ICC. This is to work around a
    // bug that generates wrong code at least up to intel 15 (see
    // tests/lac/vector-vector, tests/lac/intel-15-bug, and the discussion at
    // https://github.com/dealii/dealii/issues/598).
    template <>
    struct EnableOpenMPSimdFor<long double>
    {
      static const bool value = false;
    };
#endif

#ifdef DEAL_II_WITH_TASKFLOW
    /**
     * Internal function to do a parallel_for using taskflow
     */
    template <typename Iterator, typename Functor>
    void
    taskflow_parallel_for(Iterator           x_begin,
                          Iterator           x_end,
                          const Functor     &functor,
                          const unsigned int grainsize)
    {
      tf::Executor &executor = MultithreadInfo::get_taskflow_executor();
      tf::Taskflow  taskflow;

      // TODO: We have several choices for the Partitioner and we should spend
      // some time benchmarking them:
      // 1. StaticPartitioner(grainsize): all work items have grainsize number
      // of items
      // 2. GuidedPartitioner(grainsize):
      // "The size of a partition is proportional to the number of unassigned
      // iterations divided by the number of workers, and the size will
      // gradually decrease to the given chunk size."
      // 3. GuidedPartitioner(0): The default.
      taskflow.for_each(
        x_begin,
        x_end,
        [&](const auto &item) { functor(item); },
        tf::StaticPartitioner(grainsize));
      executor.run(taskflow).wait();
    }
#endif

#ifdef DEAL_II_WITH_TBB
    /**
     * Encapsulate tbb::parallel_for.
     */
    template <typename Iterator, typename Functor>
    void
    parallel_for(Iterator           x_begin,
                 Iterator           x_end,
                 const Functor     &functor,
                 const unsigned int grainsize)
    {
      tbb::parallel_for(tbb::blocked_range<Iterator>(x_begin, x_end, grainsize),
                        functor,
                        tbb::auto_partitioner());
    }



    /**
     * Encapsulate tbb::parallel_for when an affinite_partitioner is provided.
     */
    template <typename Iterator, typename Functor>
    void
    parallel_for(Iterator                                          x_begin,
                 Iterator                                          x_end,
                 const Functor                                    &functor,
                 const unsigned int                                grainsize,
                 const std::shared_ptr<tbb::affinity_partitioner> &partitioner)
    {
      tbb::parallel_for(tbb::blocked_range<Iterator>(x_begin, x_end, grainsize),
                        functor,
                        *partitioner);
    }

#else

    /**
     * Just execute things sequentially.
     */
    template <typename Iterator, typename Functor>
    void
    parallel_for(Iterator       x_begin,
                 Iterator       x_end,
                 const Functor &functor,
                 const unsigned int)
    {
      functor(boost::iterator_range<Iterator>(x_begin, x_end));
    }

#endif
  } // namespace internal

  /**
   * An algorithm that performs the action <code>*out++ =
   * function(*in++)</code> where the <code>in</code> iterator ranges over
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
   * topic.
   *
   * @dealiiConceptRequires{(std::invocable<Function,
   *    decltype(*std::declval<InputIterator>())> &&
   *    std::assignable_from<decltype(*std::declval<OutputIterator>()),
   *    std::invoke_result_t<Function,
   * decltype(*std::declval<InputIterator>())>>)}
   */
  template <typename InputIterator, typename OutputIterator, typename Function>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Function, decltype(*std::declval<InputIterator>())> &&
     std::assignable_from<
       decltype(*std::declval<OutputIterator>()),
       std::invoke_result_t<Function,
                            decltype(*std::declval<InputIterator>())>>))
  DEAL_II_DEPRECATED_EARLY void transform(const InputIterator &begin_in,
                                          const InputIterator &end_in,
                                          OutputIterator       out,
                                          const Function      &function,
                                          const unsigned int   grainsize)
  {
#ifdef DEAL_II_WITH_TASKFLOW
    using Iterators     = std::tuple<InputIterator, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in, out);
    Iterators x_end(end_in, OutputIterator());

    internal::taskflow_parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &it) {
        *std::get<1>(it) = function(*std::get<0>(it));
      },
      grainsize);

#elif defined(DEAL_II_WITH_TBB)
    using Iterators     = std::tuple<InputIterator, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in, out);
    Iterators x_end(end_in, OutputIterator());
    internal::parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &range) {
        for (const auto &p : range)
          *std::get<1>(p) = function(*std::get<0>(p));
      },
      grainsize);
#else
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (InputIterator in = begin_in; in != end_in;)
      *out++ = function(*in++);
#endif
  }



  /**
   * An algorithm that performs the action <code>*out++ = function(*in1++,
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
   * topic.
   *
   * @dealiiConceptRequires{(std::invocable<Function,
   *    decltype(*std::declval<InputIterator1>()),
   *    decltype(*std::declval<InputIterator2>())> &&
   *    std::assignable_from<decltype(*std::declval<OutputIterator>()),
   *    std::invoke_result_t<Function,
   * decltype(*std::declval<InputIterator1>()),
   *    decltype(*std::declval<InputIterator2>())>>)}
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename OutputIterator,
            typename Function>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Function,
                    decltype(*std::declval<InputIterator1>()),
                    decltype(*std::declval<InputIterator2>())> &&
     std::assignable_from<
       decltype(*std::declval<OutputIterator>()),
       std::invoke_result_t<Function,
                            decltype(*std::declval<InputIterator1>()),
                            decltype(*std::declval<InputIterator2>())>>))
  DEAL_II_DEPRECATED_EARLY void transform(const InputIterator1 &begin_in1,
                                          const InputIterator1 &end_in1,
                                          InputIterator2        in2,
                                          OutputIterator        out,
                                          const Function       &function,
                                          const unsigned int    grainsize)
  {
#ifdef DEAL_II_WITH_TASKFLOW
    using Iterators =
      std::tuple<InputIterator1, InputIterator2, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, out);
    Iterators x_end(end_in1, InputIterator2(), OutputIterator());

    internal::taskflow_parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &it) {
        *std::get<2>(it) = function(*std::get<0>(it), *std::get<1>(it));
      },
      grainsize);

#elif defined(DEAL_II_WITH_TBB)
    using Iterators =
      std::tuple<InputIterator1, InputIterator2, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, out);
    Iterators x_end(end_in1, InputIterator2(), OutputIterator());
    internal::parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &range) {
        for (const auto &p : range)
          *std::get<2>(p) = function(*std::get<0>(p), *std::get<1>(p));
      },
      grainsize);

#else
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (InputIterator1 in1 = begin_in1; in1 != end_in1;)
      *out++ = function(*in1++, *in2++);
#endif
  }



  /**
   * An algorithm that performs the action <code>*out++ = function(*in1++,
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
   * topic.
   *
   * @dealiiConceptRequires{(std::invocable<Function,
   *    decltype(*std::declval<InputIterator1>()),
   *    decltype(*std::declval<InputIterator2>()),
   *    decltype(*std::declval<InputIterator3>())> &&
   *    std::assignable_from<decltype(*std::declval<OutputIterator>()),
   *    std::invoke_result_t<Function,
   * decltype(*std::declval<InputIterator1>()),
   *    decltype(*std::declval<InputIterator2>()),
   *    decltype(*std::declval<InputIterator3>())>>)}
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename InputIterator3,
            typename OutputIterator,
            typename Function>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Function,
                    decltype(*std::declval<InputIterator1>()),
                    decltype(*std::declval<InputIterator2>()),
                    decltype(*std::declval<InputIterator3>())> &&
     std::assignable_from<
       decltype(*std::declval<OutputIterator>()),
       std::invoke_result_t<Function,
                            decltype(*std::declval<InputIterator1>()),
                            decltype(*std::declval<InputIterator2>()),
                            decltype(*std::declval<InputIterator3>())>>))
  DEAL_II_DEPRECATED_EARLY void transform(const InputIterator1 &begin_in1,
                                          const InputIterator1 &end_in1,
                                          InputIterator2        in2,
                                          InputIterator3        in3,
                                          OutputIterator        out,
                                          const Function       &function,
                                          const unsigned int    grainsize)
  {
#ifdef DEAL_II_WITH_TASKFLOW
    using Iterators = std::
      tuple<InputIterator1, InputIterator2, InputIterator3, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, in3, out);
    Iterators x_end(end_in1,
                    InputIterator2(),
                    InputIterator3(),
                    OutputIterator());

    internal::taskflow_parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &it) {
        *std::get<3>(it) =
          function(*std::get<0>(it), *std::get<1>(it), *std::get<2>(it));
      },
      grainsize);

#elif defined(DEAL_II_WITH_TBB)
    using Iterators = std::
      tuple<InputIterator1, InputIterator2, InputIterator3, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, in3, out);
    Iterators x_end(end_in1,
                    InputIterator2(),
                    InputIterator3(),
                    OutputIterator());
    internal::parallel_for(
      SyncIterators(x_begin),
      SyncIterators(x_end),
      [function](const auto &range) {
        for (const auto &p : range)
          *std::get<3>(p) =
            function(*std::get<0>(p), *std::get<1>(p), *std::get<2>(p));
      },
      grainsize);
#else
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = function(*in1++, *in2++, *in3++);
#endif
  }


  namespace internal
  {
#ifdef DEAL_II_WITH_TBB
    /**
     * Take a range argument and call the given function with its begin and
     * end.
     *
     * @dealiiConceptRequires{(std::invocable<Function, Iterator, Iterator>)}
     */
    template <typename Iterator, typename Function>
    DEAL_II_CXX20_REQUIRES((std::invocable<Function, Iterator, Iterator>))
    void apply_to_subranges(const tbb::blocked_range<Iterator> &range,
                            const Function                     &f)
    {
      f(range.begin(), range.end());
    }
#endif
  } // namespace internal


  /**
   * This function applies the given function argument @p f to all elements in
   * the range <code>[begin,end)</code> and may do so in parallel. An example
   * of its use is given in step-69.
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
   *         [&](const unsigned int begin_row,
   *             const unsigned int end_row)
   *         {
   *           mat_vec_on_subranges(begin_row, end_row, A, x, y);
   *         },
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
   * Note how we use the lambda function to convert
   * <code>mat_vec_on_subranges</code> from a function that takes 5 arguments
   * to one taking 2 by binding the remaining arguments. The resulting function
   * object requires only two arguments, `begin_row` and `end_row`, with all
   * other arguments fixed.
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
   * topic.
   *
   * @dealiiConceptRequires{(std::invocable<Function, Iterator, Iterator>)}
   */
  template <typename Iterator, typename Function>
  DEAL_II_CXX20_REQUIRES((std::invocable<Function, Iterator, Iterator>))
  void apply_to_subranges(const Iterator                             &begin,
                          const std_cxx20::type_identity_t<Iterator> &end,
                          const Function                             &f,
                          const unsigned int                          grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    f(begin, end);
#else
    internal::parallel_for(
      begin,
      end,
      [&f](const tbb::blocked_range<Iterator> &range) {
        internal::apply_to_subranges<Iterator, Function>(range, f);
      },
      grainsize);
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
    virtual ~ParallelForInteger() = default;

    /**
     * This function runs the for loop over the given range
     * <tt>[lower,upper)</tt>, possibly in parallel when end-begin is larger
     * than the minimum parallel grain size. This function is marked const
     * because it any operation that changes the data of a derived class will
     * inherently not be thread-safe when several threads work with the same
     * data simultaneously.
     */
    void
    apply_parallel(const std::size_t begin,
                   const std::size_t end,
                   const std::size_t minimum_parallel_grain_size) const;

    /**
     * Virtual function for working on subrange to be defined in a derived
     * class.  This function is marked const because it any operation that
     * changes the data of a derived class will inherently not be thread-safe
     * when several threads work with the same data simultaneously.
     */
    virtual void
    apply_to_subrange(const std::size_t, const std::size_t) const = 0;
  };



  /**
   * This function works a lot like the apply_to_subranges() function, but it
   * allows to accumulate numerical results computed on each subrange into one
   * number. The type of this number is given by the `ResultType` template
   * argument that needs to be explicitly specified, and results are added
   * up (i.e., the reduction of results from subranges happens by adding up
   * these results).
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
   *        ([&](const unsigned int begin_row,
   *             const unsigned int end_row)
   *         {
   *           mat_norm_sqr_on_subranges(begin_row, end_row, A, x);
   *         },
   *         0, A.n_rows(),
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
   * topic.
   *
   * @dealiiConceptRequires{(std::invocable<Function, Iterator, Iterator> &&
   *    std::convertible_to<std::invoke_result_t<Function, Iterator, Iterator>,
   *    ResultType>)}
   */
  template <typename ResultType, typename Iterator, typename Function>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Function, Iterator, Iterator> &&
     std::convertible_to<std::invoke_result_t<Function, Iterator, Iterator>,
                         ResultType>))
  ResultType
    accumulate_from_subranges(const Function                             &f,
                              const Iterator                             &begin,
                              const std_cxx20::type_identity_t<Iterator> &end,
                              const unsigned int grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    return f(begin, end);
#else
    return tbb::parallel_reduce(
      tbb::blocked_range<Iterator>(begin, end, grainsize),
      ResultType(0),
      [f](const auto &range, const ResultType &starting_value) {
        ResultType value = starting_value;
        value += f(range.begin(), range.end());
        return value;
      },
      std::plus<ResultType>(),
      tbb::auto_partitioner());
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
      TBBPartitioner();

#ifdef DEAL_II_WITH_TBB
      /**
       * Destructor. Check that the object is not in use any more, i.e., all
       * loops have been completed.
       */
      ~TBBPartitioner();

      /**
       * Return an affinity partitioner. In case the partitioner owned by the
       * class is free, it is returned here. In case another thread has not
       * released it yet, a new object is created. To free the partitioner
       * again, return it by the release_one_partitioner() call.
       */
      std::shared_ptr<tbb::affinity_partitioner>
      acquire_one_partitioner();

      /**
       * After using the partitioner in a tbb loop through
       * acquire_one_partitioner(), this call makes the partitioner available
       * again.
       */
      void
      release_one_partitioner(
        const std::shared_ptr<tbb::affinity_partitioner> &p);

    private:
      /**
       * The stored partitioner that can accumulate knowledge over several
       * runs of tbb::parallel_for
       */
      std::shared_ptr<tbb::affinity_partitioner> my_partitioner;

      /**
       * A flag to indicate whether the partitioner has been acquired but not
       * released yet, i.e., it is in use somewhere else.
       */
      bool in_use;

      /**
       * A mutex to guard the access to the in_use flag.
       */
      Threads::Mutex mutex;
#endif
    };
  } // namespace internal
} // namespace parallel


namespace internal
{
  namespace VectorImplementation
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
  } // namespace VectorImplementation


  namespace SparseMatrixImplementation
  {
    /**
     * Like internal::VectorImplementation::minimum_parallel_grain_size, but now
     * denoting the number of rows of a matrix that should be worked on as a
     * minimum.
     */
    extern unsigned int minimum_parallel_grain_size;
  } // namespace SparseMatrixImplementation

} // end of namespace internal


/* --------------------------- inline functions ------------------------- */

namespace parallel
{
  inline void
  ParallelForInteger::apply_parallel(
    const std::size_t begin,
    const std::size_t end,
    const std::size_t minimum_parallel_grain_size) const
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)minimum_parallel_grain_size;

    apply_to_subrange(begin, end);
#else
    internal::parallel_for(
      begin,
      end,
      [this](const tbb::blocked_range<std::size_t> &range) {
        apply_to_subrange(range.begin(), range.end());
      },
      minimum_parallel_grain_size);
#endif
  }

} // end of namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
