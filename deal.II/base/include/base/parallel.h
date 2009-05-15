//---------------------------------------------------------------------------
//    $Id: parallel.h 14038 2006-10-23 02:46:34Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__parallel_h
#define __deal2__parallel_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/template_constraints.h>

#include <base/std_cxx1x/tuple.h>
#include <base/std_cxx1x/bind.h>
#include <base/std_cxx1x/function.h>

#include <iterator>
#include <cstddef>

#if DEAL_II_USE_MT == 1
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/partitioner.h>
#  include <tbb/blocked_range.h>
#endif


//TODO[WB]: allow calling functions to pass along a tbb::affinity_partitioner object to ensure that subsequent calls use the same cache lines

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which we define a few algorithms that can run in parallel
 * when deal.II is configured to use multiple threads.
 *
 * @ingroup threads
 * @author Wolfgang Bangerth, 2008, 2009
 */
namespace parallel
{
  namespace internal
  {
				     /**
				      * A class that represents a set of
				      * iterators each of which are
				      * incremented by one at the same
				      * time. This is typically used in calls
				      * like <code>std::transform(a.begin(),
				      * a.end(), b.begin(), functor);</code>
				      * where we have synchronous iterators
				      * marching through the containers
				      * <code>a,b</code>. If an object of this
				      * type represents the end of a range,
				      * only the first element is considered
				      * (we only have <code>a.end()</code>,
				      * not <code>b.end()</code>)
				      *
				      * The template argument of the current
				      * class shall be of type
				      * <code>std_cxx1x::tuple</code> with
				      * arguments equal to the iterator types.
				      *
				      * This type, and the helper functions
				      * associated with it, are used as the
				      * Value concept for the blocked_range
				      * type of the Threading Building Blocks.
				      */
    template <typename Iterators>
    struct SynchronousIterators
    {
					 /**
					  * Constructor.
					  */
	SynchronousIterators (const Iterators &i);

					 /**
					  * Copy constructor.
					  */
	SynchronousIterators (const SynchronousIterators &i);

					 /**
					  * Storage for the iterators
					  * represented by the current class.
					  */
	Iterators iterators;
    };



    template <typename Iterators>
    inline
    SynchronousIterators<Iterators>::
    SynchronousIterators (const Iterators &i)
		    :
		    iterators (i)
    {}
        

    template <typename Iterators>
    inline
    SynchronousIterators<Iterators>::
    SynchronousIterators (const SynchronousIterators &i)
		    :
		    iterators (i.iterators)
    {}



				     /**
				      * Return whether the first element of
				      * the first argument is less than the
				      * first element of the second
				      * argument. Since the objects compared
				      * march forward all elements at the same
				      * time, comparing the first element is
				      * sufficient.
				      */
    template <typename Iterators>
    inline
    bool
    operator< (const SynchronousIterators<Iterators> &a,
	       const SynchronousIterators<Iterators> &b)
    {
      return std_cxx1x::get<0>(a.iterators) < std_cxx1x::get<0>(b.iterators);
    }



				     /**
				      * Return the distance between the first
				      * and the second argument. Since the
				      * objects compared march forward all
				      * elements at the same time,
				      * differencing the first element is
				      * sufficient.
				      */
    template <typename Iterators>
    inline
    std::size_t
    operator- (const SynchronousIterators<Iterators> &a,
	       const SynchronousIterators<Iterators> &b)
    {
      Assert (std::distance (std_cxx1x::get<0>(b.iterators),
			     std_cxx1x::get<0>(a.iterators)) >= 0,
	      ExcInternalError());
      return std::distance (std_cxx1x::get<0>(b.iterators),
			    std_cxx1x::get<0>(a.iterators));
    }


				     /**
				      * Advance a tuple of iterators by $n$.
				      */
    template <typename I1, typename I2>
    inline
    void advance (std_cxx1x::tuple<I1,I2> &t,
		  const unsigned int       n)
    {
      std::advance (std_cxx1x::get<0>(t), n);
      std::advance (std_cxx1x::get<1>(t), n);
    }

				     /**
				      * Advance a tuple of iterators by $n$.
				      */
    template <typename I1, typename I2, typename I3>
    inline
    void advance (std_cxx1x::tuple<I1,I2,I3> &t,
		  const unsigned int          n)
    {
      std::advance (std_cxx1x::get<0>(t), n);
      std::advance (std_cxx1x::get<1>(t), n);
      std::advance (std_cxx1x::get<2>(t), n);
    }

				     /**
				      * Advance a tuple of iterators by $n$.
				      */
    template <typename I1, typename I2,
	      typename I3, typename I4>
    inline
    void advance (std_cxx1x::tuple<I1,I2,I3, I4> &t,
		  const unsigned int              n)
    {
      std::advance (std_cxx1x::get<0>(t), n);
      std::advance (std_cxx1x::get<1>(t), n);
      std::advance (std_cxx1x::get<2>(t), n);
      std::advance (std_cxx1x::get<3>(t), n);
    }
    


				     /**
				      * Advance a tuple of iterators by 1.
				      */
    template <typename I1, typename I2>
    inline
    void advance_by_one (std_cxx1x::tuple<I1,I2> &t)
    {
      ++std_cxx1x::get<0>(t);
      ++std_cxx1x::get<1>(t);
    }

				     /**
				      * Advance a tuple of iterators by 1.
				      */
    template <typename I1, typename I2, typename I3>
    inline
    void advance_by_one (std_cxx1x::tuple<I1,I2,I3> &t)
    {
      ++std_cxx1x::get<0>(t);
      ++std_cxx1x::get<1>(t);
      ++std_cxx1x::get<2>(t);
    }
    
				     /**
				      * Advance a tuple of iterators by 1.
				      */
    template <typename I1, typename I2,
	      typename I3, typename I4>
    inline
    void advance_by_one (std_cxx1x::tuple<I1,I2,I3,I4> &t)
    {
      ++std_cxx1x::get<0>(t);
      ++std_cxx1x::get<1>(t);
      ++std_cxx1x::get<2>(t);
      ++std_cxx1x::get<3>(t);
    }
    


				     /**
				      * Advance the elements of this iterator
				      * by $n$.
				      */
    template <typename Iterators>
    inline
    SynchronousIterators<Iterators>
    operator + (const SynchronousIterators<Iterators> &a,
		const std::size_t                      n)
    {
      SynchronousIterators<Iterators> x (a);
      parallel::internal::advance (x.iterators, n);
      return x;
    }

				     /**
				      * Advance the elements of this iterator
				      * by 1.
				      */
    template <typename Iterators>
    inline
    SynchronousIterators<Iterators>
    operator ++ (SynchronousIterators<Iterators> &a)
    {
      parallel::internal::advance_by_one (a.iterators);
      return a;
    }


				     /**
				      * Compare synch iterators for
				      * inequality. Since they march in synch,
				      * comparing only the first element is
				      * sufficient.
				      */
    template <typename Iterators>
    inline
    bool
    operator != (const SynchronousIterators<Iterators> &a,
		 const SynchronousIterators<Iterators> &b)
    {
      return (std_cxx1x::get<0>(a.iterators) !=
	      std_cxx1x::get<0>(b.iterators));
    }
    

				     /**
				      * Convert a function object of type F
				      * into an object that can be applied to
				      * all elements of a range of synchronous
				      * iterators.
				      */
    template <typename F>
    struct Body
    {
					 /**
					  * Constructor. Take and package the
					  * given function object.
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
					  * Apply F to a set of iterators with
					  * two elements.
					  */
	template <typename I1, typename I2>
	static
	void
	apply (const F &f,
	       const std_cxx1x::tuple<I1,I2> &p)
	  {
	    *std_cxx1x::get<1>(p) = f (*std_cxx1x::get<0>(p));
	  }

					 /**
					  * Apply F to a set of iterators with
					  * three elements.
					  */
	template <typename I1, typename I2, typename I3>
	static
	void
	apply (const F &f,
	       const std_cxx1x::tuple<I1,I2,I3> &p)
	  {
	    *std_cxx1x::get<2>(p) = f (*std_cxx1x::get<0>(p),
				       *std_cxx1x::get<1>(p));
	  }

					 /**
					  * Apply F to a set of iterators with
					  * three elements.
					  */
	template <typename I1, typename I2,
		  typename I3, typename I4>
	static
	void
	apply (const F &f,
	       const std_cxx1x::tuple<I1,I2,I3,I4> &p)
	  {
	    *std_cxx1x::get<3>(p) = f (*std_cxx1x::get<0>(p),
				       *std_cxx1x::get<1>(p),
				       *std_cxx1x::get<2>(p));
	  }
    };


				     /**
				      * Take a function object and create a
				      * Body object from it. We do this in
				      * this helper function since
				      * alternatively we would have to specify
				      * the actual data type of F -- which for
				      * function objects is often
				      * extraordinarily complicated.
				      */
    template <typename F>
    Body<F> make_body(const F &f)
    {
      return Body<F>(f);
    }
  }
  
				   /**
				    * An algorithm that performs the action
				    * <code>*out++ = predicate(*in++)</code>
				    * where the <code>in</code> iterator
				    * ranges over the given input range.
				    *
				    * This algorithm does pretty much what
				    * std::transform does. The difference is
				    * that the function can run in parallel
				    * when deal.II is configured to use
				    * multiple threads. In that case, the last
				    * argument denotes the minimum number each
				    * thread can work on; the number must be
				    * large enough to amortize the startup
				    * cost of new threads, and small enough to
				    * ensure that (i) threads will be started
				    * at all, and (ii) threads can be
				    * reasonably load balanced.
				    *
				    * For a discussion of the kind of
				    * problems to which this function
				    * is applicable, see the
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
#if DEAL_II_USE_MT == 0
				     // make sure we don't get compiler
				     // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in = begin_in; in != end_in;)
      *out++ = predicate (*in++);
#else
    typedef std_cxx1x::tuple<InputIterator,OutputIterator> Iterators;
    typedef internal::SynchronousIterators<Iterators> SyncIterators;
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
				    * An algorithm that performs the action
				    * <code>*out++ = predicate(*in1++, *in2++)</code>
				    * where the <code>in1</code> iterator
				    * ranges over the given input range.
				    *
				    * This algorithm does pretty much what
				    * std::transform does. The difference is
				    * that the function can run in parallel
				    * when deal.II is configured to use
				    * multiple threads. In that case, the last
				    * argument denotes the minimum number each
				    * thread can work on; the number must be
				    * large enough to amortize the startup
				    * cost of new threads, and small enough to
				    * ensure that (i) threads will be started
				    * at all, and (ii) threads can be
				    * reasonably load balanced.
				    *
				    * For a discussion of the kind of
				    * problems to which this function
				    * is applicable, see the
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
#if DEAL_II_USE_MT == 0
				     // make sure we don't get compiler
				     // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate (*in1++, *in2++);
#else
    typedef
      std_cxx1x::tuple<InputIterator1,InputIterator2,OutputIterator>
      Iterators;
    typedef internal::SynchronousIterators<Iterators> SyncIterators;
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
				    * An algorithm that performs the action
				    * <code>*out++ = predicate(*in1++, *in2++, *in3++)</code>
				    * where the <code>in1</code> iterator
				    * ranges over the given input range.
				    *
				    * This algorithm does pretty much what
				    * std::transform does. The difference is
				    * that the function can run in parallel
				    * when deal.II is configured to use
				    * multiple threads. In that case, the last
				    * argument denotes the minimum number each
				    * thread can work on; the number must be
				    * large enough to amortize the startup
				    * cost of new threads, and small enough to
				    * ensure that (i) threads will be started
				    * at all, and (ii) threads can be
				    * reasonably load balanced.
				    *
				    * For a discussion of the kind of
				    * problems to which this function
				    * is applicable, see the
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
#if DEAL_II_USE_MT == 0
				     // make sure we don't get compiler
				     // warnings about unused arguments
    (void) grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate (*in1++, *in2++, *in3++);
#else
    typedef
      std_cxx1x::tuple<InputIterator1,InputIterator2,InputIterator3,OutputIterator>
      Iterators;
    typedef internal::SynchronousIterators<Iterators> SyncIterators;
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
#if DEAL_II_USE_MT
				     /**
				      * Take a range argument and call the
				      * given function with its begin and end.
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
				    * This function applies the given function
				    * argument @p f to all elements in the range
				    * <code>[begin,end)</code> and may do so
				    * in parallel.
				    *
				    * However, in many cases it is not
				    * efficient to call a function on each
				    * element, so this function calls the
				    * given function object on sub-ranges. In
				    * other words: if the given range
				    * <code>[begin,end)</code> is smaller than
				    * grainsize or if multithreading is not
				    * enabled, then we call
				    * <code>f(begin,end)</code>; otherwise, we
				    * may execute, possibly in %parallel, a
				    * sequence of calls <code>f(b,e)</code>
				    * where <code>[b,e)</code> are
				    * subintervals of <code>[begin,end)</code>
				    * and the collection of calls we do to
				    * <code>f(.,.)</code> will happen on
				    * disjoint subintervals that collectively
				    * cover the original interval
				    * <code>[begin,end)</code>.
				    *
				    * Oftentimes, the called function will of
				    * course have to get additional
				    * information, such as the object to work
				    * on for a given value of the iterator
				    * argument. This can be achieved by
				    * <i>binding</i> certain arguments. For
				    * example, here is an implementation of a
				    * matrix-vector multiplication $y=Ax$ for
				    * a full matrix $A$ and vectors $x,y$:
				    * @code
				    *   void matrix_vector_product (const FullMatrix &A,
				    *                               const Vector     &x,
				    *                               Vector           &y)
				    *   {
				    *     parallel::apply_to_subranges
				    *        (0, A.n_rows(),
				    *         std_cxx1x::bind (&mat_vec_on_subranges,
				    *                          _1, _2,
				    *                          std_cxx1x::cref(A),
				    *                          std_cxx1x::cref(x),
				    *                          std_cxx1x::ref(y)),
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
				    * Note how we use the
				    * <code>std_cxx1x::bind</code> function to
				    * convert
				    * <code>mat_vec_on_subranged</code> from a
				    * function that takes 5 arguments to one
				    * taking 2 by binding the remaining
				    * arguments (the modifiers
				    * <code>std_cxx1x::ref</code> and
				    * <code>std_cxx1x::cref</code> make sure
				    * that the enclosed variables are actually
				    * passed by reference and constant
				    * reference, rather than by value). The
				    * resulting function object requires only
				    * two arguments, begin_row and end_row,
				    * with all other arguments fixed.
				    *
				    * The code, if in single-thread mode, will
				    * call <code>mat_vec_on_subranges</code>
				    * on the entire range
				    * <code>[0,n_rows)</code> exactly once. In
				    * multi-threaded mode, however, it may be
				    * called multiple times on subranges of
				    * this interval, possibly allowing more
				    * than one CPU core to take care of part
				    * of the work.
				    *
				    * The @p grainsize argument (50 in the
				    * example above) makes sure that subranges
				    * do not become too small, to avoid
				    * spending more time on scheduling
				    * subranges to CPU resources than on doing
				    * actual work.
				    *
				    * For a discussion of the kind of
				    * problems to which this function
				    * is applicable, see also the
				    * @ref threads "Parallel computing with multiple processors"
				    * module.
				    */
  template <typename RangeType, typename Function>
  void apply_to_subranges (const RangeType                          &begin,
			   const typename identity<RangeType>::type &end,
			   const Function                           &f,
			   const unsigned int                        grainsize)
  {
#if DEAL_II_USE_MT == 0
				     // make sure we don't get compiler
				     // warnings about unused arguments
    (void) grainsize;

    f (begin, end);
#else
    tbb::parallel_for (tbb::blocked_range<RangeType>
		       (begin, end, grainsize),
		       std_cxx1x::bind (&internal::apply_to_subranges<RangeType,Function>,
					_1,
					std_cxx1x::cref(f)),
		       tbb::auto_partitioner());
#endif
  }


  
  namespace internal
  {  
#if DEAL_II_USE_MT == 1
				     /**
				      * A class that conforms to the Body
				      * requirements of the TBB
				      * parallel_reduce function. The first
				      * template argument denotes the type on
				      * which the reduction is to be done. The
				      * second denotes the type of the
				      * function object that shall be called
				      * for each subrange.
				      */
    template <typename ResultType,
	      typename Function>
    struct ReductionOnSubranges
    {
					 /**
					  * A variable that will hold the
					  * result of the reduction.
					  */
	ResultType result;

					 /**
					  * Constructor. Take the function
					  * object to call on each sub-range
					  * as well as the neutral element
					  * with respect to the reduction
					  * operation.
					  *
					  * The second argument denotes a
					  * function object that will be use
					  * to reduce the result of two
					  * computations into one number. An
					  * example if we want to simply
					  * accumulate integer results would
					  * be std::plus<int>().
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
					  * Splitting constructor. See the TBB
					  * book for more details about this.
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
					  * Join operation: merge the results
					  * from computations on different
					  * sub-intervals.
					  */
	void join (const ReductionOnSubranges &r)
	  {
	    result = reductor(result, r.result);
	  }
	
					 /**
					  * Execute the given function on the
					  * specified range.
					  */
	template <typename RangeType>
	void operator () (const tbb::blocked_range<RangeType> &range)
	  {
	    result = reductor(result,
			      f (range.begin(), range.end()));
	  }

      private:
					 /**
					  * The function object to call on
					  * every sub-range.
					  */
	const Function f;
	
					 /**
					  * The neutral element with respect
					  * to the reduction operation. This
					  * is needed when calling the
					  * splitting constructor since we
					  * have to re-set the result variable
					  * in this case.
					  */
	const ResultType neutral_element;

					 /**
					  * The function object to be used to
					  * reduce the result of two calls
					  * into one number.
					  */
	const std_cxx1x::function<ResultType (ResultType, ResultType)> reductor;
    };
#endif  
  }
  

				   /**
				    * This function works a lot like the
				    * apply_to_subranges(), but it allows to
				    * accumulate numerical results computed on
				    * each subrange into one number. The type
				    * of this number is given by the
				    * ResultType template argument that needs
				    * to be explicitly specified.
				    *
				    * An example of use of this function is to
				    * compute the value of the expression $x^T
				    * A x$ for a square matrix $A$ and a
				    * vector $x$. The sum over rows can be
				    * parallelized and the whole code might
				    * look like this:
				    * @code
				    *   void matrix_norm (const FullMatrix &A,
				    *                     const Vector     &x)
				    *   {
				    *     return
				    *      std::sqrt
				    *       (parallel::accumulate_from_subranges<double>
				    *        (0, A.n_rows(),
				    *         std_cxx1x::bind (&mat_norm_sqr_on_subranges,
				    *                          _1, _2,
				    *                          std_cxx1x::cref(A),
				    *                          std_cxx1x::cref(x)),
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
				    * Here,
				    * <code>mat_norm_sqr_on_subranges</code>
				    * is called on the entire range
				    * <code>[0,A.n_rows())</code> if this
				    * range is less than the minimum grainsize
				    * (above chosen as 50) or if deal.II is
				    * configured to not use
				    * multithreading. Otherwise, it may be
				    * called on subsets of the given range,
				    * with results from the individual
				    * subranges accumulated internally.
				    *
				    * @warning If ResultType is a floating point
				    * type, then accumulation is not a
				    * commutative operation. In other words,
				    * if the given function object is called
				    * three times on three subranges,
				    * returning values $a,b,c$, then the
				    * returned result of this function is
				    * $a+b+c$. However, depending on how the
				    * three sub-tasks are distributed on
				    * available CPU resources, the result may
				    * also be $a+c+b$ or any other
				    * permutation; because floating point
				    * addition does not commute (as oppose, of
				    * course, to addition of real %numbers),
				    * the result of invoking this function
				    * several times may differ on the order of
				    * round-off.
				    *
				    * For a discussion of the kind of
				    * problems to which this function
				    * is applicable, see also the
				    * @ref threads "Parallel computing with multiple processors"
				    * module.
				    */
  template <typename ResultType, typename RangeType, typename Function>
  ResultType accumulate_from_subranges (const Function &f,
					const RangeType                          &begin,
					const typename identity<RangeType>::type &end,
					const unsigned int grainsize)
  {
#if DEAL_II_USE_MT == 0
				     // make sure we don't get compiler
				     // warnings about unused arguments
    (void) grainsize;

    return f(begin,end);
#else
    internal::ReductionOnSubranges<ResultType,Function>
      reductor (f, std::plus<ResultType>(), 0);
    tbb::parallel_reduce (tbb::blocked_range<RangeType>(begin, end, grainsize),
			  reductor,
			  tbb::auto_partitioner());
    return reductor.result;
#endif
  }
  
}


DEAL_II_NAMESPACE_CLOSE

#endif
