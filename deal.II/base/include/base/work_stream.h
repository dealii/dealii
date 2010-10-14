//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__work_stream_h
#define __deal2__work_stream_h


#include <base/config.h>
#include <base/multithread_info.h>
#include <base/template_constraints.h>
#include <base/std_cxx1x/function.h>
#include <base/std_cxx1x/bind.h>

#if DEAL_II_USE_MT == 1
#  include <base/thread_management.h>
#  include <tbb/pipeline.h>
#endif

#include <vector>
#include <utility>


DEAL_II_NAMESPACE_OPEN



/**
 * A class whose main template function supports running multiple
 * threads each of which operates on a subset of the given range of
 * objects. The class uses the Intel Threading Building Blocks (TBB)
 * to load balance the individual subranges onto the available
 * threads. For a lengthy discussion of the rationale of this class,
 * see the @ref threads "Parallel computing with multiple processors"
 * module.
 *
 * The class is built on the following premise: One frequently has some work
 * that needs to be done on a sequence of objects; a prototypical example is
 * assembling cell contributions to a system matrix or right hand side. In
 * many such examples, part of the work can be done entirely independently and
 * in parallel, possibly using several processor cores on a machine with
 * shared memory. However, some other part of this work may need to be
 * synchronised and be done in order. In the example of assembling a matrix,
 * the computation of local contributions can be done entirely in parallel,
 * but copying the the local contributions into the global matrix requires
 * some care: First, several threads can't write at the same time, but need to
 * synchronise writing using a mutex; secondly, we want the order in which
 * local contributions are added to the global matrix to be always the same
 * because floating point addition is not commutative and adding local
 * contributions to the global matrix in different orders leads to subtly
 * different results that can affect the number of iterations for iterative
 * solvers as well as the round-off error in the solution in random
 * ways. Consequently, we want to ensure that only one thread at a time writes
 * into the global matrix, and that results are copied in a stable and
 * reproducible order.
 *
 * This class implements a framework for this work model. It works with a
 * stream of objects given by an iterator range, runs a worker function in
 * parallel on all of these objects and then passes each object to a
 * postprocessor function that runs sequentially and gets objects in exactly
 * the order in which they appear in the input iterator range. None of the
 * synchronisation work is exposed to the user of this class.
 *
 * Internally, the range given to the run() function of this class is split
 * into a sequence of "items", which are then distributed according to some
 * %internal algorithm onto the number of available threads. An item is an
 * element of the range of iterators on which we are to operate; for example,
 * for the purpose of assembling matrices or evaluating error indicators, an
 * item could be a cell. The TBB library determines how many threads are
 * created (typically as many as there are processor cores), but the number of
 * items that may be active at any given time is specified by the argument to
 * the constructor. It should be bigger or equal to the number of processor
 * cores - the default is four times the number of cores on the current system.
 *
 * Items are created upon request by the TBB whenever one of the worker
 * threads is idle or is expected to become idle. It is then handed off to a
 * worker function, typically a member function of a main class. These worker
 * functions are run in parallel on a number of threads, and there is no
 * guarantee that they are asked to work on items in any particular order, in
 * particular not necessarily in the order in which items are generated from
 * the iterator range.
 *
 * Typically, worker functions need additional data, for example FEValues
 * objects, input data vectors, etc, some of which can not be shared among
 * threads. To this end, the run() function takes another template argument,
 * ScratchData, which designates a type objects of which are stored with
 * each item and which threads can use as private data without having to
 * share them with other threads. The run() function takes an additional
 * argument with an object of type ScratchData that is going to be copied
 * for the arguments passed to each of the worker functions.
 *
 * In addition, worker functions store their results in objects of template type
 * CopyData. These are then handed off to a separate function, called copier,
 * that may use the stored results to transfer them into permanent
 * storage. For example, it may copy the results of local contributions to a
 * matrix computed by a worker function into the global matrix. In contrast to
 * the worker function, however, only one instance of the copier is run at any
 * given time; it can therefore safely copy local contributions into the
 * global matrix without the need to lock the global object using a mutex or
 * similar means. Furthermore, it is guaranteed that the copier is run with
 * CopyData objects in the same order in which their associated items
 * were created; consequently, even if worker threads may compute results in
 * unspecified order, the copier always receives the results in exactly the
 * same order as the items were created.
 *
 * Once an item is processed by the copier, it is deleted and the
 * ScratchData and CopyData objects that were used in its computation
 * are considered unused and may be re-used for the next invokation of
 * the worker function, on this or another thread.
 *
 * This class only really works in parallel when multithread mode was selected
 * during deal.II configuration. Otherwise it simply works on each item
 * sequentially.
 *
 * @ingroup threads
 * @author Wolfgang Bangerth, 2007, 2008, 2009
 */
namespace WorkStream
{

#if DEAL_II_USE_MT == 1


  namespace internal
  {
				     /**
				      * A class that creates a sequence of
				      * items from a range of iterators.
				      */
    template <typename Iterator,
	      typename ScratchData,
	      typename CopyData>
    class IteratorRangeToItemStream : public tbb::filter
    {
      public:
					 /**
					  * A data type that we use to identify
					  * items to be worked on.
					  *
					  * The first element indicates an array
					  * of iterators to work on; the second
					  * the scratch space; the third an
					  * array of copy data spaces; and the
					  * last the number of elements to work
					  * on. This last argument is an integer
					  * between one and chunk_size. The
					  * arrays have a length of chunk_size.
					  */
	typedef
	std_cxx1x::tuple<std::vector<Iterator>,
			 ScratchData*,
			 std::vector<CopyData>,
			 unsigned int>
	ItemType;
	

					 /**
					  * Constructor. Take an iterator
					  * range, the size of a buffer that
					  * can hold items, and the sample
					  * additional data object that will
					  * be passed to each worker and
					  * copier function invokation.
					  */
	IteratorRangeToItemStream (const Iterator       &begin,
				   const Iterator       &end,
				   const unsigned int    buffer_size,
				   const unsigned int    chunk_size,
				   const ScratchData    &sample_scratch_data,
				   const CopyData       &sample_copy_data)
			:
			tbb::filter (/*is_serial=*/true),
			remaining_iterator_range (begin, end),
			ring_buffer (buffer_size),
			n_emitted_items (0),
			chunk_size (chunk_size)
	  {
					     // initialize copies of
					     // additional_data. since
					     // this is frequently
					     // expensive (creating
					     // FEValues objects etc) do
					     // that in parallel
	    Threads::TaskGroup<> tasks;
	    for (unsigned int i=0; i<ring_buffer.size(); ++i)
	      tasks += Threads::new_task (&IteratorRangeToItemStream::init_buffer_elements,
					  *this,
					  i,
					  std_cxx1x::cref(sample_scratch_data),
					  std_cxx1x::cref(sample_copy_data));
	    tasks.join_all ();
	  }

					 /**
					  * Destructor.
					  */
	~IteratorRangeToItemStream ()
	  {
	    for (unsigned int i=0; i<ring_buffer.size(); ++i)
	      delete std_cxx1x::get<1>(ring_buffer[i]);
	  }
      
					 /**
					  * Create a item and return a
					  * pointer to it.
					  */
	virtual void * operator () (void *)
	  {
					     // store the current
					     // position of the pointer
	    ItemType *current_item
	      = &ring_buffer[n_emitted_items % ring_buffer.size()];

					     // initialize the next item. it may
					     // consist of at most chunk_size
					     // elements
	    std_cxx1x::get<3>(*current_item) = 0;
	    while ((remaining_iterator_range.first !=
		    remaining_iterator_range.second)
		   &&
		   (std_cxx1x::get<3>(*current_item) < chunk_size))
	      {
		std_cxx1x::get<0>(*current_item)[std_cxx1x::get<3>(*current_item)]
		  = remaining_iterator_range.first;

		++remaining_iterator_range.first;
		++std_cxx1x::get<3>(*current_item);
	      }

	    if (std_cxx1x::get<3>(*current_item) == 0)
					       // there were no items
					       // left. terminate the pipeline
	      return 0;
	    else
	      {
		++n_emitted_items;
		return current_item;
	      }
	  }
	
      private:
					 /**
					  * The interval of iterators still to
					  * be worked on. This range will shrink
					  * over time.
					  */
	std::pair<Iterator,Iterator> remaining_iterator_range;
      
					 /**
					  * A ring buffer that will store items.
					  */
	std::vector<ItemType>        ring_buffer;

					 /**
					  * Counter for the number of emitted
					  * items. Each item may consist of up
					  * to chunk_size iterator elements.
					  */
	unsigned int                 n_emitted_items;

					 /**
					  * Number of elements of the
					  * iterator range that each
					  * thread should work on
					  * sequentially; a large number
					  * makes sure that each thread
					  * gets a significant amount of
					  * work before the next task
					  * switch happens, whereas a
					  * small number is better for
					  * load balancing.
					  */
	const unsigned int           chunk_size;

					 /**
					  * Initialize the pointers and vector
					  * elements in the specified entry of
					  * the ring buffer.
					  */
	void init_buffer_elements (const unsigned int element,
				   const ScratchData &sample_scratch_data,
				   const CopyData    &sample_copy_data)
	  {
	    Assert (std_cxx1x::get<1>(ring_buffer[element]) == 0,
		    ExcInternalError());
	      
	    std_cxx1x::get<0>(ring_buffer[element])
	      .resize (chunk_size, remaining_iterator_range.second);
	    std_cxx1x::get<1>(ring_buffer[element])
	      = new ScratchData(sample_scratch_data);
	    std_cxx1x::get<2>(ring_buffer[element])
	      .resize (chunk_size, sample_copy_data);
	  }
    };



				     /**
				      * A class that manages calling the
				      * worker function on a number of
				      * parallel threads. Note that it is, in
				      * the TBB notation, a filter that can
				      * run in parallel.
				      */
    template <typename Iterator,
	      typename ScratchData,
	      typename CopyData>
    class Worker : public tbb::filter
    {
      public:
					 /**
					  * Constructor. Takes a
					  * reference to the object on
					  * which we will operate as
					  * well as a pointer to the
					  * function that will do the
					  * assembly.
					  */
	Worker (const std_cxx1x::function<void (const Iterator &,
						ScratchData &,
						CopyData &)> &worker)
	:
	tbb::filter (/* is_serial= */ false),
	worker (worker)
	{}


					 /**
					  * Work on an item.
					  */
	void * operator () (void *item)
	  {
					     // first unpack the current item
	    typedef
	      typename IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
	      ItemType;
	    
	    ItemType *current_item = reinterpret_cast<ItemType*> (item);
	    
					     // then call the worker function on
					     // each element of the chunk we
					     // were given
	    for (unsigned int i=0; i<std_cxx1x::get<3>(*current_item); ++i)
	      worker (std_cxx1x::get<0>(*current_item)[i],
		      *std_cxx1x::get<1>(*current_item),
		      std_cxx1x::get<2>(*current_item)[i]);
	    
					     // then return the original pointer
					     // to the now modified object
	    return item;
	  }
	

      private:
					 /**
					  * Pointer to the function
					  * that does the assembling
					  * on the sequence of cells.
					  */
	const std_cxx1x::function<void (const Iterator &,
					ScratchData &,
					CopyData &)> worker;
    };



				     /**
				      * A class that manages calling the
				      * copier function. Note that it is, in
				      * the TBB notation, a filter that runs
				      * sequentially, ensuring that all items
				      * are copied in the same order in which
				      * they are created.
				      */
    template <typename Iterator,
	      typename ScratchData,
	      typename CopyData>
    class Copier : public tbb::filter
    {
      public:
					 /**
					  * Constructor. Takes a
					  * reference to the object on
					  * which we will operate as
					  * well as a pointer to the
					  * function that will do the
					  * copying from the
					  * additional data object to
					  * the global matrix or
					  * similar.
					  */
	Copier (const std_cxx1x::function<void (const CopyData &)> &copier)
			:
			tbb::filter (/* is_serial= */ true),
			copier (copier)
	  {}


					 /**
					  * Work on a single item.
					  */
	void * operator () (void *item)
	  {
					     // first unpack the current item
	    typedef
	      typename IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
	      ItemType;
	    
	    ItemType *current_item = reinterpret_cast<ItemType*> (item);

					     // initiate copying data
	    for (unsigned int i=0; i<std_cxx1x::get<3>(*current_item); ++i)
	      copier (std_cxx1x::get<2>(*current_item)[i]);

					     // return an invalid
					     // item since we are at
					     // the end of the
					     // pipeline
	    return 0;
	  }
	

      private:
					 /**
					  * Pointer to the function
					  * that does the copying of
					  * data.
					  */
	const std_cxx1x::function<void (const CopyData &)> copier;
    };
    
  }
  

#endif // DEAL_II_USE_MT



				   /**
				    * This is the main function of the
				    * WorkStream concept, doing work as
				    * described in the introduction to this
				    * namespace.
				    *
				    * This is the function that can be used
				    * for worker and copier objects that are
				    * either pointers to non-member
				    * functions or objects that allow to be
				    * called with an operator(), for example
				    * objects created by std::bind.
				    *
				    * The argument passed as @p end must be
				    * convertible to the same type as
				    * @p begin, but doesn't have to be of the
				    * same type itself. This allows to write
				    * code like
				    * <code>WorkStream().run(dof_handler.begin_active(),
				    * dof_handler.end(), ...</code> where
				    * the first is of type
				    * DoFHandler::active_cell_iterator
				    * whereas the second is of type
				    * DoFHandler::raw_cell_iterator.
				    *
				    * The two data types
				    * <tt>ScratchData</tt> and
				    * <tt>CopyData</tt> need to have a
				    * working copy
				    * constructor. <tt>ScratchData</tt>
				    * is only used in the
				    * <tt>worker</tt> function, while
				    * <tt>CopyData</tt> is the object
				    * passed from the <tt>worker</tt>
				    * to the <tt>copier</tt>.
				    *
				    * The @p queue_length argument indicates
				    * the number of items that can be live
				    * at any given time. Each item consists
				    * of @p chunk_size elements of the input
				    * stream that will be worked on by the
				    * worker and copier functions one after
				    * the other on the same thread.
				    */
  template <typename Worker,
	    typename Copier,
	    typename Iterator,
	    typename ScratchData,
	    typename CopyData>
  void
  run (const Iterator                          &begin,
       const typename identity<Iterator>::type &end,
       Worker                                   worker,
       Copier                                   copier,
       const ScratchData                       &sample_scratch_data,
       const CopyData                          &sample_copy_data,
       const unsigned int queue_length = 2*multithread_info.n_default_threads,
       const unsigned int chunk_size = 8)
  {
    Assert (queue_length > 0,
	    ExcMessage ("The queue length must be at least one, and preferably "
			"larger than the number of processors on this system."));
    Assert (chunk_size > 0,
	    ExcMessage ("The chunk_size must be at least one."));  

				     // if no work then skip. (only use
				     // operator!= for iterators since we may
				     // not have an equality comparison
				     // operator)
    if (!(begin != end))
      return;
    
#if DEAL_II_USE_MT == 1  
				     // create the three stages of the
				     // pipeline
    internal::IteratorRangeToItemStream<Iterator,ScratchData,CopyData>
      iterator_range_to_item_stream (begin, end,
				     queue_length,
				     chunk_size,
				     sample_scratch_data,
				     sample_copy_data);

    internal::Worker<Iterator, ScratchData, CopyData> worker_filter (worker);
    internal::Copier<Iterator, ScratchData, CopyData> copier_filter (copier);

				     // now create a pipeline from
				     // these stages
    tbb::pipeline assembly_line;
    assembly_line.add_filter (iterator_range_to_item_stream);
    assembly_line.add_filter (worker_filter);
    assembly_line.add_filter (copier_filter);

				     // and run it
    assembly_line.run (queue_length);    

    assembly_line.clear ();

#else

				     // need to copy the sample since it is
				     // marked const
    ScratchData scratch_data = sample_scratch_data;
    CopyData    copy_data    = sample_copy_data;
    
    for (Iterator i=begin; i!=end; ++i)
      {
	worker (i, scratch_data, copy_data);
	copier (copy_data);
      }  
#endif
  }



				   /**
				    * This is the main function of the
				    * WorkStream concept, doing work as
				    * described in the introduction to this
				    * namespace.
				    *
				    * This is the function that can be
				    * used for worker and copier functions
				    * that are member functions of a class.
				    *
				    * The argument passed as @p end must be
				    * convertible to the same type as
				    * @p begin, but doesn't have to be of the
				    * same type itself. This allows to write
				    * code like
				    * <code>WorkStream().run(dof_handler.begin_active(),
				    * dof_handler.end(), ...</code> where
				    * the first is of type
				    * DoFHandler::active_cell_iterator
				    * whereas the second is of type
				    * DoFHandler::raw_cell_iterator.
				    *
				    * The @p queue_length argument indicates
				    * the number of items that can be live
				    * at any given time. Each item consists
				    * of @p chunk_size elements of the input
				    * stream that will be worked on by the
				    * worker and copier functions one after
				    * the other on the same thread.
				    */
  template <typename MainClass,
	    typename Iterator,
	    typename ScratchData,
	    typename CopyData>
  void
  run (const Iterator                          &begin,
       const typename identity<Iterator>::type &end,
       MainClass                               &main_object,
       void (MainClass::*worker) (const Iterator &,
				  ScratchData &,
				  CopyData &),
       void (MainClass::*copier) (const CopyData &),
       const ScratchData                    &sample_scratch_data,
       const CopyData                       &sample_copy_data,
       const unsigned int queue_length = 2*multithread_info.n_default_threads,
       const unsigned int chunk_size = 8)
  {
				     // forward to the other function
    run (begin, end,
	 std_cxx1x::bind (worker,
			  std_cxx1x::ref (main_object),
			  _1, _2, _3),
	 std_cxx1x::bind (copier,
			  std_cxx1x::ref (main_object),
			  _1),
	 sample_scratch_data,
	 sample_copy_data,
	 queue_length,
	 chunk_size);
  }

}




DEAL_II_NAMESPACE_CLOSE




//----------------------------   work_stream.h     ---------------------------
// end of #ifndef __deal2__work_stream_h
#endif
//----------------------------   work_stream.h     ---------------------------
