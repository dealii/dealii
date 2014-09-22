// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__work_stream_h
#define __deal2__work_stream_h


#include <deal.II/base/config.h>
#include <deal.II/base/graph_coloring.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/bind.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/parallel.h>

#ifdef DEAL_II_WITH_THREADS
#  include <deal.II/base/thread_management.h>
#  include <tbb/pipeline.h>
#endif

#include <vector>
#include <utility>
#include <memory>


DEAL_II_NAMESPACE_OPEN



/**
 * A namespace whose main template function supports running multiple
 * threads each of which operates on a subset of the given range of
 * objects. The class uses the Intel Threading Building Blocks (TBB)
 * to load balance the individual subranges onto the available
 * threads. For a lengthy discussion of the rationale of this class,
 * see the @ref threads "Parallel computing with multiple processors"
 * module. It is used in the tutorial first in step-9, and again in
 * step-13, step-14, step-32 and others.
 *
 * The class is built on the following premise: One frequently has some work
 * that needs to be done on a sequence of objects; a prototypical example is
 * assembling cell contributions to a system matrix or right hand side. In
 * many such examples, part of the work can be done entirely independently and
 * in parallel, possibly using several processor cores on a machine with
 * shared memory. However, some other part of this work may need to be
 * synchronised and be done in order. In the example of assembling a matrix,
 * the computation of local contributions can be done entirely in parallel,
 * but copying the local contributions into the global matrix requires
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
 * The functions in this namespace only really work in parallel when
 * multithread mode was selected during deal.II configuration. Otherwise they
 * simply work on each item sequentially.
 *
 * @ingroup threads
 * @author Wolfgang Bangerth, 2007, 2008, 2009, 2013. Bruno Turcksin, 2013.
 */
namespace WorkStream
{

#ifdef DEAL_II_WITH_THREADS

  namespace internal
  {

//TODO: The following classes all use std_cxx11::shared_ptr, but the
//  correct pointer class would actually be std::unique_ptr. make this
//  replacement whenever we have a class that provides these semantics
//  and that is available also as a fall-back whenever via boost or similar

    /**
     * A namespace for the implementation of details of the WorkStream pattern
     * and function. This namespace holds classes that deal with the second
     * implementation described in the paper by Turcksin, Kronbichler and
     * Bangerth (see @ref workstream_paper).
     *
     * Even though this implementation is slower than the third implementation
     * discussed in that paper, we need to keep it around for two reasons:
     * (i) a user may not give us a graph coloring, (ii) we want to use
     * this implementation for colors that are just too small.
     */
    namespace Implementation2
    {
      /**
       * A class that creates a sequence of items from a range of iterators.
       */
      template <typename Iterator,
                typename ScratchData,
                typename CopyData>
      class IteratorRangeToItemStream : public tbb::filter
      {
      public:
        /**
         * A data type that we use to identify items to be worked on. This is
         * the structure that is passed around between the different parts of
         * the WorkStream implementation to identify what needs to be done by
         * the various stages of the pipeline.
         */
        struct ItemType
        {
          /**
           * A structure that contains a pointer to a scratch data object along
           * with a flag that indicates whether this object is currently in use.
           */
          struct ScratchDataObject
          {
            std_cxx11::shared_ptr<ScratchData> scratch_data;
            bool                               currently_in_use;

            /**
             * Default constructor.
             */
            ScratchDataObject ()
              :
              currently_in_use (false)
            {}

            ScratchDataObject (ScratchData *p,
                               const bool in_use)
              :
              scratch_data (p),
              currently_in_use (in_use)
            {}

//TODO: when we push back an object to the list of scratch objects, in
//  Worker::operator(), we first create an object and then copy
//  it to the end of this list. this involves having two objects
//      of the current type having pointers to it, each with their own
//      currently_in_use flag. there is probably little harm in this because
//      the original one goes out of scope right away again, but it's
//      certainly awkward. one way to avoid this would be to use unique_ptr
//      but we'd need to figure out a way to use it in non-C++11 mode
            ScratchDataObject (const ScratchDataObject &o)
              :
              scratch_data (o.scratch_data),
              currently_in_use (o.currently_in_use)
            {}
          };


          /**
           * Typedef to a list of scratch data objects. The rationale for this
           * list is provided in the variables that use these lists.
           */
          typedef std::list<ScratchDataObject> ScratchDataList;

          /**
           * A list of iterators that need to be worked on. Only the first
           * n_items are relevant.
           */
          std::vector<Iterator> work_items;

          /**
           * The CopyData objects that the Worker part of the pipeline
           * fills for each work item. Again, only the first n_items
           * elements are what we care about.
           */
          std::vector<CopyData> copy_datas;

          /**
           * Number of items identified by the work_items array that the
           * Worker and Copier pipeline stage need to work on. The maximum
           * value of this variable will be chunk_size.
           */
          unsigned int          n_items;

          /**
           * Pointer to a thread local variable identifying the scratch data objects
           * this thread will use. The initial implementation of this
           * class using thread local variables provided only a single
           * scratch object per thread. This doesn't work, because
           * the worker functions may start tasks itself and then call
           * Threads::TaskGroup::join_all() or a similar function, which the
           * TBB scheduler may use to run something else on the current
           * thread -- for example another instance of the worker function.
           * Consequently, there would be two instances of the worker
           * function that use the same scratch object if we only
           * provided a single scratch object per thread. The solution is
           * to provide a list of scratch objects for each thread, together
           * with a flag indicating whether this scratch object is currently
           * used. If a thread needs a scratch object, it walks this list
           * until it finds an unused object, or, if there is none, creates one
           * itself. Note that we need not use synchronization primitives
           * for this process since the lists are thread-local and
           * we are guaranteed that only a single thread accesses them as long
           * as we have no yield point in between the accesses to the list.
           *
           * The pointers to scratch objects stored in each of these lists must
           * be so that they are deleted on all threads when the thread
           * local object is destroyed. This is achieved by using shared_ptr.
           *
           * Note that when a worker needs to create a scratch object, it allocates
           * it using sample_scratch_data to copy from. This has
           * the advantage of a first-touch initialization, i.e., the
           * memory for the scratch data object is allocated and initialized
           * by the same thread that will later use it.
           */
          Threads::ThreadLocalStorage<ScratchDataList> *scratch_data;

          /**
           * Pointer to a sample scratch data object, to be used to initialize
           * the scratch data objects created for each individual thread.
           */
          const ScratchData *sample_scratch_data;

          /**
           * Flag is true if the buffer is used and false if the buffer can be
           * used.
           */
          bool currently_in_use;


          /**
           * Default constructor.
           * Initialize everything that doesn't have a default constructor
           * itself.
           */
          ItemType ()
            :
            n_items (0),
            scratch_data (0),
            sample_scratch_data (0),
            currently_in_use (false)
          {}
        };


        /**
         * Constructor. Take an iterator range, the size of a buffer that can
         * hold items, and the sample additional data object that will be passed
         * to each worker and copier function invokation.
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
          item_buffer (buffer_size),
          sample_scratch_data (sample_scratch_data),
          chunk_size (chunk_size)
        {
          // initialize the elements of the ring buffer
          for (unsigned int element=0; element<item_buffer.size(); ++element)
            {
              Assert (item_buffer[element].n_items == 0,
                      ExcInternalError());

              item_buffer[element].work_items.resize (chunk_size,
                                                      remaining_iterator_range.second);
              item_buffer[element].scratch_data = &thread_local_scratch;
              item_buffer[element].sample_scratch_data = &sample_scratch_data;
              item_buffer[element].copy_datas.resize (chunk_size,
                                                      sample_copy_data);
              item_buffer[element].currently_in_use = false;
            }
        }


        /**
         * Create an item and return a
         * pointer to it.
         */
        virtual void *operator () (void *)
        {
          // find first unused item. we know that there must be one
          // because we have set the maximal number of tokens in flight
          // and have set the ring buffer to have exactly this size. so
          // if this function is called, we know that less than the
          // maximal number of items in currently in flight
          //
          // note that we need not lock access to this array since
          // the current stage is run sequentially and we can therefore
          // enter the following block only once at any given time.
          // thus, there can be no race condition between checking that
          // a flag is false and setting it to true. (there may be
          // another thread where we release items and set 'false'
          // flags to 'true', but that too does not produce any
          // problems)
          ItemType *current_item = 0;
          for (unsigned int i=0; i<item_buffer.size(); ++i)
            if (item_buffer[i].currently_in_use == false)
              {
                item_buffer[i].currently_in_use = true;
                current_item = &item_buffer[i];
                break;
              }
          Assert (current_item != 0, ExcMessage ("This can't be. There must be a free item!"));

          // initialize the next item. it may
          // consist of at most chunk_size
          // elements
          current_item->n_items = 0;
          while ((remaining_iterator_range.first !=
                  remaining_iterator_range.second)
                 &&
                 (current_item->n_items < chunk_size))
            {
              current_item->work_items[current_item->n_items]
                = remaining_iterator_range.first;

              ++remaining_iterator_range.first;
              ++current_item->n_items;
            }

          if (current_item->n_items == 0)
            // there were no items
            // left. terminate the pipeline
            return 0;
          else
            return current_item;
        }

      private:
        /**
         * The interval of iterators still to
         * be worked on. This range will shrink
         * over time.
         */
        std::pair<Iterator,Iterator> remaining_iterator_range;

        /**
         * A buffer that will store items.
         */
        std::vector<ItemType>        item_buffer;

        /**
         * Pointer to a thread local variable identifying the scratch data objects
         * this thread will use. The initial implementation of this
         * class using thread local variables provided only a single
         * scratch object per thread. This doesn't work, because
         * the worker functions may start tasks itself and then call
         * Threads::TaskGroup::join_all() or a similar function, which the
         * TBB scheduler may use to run something else on the current
         * thread -- for example another instance of the worker function.
         * Consequently, there would be two instances of the worker
         * function that use the same scratch object if we only
         * provided a single scratch object per thread. The solution is
         * to provide a list of scratch objects for each thread, together
         * with a flag indicating whether this scratch object is currently
         * used. If a thread needs a scratch object, it walks this list
         * until it finds an unused object, or, if there is none, creates one
         * itself. Note that we need not use synchronization primitives
         * for this process since the lists are thread-local and
         * we are guaranteed that only a single thread accesses them as long
         * as we have no yield point in between the accesses to the list.
         *
         * The pointers to scratch objects stored in each of these lists must
         * be so that they are deleted on all threads when the thread
         * local object is destroyed. This is achieved by using shared_ptr.
         *
         * Note that when a worker needs to create a scratch object, it allocates
         * it using sample_scratch_data to copy from. This has
         * the advantage of a first-touch initialization, i.e., the
         * memory for the scratch data object is allocated and initialized
         * by the same thread that will later use it.
         */
        Threads::ThreadLocalStorage<typename ItemType::ScratchDataList> thread_local_scratch;

        /**
         * A reference to a sample scratch data that will be used to
         * initialize the thread-local pointers to a scratch data object
         * each of the worker tasks uses.
         */
        const ScratchData &sample_scratch_data;

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
         * the item_buffer.
         */
        void init_buffer_elements (const unsigned int element,
                                   const CopyData    &sample_copy_data)
        {
          Assert (item_buffer[element].n_items == 0,
                  ExcInternalError());

          item_buffer[element].work_items
          .resize (chunk_size, remaining_iterator_range.second);
          item_buffer[element].scratch_data
            = &thread_local_scratch;
          item_buffer[element].sample_scratch_data
            = &sample_scratch_data;
          item_buffer[element].copy_datas
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
        Worker (const std_cxx11::function<void (const Iterator &,
                                                ScratchData &,
                                                CopyData &)> &worker,
                bool copier_exist=true)
          :
          tbb::filter (/* is_serial= */ false),
          worker (worker),
          copier_exist (copier_exist)
        {}


        /**
         * Work on an item.
         */
        void *operator () (void *item)
        {
          // first unpack the current item
          typedef
          typename IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
          ItemType;

          ItemType *current_item = static_cast<ItemType *> (item);

          // we need to find an unused scratch data object in the list that
          // corresponds to the current thread and then mark it as used. if
          // we can't find one, create one
          //
          // as discussed in the discussion of the documentation of the
          // IteratorRangeToItemStream::scratch_data variable, there is no
          // need to synchronize access to this variable using a mutex
          // as long as we have no yield-point in between. this means that
          // we can't take an iterator into the list now and expect it to
          // still be valid after calling the worker, but we at least do
          // not have to lock the following section
          ScratchData *scratch_data = 0;
          {
            typename ItemType::ScratchDataList &
            scratch_data_list = current_item->scratch_data->get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ItemType::ScratchDataList::iterator
                 p = scratch_data_list.begin();
                 p != scratch_data_list.end(); ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data = p->scratch_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no object was found, create one and mark it as used
            if (scratch_data == 0)
              {
                scratch_data = new ScratchData(*current_item->sample_scratch_data);

                typename ItemType::ScratchDataList::value_type
                new_scratch_object (scratch_data, true);
                scratch_data_list.push_back (new_scratch_object);
              }
          }

          // then call the worker function on each element of the chunk we were
          // given. since these worker functions are called on separate threads,
          // nothing good can happen if they throw an exception and we are best
          // off catching it and showing an error message
          for (unsigned int i=0; i<current_item->n_items; ++i)
            {
              try
                {
                  if (worker)
                    worker (current_item->work_items[i],
                            *scratch_data,
                            current_item->copy_datas[i]);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception (exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception ();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            typename ItemType::ScratchDataList &
            scratch_data_list = current_item->scratch_data->get();

            for (typename ItemType::ScratchDataList::iterator p =
                   scratch_data_list.begin(); p != scratch_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }

          // if there is no copier, mark current item as usable again
          if (copier_exist==false)
            current_item->currently_in_use = false;


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
        const std_cxx11::function<void (const Iterator &,
                                        ScratchData &,
                                        CopyData &)> worker;

        /**
         * This flag is true if the copier stage exist. If it does not,
         * the worker has to free the buffer. Otherwise the copier will do it.
         */
        bool copier_exist;
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
        Copier (const std_cxx11::function<void (const CopyData &)> &copier)
          :
          tbb::filter (/*is_serial=*/true),
          copier (copier)
        {}


        /**
         * Work on a single item.
         */
        void *operator () (void *item)
        {
          // first unpack the current item
          typedef
          typename IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
          ItemType;

          ItemType *current_item = static_cast<ItemType *> (item);

          // initiate copying data. for the same reasons as in the worker class
          // above, catch exceptions rather than letting it propagate into
          // unknown territories
          for (unsigned int i=0; i<current_item->n_items; ++i)
            {
              try
                {
                  if (copier)
                    copier (current_item->copy_datas[i]);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception (exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception ();
                }
            }

          // mark current item as useable again
          current_item->currently_in_use = false;


          // return an invalid item since we are at the end of the
          // pipeline
          return 0;
        }


      private:
        /**
         * Pointer to the function that does the copying of data.
         */
        const std_cxx11::function<void (const CopyData &)> copier;
      };

    }


    /**
     * A namespace for the implementation of details of the WorkStream pattern
     * and function. This namespace holds classes that deal with the third
     * implementation described in the paper by Turcksin, Kronbichler and
     * Bangerth (see @ref workstream_paper).
     */
    namespace Implementation3
    {
      /**
       * A class that creates a sequence of items from a range of iterators.
       */
      template <typename Iterator,
                typename ScratchData,
                typename CopyData>
      class IteratorRangeToItemStream : public tbb::filter
      {
      public:
        /**
         * A data type that we use to identify items to be worked on. This is
         * the structure that is passed around between the different parts of
         * the WorkStream implementation to identify what needs to be done by
         * the various stages of the pipeline.
         */
        struct ItemType
        {
          /**
           * A structure that contains a pointer to scratch and copy data objects along
           * with a flag that indicates whether this object is currently in use.
           */
          struct ScratchAndCopyDataObjects
          {
            std_cxx11::shared_ptr<ScratchData> scratch_data;
            std_cxx11::shared_ptr<CopyData>    copy_data;
            bool                               currently_in_use;

            /**
             * Default constructor.
             */
            ScratchAndCopyDataObjects ()
              :
              currently_in_use (false)
            {}

            ScratchAndCopyDataObjects (ScratchData *p,
                                       CopyData *q,
                                       const bool in_use)
              :
              scratch_data (p),
              copy_data (q),
              currently_in_use (in_use)
            {}

//TODO: when we push back an object to the list of scratch objects, in
//      Worker::operator(), we first create an object and then copy
//      it to the end of this list. this involves having two objects
//      of the current type having pointers to it, each with their own
//      currently_in_use flag. there is probably little harm in this because
//      the original one goes out of scope right away again, but it's
//      certainly awkward. one way to avoid this would be to use unique_ptr
//      but we'd need to figure out a way to use it in non-C++11 mode
            ScratchAndCopyDataObjects (const ScratchAndCopyDataObjects &o)
              :
              scratch_data (o.scratch_data),
              copy_data (o.copy_data),
              currently_in_use (o.currently_in_use)
            {}
          };


          /**
           * Typedef to a list of scratch data objects. The rationale for this
           * list is provided in the variables that use these lists.
           */
          typedef std::list<ScratchAndCopyDataObjects> ScratchAndCopyDataList;

          /**
           * A list of iterators that need to be worked on. Only the first
           * n_items are relevant.
           */
          std::vector<Iterator> work_items;

          /**
           * Number of items identified by the work_items array that the
           * Worker and Copier pipeline stage need to work on. The maximum
           * value of this variable will be chunk_size.
           */
          unsigned int          n_items;

          /**
           * Pointer to a thread local variable identifying the scratch
           * and copy data objects this thread will use. The same considerations
           * apply as documented in the Implementation2::IteratorRangeToItemStream
           * class as well as in the paper by Turcksin, Kronbichler and Bangerth
          * (see @ref workstream_paper).
           */
          Threads::ThreadLocalStorage<ScratchAndCopyDataList> *scratch_and_copy_data;

          /**
           * Pointer to a sample scratch data object, to be used to initialize
           * the scratch data objects created for each individual thread.
           */
          const ScratchData *sample_scratch_data;

          /**
           * Pointer to a sample copy data object.
           */
          const CopyData *sample_copy_data;

          /**
           * Flag is true if the buffer is used and false if the buffer can be
           * used.
           */
          bool currently_in_use;

          /**
           * Default constructor.
           * Initialize everything that doesn't have a default constructor
           * itself.
           */
          ItemType ()
            :
            n_items (0),
            sample_scratch_data (0),
            sample_copy_data (0),
            currently_in_use (false)
          {}
        };


        /**
         * Constructor. Take range of iterators into an array of the kind of object we
         * want to loop over, the size of a buffer that can
         * hold items, and the sample additional data object that will be passed
         * to each worker and copier function invokation.
         */
        IteratorRangeToItemStream (const typename std::vector<Iterator>::const_iterator &begin,
                                   const typename std::vector<Iterator>::const_iterator &end,
                                   const unsigned int    buffer_size,
                                   const unsigned int    chunk_size,
                                   const ScratchData    &sample_scratch_data,
                                   const CopyData       &sample_copy_data)
          :
          tbb::filter (/*is_serial=*/true),
          remaining_iterator_range (begin, end),
          item_buffer (buffer_size),
          sample_scratch_data (sample_scratch_data),
          sample_copy_data (sample_copy_data),
          chunk_size (chunk_size)
        {
          Assert (begin != end, ExcMessage ("This class is not prepared to deal with empty ranges!"));
          // initialize the elements of the item_buffer
          for (unsigned int element=0; element<item_buffer.size(); ++element)
            {
              Assert (item_buffer[element].n_items == 0,
                      ExcInternalError());

              // resize the item_buffer. we have to initialize the new
              // elements with something and because we can't rely on there
              // being a default constructor for 'Iterator', we take the first
              // element in the range [begin,end) pointed to.
              item_buffer[element].work_items.resize (chunk_size,*begin);
              item_buffer[element].scratch_and_copy_data = &thread_local_scratch_and_copy;
              item_buffer[element].sample_scratch_data = &sample_scratch_data;
              item_buffer[element].sample_copy_data = &sample_copy_data;
              item_buffer[element].currently_in_use = false;
            }
        }


        /**
         * Create an item and return a
         * pointer to it.
         */
        virtual void *operator () (void *)
        {
          // find first unused item. we know that there must be one
          // because we have set the maximal number of tokens in flight
          // and have set the ring buffer to have exactly this size. so
          // if this function is called, we know that less than the
          // maximal number of items in currently in flight
          ItemType *current_item = 0;
          for (unsigned int i=0; i<item_buffer.size(); ++i)
            if (item_buffer[i].currently_in_use == false)
              {
                item_buffer[i].currently_in_use = true;
                current_item = &item_buffer[i];
                break;
              }
          Assert (current_item != 0, ExcMessage ("This can't be. There must be a free item!"));


          // initialize the next item. it may
          // consist of at most chunk_size
          // elements
          current_item->n_items = 0;
          while ((remaining_iterator_range.first !=
                  remaining_iterator_range.second)
                 &&
                 (current_item->n_items < chunk_size))
            {
              // initialize the iterators to work on with the elements
              // of the vector that remaining_iterator_range
              // points into
              current_item->work_items[current_item->n_items]
                = *remaining_iterator_range.first;

              ++remaining_iterator_range.first;
              ++current_item->n_items;
            }

          if (current_item->n_items == 0)
            // there were no items
            // left. terminate the pipeline
            return 0;
          else
            return current_item;
        }

      private:
        /**
         * The interval of iterators still to
         * be worked on. This range will shrink
         * over time.
         */
        std::pair<typename std::vector<Iterator>::const_iterator,typename std::vector<Iterator>::const_iterator> remaining_iterator_range;

        /**
         * A buffer that will store items.
         */
        std::vector<ItemType>        item_buffer;

        /**
         * Pointer to a thread local variable identifying the scratch and
         * copy data objects each thread will use. The same is true as
         * discussed for the implementation in the
         * Implementation2::IteratorRangeToItemStream class and the paper
         * by Turcksin, Kronbichler and Bangerth (see @ref workstream_paper).
         */
        Threads::ThreadLocalStorage<typename ItemType::ScratchAndCopyDataList> thread_local_scratch_and_copy;

        /**
         * A reference to a sample scratch data that will be used to
         * initialize the thread-local pointers to a scratch data object
         * each of the worker tasks uses.
         */
        const ScratchData &sample_scratch_data;

        /**
         * A reference to a sample scratch data that will be used to
         * initialize the thread-local pointers to a scratch data object
         * each of the worker tasks uses.
         */
        const CopyData &sample_copy_data;

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
         * the item_buffer.
         */
        void init_buffer_elements (const unsigned int element)
        {
          Assert (item_buffer[element].n_items == 0,
                  ExcInternalError());

          item_buffer[element].work_items
          .resize (chunk_size, remaining_iterator_range.second);
          item_buffer[element].scratch_and_copy_data
            = &thread_local_scratch_and_copy;
          item_buffer[element].sample_scratch_data
            = &sample_scratch_data;
          item_buffer[element].sample_copy_data
            = &sample_copy_data;
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
      class WorkerAndCopier : public tbb::filter
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
        WorkerAndCopier (const std_cxx11::function<void (const Iterator &,
                                                         ScratchData &,
                                                         CopyData &)> &worker,
                         const std_cxx11::function<void (const CopyData &)> &copier)
          :
          tbb::filter (/* is_serial= */ false),
          worker (worker),
          copier (copier)
        {}


        /**
         * Work on an item.
         */
        void *operator () (void *item)
        {
          // first unpack the current item
          typedef
          typename IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
          ItemType;

          ItemType *current_item = static_cast<ItemType *> (item);

          // we need to find an unused scratch and corresponding copy
          // data object in the list that
          // corresponds to the current thread and then mark it as used. if
          // we can't find one, create one
          //
          // as discussed in the discussion of the documentation of the
          // IteratorRangeToItemStream::scratch_data variable, there is no
          // need to synchronize access to this variable using a mutex
          // as long as we have no yield-point in between. this means that
          // we can't take an iterator into the list now and expect it to
          // still be valid after calling the worker, but we at least do
          // not have to lock the following section
          ScratchData *scratch_data = 0;
          CopyData    *copy_data    = 0;
          {
            typename ItemType::ScratchAndCopyDataList &
            scratch_and_copy_data_list = current_item->scratch_and_copy_data->get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ItemType::ScratchAndCopyDataList::iterator
                 p = scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end(); ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data = p->scratch_data.get();
                  copy_data    = p->copy_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no element in the list was found, create one and mark it as used
            if (scratch_data == 0)
              {
                Assert (copy_data==0, ExcInternalError());
                scratch_data = new ScratchData(*current_item->sample_scratch_data);
                copy_data    = new CopyData(*current_item->sample_copy_data);

                typename ItemType::ScratchAndCopyDataList::value_type
                new_scratch_object (scratch_data, copy_data, true);
                scratch_and_copy_data_list.push_back (new_scratch_object);
              }
          }

          // then call the worker and copier function on each element of the chunk we were
          // given. since these functions are called on separate threads,
          // nothing good can happen if they throw an exception and we are best
          // off catching it and showing an error message
          for (unsigned int i=0; i<current_item->n_items; ++i)
            {
              try
                {
                  if (worker)
                    worker (current_item->work_items[i],
                            *scratch_data,
                            *copy_data);
                  if (copier)
                    copier (*copy_data);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception (exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception ();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            typename ItemType::ScratchAndCopyDataList &
            scratch_and_copy_data_list = current_item->scratch_and_copy_data->get();

            for (typename ItemType::ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin(); p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }

          // mark current item as usable again
          current_item->currently_in_use = false;

          // return an invalid item since we are at the end of the
          // pipeline
          return 0;
        }


      private:
        /**
         * Pointer to the function
         * that does the assembling
         * on the sequence of cells.
         */
        const std_cxx11::function<void (const Iterator &,
                                        ScratchData &,
                                        CopyData &)> worker;

        /**
         * Pointer to the function that does the copying from
         * local contribution to global object.
         */
        const std_cxx11::function<void (const CopyData &)> copier;
      };
    }


    /**
     * A namespace for functions used in the implementation of
     * Implementation3 whenever we don't actually have a copier
     * function. In that case, we can simply use parallel_for, but we
     * have to be careful with the use of thread-local objects as we
     * had to be with Implementation3 as well.
     */
    namespace ParallelFor
    {
      template <typename Iterator,
                typename ScratchData,
                typename CopyData>
      class Worker
      {
      public:
        /**
         * Constructor.
         */
        Worker (const std_cxx11::function<void (const Iterator &,
                                                ScratchData &,
                                                CopyData &)> &worker,
                const ScratchData    &sample_scratch_data,
                const CopyData       &sample_copy_data)
          :
          worker (worker),
          sample_scratch_data (sample_scratch_data),
          sample_copy_data (sample_copy_data)
        {}


        /**
         * The function that calls the worker function on a
         * range of items denoted by the two arguments.
         */
        void operator() (const tbb::blocked_range<typename std::vector<Iterator>::const_iterator> &range)
        {
          // we need to find an unused scratch and corresponding copy
          // data object in the list that
          // corresponds to the current thread and then mark it as used. if
          // we can't find one, create one
          //
          // as discussed in the discussion of the documentation of the
          // IteratorRangeToItemStream::scratch_data variable, there is no
          // need to synchronize access to this variable using a mutex
          // as long as we have no yield-point in between. this means that
          // we can't take an iterator into the list now and expect it to
          // still be valid after calling the worker, but we at least do
          // not have to lock the following section
          ScratchData *scratch_data = 0;
          CopyData    *copy_data    = 0;
          {
            typename ItemType::ScratchAndCopyDataList &
            scratch_and_copy_data_list = data.get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ItemType::ScratchAndCopyDataList::iterator
                 p = scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end(); ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data = p->scratch_data.get();
                  copy_data    = p->copy_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no element in the list was found, create one and mark it as used
            if (scratch_data == 0)
              {
                Assert (copy_data==0, ExcInternalError());
                scratch_data = new ScratchData(sample_scratch_data);
                copy_data    = new CopyData(sample_copy_data);

                typename ItemType::ScratchAndCopyDataList::value_type
                new_scratch_object (scratch_data, copy_data, true);
                scratch_and_copy_data_list.push_back (new_scratch_object);
              }
          }

          // then call the worker and copier function on each
          // element of the chunk we were given. since these
          // functions are called on separate threads, nothing good
          // can happen if they throw an exception and we are best
          // off catching it and showing an error message
          for (typename std::vector<Iterator>::const_iterator p=range.begin();
               p != range.end(); ++p)
            {
              try
                {
                  worker (*p,
                          *scratch_data,
                          *copy_data);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception (exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception ();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            typename ItemType::ScratchAndCopyDataList &
            scratch_and_copy_data_list = data.get();

            for (typename ItemType::ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin(); p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }

        }

      private:
        typedef
        typename Implementation3::IteratorRangeToItemStream<Iterator,ScratchData,CopyData>::ItemType
        ItemType;

        typedef
        typename ItemType::ScratchAndCopyDataList
        ScratchAndCopyDataList;

        Threads::ThreadLocalStorage<ScratchAndCopyDataList> data;

        /**
         * Pointer to the function
         * that does the assembling
         * on the sequence of cells.
         */
        const std_cxx11::function<void (const Iterator &,
                                        ScratchData &,
                                        CopyData &)> worker;

        /**
         * References to sample scratch and copy data for
         * when we need them.
         */
        const ScratchData    &sample_scratch_data;
        const CopyData       &sample_copy_data;
      };
    }

  }


#endif // DEAL_II_WITH_THREADS


  /**
   * This is one of two main functions of the WorkStream concept, doing work as
   * described in the introduction to this namespace. It corresponds to
   * implementation 3 of the paper by Turcksin, Kronbichler and Bangerth,
   * see @ref workstream_paper .
   * As such, it takes not a range of iterators described by a begin
   * and end iterator, but a "colored" graph of iterators where each
   * color represents cells for which writing the cell contributions into
   * the global object does not conflict (in other words, these cells
   * are not neighbors). Each "color" is represented by std::vectors of cells.
   * The first argument to this function, a set of sets of cells (which are
   * represent as a vector of vectors, for efficiency), is typically
   * constructed by calling GraphColoring::make_graph_coloring(). See there
   * for more information.
   *
   * This function that can be used for worker and copier objects that
   * are either pointers to non-member functions or objects that allow to be
   * called with an operator(), for example objects created by std::bind.
   *
   * The two data types <tt>ScratchData</tt> and <tt>CopyData</tt> need to
   * have a working copy constructor. <tt>ScratchData</tt> is only used in the
   * <tt>worker</tt> function, while <tt>CopyData</tt> is the object passed
   * from the <tt>worker</tt> to the <tt>copier</tt>.
   *
   * The @p queue_length argument indicates the number of items that can be
   * live at any given time. Each item consists of @p chunk_size elements of
   * the input stream that will be worked on by the worker and copier
   * functions one after the other on the same thread.
   *
   * @note If your data objects are large, or their constructors are
   * expensive, it is helpful to keep in mind that <tt>queue_length</tt>
   * copies of the <tt>ScratchData</tt> object and
   * <tt>queue_length*chunk_size</tt> copies of the <tt>CopyData</tt> object
   * are generated.
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run (const std::vector<std::vector<Iterator> > &colored_iterators,
       Worker                                     worker,
       Copier                                     copier,
       const ScratchData                         &sample_scratch_data,
       const CopyData                            &sample_copy_data,
       const unsigned int queue_length = 2*multithread_info.n_threads(),
       const unsigned int                         chunk_size = 8);


  /**
   * This is one of two main functions of the WorkStream concept, doing work as
   * described in the introduction to this namespace. It corresponds to
   * implementation 2 of the paper by Turcksin, Kronbichler and Bangerth
   * (see @ref workstream_paper).
   *
   * This function that can be used for worker and copier objects that
   * are either pointers to non-member functions or objects that allow to be
   * called with an operator(), for example objects created by std::bind. If
   * the copier is an empty function, it is ignored in the pipeline.
   *
   * The argument passed as @p end must be convertible to the same type as @p
   * begin, but doesn't have to be of the same type itself. This allows to
   * write code like <code>WorkStream().run(dof_handler.begin_active(),
   * dof_handler.end(), ...</code> where the first is of type
   * DoFHandler::active_cell_iterator whereas the second is of type
   * DoFHandler::raw_cell_iterator.
   *
   * The two data types <tt>ScratchData</tt> and <tt>CopyData</tt> need to
   * have a working copy constructor. <tt>ScratchData</tt> is only used in the
   * <tt>worker</tt> function, while <tt>CopyData</tt> is the object passed
   * from the <tt>worker</tt> to the <tt>copier</tt>.
   *
   * The @p queue_length argument indicates the number of items that can be
   * live at any given time. Each item consists of @p chunk_size elements of
   * the input stream that will be worked on by the worker and copier
   * functions one after the other on the same thread.
   *
   * @note If your data objects are large, or their constructors are
   * expensive, it is helpful to keep in mind that <tt>queue_length</tt>
   * copies of the <tt>ScratchData</tt> object and
   * <tt>queue_length*chunk_size</tt> copies of the <tt>CopyData</tt> object
   * are generated.
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
       const unsigned int queue_length = 2*multithread_info.n_threads(),
       const unsigned int                       chunk_size = 8)
  {
    Assert (queue_length > 0,
            ExcMessage ("The queue length must be at least one, and preferably "
                        "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert (chunk_size > 0,
            ExcMessage ("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode

    // if no work then skip. (only use operator!= for iterators since we may
    // not have an equality comparison operator)
    if (!(begin != end))
      return;

    // we want to use TBB if we have support and if it is not disabled at
    // runtime:
#ifdef DEAL_II_WITH_THREADS
    if (multithread_info.n_threads()==1)
#endif
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data;

        for (Iterator i=begin; i!=end; ++i)
          {
            // need to check if the function is not the zero function. To
            // check zero-ness, create a C++ function out of it and check that
            if (static_cast<const std_cxx11::function<void (const Iterator &,
                                                            ScratchData &,
                                                            CopyData &)>& >(worker))
              worker (i, scratch_data, copy_data);
            if (static_cast<const std_cxx11::function<void (const CopyData &)>& >
                (copier))
              copier (copy_data);
          }
      }
#ifdef DEAL_II_WITH_THREADS
    else // have TBB and use more than one thread
      {
        // Check that the copier exist
        if (static_cast<const std_cxx11::function<void (const CopyData &)>& >(copier))
          {
            // create the three stages of the pipeline
            internal::Implementation2::IteratorRangeToItemStream<Iterator,ScratchData,CopyData>
            iterator_range_to_item_stream (begin, end,
                                           queue_length,
                                           chunk_size,
                                           sample_scratch_data,
                                           sample_copy_data);

            internal::Implementation2::Worker<Iterator, ScratchData, CopyData> worker_filter (worker);
            internal::Implementation2::Copier<Iterator, ScratchData, CopyData> copier_filter (copier);

            // now create a pipeline from these stages
            tbb::pipeline assembly_line;
            assembly_line.add_filter (iterator_range_to_item_stream);
            assembly_line.add_filter (worker_filter);
            assembly_line.add_filter (copier_filter);

            // and run it
            assembly_line.run (queue_length);

            assembly_line.clear ();
          }
        else
          {
            // there is no copier function. in this case, we have an
            // embarrassingly parallel problem where we can
            // essentially apply parallel_for. because parallel_for
            // requires subdividing the range for which operator- is
            // necessary between iterators, it is often inefficient to
            // apply it directory to cell ranges and similar iterator
            // types for which operator- is expensive or, in fact,
            // nonexistent. rather, in that case, we simply copy the
            // iterators into a large array and use operator- on
            // iterators to this array of iterators.
            //
            // instead of duplicating code, this is essentially the
            // same situation we have in Implementation3 below, so we
            // just defer to that place
            std::vector<std::vector<Iterator> > all_iterators (1);
            for (Iterator p=begin; p!=end; ++p)
              all_iterators[0].push_back (p);

            run (all_iterators,
                 worker, copier,
                 sample_scratch_data,
                 sample_copy_data,
                 queue_length,
                 chunk_size);
          }
      }
#endif
  }


  // Implementation 3:
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run (const std::vector<std::vector<Iterator> > &colored_iterators,
       Worker                                     worker,
       Copier                                     copier,
       const ScratchData                         &sample_scratch_data,
       const CopyData                            &sample_copy_data,
       const unsigned int                         queue_length,
       const unsigned int                         chunk_size)
  {
    Assert (queue_length > 0,
            ExcMessage ("The queue length must be at least one, and preferably "
                        "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert (chunk_size > 0,
            ExcMessage ("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode

    // we want to use TBB if we have support and if it is not disabled at
    // runtime:
#ifdef DEAL_II_WITH_THREADS
    if (multithread_info.n_threads()==1)
#endif
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data;

        for (unsigned int color=0; color<colored_iterators.size(); ++color)
          for (typename std::vector<Iterator>::const_iterator p = colored_iterators[color].begin();
               p != colored_iterators[color].end(); ++p)
            {
              // need to check if the function is not the zero function. To
              // check zero-ness, create a C++ function out of it and check that
              if (static_cast<const std_cxx11::function<void (const Iterator &,
                                                              ScratchData &,
                                                              CopyData &)>& >(worker))
                worker (*p, scratch_data, copy_data);
              if (static_cast<const std_cxx11::function<void (const CopyData &)>& >(copier))
                copier (copy_data);
            }
      }
#ifdef DEAL_II_WITH_THREADS
    else // have TBB and use more than one thread
      {
        // loop over the various colors of what we're given
        for (unsigned int color=0; color<colored_iterators.size(); ++color)
          if (colored_iterators[color].size() > 0)
            {
              if (static_cast<const std_cxx11::function<void (const CopyData &)>& >(copier))
                {
                  // there is a copier function, so we have to go with
                  // the full three-stage design of the pipeline
                  internal::Implementation3::IteratorRangeToItemStream<Iterator,ScratchData,CopyData>
                  iterator_range_to_item_stream (colored_iterators[color].begin(),
                                                 colored_iterators[color].end(),
                                                 queue_length,
                                                 chunk_size,
                                                 sample_scratch_data,
                                                 sample_copy_data);


                  internal::Implementation3::WorkerAndCopier<Iterator, ScratchData, CopyData>
                  worker_and_copier_filter (worker, copier);

                  // now create a pipeline from these stages
                  tbb::pipeline assembly_line;
                  assembly_line.add_filter (iterator_range_to_item_stream);
                  assembly_line.add_filter (worker_and_copier_filter);

                  // and run it
                  assembly_line.run (queue_length);

                  assembly_line.clear ();
                }
              else
                {
                  // no copier function, we can implement things as a parallel for
                  Assert (static_cast<const std_cxx11::function<void (const Iterator &,
                                                                      ScratchData &,
                                                                      CopyData &)>& >(worker),
                          ExcMessage ("It makes no sense to call this function with "
                                      "empty functions for both the worker and the "
                                      "copier!"));

                  typedef
                  internal::ParallelFor::Worker<Iterator,ScratchData,CopyData>
                  ParallelForWorker;

                  typedef
                  typename std::vector<Iterator>::const_iterator
                  RangeType;

                  ParallelForWorker parallel_for_worker (worker,
                                                         sample_scratch_data,
                                                         sample_copy_data);

                  tbb::parallel_for (tbb::blocked_range<RangeType>
                                     (colored_iterators[color].begin(),
                                      colored_iterators[color].end(),
                                      /*grain_size=*/chunk_size),
                                     std_cxx11::bind (&ParallelForWorker::operator(),
                                                      std_cxx11::ref(parallel_for_worker),
                                                      std_cxx11::_1),
                                     tbb::auto_partitioner());
                }
            }
      }
#endif
  }





  /**
   * This is a variant of one of the two main functions of the WorkStream
   * concept, doing work as described in the introduction to this namespace.
   * It corresponds to implementation 2 of the paper by Turcksin, Kronbichler
   * and Bangerth (see @ref workstream_paper).
   *
   * This is the function that can be used for worker and copier functions
   * that are member functions of a class. If the copier is an empty function,
   * it is ignored in the pipeline.
   *
   * The argument passed as @p end must be convertible to the same type as @p
   * begin, but doesn't have to be of the same type itself. This allows to
   * write code like <code>WorkStream().run(dof_handler.begin_active(),
   * dof_handler.end(), ...</code> where the first is of type
   * DoFHandler::active_cell_iterator whereas the second is of type
   * DoFHandler::raw_cell_iterator.
   *
   * The @p queue_length argument indicates the number of items that can be
   * live at any given time. Each item consists of @p chunk_size elements of
   * the input stream that will be worked on by the worker and copier
   * functions one after the other on the same thread.
   *
   * @note If your data objects are large, or their constructors are
   * expensive, it is helpful to keep in mind that <tt>queue_length</tt>
   * copies of the <tt>ScratchData</tt> object and
   * <tt>queue_length*chunk_size</tt> copies of the <tt>CopyData</tt> object
   * are generated.
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
       const ScratchData                       &sample_scratch_data,
       const CopyData                          &sample_copy_data,
       const unsigned int queue_length =        2*multithread_info.n_threads(),
       const unsigned int chunk_size =          8)
  {
    // forward to the other function
    run (begin, end,
         std_cxx11::bind (worker,
                          std_cxx11::ref (main_object),
                          std_cxx11::_1, std_cxx11::_2, std_cxx11::_3),
         std_cxx11::bind (copier,
                          std_cxx11::ref (main_object),
                          std_cxx11::_1),
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
