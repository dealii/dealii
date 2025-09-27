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

#ifndef dealii_work_stream_h
#  define dealii_work_stream_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/graph_coloring.h>
#  include <deal.II/base/iterator_range.h>
#  include <deal.II/base/multithread_info.h>
#  include <deal.II/base/parallel.h>
#  include <deal.II/base/std_cxx20/type_traits.h>
#  include <deal.II/base/template_constraints.h>
#  include <deal.II/base/thread_local_storage.h>
#  include <deal.II/base/thread_management.h>

#  ifdef DEAL_II_WITH_TBB
#    ifdef DEAL_II_TBB_WITH_ONEAPI
#      include <tbb/parallel_pipeline.h>
#    else
#      include <tbb/pipeline.h>
#    endif
#    include <tbb/blocked_range.h>
#  endif

#  ifdef DEAL_II_WITH_TASKFLOW
#    include <taskflow/taskflow.hpp>
#  endif

#  include <functional>
#  include <iterator>
#  include <list>
#  include <memory>
#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN



/**
 * A namespace whose main template function supports running multiple threads
 * each of which operates on a subset of the given range of objects. The class
 * uses the Intel Threading Building Blocks (TBB) to load balance the
 * individual subranges onto the available threads. For a lengthy discussion
 * of the rationale of this class, see the
 * @ref threads "Parallel computing with multiple processors"
 * topic. It is used in the tutorial first in step-9, and again in step-13,
 * step-14, step-32 and others.
 *
 * The class is built on the following premise: One frequently has some work
 * that needs to be done on a sequence of objects; a prototypical example is
 * assembling cell contributions to a system matrix or right hand side. In
 * many such examples, part of the work can be done entirely independently and
 * in parallel, possibly using several processor cores on a machine with
 * shared memory. However, some other part of this work may need to be
 * synchronized and be done in order. In the example of assembling a matrix,
 * the computation of local contributions can be done entirely in parallel,
 * but copying the local contributions into the global matrix requires some
 * care: First, several threads can't write at the same time, but need to
 * synchronize writing using a mutex; secondly, we want the order in which
 * local contributions are added to the global matrix to be always the same
 * because floating point addition is not commutative and adding local
 * contributions to the global matrix in different orders leads to subtly
 * different results that can affect the number of iterations for iterative
 * solvers as well as the round-off error in the solution in random ways.
 * Consequently, we want to ensure that only one thread at a time writes into
 * the global matrix, and that results are copied in a stable and reproducible
 * order.
 *
 * This class implements a framework for this work model. It works with a
 * stream of objects given by an iterator range, runs a worker function in
 * parallel on all of these objects and then passes each object to a
 * postprocessor function that runs sequentially and gets objects in exactly
 * the order in which they appear in the input iterator range. None of the
 * synchronization work is exposed to the user of this class.
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
 * cores - the default is four times the number of cores on the current
 * system.
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
 * ScratchData, which designates a type objects of which are stored with each
 * item and which threads can use as private data without having to share them
 * with other threads. The run() function takes an additional argument with an
 * object of type ScratchData that is going to be copied for the arguments
 * passed to each of the worker functions.
 *
 * In addition, worker functions store their results in objects of template
 * type CopyData. These are then handed off to a separate function, called
 * copier, that may use the stored results to transfer them into permanent
 * storage. For example, it may copy the results of local contributions to a
 * matrix computed by a worker function into the global matrix. In contrast to
 * the worker function, however, only one instance of the copier is run at any
 * given time; it can therefore safely copy local contributions into the
 * global matrix without the need to lock the global object using a mutex or
 * similar means. Furthermore, it is guaranteed that the copier is run with
 * CopyData objects in the same order in which their associated items were
 * created; consequently, even if worker threads may compute results in
 * unspecified order, the copier always receives the results in exactly the
 * same order as the items were created.
 *
 * Once an item is processed by the copier, it is deleted and the ScratchData
 * and CopyData objects that were used in its computation are considered
 * unused and may be re-used for the next invocation of the worker function,
 * on this or another thread. However, the WorkStream functions make no
 * attempt to reset these objects to any kind of pristine state -- a worker
 * should assume that the CopyData object it gets handed has prior content
 * and clear it first in whatever manner seems appropriate, before putting
 * content into it that can later be processed again by the copier.
 *
 * The member variables in ScratchData and CopyData can be accessed
 * independently of other concurrent uses of copies of these data structures.
 * Therefore, it is perfectly fine to resize auxiliary data structures
 * associated with ScratchData and CopyData to different lengths on each cell.
 * For example, a vector holding densities at each quadrature point which is
 * used with LocalIntegrators::L2::weighted_mass_matrix() to assemble the local
 * matrix could be resized to the corresponding number of quadrature points of
 * the current cell in DoFHandlers with hp-capabilities. Similarly, local
 * stiffness matrix in CopyData can be resized in accordance with the number of
 * local DoFs on the current cell.
 *
 * @note For integration over cells and faces, it is often useful to use
 * methods more specific to the task than the current function (which doesn't
 * care whether the iterators are over cells, vector elements, or any other
 * kind of range). An implementation of an interface specifically suited to
 * integration is the MeshWorker::mesh_loop() function.
 *
 * @note The functions in this namespace only really work in parallel when
 * multithread mode was selected during deal.II configuration. Otherwise they
 * simply work on each item sequentially.
 *
 * @ingroup threads
 */
namespace WorkStream
{
  /**
   * The nested namespaces contain various implementations of the workstream
   * algorithms.
   */
  namespace internal
  {
    /**
     * A structure that contains a pointer to a scratch data object and a copy
     * data object along with a flag that indicates whether this object is
     * currently in use.
     */
    template <typename Iterator, typename ScratchData, typename CopyData>
    struct ScratchAndCopyDataObjects
    {
      std::unique_ptr<ScratchData> scratch_data;
      std::unique_ptr<CopyData>    copy_data;
      bool                         currently_in_use;

      /**
       * Default constructor.
       */
      ScratchAndCopyDataObjects()
        : currently_in_use(false)
      {}

      ScratchAndCopyDataObjects(std::unique_ptr<ScratchData> &&p,
                                std::unique_ptr<CopyData>    &&q,
                                const bool                     in_use)
        : scratch_data(std::move(p))
        , copy_data(std::move(q))
        , currently_in_use(in_use)
      {}

      // Provide a copy constructor that actually doesn't copy the
      // internal state. This makes handling ScratchAndCopyDataObjects
      // easier to handle with STL containers.
      ScratchAndCopyDataObjects(const ScratchAndCopyDataObjects &)
        : currently_in_use(false)
      {}
    };

    /**
     * A structure that contains a pointer to a scratch data object
     * along with a flag that indicates whether this object is currently
     * in use.
     */
    template <typename ScratchData>
    struct ScratchDataObject
    {
      std::unique_ptr<ScratchData> scratch_data;
      bool                         currently_in_use;

      /**
       * Default constructor.
       */
      ScratchDataObject()
        : currently_in_use(false)
      {}

      ScratchDataObject(std::unique_ptr<ScratchData> &&p, const bool in_use)
        : scratch_data(std::move(p))
        , currently_in_use(in_use)
      {}

      ScratchDataObject(ScratchData *p, const bool in_use)
        : scratch_data(p)
        , currently_in_use(in_use)
      {}

      // Provide a copy constructor that actually doesn't copy the
      // internal state. This makes handling ScratchAndCopyDataObjects
      // easier to handle with STL containers.
      ScratchDataObject(const ScratchDataObject &)
        : currently_in_use(false)
      {}

      ScratchDataObject(ScratchDataObject &&o) noexcept = default;
    };

#  ifdef DEAL_II_WITH_TBB
    /**
     * A namespace for the implementation of details of the WorkStream pattern
     * and function. This namespace holds classes that deal with the second
     * implementation described in the paper by Turcksin, Kronbichler and
     * Bangerth (see
     * @ref workstream_paper).
     * Here, no coloring is provided, so copying is done sequentially using a
     * TBB filter.
     *
     * Even though this implementation is slower than the third implementation
     * discussed in that paper, we need to keep it around for two reasons: (i)
     * a user may not give us a graph coloring, (ii) we want to use this
     * implementation for colors that are just too small.
     */
    namespace tbb_no_coloring
    {
      /**
       * A class that creates a sequence of items from a range of iterators.
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class IteratorRangeToItemStream
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
           * Typedef to a list of scratch data objects. The rationale for this
           * list is provided in the variables that use these lists.
           */
          using ScratchDataList = std::list<ScratchDataObject<ScratchData>>;

          /**
           * A list of iterators that need to be worked on. Only the first
           * n_iterators are relevant.
           */
          std::vector<Iterator> iterators;

          /**
           * The CopyData objects that the Worker part of the pipeline fills
           * for each work item. Again, only the first n_iterators elements are
           * what we care about.
           */
          std::vector<CopyData> copy_datas;

          /**
           * Number of items identified by the iterators array that the
           * Worker and Copier pipeline stage need to work on. The maximum
           * value of this variable will be chunk_size.
           */
          unsigned int n_iterators;

          /**
           * Pointer to a thread local variable identifying the scratch data
           * objects this thread will use. The initial implementation of this
           * class using thread local variables provided only a single scratch
           * object per thread. This doesn't work, because the worker
           * functions may start tasks itself and then call
           * Threads::TaskGroup::join_all() or a similar function, which the
           * TBB scheduler may use to run something else on the current thread
           * -- for example another instance of the worker function.
           * Consequently, there would be two instances of the worker function
           * that use the same scratch object if we only provided a single
           * scratch object per thread. The solution is to provide a list of
           * scratch objects for each thread, together with a flag indicating
           * whether this scratch object is currently used. If a thread needs
           * a scratch object, it walks this list until it finds an unused
           * object, or, if there is none, creates one itself. Note that we
           * need not use synchronization primitives for this process since
           * the lists are thread-local and we are guaranteed that only a
           * single thread accesses them as long as we have no yield point in
           * between the accesses to the list.
           *
           * The pointers to scratch objects stored in each of these lists
           * must be so that they are deleted on all threads when the thread
           * local object is destroyed. This is achieved by using unique_ptr.
           *
           * Note that when a worker needs to create a scratch object, it
           * allocates it using sample_scratch_data to copy from. This has the
           * advantage of a first-touch initialization, i.e., the memory for
           * the scratch data object is allocated and initialized by the same
           * thread that will later use it.
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
           * Default constructor. Initialize everything that doesn't have a
           * default constructor itself.
           */
          ItemType()
            : n_iterators(0)
            , scratch_data(nullptr)
            , sample_scratch_data(nullptr)
            , currently_in_use(false)
          {}
        };


        /**
         * Constructor. Take an iterator range, the size of a buffer that can
         * hold items, and the sample additional data object that will be
         * passed to each worker and copier function invocation.
         */
        IteratorRangeToItemStream(const Iterator    &begin,
                                  const Iterator    &end,
                                  const unsigned int buffer_size,
                                  const unsigned int chunk_size,
                                  const ScratchData &sample_scratch_data,
                                  const CopyData    &sample_copy_data)
          : remaining_iterator_range(begin, end)
          , item_buffer(buffer_size)
          , sample_scratch_data(sample_scratch_data)
          , chunk_size(chunk_size)
        {
          // initialize the elements of the ring buffer
          for (auto &item : item_buffer)
            {
              Assert(item.n_iterators == 0, ExcInternalError());

              item.iterators.resize(chunk_size,
                                    remaining_iterator_range.second);
              item.scratch_data        = &thread_local_scratch;
              item.sample_scratch_data = &sample_scratch_data;
              item.copy_datas.resize(chunk_size, sample_copy_data);
              item.currently_in_use = false;
            }
        }


        /**
         * Create an item and return a pointer to it.
         */
        ItemType *
        get_item()
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
          ItemType *current_item = nullptr;
          for (unsigned int i = 0; i < item_buffer.size(); ++i)
            if (item_buffer[i].currently_in_use == false)
              {
                item_buffer[i].currently_in_use = true;
                current_item                    = &item_buffer[i];
                break;
              }
          Assert(current_item != nullptr,
                 ExcMessage("This can't be. There must be a free item!"));

          // initialize the next item. it may
          // consist of at most chunk_size
          // elements
          current_item->n_iterators = 0;
          while ((remaining_iterator_range.first !=
                  remaining_iterator_range.second) &&
                 (current_item->n_iterators < chunk_size))
            {
              current_item->iterators[current_item->n_iterators] =
                remaining_iterator_range.first;

              ++remaining_iterator_range.first;
              ++current_item->n_iterators;
            }

          if (current_item->n_iterators == 0)
            // there were no items
            // left. terminate the pipeline
            return nullptr;
          else
            return current_item;
        }

      private:
        /**
         * The interval of iterators still to be worked on. This range will
         * shrink over time.
         */
        std::pair<Iterator, Iterator> remaining_iterator_range;

        /**
         * A buffer that will store items.
         */
        std::vector<ItemType> item_buffer;

        /**
         * Pointer to a thread local variable identifying the scratch data
         * objects this thread will use. The initial implementation of this
         * class using thread local variables provided only a single scratch
         * object per thread. This doesn't work, because the worker functions
         * may start tasks itself and then call Threads::TaskGroup::join_all()
         * or a similar function, which the TBB scheduler may use to run
         * something else on the current thread -- for example another
         * instance of the worker function. Consequently, there would be two
         * instances of the worker function that use the same scratch object
         * if we only provided a single scratch object per thread. The
         * solution is to provide a list of scratch objects for each thread,
         * together with a flag indicating whether this scratch object is
         * currently used. If a thread needs a scratch object, it walks this
         * list until it finds an unused object, or, if there is none, creates
         * one itself. Note that we need not use synchronization primitives
         * for this process since the lists are thread-local and we are
         * guaranteed that only a single thread accesses them as long as we
         * have no yield point in between the accesses to the list.
         *
         * The pointers to scratch objects stored in each of these lists must
         * be so that they are deleted on all threads when the thread local
         * object is destroyed. This is achieved by using unique_ptr.
         *
         * Note that when a worker needs to create a scratch object, it
         * allocates it using sample_scratch_data to copy from. This has the
         * advantage of a first-touch initialization, i.e., the memory for the
         * scratch data object is allocated and initialized by the same thread
         * that will later use it.
         */
        Threads::ThreadLocalStorage<typename ItemType::ScratchDataList>
          thread_local_scratch;

        /**
         * A reference to a sample scratch data that will be used to
         * initialize the thread-local pointers to a scratch data object each
         * of the worker tasks uses.
         */
        const ScratchData &sample_scratch_data;

        /**
         * Number of elements of the iterator range that each thread should
         * work on sequentially; a large number makes sure that each thread
         * gets a significant amount of work before the next task switch
         * happens, whereas a small number is better for load balancing.
         */
        const unsigned int chunk_size;
      };



      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const Iterator                             &begin,
          const std_cxx20::type_identity_t<Iterator> &end,
          Worker                                      worker,
          Copier                                      copier,
          const ScratchData                          &sample_scratch_data,
          const CopyData                             &sample_copy_data,
          const unsigned int                          queue_length,
          const unsigned int                          chunk_size)
      {
        using ItemType = typename IteratorRangeToItemStream<Iterator,
                                                            ScratchData,
                                                            CopyData>::ItemType;

        // Define the three stages of the pipeline:

        //
        // ----- Stage 1 -----
        //
        // The first stage is the one that provides us with chunks of data
        // to work on (the stream of "items"). This stage will run sequentially.
        IteratorRangeToItemStream<Iterator, ScratchData, CopyData>
             iterator_range_to_item_stream(begin,
                                        end,
                                        queue_length,
                                        chunk_size,
                                        sample_scratch_data,
                                        sample_copy_data);
        auto item_generator = [&](tbb::flow_control &fc) -> ItemType * {
          if (const auto item = iterator_range_to_item_stream.get_item())
            return item;
          else
            {
              fc.stop();
              return nullptr;
            }
        };

        //
        // ----- Stage 2 -----
        //
        // The second stage is the one that does the actual work. This is the
        // stage that runs in parallel
        auto item_worker =
          [worker =
             std::function<void(const Iterator &, ScratchData &, CopyData &)>(
               worker),
           copier_exists =
             static_cast<bool>(std::function<void(const CopyData &)>(copier))](
            ItemType *current_item) {
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
            ScratchData *scratch_data = nullptr;
            {
              // see if there is an unused object. if so, grab it and mark
              // it as used
              for (auto &p : current_item->scratch_data->get())
                if (p.currently_in_use == false)
                  {
                    scratch_data       = p.scratch_data.get();
                    p.currently_in_use = true;

                    break;
                  }

              // if no object was found, create one and mark it as used
              if (scratch_data == nullptr)
                {
                  scratch_data =
                    new ScratchData(*current_item->sample_scratch_data);
                  current_item->scratch_data->get().emplace_back(scratch_data,
                                                                 true);
                }
            };

            // then call the worker function on each element of the chunk we
            // were given. since these worker functions are called on separate
            // threads, nothing good can happen if they throw an exception and
            // we are best off catching it and showing an error message
            for (unsigned int i = 0; i < current_item->n_iterators; ++i)
              {
                try
                  {
                    if (worker)
                      worker(current_item->iterators[i],
                             *scratch_data,
                             current_item->copy_datas[i]);
                  }
                catch (const std::exception &exc)
                  {
                    Threads::internal::handle_std_exception(exc);
                  }
                catch (...)
                  {
                    Threads::internal::handle_unknown_exception();
                  }
              }

            // finally mark the scratch object as unused again. as above, there
            // is no need to lock anything here since the object we work on
            // is thread-local
            for (auto &p : current_item->scratch_data->get())
              if (p.scratch_data.get() == scratch_data)
                {
                  Assert(p.currently_in_use == true, ExcInternalError());
                  p.currently_in_use = false;

                  break;
                }

            // if there is no copier, mark current item as usable again
            if (copier_exists == false)
              current_item->currently_in_use = false;


            // Then return the original pointer
            // to the now modified object. The copier will work on it next.
            return current_item;
          };

        //
        // ----- Stage 3 -----
        //
        // The last stage is the one that copies data from the CopyData objects
        // to the final destination. This stage runs sequentially again.
        auto item_copier = [copier = std::function<void(const CopyData &)>(
                              copier)](ItemType *current_item) {
          if (copier)
            {
              // Initiate copying data. For the same reasons as in the worker
              // class above, catch exceptions rather than letting them
              // propagate into unknown territories:
              for (unsigned int i = 0; i < current_item->n_iterators; ++i)
                {
                  try
                    {
                      copier(current_item->copy_datas[i]);
                    }
                  catch (const std::exception &exc)
                    {
                      Threads::internal::handle_std_exception(exc);
                    }
                  catch (...)
                    {
                      Threads::internal::handle_unknown_exception();
                    }
                }
            }
          // mark current item as usable again
          current_item->currently_in_use = false;
        };


        // Now we just have to set up the pipeline and run it:
        auto tbb_item_stream_filter = tbb::make_filter<void, ItemType *>(
#    ifdef DEAL_II_TBB_WITH_ONEAPI
          tbb::filter_mode::serial_in_order,
#    else
          tbb::filter::serial,
#    endif
          item_generator);

        auto tbb_worker_filter = tbb::make_filter<ItemType *, ItemType *>(
#    ifdef DEAL_II_TBB_WITH_ONEAPI
          tbb::filter_mode::parallel,
#    else
          tbb::filter::parallel,
#    endif
          item_worker);

        auto tbb_copier_filter = tbb::make_filter<ItemType *, void>(
#    ifdef DEAL_II_TBB_WITH_ONEAPI
          tbb::filter_mode::serial_in_order,
#    else
          tbb::filter::serial,
#    endif
          item_copier);

        tbb::parallel_pipeline(queue_length,
                               tbb_item_stream_filter & tbb_worker_filter &
                                 tbb_copier_filter);
      }

    }    // namespace tbb_no_coloring
#  endif // DEAL_II_WITH_TBB



#  ifdef DEAL_II_WITH_TASKFLOW
    /**
     * Mostly a copy of the 2nd implementation of the Workstream paper taking
     * advantage of thread local lists for re-use. Uses taskflow for task
     * scheduling rather than TBB. Currently does not support chunking.
     */


    namespace taskflow_no_coloring
    {
      /**
       * The main run function for the taskflow colorless implementation. The
       * last two arguments in this function are for chunking support which
       * currently does not exist but ideally will later. For now they are
       * ignored but still here to permit existing programs to function
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const Iterator                             &begin,
          const std_cxx20::type_identity_t<Iterator> &end,
          Worker                                      worker,
          Copier                                      copier,
          const ScratchData                          &sample_scratch_data,
          const CopyData                             &sample_copy_data,
          const unsigned int /*queue_length*/ = 2 *
                                                MultithreadInfo::n_threads(),
          const unsigned int chunk_size = 8)

      {
        tf::Executor &executor = MultithreadInfo::get_taskflow_executor();
        tf::Taskflow  taskflow;

        using ScratchDataList = std::list<ScratchDataObject<ScratchData>>;

        Threads::ThreadLocalStorage<ScratchDataList>
          thread_safe_scratch_data_list;

        tf::Task last_copier;


        // A collection of chunk_size copy objects which each represent the
        // contribution of a single worker. Chunk_size workers are thus chunked
        // together by having their outputs stored in this object in a single
        // task.
        struct CopyChunk
        {
          std::vector<CopyData> copy_datas;

          CopyChunk(const unsigned int chunk_size,
                    const CopyData    &sample_copy_data)
            : copy_datas(chunk_size, sample_copy_data)
          {}
        };

        // idx is used to connect each worker to its copier as communication
        // between tasks is not supported. It does this by providing a unique
        // index in the vector of pointers copy_datas at which the copy data
        // object where the work done by work task #idx is stored
        unsigned int idx = 0;

        // A collection of handles to a CopyChunk object for each chunk. The
        // actual data will be allocated only when a worker arrives at this
        // chunk and will be freed as soon as the result has been copied.
        std::vector<std::unique_ptr<CopyChunk>> copy_chunks;

        std::vector<Iterator> chunk;
        chunk.reserve(chunk_size);

        unsigned int total_chunks  = 0;
        unsigned int chunk_counter = 0;

        // Generate a static task graph. Here we generate a task for each cell
        // that will be worked on. The tasks are not executed until all of them
        // are created, this code runs sequentially.
        for (Iterator it = begin; it != end;)
          {
            for (; (chunk_counter < chunk_size) && (it != end);
                 ++it, ++chunk_counter)
              chunk.emplace_back(it);
            ++total_chunks;
            // Create a worker task.
            auto worker_task =
              taskflow
                .emplace([chunk,
                          idx,
                          chunk_counter,
                          &thread_safe_scratch_data_list,
                          &sample_scratch_data,
                          &sample_copy_data,
                          &copy_chunks,
                          &worker]() {
                  ScratchData *scratch_data = nullptr;

                  ScratchDataList &scratch_data_list =
                    thread_safe_scratch_data_list.get();
                  // We need to find an unused scratch data object in the list
                  // that corresponds to the current thread and then mark it as
                  // used. if we can't find one, create one. There is no need to
                  // synchronize access to this variable using a mutex since
                  // each object is local to its own thread.
                  for (auto &p : scratch_data_list)
                    {
                      if (p.currently_in_use == false)
                        {
                          scratch_data       = p.scratch_data.get();
                          p.currently_in_use = true;
                          break;
                        }
                    }
                  // If no element in the list was found, create
                  // one and mark it as used.
                  if (scratch_data == nullptr)
                    {
                      scratch_data_list.emplace_back(
                        std::make_unique<ScratchData>(sample_scratch_data),
                        true);
                      scratch_data =
                        scratch_data_list.back().scratch_data.get();
                    }

                  // Create a unique copy chunk object where this
                  // worker's work will be stored.
                  copy_chunks[idx] =
                    std::make_unique<CopyChunk>(chunk_counter,
                                                sample_copy_data);
                  auto         copy_chunk = copy_chunks[idx].get();
                  unsigned int i          = 0;
                  for (auto &it : chunk)
                    {
                      worker(it, *scratch_data, copy_chunk->copy_datas[i]);
                      ++i;
                    }

                  // Find our currently used scratch data and
                  // mark it as unused.
                  for (auto &p : scratch_data_list)
                    {
                      if (p.scratch_data.get() == scratch_data)
                        {
                          Assert(p.currently_in_use == true,
                                 ExcInternalError());
                          p.currently_in_use = false;
                        }
                    }
                })
                .name("worker");

            // Create a copier task. This task is a separate object from the
            // worker task.
            tf::Task copier_task =
              taskflow
                .emplace([idx, &copy_chunks, &copier]() {
                  auto copy_chunk = copy_chunks[idx].get();
                  for (auto &copy_data : copy_chunk->copy_datas)
                    {
                      copier(copy_data);
                    }
                  // Finally free the memory.
                  copy_chunks[idx].reset();
                })
                .name("copy");

            // Ensure the copy task runs after the worker task.
            worker_task.precede(copier_task);

            // Ensure that only one copy task can run at a time. The code
            // below makes each copy task wait until the previous one has
            // finished before it can start
            if (!last_copier.empty())
              last_copier.precede(copier_task);

            // Keep a handle to the last copier. Tasks in taskflow are
            // basically handles to internally stored data, so this does not
            // perform a copy:
            last_copier = copier_task;
            ++idx;
            chunk_counter = 0;
            chunk.clear();
          }
        copy_chunks.resize(total_chunks);
        // Now we run all the tasks in the task graph. They will be run in
        // parallel and are eligible to run when their dependencies established
        // above are met.
        executor.run(taskflow).wait();
      }
    } // namespace taskflow_no_coloring
#  endif

    /**
     * A reference implementation without using multithreading to be used if we
     * don't have multithreading support or if the user requests to run things
     * sequentially. This is more efficient than using TBB or taskflow if we
     * only have a single thread.
     */
    namespace sequential
    {
      /**
       * Sequential version without colors.
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const Iterator                             &begin,
          const std_cxx20::type_identity_t<Iterator> &end,
          Worker                                      worker,
          Copier                                      copier,
          const ScratchData                          &sample_scratch_data,
          const CopyData                             &sample_copy_data)
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data; // NOLINT

        // Optimization: Check if the functions are not the zero function. To
        // check zero-ness, create a C++ function out of it:
        const bool have_worker =
          (static_cast<const std::function<
             void(const Iterator &, ScratchData &, CopyData &)> &>(worker)) !=
          nullptr;
        const bool have_copier =
          (static_cast<const std::function<void(const CopyData &)> &>(
            copier)) != nullptr;

        // Finally loop over all items and perform the necessary work:
        for (Iterator i = begin; i != end; ++i)
          {
            if (have_worker)
              worker(i, scratch_data, copy_data);
            if (have_copier)
              copier(copy_data);
          }
      }



      /**
       * Sequential version with colors
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const std::vector<std::vector<Iterator>> &colored_iterators,
          Worker                                    worker,
          Copier                                    copier,
          const ScratchData                        &sample_scratch_data,
          const CopyData                           &sample_copy_data)
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data; // NOLINT

        // Optimization: Check if the functions are not the zero function. To
        // check zero-ness, create a C++ function out of it:
        const bool have_worker =
          (static_cast<const std::function<
             void(const Iterator &, ScratchData &, CopyData &)> &>(worker)) !=
          nullptr;
        const bool have_copier =
          (static_cast<const std::function<void(const CopyData &)> &>(
            copier)) != nullptr;

        // Finally loop over all items and perform the necessary work:
        for (unsigned int color = 0; color < colored_iterators.size(); ++color)
          if (colored_iterators[color].size() > 0)
            for (auto &it : colored_iterators[color])
              {
                if (have_worker)
                  worker(it, scratch_data, copy_data);
                if (have_copier)
                  copier(copy_data);
              }
      }

    } // namespace sequential



#  ifdef DEAL_II_WITH_TBB
    /**
     * A namespace for the implementation of details of the WorkStream pattern
     * and function. This namespace holds classes that deal with the third
     * implementation described in the paper by Turcksin, Kronbichler and
     * Bangerth (see
     * @ref workstream_paper).
     */
    namespace tbb_colored
    {
      /**
       * A class that manages calling the worker and copier functions. Unlike
       * the other implementations, parallel_for is used instead of a
       * pipeline.
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class WorkerAndCopier
      {
      public:
        /**
         * Constructor.
         */
        WorkerAndCopier(
          const std::function<void(const Iterator &, ScratchData &, CopyData &)>
                                                      &worker,
          const std::function<void(const CopyData &)> &copier,
          const ScratchData                           &sample_scratch_data,
          const CopyData                              &sample_copy_data)
          : worker(worker)
          , copier(copier)
          , sample_scratch_data(sample_scratch_data)
          , sample_copy_data(sample_copy_data)
        {}


        /**
         * The function that calls the worker and the copier functions on a
         * range of items denoted by the two arguments.
         */
        void
        operator()(const tbb::blocked_range<
                   typename std::vector<Iterator>::const_iterator> &range)
        {
          // we need to find an unused scratch and corresponding copy
          // data object in the list that corresponds to the current
          // thread and then mark it as used. If we can't find one,
          // create one as discussed in the discussion of the documentation
          // of the IteratorRangeToItemStream::scratch_data variable,
          // there is no need to synchronize access to this variable
          // using a mutex as long as we have no yield-point in between.
          // This means that we can't take an iterator into the list
          // now and expect it to still be valid after calling the worker,
          // but we at least do not have to lock the following section.
          ScratchData *scratch_data = nullptr;
          CopyData    *copy_data    = nullptr;
          {
            ScratchAndCopyDataList &scratch_and_copy_data_list = data.get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data        = p->scratch_data.get();
                  copy_data           = p->copy_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no element in the list was found, create one and mark it as
            // used
            if (scratch_data == nullptr)
              {
                Assert(copy_data == nullptr, ExcInternalError());

                scratch_and_copy_data_list.emplace_back(
                  std::make_unique<ScratchData>(sample_scratch_data),
                  std::make_unique<CopyData>(sample_copy_data),
                  true);
                scratch_data =
                  scratch_and_copy_data_list.back().scratch_data.get();
                copy_data = scratch_and_copy_data_list.back().copy_data.get();
              }
          }

          // then call the worker and copier functions on each
          // element of the chunk we were given.
          for (typename std::vector<Iterator>::const_iterator p = range.begin();
               p != range.end();
               ++p)
            {
              try
                {
                  if (worker)
                    worker(*p, *scratch_data, *copy_data);
                  if (copier)
                    copier(*copy_data);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception(exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            ScratchAndCopyDataList &scratch_and_copy_data_list = data.get();

            for (typename ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }
        }

      private:
        using ScratchAndCopyDataObjects = typename internal::
          ScratchAndCopyDataObjects<Iterator, ScratchData, CopyData>;

        /**
         * Typedef to a list of scratch data objects. The rationale for this
         * list is provided in the variables that use these lists.
         */
        using ScratchAndCopyDataList = std::list<ScratchAndCopyDataObjects>;

        Threads::ThreadLocalStorage<ScratchAndCopyDataList> data;

        /**
         * Pointer to the function that does the assembling on the sequence of
         * cells.
         */
        const std::function<void(const Iterator &, ScratchData &, CopyData &)>
          worker;

        /**
         * Pointer to the function that does the copying from local
         * contribution to global object.
         */
        const std::function<void(const CopyData &)> copier;

        /**
         * References to sample scratch and copy data for when we need them.
         */
        const ScratchData &sample_scratch_data;
        const CopyData    &sample_copy_data;
      };

      /**
       * The colored run function using TBB.
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const std::vector<std::vector<Iterator>> &colored_iterators,
          Worker                                    worker,
          Copier                                    copier,
          const ScratchData                        &sample_scratch_data,
          const CopyData                           &sample_copy_data,
          const unsigned int                        chunk_size)
      {
        // loop over the various colors of what we're given
        for (unsigned int color = 0; color < colored_iterators.size(); ++color)
          if (colored_iterators[color].size() > 0)
            {
              using WorkerAndCopier = internal::tbb_colored::
                WorkerAndCopier<Iterator, ScratchData, CopyData>;

              WorkerAndCopier worker_and_copier(worker,
                                                copier,
                                                sample_scratch_data,
                                                sample_copy_data);

              parallel::internal::parallel_for(
                colored_iterators[color].begin(),
                colored_iterators[color].end(),
                [&worker_and_copier](
                  const tbb::blocked_range<
                    typename std::vector<Iterator>::const_iterator> &range) {
                  worker_and_copier(range);
                },
                chunk_size);
            }
      }

    }    // namespace tbb_colored
#  endif // DEAL_II_WITH_TBB



#  ifdef DEAL_II_WITH_TASKFLOW
    /**
     * Mostly a copy of the 3rd implementation of the Workstream paper taking
     * advantage of thread local lists for re-use. Uses taskflow for task
     * scheduling rather than TBB. Currently does not support chunking.
     */
    namespace taskflow_colored
    {
      /**
       * The main run function for the taskflow colored implementation. The
       * last two arguments in this function are for chunking support which
       * currently does not exist but ideally will later. For now they are
       * ignored but still here to permit existing programs to function.
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const std::vector<std::vector<Iterator>> &colored_iterators,
          Worker                                    worker,
          Copier                                    copier,
          const ScratchData                        &sample_scratch_data,
          const CopyData                           &sample_copy_data,
          const unsigned int /*queue_length*/ = 2 *
                                                MultithreadInfo::n_threads(),
          const unsigned int chunk_size = 8)

      {
        tf::Executor &executor = MultithreadInfo::get_taskflow_executor();
        using ScratchAndCopyDataObjects = typename internal::
          ScratchAndCopyDataObjects<Iterator, ScratchData, CopyData>;

        using ScratchAndCopyDataList = std::list<ScratchAndCopyDataObjects>;

        Threads::ThreadLocalStorage<ScratchAndCopyDataList>
          thread_safe_scratch_and_copy_data_list;

        tf::Taskflow taskflow;

        // Create a "future" object which eventually contains the execution
        // result of a taskflow graph and can be used to yield execution
        tf::Future<void> execution_future;

        const bool have_worker =
          (static_cast<const std::function<
             void(const Iterator &, ScratchData &, CopyData &)> &>(worker)) !=
          nullptr;
        const bool have_copier =
          (static_cast<const std::function<void(const CopyData &)> &>(
            copier)) != nullptr;

        std::vector<Iterator> chunk;
        chunk.reserve(chunk_size);

        unsigned int chunk_counter = 0;
        // Generate a static task graph. Here we generate a task for each cell
        // that will be worked on. The tasks are not executed until all of them
        // are created, this code runs sequentially. Cells have been grouped
        // into "colors" data from cells in the same color are safe to copy in
        // parallel so copying need not be sequential.
        for (unsigned int color = 0; color < colored_iterators.size(); ++color)
          // Ignore color blocks which are empty.
          if (colored_iterators[color].size() > 0)
            {
              // For each cell queue up a combined worker and copier task. These
              // are not yet run.
              for (auto it = colored_iterators[color].begin();
                   it != colored_iterators[color].end();)
                {
                  for (; (chunk_counter < chunk_size) &&
                         (it != colored_iterators[color].end());
                       ++it, ++chunk_counter)
                    chunk.emplace_back(*it);
                  taskflow
                    .emplace([chunk,
                              have_worker,
                              have_copier,
                              &thread_safe_scratch_and_copy_data_list,
                              &sample_scratch_data,
                              &sample_copy_data,
                              &worker,
                              &copier]() {
                      ScratchData *scratch_data = nullptr;
                      CopyData    *copy_data    = nullptr;

                      ScratchAndCopyDataList &scratch_and_copy_data_list =
                        thread_safe_scratch_and_copy_data_list.get();
                      // See if there is an unused object. if so, grab it
                      // and mark it as used.
                      for (typename ScratchAndCopyDataList::iterator p =
                             scratch_and_copy_data_list.begin();
                           p != scratch_and_copy_data_list.end();
                           ++p)
                        {
                          if (p->currently_in_use == false)
                            {
                              scratch_data        = p->scratch_data.get();
                              copy_data           = p->copy_data.get();
                              p->currently_in_use = true;
                              break;
                            }
                        }
                      // If no element in the list was found, create one and
                      // mark it as used.
                      if (scratch_data == nullptr)
                        {
                          Assert(copy_data == nullptr, ExcInternalError());
                          scratch_and_copy_data_list.emplace_back(
                            std::make_unique<ScratchData>(sample_scratch_data),
                            std::make_unique<CopyData>(sample_copy_data),
                            true);
                          scratch_data = scratch_and_copy_data_list.back()
                                           .scratch_data.get();
                          copy_data =
                            scratch_and_copy_data_list.back().copy_data.get();
                        }

                      for (const Iterator &it : chunk)
                        {
                          if (have_worker)
                            worker(it, *scratch_data, *copy_data);
                          if (have_copier)
                            copier(*copy_data);
                        }

                      // Mark objects as free to be used again.
                      for (typename ScratchAndCopyDataList::iterator p =
                             scratch_and_copy_data_list.begin();
                           p != scratch_and_copy_data_list.end();
                           ++p)
                        {
                          if (p->scratch_data.get() == scratch_data)
                            {
                              Assert(p->currently_in_use == true,
                                     ExcInternalError());
                              p->currently_in_use = false;
                            }
                        }
                    })
                    .name("worker_and_copier");
                  chunk.clear();
                  chunk_counter = 0;
                }

              if (color > 0)
                // Wait for the previous color to finish executing
                execution_future.wait();
              execution_future = executor.run(std::move(taskflow));
            }
        // Wait for our final execution to finish
        if (colored_iterators.size() > 0)
          execution_future.wait();
      }
    }    // namespace taskflow_colored
#  endif // DEAL_II_WITH_TASKFLOW


  } // namespace internal



  /**
   * This is one of two main functions of the WorkStream concept, doing work
   * as described in the introduction to this namespace. It corresponds to
   * implementation 3 of the paper by Turcksin, Kronbichler and Bangerth, see
   * @ref workstream_paper.
   * As such, it takes not a range of iterators described by a begin and end
   * iterator, but a "colored" graph of iterators where each color represents
   * cells for which writing the cell contributions into the global object
   * does not conflict (in other words, these cells are not neighbors). Each
   * "color" is represented by std::vectors of cells. The first argument to
   * this function, a set of sets of cells (which are represent as a vector of
   * vectors, for efficiency), is typically constructed by calling
   * GraphColoring::make_graph_coloring(). See there for more information.
   *
   * This function that can be used for worker and copier objects that are
   * either pointers to non-member functions or objects that allow to be
   * called with an operator(), for example objects created by lambda functions
   * or std::bind.
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
   *
   * @note In case the copier does not do anything, pass
   * `std::function<void(const CopyData &)>()` as @p copier to make sure
   * a more efficient algorithm is used internally. It is important, however,
   * to recognize that the empty function object created above is *not*
   * the same as a lambda function with an empty body,
   * `[](const CopyData &) {}` -- from the perspective of this function,
   * there is no way to recognize whether a lambda function provided as
   * a copier does something or does not do something in its body,
   * and so it needs to be copied. On the other hand, a default-constructed
   * `std::function` object *can* be recognized, and is then used to select
   * a more efficient algorithm.
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const std::vector<std::vector<Iterator>> &colored_iterators,
      Worker                                    worker,
      Copier                                    copier,
      const ScratchData                        &sample_scratch_data,
      const CopyData                           &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8);


  /**
   * This is one of two main functions of the WorkStream concept, doing work
   * as described in the introduction to this namespace. It corresponds to
   * implementation 2 of the paper by Turcksin, Kronbichler and Bangerth (see
   * @ref workstream_paper).
   *
   * This function that can be used for worker and copier objects that are
   * either pointers to non-member functions or objects that allow to be
   * called with an operator(), for example lambda functions
   * or objects created by std::bind. If the copier is an empty function, it is
   * ignored in the pipeline. (However, a lambda function with an empty body is
   * *not* equivalent to an empty `std::function` object and will, consequently,
   * not be ignored.
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
   *
   * @note In case the copier does not do anything, pass
   * `std::function<void(const CopyData &)>()` as @p copier to make sure
   * a more efficient algorithm is used internally. It is important, however,
   * to recognize that the empty function object created above is *not*
   * the same as a lambda function with an empty body,
   * `[](const CopyData &) {}` -- from the perspective of this function,
   * there is no way to recognize whether a lambda function provided as
   * a copier does something or does not do something in its body,
   * and so it needs to be copied. On the other hand, a default-constructed
   * `std::function` object *can* be recognized, and is then used to select
   * a more efficient algorithm.
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator                             &begin,
      const std_cxx20::type_identity_t<Iterator> &end,
      Worker                                      worker,
      Copier                                      copier,
      const ScratchData                          &sample_scratch_data,
      const CopyData                             &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    Assert(queue_length > 0,
           ExcMessage("The queue length must be at least one, and preferably "
                      "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert(chunk_size > 0, ExcMessage("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode

    // If no work then skip. (only use operator!= for iterators since we may
    // not have an equality comparison operator)
    if (!(begin != end))
      return;

    if (MultithreadInfo::n_threads() > 1)
      {
#  if defined(DEAL_II_WITH_TBB) || defined(DEAL_II_WITH_TASKFLOW)
        if (static_cast<const std::function<void(const CopyData &)> &>(copier))
          {
            // If we have a copier, run the algorithm:
#    if defined(DEAL_II_WITH_TASKFLOW)
            internal::taskflow_no_coloring::run(begin,
                                                end,
                                                worker,
                                                copier,
                                                sample_scratch_data,
                                                sample_copy_data,
                                                queue_length,
                                                chunk_size);
#    elif defined(DEAL_II_WITH_TBB)
            internal::tbb_no_coloring::run(begin,
                                           end,
                                           worker,
                                           copier,
                                           sample_scratch_data,
                                           sample_copy_data,
                                           queue_length,
                                           chunk_size);
#    endif
          }
        else
          {
            // There is no copier function. in this case, we have an
            // embarrassingly parallel problem where we can
            // essentially apply parallel_for. because parallel_for
            // requires subdividing the range for which operator- is
            // necessary between iterators, it is often inefficient to
            // apply it directly to cell ranges and similar iterator
            // types for which operator- is expensive or, in fact,
            // nonexistent. rather, in that case, we simply copy the
            // iterators into a large array and use operator- on
            // iterators to this array of iterators.
            //
            // instead of duplicating code, this is essentially the
            // same situation we have in the colored implementation below, so we
            // just defer to that place
            std::vector<std::vector<Iterator>> all_iterators(1);
            for (Iterator p = begin; p != end; ++p)
              all_iterators[0].push_back(p);

            run(all_iterators,
                worker,
                copier,
                sample_scratch_data,
                sample_copy_data,
                queue_length,
                chunk_size);
          }

        // exit this function to not run the sequential version below:
        return;
#  endif
      }

    // no TBB or Taskflow installed or we are requested to run sequentially:
    internal::sequential::run(
      begin, end, worker, copier, sample_scratch_data, sample_copy_data);
  }



  /**
   * Same as the function above, but for iterator ranges and C-style arrays.
   * A class that fulfills the requirements of an iterator range defines the
   * functions `IteratorRangeType::begin()` and `IteratorRangeType::end()`,
   * both of which return iterators to elements that form the bounds of the
   * range.
   */
  template <
    typename Worker,
    typename Copier,
    typename IteratorRangeType,
    typename ScratchData,
    typename CopyData,
    typename = std::enable_if_t<
      has_begin_and_end<IteratorRangeType> &&
      !std::is_same_v<IteratorRangeType,
                      IteratorRange<typename IteratorRangeType::iterator>>>>
  void
  run(IteratorRangeType  iterator_range,
      Worker             worker,
      Copier             copier,
      const ScratchData &sample_scratch_data,
      const CopyData    &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(iterator_range.begin(),
        iterator_range.end(),
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  /**
   * Same as the function above, but for deal.II's IteratorRange.
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const IteratorRange<Iterator> &iterator_range,
      Worker                         worker,
      Copier                         copier,
      const ScratchData             &sample_scratch_data,
      const CopyData                &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(iterator_range.begin(),
        iterator_range.end(),
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const std::vector<std::vector<Iterator>> &colored_iterators,
      Worker                                    worker,
      Copier                                    copier,
      const ScratchData                        &sample_scratch_data,
      const CopyData                           &sample_copy_data,
      const unsigned int                        queue_length,
      const unsigned int                        chunk_size)
  {
    Assert(queue_length > 0,
           ExcMessage("The queue length must be at least one, and preferably "
                      "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert(chunk_size > 0, ExcMessage("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode


    if (MultithreadInfo::n_threads() > 1)
      {
#  ifdef DEAL_II_WITH_TASKFLOW
        internal::taskflow_colored::run(colored_iterators,
                                        worker,
                                        copier,
                                        sample_scratch_data,
                                        sample_copy_data,
                                        chunk_size);

        // exit this function to not run the sequential version below:
        return;
#  elif defined(DEAL_II_WITH_TBB)
        internal::tbb_colored::run(colored_iterators,
                                   worker,
                                   copier,
                                   sample_scratch_data,
                                   sample_copy_data,
                                   chunk_size);

        // exit this function to not run the sequential version below:
        return;
#  endif
      }

    // run all colors sequentially:
    {
      internal::sequential::run(colored_iterators,
                                worker,
                                copier,
                                sample_scratch_data,
                                sample_copy_data);
    }
  }



  /**
   * This is a variant of one of the two main functions of the WorkStream
   * concept, doing work as described in the introduction to this namespace.
   * It corresponds to implementation 2 of the paper by Turcksin, Kronbichler
   * and Bangerth (see
   * @ref workstream_paper).
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
   *
   * @note In case the copier does not do anything, pass
   * `std::function<void(const CopyData &)>()` as @p copier to make sure
   * a more efficient algorithm is used internally. It is important, however,
   * to recognize that the empty function object created above is *not*
   * the same as a lambda function with an empty body,
   * `[](const CopyData &) {}` -- from the perspective of this function,
   * there is no way to recognize whether a lambda function provided as
   * a copier does something or does not do something in its body,
   * and so it needs to be copied. On the other hand, a default-constructed
   * `std::function` object *can* be recognized, and is then used to select
   * a more efficient algorithm.
   */
  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator                             &begin,
      const std_cxx20::type_identity_t<Iterator> &end,
      MainClass                                  &main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData    &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // forward to the other function
    run(
      begin,
      end,
      [&main_object, worker](const Iterator &iterator,
                             ScratchData    &scratch_data,
                             CopyData       &copy_data) {
        (main_object.*worker)(iterator, scratch_data, copy_data);
      },
      [&main_object, copier](const CopyData &copy_data) {
        (main_object.*copier)(copy_data);
      },
      sample_scratch_data,
      sample_copy_data,
      queue_length,
      chunk_size);
  }


  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const IteratorOverIterators<Iterator>                             &begin,
      const IteratorOverIterators<std_cxx20::type_identity_t<Iterator>> &end,
      MainClass &main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData    &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // forward to the other function
    run(
      begin,
      end,
      [&main_object, worker](const Iterator &iterator,
                             ScratchData    &scratch_data,
                             CopyData       &copy_data) {
        (main_object.*worker)(iterator, scratch_data, copy_data);
      },
      [&main_object, copier](const CopyData &copy_data) {
        (main_object.*copier)(copy_data);
      },
      sample_scratch_data,
      sample_copy_data,
      queue_length,
      chunk_size);
  }



  /**
   * Same as the function above, but for iterator ranges and C-style arrays.
   * A class that fulfills the requirements of an iterator range defines the
   * functions `IteratorRangeType::begin()` and `IteratorRangeType::end()`,
   * both of which return iterators to elements that form the bounds of the
   * range.
   */
  template <
    typename MainClass,
    typename IteratorRangeType,
    typename ScratchData,
    typename CopyData,
    typename = std::enable_if_t<
      has_begin_and_end<IteratorRangeType> &&
      !std::is_same_v<IteratorRangeType,
                      IteratorRange<typename IteratorRangeType::iterator>>>>
  void
  run(
    IteratorRangeType iterator_range,
    MainClass        &main_object,
    void (MainClass::*worker)(
      const typename std_cxx20::type_identity_t<IteratorRangeType>::iterator &,
      ScratchData &,
      CopyData &),
    void (MainClass::*copier)(const CopyData &),
    const ScratchData &sample_scratch_data,
    const CopyData    &sample_copy_data,
    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(std::begin(iterator_range),
        std::end(iterator_range),
        main_object,
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  /**
   * Same as the function above, but for deal.II's IteratorRange.
   */
  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(IteratorRange<Iterator> iterator_range,
      MainClass              &main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData    &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(std::begin(iterator_range),
        std::end(iterator_range),
        main_object,
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }

} // namespace WorkStream



DEAL_II_NAMESPACE_CLOSE



//----------------------------   work_stream.h     ---------------------------
// end of #ifndef dealii_work_stream_h
#endif
//----------------------------   work_stream.h     ---------------------------
