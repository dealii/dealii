/*
    Copyright 2005-2010 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#ifndef _TBB_arena_H
#define _TBB_arena_H

#include "tbb/tbb_stddef.h"
#include "tbb/atomic.h"

#if __TBB_ARENA_PER_MASTER
#include "scheduler_common.h"
#include "market.h"
#include "intrusive_list.h"
#include "task_stream.h"
#else /* !__TBB_ARENA_PER_MASTER */
#include "../rml/include/rml_tbb.h"
#endif /* !__TBB_ARENA_PER_MASTER */

#include "mailbox.h"

namespace tbb {

#if __TBB_ARENA_PER_MASTER
class task_group_context;
class allocate_root_with_context_proxy;
#endif /* __TBB_ARENA_PER_MASTER */

namespace internal {

class governor;
class arena;
class generic_scheduler;
template<typename SchedulerTraits> class custom_scheduler;

#if !__TBB_ARENA_PER_MASTER
//------------------------------------------------------------------------
// UnpaddedArenaPrefix
//------------------------------------------------------------------------

struct WorkerDescriptor {
    //! NULL until worker is published.  -1 if worker should not be published.
    generic_scheduler* scheduler;
};

//! The useful contents of an ArenaPrefix
class UnpaddedArenaPrefix: no_copy, rml::tbb_client {
    friend class generic_scheduler;
    template<typename SchedulerTraits> friend class custom_scheduler;
    friend class arena;
    friend class governor;
    friend struct WorkerDescriptor;

    //! Arena slot to try to acquire first for the next new master.
    unsigned limit;

    //! Number of masters that own this arena.
    /** This may be smaller than the number of masters who have entered the arena. */
    unsigned number_of_masters;

    //! Total number of slots in the arena
    const unsigned number_of_slots;

    //! Number of workers that belong to this arena
    const unsigned number_of_workers;

    //! Pointer to the RML server object that services requests for this arena.
    rml::tbb_server* server;

    //! Counter used to allocate job indices
    tbb::atomic<size_t> next_job_index;

    //! Stack size of worker threads
    size_t stack_size;

    //! Array of workers.
    WorkerDescriptor* worker_list;

#if __TBB_COUNT_TASK_NODES
    //! Net number of nodes that have been allocated from heap.
    /** Updated each time a scheduler is destroyed. */
    atomic<intptr_t> task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */

    //! Estimate of number of available tasks.  
    /** The estimate is either 0 (SNAPSHOT_EMPTY), infinity (SNAPSHOT_FULL), or a special value. 
        The implementation of arena::is_busy_or_empty requires that pool_state_t be unsigned. */
    typedef uintptr_t pool_state_t;

    //! Current estimate of number of available tasks.  
    tbb::atomic<pool_state_t> pool_state;
 
protected:
    UnpaddedArenaPrefix( unsigned number_of_slots_, unsigned number_of_workers_ ) :
        number_of_masters(1),
        number_of_slots(number_of_slots_),
        number_of_workers(number_of_workers_)
    {
#if __TBB_COUNT_TASK_NODES
        task_node_count = 0;
#endif /* __TBB_COUNT_TASK_NODES */
        limit = number_of_workers_;
        server = NULL;
        stack_size = 0;
        next_job_index = 0;
    }
        
private:
    //! Return reference to corresponding arena.
    arena& Arena();

    /*override*/ version_type version() const {
        return 0;
    }

    /*override*/ unsigned max_job_count() const {
        return number_of_workers;
    }

    /*override*/ size_t min_stack_size() const {
        return stack_size;
    }

    /*override*/ policy_type policy() const {
        return throughput;
    }

    /*override*/ job* create_one_job();

    /*override*/ void cleanup( job& j );

    /*override*/ void acknowledge_close_connection();

    /*override*/ void process( job& j );
}; // class UnpaddedArenaPrefix

//------------------------------------------------------------------------
// ArenaPrefix
//------------------------------------------------------------------------

//! The prefix to arena with padding.
class ArenaPrefix: public UnpaddedArenaPrefix {
    //! Padding to fill out to multiple of cache line size.
    char pad[(sizeof(UnpaddedArenaPrefix)/NFS_MaxLineSize+1)*NFS_MaxLineSize-sizeof(UnpaddedArenaPrefix)];

public:
    ArenaPrefix( unsigned number_of_slots_, unsigned number_of_workers_ ) :
        UnpaddedArenaPrefix(number_of_slots_,number_of_workers_)
    {
    }
}; // class ArenaPrefix

#endif /* !__TBB_ARENA_PER_MASTER */

//------------------------------------------------------------------------
// arena_slot
//------------------------------------------------------------------------

struct arena_slot {
#if __TBB_ARENA_PER_MASTER
    //! Scheduler of the thread attached to the slot
    /** Marks the slot as busy, and is used to iterate through the schedulers belonging to this arena **/
    generic_scheduler* my_scheduler;
#endif /* __TBB_ARENA_PER_MASTER */

    // Task pool (the deque of task pointers) of the scheduler that owns this slot
    /** Also is used to specify if the slot is empty or locked:
         0 - empty
        -1 - locked **/
    task** task_pool;

    //! Index of the first ready task in the deque.
    /** Modified by thieves, and by the owner during compaction/reallocation **/
    size_t head;

    //! Padding to avoid false sharing caused by the thieves accessing this slot
    char pad1[NFS_MaxLineSize - sizeof(size_t) - sizeof(task**)
#if __TBB_ARENA_PER_MASTER
              - sizeof(generic_scheduler*)
#endif /* __TBB_ARENA_PER_MASTER */
             ];

    //! Index of the element following the last ready task in the deque.
    /** Modified by the owner thread. **/
    size_t tail;

#if __TBB_ARENA_PER_MASTER
    //! Hints provided for operations with the container of starvation-resistant tasks.
    /** Modified by the owner thread (during these operations). **/
    unsigned hint_for_push, hint_for_pop;

#endif /* __TBB_ARENA_PER_MASTER */

#if __TBB_STATISTICS
    //! Set of counters to accumulate internal statistics related to this arena
    statistics_counters *my_counters;
#endif /* __TBB_STATISTICS */
    //! Padding to avoid false sharing caused by the thieves accessing the next slot
    char pad2[NFS_MaxLineSize - sizeof(size_t)
#if __TBB_ARENA_PER_MASTER
              - 2*sizeof(unsigned)
#endif /* __TBB_ARENA_PER_MASTER */
#if __TBB_STATISTICS
              - sizeof(statistics_counters*)
#endif /* __TBB_STATISTICS */
             ];
}; // class arena_slot

//------------------------------------------------------------------------
// arena
//------------------------------------------------------------------------

#if __TBB_ARENA_PER_MASTER

//! arena data except the array of slots
/** Separated in order to simplify padding. 
    Intrusive list node base class is used by market to form a list of arenas. **/
struct arena_base : intrusive_list_node {
    //! Market owning this arena
    market* my_market;

    //! Maximal currently busy slot.
    atomic<unsigned> my_limit;

    //! Number of slots in the arena
    unsigned my_num_slots;

    //! Number of workers requested by the master thread owning the arena
    unsigned my_max_num_workers;

    //! Number of workers that are currently requested from the resource manager
    atomic<int> my_num_workers_requested;

    //! Number of workers that have been marked out by the resource manager to service the arena
    unsigned my_num_workers_allotted;

    //! Number of threads in the arena at the moment
    /** Consists of the workers servicing the arena and one master until it starts 
        arena shutdown and detaches from it. Plays the role of the arena's ref count. **/
    atomic<unsigned> my_num_threads_active;

    //! Number of threads that has exited the dispatch loop but has not left the arena yet
    atomic<unsigned> my_num_threads_leaving;

    //! Current task pool state and estimate of available tasks amount.
    /** The estimate is either 0 (SNAPSHOT_EMPTY) or infinity (SNAPSHOT_FULL). 
        Special state is "busy" (any other unsigned value). 
        Note that the implementation of arena::is_busy_or_empty() requires 
        pool_state to be unsigned. */
    tbb::atomic<uintptr_t> pool_state;

#if __TBB_TASK_GROUP_CONTEXT
    //! Pointer to the "default" task_group_context allocated by the arena's master.
    task_group_context* my_master_default_ctx;
#endif

    //! The task pool that guarantees eventual execution even if new tasks are constantly coming.
    task_stream my_task_stream;

    bool my_mandatory_concurrency;

#if TBB_USE_ASSERT
    uintptr_t my_guard;
#endif /* TBB_USE_ASSERT */
}; // struct arena_base

#endif /* __TBB_ARENA_PER_MASTER */

class arena
#if __TBB_ARENA_PER_MASTER
#if (__GNUC__<4 || __GNUC__==4 && __GNUC_MINOR__==0) && !__INTEL_COMPILER
    : public padded<arena_base>
#else
    : padded<arena_base>
#endif
#endif /* __TBB_ARENA_PER_MASTER */
{
    friend class generic_scheduler;
    template<typename SchedulerTraits> friend class custom_scheduler;
    friend class governor;

#if __TBB_ARENA_PER_MASTER
    friend class market;
    friend class tbb::task_group_context;
    friend class allocate_root_with_context_proxy;
    friend class intrusive_list<arena>;

    typedef padded<arena_base> base_type;

    //! Constructor
    arena ( market&, unsigned max_num_workers );

    arena& prefix() const { return const_cast<arena&>(*this); }

    //! Allocate an instance of arena.
    static arena& allocate_arena( market&, unsigned max_num_workers );

#if __TBB_TASK_GROUP_CONTEXT
    //! Propagates cancellation request to all descendants of the context.
    /** The propagation is relayed to the market because tasks created by one 
        master thread can be passed to and executed by other masters. This means 
        that context trees can span several arenas at once and thus cancellation
        propagation cannot be generally localized to one arena only. **/
    void propagate_cancellation ( task_group_context& ctx ) {
        my_market->propagate_cancellation( ctx );
    }
#endif /* __TBB_TASK_GROUP_CONTEXT */

#else /* !__TBB_ARENA_PER_MASTER */

    friend class UnpaddedArenaPrefix;
    friend struct WorkerDescriptor;

    //! Get reference to prefix portion
    ArenaPrefix& prefix() const {return ((ArenaPrefix*)(void*)this)[-1];}

    //! Allocate an instance of arena, and prepare everything to start workers.
    static arena* allocate_arena( unsigned num_slots, unsigned num_workers, size_t stack_size );
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Get reference to mailbox corresponding to given affinity_id.
    mail_outbox& mailbox( affinity_id id ) {
        __TBB_ASSERT( 0<id, "affinity id must be positive integer" );
#if __TBB_ARENA_PER_MASTER
        __TBB_ASSERT( id <= my_num_slots, "affinity id out of bounds" );
#else /* !__TBB_ARENA_PER_MASTER */
        __TBB_ASSERT( id <= prefix().number_of_slots, "id out of bounds" );
#endif /* !__TBB_ARENA_PER_MASTER */

        return ((mail_outbox*)&prefix())[-(int)id];
    }

    //! Completes arena shutdown, destructs and deallocates it.
    void free_arena ();

    typedef uintptr_t pool_state_t;

    //! No tasks to steal since last snapshot was taken
    static const pool_state_t SNAPSHOT_EMPTY = 0;

    //! At least one task has been offered for stealing since the last snapshot started
    static const pool_state_t SNAPSHOT_FULL = pool_state_t(-1);

#if __TBB_ARENA_PER_MASTER
    //! No tasks to steal or snapshot is being taken.
    static bool is_busy_or_empty( pool_state_t s ) { return s < SNAPSHOT_FULL; }

    //! The number of workers active in the arena.
    unsigned num_workers_active( ) {
        return my_num_threads_active - my_num_threads_leaving - (slot[0].my_scheduler? 1: 0);
    }

    //! If necessary, raise a flag that there is new job in arena.
    template<bool Spawned> void advertise_new_work();
#else /*__TBB_ARENA_PER_MASTER*/
    //! Server is going away and hence further calls to adjust_job_count_estimate are unsafe.
    static const pool_state_t SNAPSHOT_SERVER_GOING_AWAY = pool_state_t(-2);

    //! No tasks to steal or snapshot is being taken.
    static bool is_busy_or_empty( pool_state_t s ) { return s < SNAPSHOT_SERVER_GOING_AWAY; }

    //! If necessary, raise a flag that task was added to pool recently.
    inline void mark_pool_full();
#endif /* __TBB_ARENA_PER_MASTER */

    //! Check if there is job anywhere in arena.
    /** Return true if no job or if arena is being cleaned up. */
    bool is_out_of_work();

    //! Initiates arena shutdown.
    void close_arena ();

#if __TBB_ARENA_PER_MASTER
    //! Registers the worker with the arena and enters TBB scheduler dispatch loop
    void process( generic_scheduler& s );

#if __TBB_STATISTICS
    //! Outputs internal statistics accumulated by the arena
    void dump_arena_statistics ();
#endif /* __TBB_STATISTICS */
#endif /* __TBB_ARENA_PER_MASTER */

#if __TBB_COUNT_TASK_NODES
    //! Returns the number of task objects "living" in worker threads
    intptr_t workers_task_node_count();
#endif

    /** Must be the last data field */
    arena_slot slot[1];
}; // class arena


#if __TBB_ARENA_PER_MASTER
template<bool Spawned> void arena::advertise_new_work() {
    if( !Spawned ) { // i.e. the work was enqueued
        if( my_max_num_workers==0 ) {
            my_max_num_workers = 1;
            my_mandatory_concurrency = true;
            prefix().pool_state = SNAPSHOT_FULL;
            my_market->adjust_demand( *this, 1 );
            return;
        }
        // Local memory fence is required to avoid missed wakeups; see the comment below.
        // Starvation resistant tasks require mandatory concurrency, so missed wakeups are unacceptable.
        __TBB_full_memory_fence(); 
    }
    // Double-check idiom that, in case of spawning, is deliberately sloppy about memory fences.
    // Technically, to avoid missed wakeups, there should be a full memory fence between the point we 
    // released the task pool (i.e. spawned task) and read the arena's state.  However, adding such a 
    // fence might hurt overall performance more than it helps, because the fence would be executed 
    // on every task pool release, even when stealing does not occur.  Since TBB allows parallelism, 
    // but never promises parallelism, the missed wakeup is not a correctness problem.
    pool_state_t snapshot = prefix().pool_state;
    if( is_busy_or_empty(snapshot) ) {
        // Attempt to mark as full.  The compare_and_swap below is a little unusual because the 
        // result is compared to a value that can be different than the comparand argument.
        if( prefix().pool_state.compare_and_swap( SNAPSHOT_FULL, snapshot )==SNAPSHOT_EMPTY ) {
            if( snapshot!=SNAPSHOT_EMPTY ) {
                // This thread read "busy" into snapshot, and then another thread transitioned 
                // pool_state to "empty" in the meantime, which caused the compare_and_swap above 
                // to fail.  Attempt to transition pool_state from "empty" to "full".
                if( prefix().pool_state.compare_and_swap( SNAPSHOT_FULL, SNAPSHOT_EMPTY )!=SNAPSHOT_EMPTY ) {
                    // Some other thread transitioned pool_state from "empty", and hence became
                    // responsible for waking up workers.
                    return;
                }
            }
            // This thread transitioned pool from empty to full state, and thus is responsible for
            // telling RML that there is work to do.
            if( Spawned ) {
                if( my_mandatory_concurrency ) {
                    __TBB_ASSERT(my_max_num_workers==1, "");
                    // There was deliberate oversubscription on 1 core for sake of starvation-resistant tasks.
                    // Now a single active thread (must be the master) supposedly starts a new parallel region
                    // with relaxed sequential semantics, and oversubscription should be avoided.
                    // Demand for workers has been decreased to 0 during SNAPSHOT_EMPTY, so just keep it.
                    my_max_num_workers = 0;
                    my_mandatory_concurrency = false;
                    return;
                }
            }
            my_market->adjust_demand( *this, my_max_num_workers );
        }
    }
}
#else /* !__TBB_ARENA_PER_MASTER */
inline void arena::mark_pool_full()  {
    // Double-check idiom that is deliberately sloppy about memory fences.
    // Technically, to avoid missed wakeups, there should be a full memory fence between the point we 
    // released the task pool (i.e. spawned task) and read the arena's state.  However, adding such a 
    // fence might hurt overall performance more than it helps, because the fence would be executed 
    // on every task pool release, even when stealing does not occur.  Since TBB allows parallelism, 
    // but never promises parallelism, the missed wakeup is not a correctness problem.
    pool_state_t snapshot = prefix().pool_state;
    if( is_busy_or_empty(snapshot) ) {
        // Attempt to mark as full.  The compare_and_swap below is a little unusual because the 
        // result is compared to a value that can be different than the comparand argument.
        if( prefix().pool_state.compare_and_swap( SNAPSHOT_FULL, snapshot )==SNAPSHOT_EMPTY ) {
            if( snapshot!=SNAPSHOT_EMPTY ) {
                // This thread read "busy" into snapshot, and then another thread transitioned 
                // pool_state to "empty" in the meantime, which caused the compare_and_swap above 
                // to fail.  Attempt to transition pool_state from "empty" to "full".
                if( prefix().pool_state.compare_and_swap( SNAPSHOT_FULL, SNAPSHOT_EMPTY )!=SNAPSHOT_EMPTY ) {
                    // Some other thread transitioned pool_state from "empty", and hence became
                    // responsible for waking up workers.
                    return;
                }
            }
            // This thread transitioned pool from empty to full state, and thus is responsible for
            // telling RML that there is work to do.
            prefix().server->adjust_job_count_estimate( int(prefix().number_of_workers) );
        }
    }
}
#endif /* !__TBB_ARENA_PER_MASTER */

} // namespace internal
} // namespace tbb

#endif /* _TBB_arena_H */
