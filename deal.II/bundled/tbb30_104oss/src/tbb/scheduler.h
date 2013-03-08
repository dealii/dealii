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

#ifndef _TBB_scheduler_H
#define _TBB_scheduler_H

#include "scheduler_common.h"
#include "arena.h"
#include "mailbox.h"
#include "tbb_misc.h" // for FastRandom

#if __TBB_TASK_GROUP_CONTEXT
#include "tbb/spin_mutex.h"
#endif /* __TBB_TASK_GROUP_CONTEXT */

#if __TBB_SURVIVE_THREAD_SWITCH
#include "cilk-tbb-interop.h"
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

namespace tbb {
namespace internal {

template<typename SchedulerTraits> class custom_scheduler;

//------------------------------------------------------------------------
// generic_scheduler
//------------------------------------------------------------------------

#if __TBB_TASK_GROUP_CONTEXT
struct scheduler_list_node_t {
    scheduler_list_node_t *my_prev,
                          *my_next;
};
#endif /* __TBB_TASK_GROUP_CONTEXT */

#define EmptyTaskPool ((task**)0)
#define LockedTaskPool ((task**)~(intptr_t)0)

class governor;

#if __TBB_SCHEDULER_OBSERVER
class task_scheduler_observer_v3;
class observer_proxy;
#endif /* __TBB_SCHEDULER_OBSERVER */

#if __TBB_ARENA_PER_MASTER
class market;
#endif

//! Cilk-style task scheduler.
/** None of the fields here are every read or written by threads other than
    the thread that creates the instance.

    Class generic_scheduler is an abstract base class that contains most of the scheduler,
    except for tweaks specific to processors and tools (e.g. VTune).
    The derived template class custom_scheduler<SchedulerTraits> fills in the tweaks. */
class generic_scheduler: public scheduler, public ::rml::job {
    friend class tbb::task;
#if __TBB_ARENA_PER_MASTER
    friend class market;
#else
    friend class UnpaddedArenaPrefix;
#endif /* !__TBB_ARENA_PER_MASTER */
    friend class arena;
    friend class allocate_root_proxy;
    friend class governor;
#if __TBB_TASK_GROUP_CONTEXT
    friend class allocate_root_with_context_proxy;
    friend class tbb::task_group_context;
#endif /* __TBB_TASK_GROUP_CONTEXT */
#if __TBB_SCHEDULER_OBSERVER
    friend class task_scheduler_observer_v3;
#endif /* __TBB_SCHEDULER_OBSERVER */
    friend class scheduler;
    template<typename SchedulerTraits> friend class custom_scheduler;

    //! If sizeof(task) is <=quick_task_size, it is handled on a free list instead of malloc'd.
    static const size_t quick_task_size = 256-task_prefix_reservation_size;

    static bool is_version_3_task( task& t ) {
        return (t.prefix().extra_state & 0x0F)>=0x1;
    }

    //! Position in the call stack specifying its maximal filling when stealing is still allowed
    uintptr_t my_stealing_threshold;
#if __TBB_ipf
    //! Position in the RSE backup area specifying its maximal filling when stealing is still allowed
    uintptr_t my_rsb_stealing_threshold;
#endif

    static const size_t null_arena_index = ~size_t(0);

    //! Index of the arena slot the scheduler occupies now, or occupied last time.
    size_t arena_index;

    //! Capacity of ready tasks deque (number of elements - pointers to task).
    size_t task_pool_size;

    //! Pointer to the slot in the arena we own at the moment.
    /** When out of arena it points to this scheduler's dummy_slot. **/
    mutable arena_slot* my_arena_slot;

    bool in_arena () const { return my_arena_slot != &dummy_slot; }

    bool local_task_pool_empty () {
        return my_arena_slot->task_pool == EmptyTaskPool || my_arena_slot->head >= my_arena_slot->tail;
    }

#if __TBB_ARENA_PER_MASTER
    //! The market I am in
    market* my_market;

    //! The arena that I own (if master) or am servicing at the moment (if worker)
    arena* my_arena;
#else /* !__TBB_ARENA_PER_MASTER */
    //! The arena that I own (if master) or am servicing at the moment (if worker)
    arena* const my_arena;
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Random number generator used for picking a random victim from which to steal.
    FastRandom random;

    //! Free list of small tasks that can be reused.
    task* free_list;

    //! Innermost task whose task::execute() is running.
    task* innermost_running_task;

    //! Fake root task created by slave threads.
    /** The task is used as the "parent" argument to method wait_for_all. */
    task* dummy_task;

    //! Reference count for scheduler
    /** Number of task_scheduler_init objects that point to this scheduler */
    long ref_count;

    mail_inbox inbox;

    void attach_mailbox( affinity_id id ) {
        __TBB_ASSERT(id>0,NULL);
        inbox.attach( my_arena->mailbox(id) );
        my_affinity_id = id;
    }

    //! The mailbox id assigned to this scheduler.
    /** The id is assigned upon first entry into the arena.
        TODO: how are id's being garbage collected? 
        TODO: master thread may enter arena and leave and then reenter.
                We want to give it the same affinity_id upon reentry, if practical.
      */
    affinity_id my_affinity_id;

    /* A couple of bools can be located here because space is otherwise just padding after my_affinity_id. */

    //! True if this is assigned to thread local storage by registering with governor.
    bool is_registered;

    //! True if *this was created by automatic TBB initialization
    bool is_auto_initialized;

#if __TBB_SCHEDULER_OBSERVER
    //! Last observer_proxy processed by this scheduler
    observer_proxy* local_last_observer_proxy;

    //! Notify any entry observers that have been created since the last call by this thread.
    void notify_entry_observers();
 
    //! Notify all exit observers that this thread is no longer participating in task scheduling.
    void notify_exit_observers( bool is_worker );
#endif /* __TBB_SCHEDULER_OBSERVER */

#if __TBB_COUNT_TASK_NODES
    //! Net number of big task objects that have been allocated but not yet freed.
    intptr_t task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */

    //! Sets up the data necessary for the stealing limiting heuristics
    void init_stack_info ();

    //! Returns true if stealing is allowed
    bool can_steal () {
        int anchor;
#if __TBB_ipf
        return my_stealing_threshold < (uintptr_t)&anchor && (uintptr_t)__TBB_get_bsp() < my_rsb_stealing_threshold;
#else
        return my_stealing_threshold < (uintptr_t)&anchor;
#endif
    }

    //! Actions common to enter_arena and try_enter_arena
    void do_enter_arena();

    //! Used by workers to enter the arena 
    /** Does not lock the task pool in case if arena slot has been successfully grabbed. **/
    void enter_arena();

#if !__TBB_ARENA_PER_MASTER
    //! Used by masters to try to enter the arena
    /** Does not lock the task pool in case if arena slot has been successfully grabbed. **/
    void try_enter_arena();
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Leave the arena
    void leave_arena();

    //! Locks victim's task pool, and returns pointer to it. The pointer can be NULL.
    task** lock_task_pool( arena_slot* victim_arena_slot ) const;

    //! Unlocks victim's task pool
    void unlock_task_pool( arena_slot* victim_arena_slot, task** victim_task_pool ) const;


    //! Locks the local task pool
    void acquire_task_pool() const;

    //! Unlocks the local task pool
    void release_task_pool() const;

    //! Checks if t is affinitized to another thread, and if so, bundles it as proxy.
    /** Returns either t or proxy containing t. **/
    task* prepare_for_spawning( task* t );

    //! Get a task from the local pool.
    /** Called only by the pool owner.
        Returns the pointer to the task or NULL if the pool is empty. 
        In the latter case compacts the pool. **/
    task* get_task();

    //! Attempt to get a task from the mailbox.
    /** Called only by the thread that owns *this.
        Gets a task only if there is one not yet executed by another thread.
        If successful, unlinks the task and returns a pointer to it.
        Otherwise returns NULL. */
    task* get_mailbox_task();

    //! True if t is a task_proxy
    static bool is_proxy( const task& t ) {
        return t.prefix().extra_state==es_task_proxy;
    }

    //! Extracts task pointer from task_proxy, and frees the proxy.
    /** Return NULL if underlying task was claimed by mailbox. */
    task* strip_proxy( task_proxy* result );

#if __TBB_ARENA_PER_MASTER
    //! Get a task from the starvation-resistant task stream of the current arena.
    /** Returns the pointer to the task, or NULL if the attempt was unsuccessful. 
        The latter case does not mean that the stream is drained, however. **/
    task* dequeue_task();

#endif /* __TBB_ARENA_PER_MASTER */
    //! Steal task from another scheduler's ready pool.
    task* steal_task( arena_slot& victim_arena_slot );

    /** Initial size of the task deque sufficient to serve without reallocation
        4 nested parallel_for calls with iteration space of 65535 grains each. **/
    static const size_t min_task_pool_size = 64;

    //! Allocate task pool containing at least n elements.
    task** allocate_task_pool( size_t n );

    //! Deallocate task pool that was allocated by means of allocate_task_pool.
    static void free_task_pool( task** pool ) {
        __TBB_ASSERT( pool, "attempt to free NULL TaskPool" );
        NFS_Free( pool );
    }

    //! Grow ready task deque to at least n elements.
    void grow_task_pool( size_t n );

    //! Initialize a scheduler for a master thread.
    static generic_scheduler* create_master( arena& a );

    //! Perform necessary cleanup when a master thread stops using TBB.
    void cleanup_master();

    //! Initialize a scheduler for a worker thread.
#if __TBB_ARENA_PER_MASTER
    static generic_scheduler* create_worker( market& m, size_t index );
#else /* !__TBB_ARENA_PER_MASTER */
    static generic_scheduler* create_worker( arena& a, size_t index );
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Perform necessary cleanup when a worker thread finishes.
    static void cleanup_worker( void* arg, bool is_worker );

protected:
    generic_scheduler( arena*, size_t index );

#if TBB_USE_ASSERT > 1
    //! Check that internal data structures are in consistent state.
    /** Raises __TBB_ASSERT failure if inconsistency is found. */
    void assert_task_pool_valid() const;
#else
    void assert_task_pool_valid() const {}
#endif /* TBB_USE_ASSERT <= 1 */

public:
    /*override*/ 
    void spawn( task& first, task*& next );

    /*override*/ 
    void spawn_root_and_wait( task& first, task*& next );

#if __TBB_ARENA_PER_MASTER
    /*override*/ 
    void enqueue( task& task_, void* reserved );

    void local_enqueue( task& task_ );
#endif /* __TBB_ARENA_PER_MASTER */

    void local_spawn( task& first, task*& next );
    void local_spawn_root_and_wait( task& first, task*& next );
    virtual void local_wait_for_all( task& parent, task* child ) = 0;

    //! Destroy and deallocate this scheduler object
    void free_scheduler();

    //! Allocate task object, either from the heap or a free list.
    /** Returns uninitialized task object with initialized prefix. */
    task& allocate_task( size_t number_of_bytes, 
                       __TBB_CONTEXT_ARG(task* parent, task_group_context* context) );

    //! Put task on free list.
    /** Does not call destructor. */
    template<free_task_hint h>
    void free_task( task& t );

    void free_task_proxy( task_proxy& tp ) {
#if TBB_USE_ASSERT
        poison_pointer( tp.outbox );
        poison_pointer( tp.next_in_mailbox );
        tp.task_and_tag = 0xDEADBEEF;
#endif /* TBB_USE_ASSERT */
        free_task<small_task>(tp);
    }

    //! Return task object to the memory allocator.
    void deallocate_task( task& t ) {
#if TBB_USE_ASSERT
        task_prefix& p = t.prefix();
        p.state = 0xFF;
        p.extra_state = 0xFF; 
        poison_pointer(p.next);
#endif /* TBB_USE_ASSERT */
        NFS_Free((char*)&t-task_prefix_reservation_size);
#if __TBB_COUNT_TASK_NODES
        --task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */
    }

    //! True if running on a worker thread, false otherwise.
    bool is_worker() {
#if __TBB_ARENA_PER_MASTER
        return arena_index != 0;
#else /* !__TBB_ARENA_PER_MASTER */
        return arena_index < my_arena->prefix().number_of_workers;
#endif /* !__TBB_ARENA_PER_MASTER */
    }

#if __TBB_ARENA_PER_MASTER
    //! Returns number of worker threads in the arena this thread belongs to.
    unsigned number_of_workers_in_my_arena() {
        return my_arena->my_max_num_workers;
    }
#endif /* __TBB_ARENA_PER_MASTER */

#if __TBB_COUNT_TASK_NODES
    intptr_t get_task_node_count( bool count_arena_workers = false ) {
        return task_node_count + (count_arena_workers? my_arena->workers_task_node_count(): 0);
    }
#endif /* __TBB_COUNT_TASK_NODES */

    //! Special value used to mark return_list as not taking any more entries.
    static task* plugged_return_list() {return (task*)(intptr_t)(-1);}

    //! Number of small tasks that have been allocated by this scheduler. 
    intptr_t small_task_count;

    //! List of small tasks that have been returned to this scheduler by other schedulers.
    task* return_list;

    //! Try getting a task from the mailbox or stealing from another scheduler.
    /** Redirects to a customization. */
    virtual task* receive_or_steal_task( reference_count&, bool ) = 0; 

    //! Free a small task t that that was allocated by a different scheduler 
    void free_nonlocal_small_task( task& t ); 

#if __TBB_TASK_GROUP_CONTEXT
    //! Padding isolating thread local members from members that can be written to by other threads.
    char _padding1[NFS_MaxLineSize - sizeof(context_list_node_t)];

    //! Head of the thread specific list of task group contexts.
    context_list_node_t context_list_head;

    //! Mutex protecting access to the list of task group contexts.
    spin_mutex context_list_mutex;

#if !__TBB_ARENA_PER_MASTER
    //! Used to form the list of master thread schedulers.
    scheduler_list_node_t my_node;
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Thread local cancellation epoch.
    /** When local epoch equals the global one, the cancellation state known
        to this thread is synchronized with the global cancellation state. **/
    uintptr_t local_cancel_count;

    //! Flag indicating that a context is being destructed by its owner thread 
    /** Together with nonlocal_ctx_list_update constitue a synchronization protocol
        that keeps hot path of context destruction (by the owner thread) mostly 
        lock-free. **/
    uintptr_t local_ctx_list_update;

    //! Detaches abandoned contexts
    /** These contexts must be destroyed by other threads. **/
    void cleanup_local_context_list ();

#if !__TBB_ARENA_PER_MASTER
    //! Propagates cancellation request to all descendants of the context.
    void propagate_cancellation ( task_group_context& ctx );
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Propagates cancellation request to contexts registered by this scheduler.
    void propagate_cancellation ();
#endif /* __TBB_TASK_GROUP_CONTEXT */

#if _WIN32||_WIN64
private:
    //! Handle returned by RML when registering a master with RML
    ::rml::server::execution_resource_t master_exec_resource;

#if !__TBB_ARENA_PER_MASTER
    //! register master with the resource manager
    void register_master() {
        __TBB_ASSERT( my_arena->prefix().server, "RML server not defined?" );
        // the server may ignore registration and set master_exec_resource to NULL.
        my_arena->prefix().server->register_master( master_exec_resource );
    }

    //! unregister master with the resource manager
    void unregister_master() const {
        my_arena->prefix().server->unregister_master( master_exec_resource );
    }
#endif /* !__TBB_ARENA_PER_MASTER && ( _WIN32||_WIN64 ) */
#endif /* _WIN32||_WIN64 */

    //! Dummy slot used when scheduler is not in arena
    /** The data structure is heavily padded, therefore it should be placed after 
        other data fields used by the owner thread only to allow compiler using 
        instructions with short offsets when accessing the majority of data members. **/
    arena_slot dummy_slot;

#if __TBB_TASK_GROUP_CONTEXT
    //! Flag indicating that a context is being destructed by non-owner thread.
    /** See also local_ctx_list_update. **/
    uintptr_t nonlocal_ctx_list_update;
#endif /* __TBB_TASK_GROUP_CONTEXT */

#if __TBB_SURVIVE_THREAD_SWITCH
    __cilk_tbb_unwatch_thunk my_cilk_unwatch_thunk;
#if TBB_USE_ASSERT
    //! State values used to check interface contract with Cilk runtime.
    /** Names of cs_running...cs_freed derived from state machine diagram in cilk-tbb-interop.h */
    enum cilk_state_t {
        cs_none=0xF000, // Start at nonzero value so that we can detect use of zeroed memory.
        cs_running,
        cs_limbo,
        cs_freed
    };
    cilk_state_t my_cilk_state;
#endif /* TBB_USE_ASSERT */
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

#if __TBB_STATISTICS
    //! Set of counters to track internal statistics on per thread basis
    /** Placed at the end of the class definition to minimize the disturbance of
        the core logic memory operations. **/
    mutable statistics_counters my_counters;
#endif /* __TBB_STATISTICS */

}; // class generic_scheduler


template<free_task_hint h>
void generic_scheduler::free_task( task& t ) {
    GATHER_STATISTIC(--my_counters.active_tasks);
    task_prefix& p = t.prefix();
    // Verify that optimization hints are correct.
    __TBB_ASSERT( h!=small_local_task || p.origin==this, NULL );
    __TBB_ASSERT( !(h&small_task) || p.origin, NULL );
#if TBB_USE_ASSERT
    p.depth = 0xDEADBEEF;
    p.ref_count = 0xDEADBEEF;
    poison_pointer(p.owner);
#endif /* TBB_USE_ASSERT */
    __TBB_ASSERT( 1L<<t.state() & (1L<<task::executing|1L<<task::allocated), NULL );
    p.state = task::freed;
    if( h==small_local_task || p.origin==this ) {
        GATHER_STATISTIC(++my_counters.free_list_length);
        p.next = free_list;
        free_list = &t;
    } else if( !(h&local_task) && p.origin ) {
        free_nonlocal_small_task(t);
    } else {
        GATHER_STATISTIC(--my_counters.big_tasks);
        deallocate_task(t);
    }
}

} // namespace internal
} // namespace tbb

#include "governor.h"

inline void tbb::internal::generic_scheduler::spawn( task& first, task*& next ) {
    governor::local_scheduler()->local_spawn( first, next );
}

inline void tbb::internal::generic_scheduler::spawn_root_and_wait( task& first, task*& next ) {
    governor::local_scheduler()->local_spawn_root_and_wait( first, next );
}

#if __TBB_ARENA_PER_MASTER
inline void tbb::internal::generic_scheduler::enqueue( task& task_, void* /*reserved*/ ) {
    governor::local_scheduler()->local_enqueue( task_ );
}

#endif /* __TBB_ARENA_PER_MASTER */
#endif /* _TBB_scheduler_H */
