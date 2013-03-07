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

#include "tbb/tbb_machine.h"

#include "custom_scheduler.h"
#include "scheduler_utility.h"
#include "governor.h"
#include "market.h"
#include "arena.h"
#include "mailbox.h"
#include "observer_proxy.h"
#include "itt_notify.h"

namespace tbb {
namespace internal {

/** Defined in tbb_main.cpp **/
extern generic_scheduler* (*AllocateSchedulerPtr)( arena*, size_t index );

inline generic_scheduler* allocate_scheduler ( arena* a, size_t index ) { 
    return AllocateSchedulerPtr(a, index);
}

#if __TBB_TASK_GROUP_CONTEXT
#if !__TBB_ARENA_PER_MASTER
//! Head of the list of master thread schedulers.
static scheduler_list_node_t the_scheduler_list_head;

//! Mutex protecting access to the list of schedulers.
static mutex the_scheduler_list_mutex;
#endif /* !__TBB_ARENA_PER_MASTER */

//! Counter that is incremented whenever new cancellation signal is sent to a task group.
/** Together with generic_scheduler::local_cancel_count forms cross-thread signaling
    mechanism that allows to avoid locking at the hot path of normal execution flow.

    When a descendant task group context is being registered or unregistered,
    the global and local counters are compared. If they differ, it means that 
    a cancellation signal is being propagated, and registration/deregistration
    routines take slower branch that may block (at most one thread of the pool
    can be blocked at any moment). Otherwise the control path is lock-free and fast. **/
uintptr_t global_cancel_count = 0;

//! Context to be associated with dummy tasks of worker threads schedulers.
/** It is never used for its direct purpose, and is introduced solely for the sake 
    of avoiding one extra conditional branch in the end of wait_for_all method. **/
static task_group_context dummy_context(task_group_context::isolated);
#endif /* __TBB_TASK_GROUP_CONTEXT */

void Scheduler_OneTimeInitialization ( bool itt_present ) {
    AllocateSchedulerPtr = itt_present ? &custom_scheduler<DefaultSchedulerTraits>::allocate_scheduler :
                                      &custom_scheduler<IntelSchedulerTraits>::allocate_scheduler;
#if __TBB_TASK_GROUP_CONTEXT && !__TBB_ARENA_PER_MASTER
    ITT_SYNC_CREATE(&the_scheduler_list_mutex, SyncType_GlobalLock, SyncObj_SchedulersList);
    the_scheduler_list_head.my_next = &the_scheduler_list_head;
    the_scheduler_list_head.my_prev = &the_scheduler_list_head;
#endif /* __TBB_TASK_GROUP_CONTEXT && !__TBB_ARENA_PER_MASTER */
}

//------------------------------------------------------------------------
// scheduler interface
//------------------------------------------------------------------------

//  A pure virtual destructor should still have a body
//  so the one for tbb::internal::scheduler::~scheduler() is provided here
scheduler::~scheduler( ) {}

//------------------------------------------------------------------------
// generic_scheduler
//------------------------------------------------------------------------

#if _MSC_VER && !defined(__INTEL_COMPILER)
    // Suppress overzealous compiler warning about using 'this' in base initializer list. 
    #pragma warning(push)
    #pragma warning(disable:4355)
#endif

generic_scheduler::generic_scheduler( arena* a, size_t index ) :
    my_stealing_threshold(0),
    arena_index(index),
    task_pool_size(0),
    my_arena_slot(&dummy_slot),
#if __TBB_ARENA_PER_MASTER
    my_market(NULL),
#endif /* __TBB_ARENA_PER_MASTER */
    my_arena(a),
    random( unsigned(this-(generic_scheduler*)NULL) ),
    free_list(NULL),
    innermost_running_task(NULL),
    dummy_task(NULL),
    ref_count(1),
    my_affinity_id(0),
    is_registered(false),
    is_auto_initialized(false),
#if __TBB_SCHEDULER_OBSERVER
    local_last_observer_proxy(NULL),
#endif /* __TBB_SCHEDULER_OBSERVER */
#if __TBB_COUNT_TASK_NODES
    task_node_count(0),
#endif /* __TBB_COUNT_TASK_NODES */
    small_task_count(1),   // Extra 1 is a guard reference
    return_list(NULL),
#if __TBB_TASK_GROUP_CONTEXT
    local_ctx_list_update(0),
    nonlocal_ctx_list_update(0)
#endif /* __TBB_TASK_GROUP_CONTEXT */
#if __TBB_SURVIVE_THREAD_SWITCH && TBB_USE_ASSERT
   ,my_cilk_state(cs_none)
#endif /* __TBB_SURVIVE_THREAD_SWITCH && TBB_USE_ASSERT */
{
    dummy_slot.task_pool = allocate_task_pool( min_task_pool_size );
    dummy_slot.head = dummy_slot.tail = 0;
    dummy_task = &allocate_task( sizeof(task), __TBB_CONTEXT_ARG(NULL, NULL) );
#if __TBB_TASK_GROUP_CONTEXT
    context_list_head.my_prev = &context_list_head;
    context_list_head.my_next = &context_list_head;
    ITT_SYNC_CREATE(&context_list_mutex, SyncType_Scheduler, SyncObj_ContextsList);
#endif /* __TBB_TASK_GROUP_CONTEXT */
    dummy_task->prefix().ref_count = 2;
    ITT_SYNC_CREATE(&dummy_task->prefix().ref_count, SyncType_Scheduler, SyncObj_WorkerLifeCycleMgmt);
    ITT_SYNC_CREATE(&return_list, SyncType_Scheduler, SyncObj_TaskReturnList);
    assert_task_pool_valid();
#if __TBB_SURVIVE_THREAD_SWITCH
    my_cilk_unwatch_thunk.routine = NULL;
#endif /* __TBB_SURVIVE_THREAD_SWITCH */
}

#if _MSC_VER && !defined(__INTEL_COMPILER)
    #pragma warning(pop)
#endif // warning 4355 is back

#if TBB_USE_ASSERT > 1
bool generic_scheduler::assert_task_pool_valid() const {
    acquire_task_pool();
    task** tp = dummy_slot.task_pool;
    __TBB_ASSERT( task_pool_size >= min_task_pool_size, NULL );
    __TBB_ASSERT( my_arena_slot->head <= my_arena_slot->tail, NULL );
    for ( size_t i = 0; i < my_arena_slot->head; ++i )
        __TBB_ASSERT( tp[i] == poisoned_ptr, "Task pool corrupted" );
    for ( size_t i = my_arena_slot->head; i < my_arena_slot->tail; ++i ) {
        __TBB_ASSERT( (uintptr_t)tp[i] + 1 > 1u, "nil or invalid task pointer in the deque" );
        __TBB_ASSERT( tp[i]->prefix().state == task::ready ||
                      tp[i]->prefix().extra_state == es_task_proxy, "task in the deque has invalid state" );
    }
    for ( size_t i = my_arena_slot->tail; i < task_pool_size; ++i )
        __TBB_ASSERT( tp[i] == poisoned_ptr, "Task pool corrupted" );
    release_task_pool();
}
#endif /* TBB_USE_ASSERT > 1 */

#if __TBB_TASK_GROUP_CONTEXT
void generic_scheduler::propagate_cancellation () {
    spin_mutex::scoped_lock lock(context_list_mutex);
    // Acquire fence is necessary to ensure that the subsequent node->my_next load 
    // returned the correct value in case it was just inserted in another thread.
    // The fence also ensures visibility of the correct my_parent value.
    context_list_node_t *node = __TBB_load_with_acquire(context_list_head.my_next);
    while ( node != &context_list_head ) {
        task_group_context &ctx = __TBB_get_object_ref(task_group_context, my_node, node);
            // The absence of acquire fence while reading my_cancellation_requested may result 
            // in repeated traversals of the same parents chain if another group (precedent or 
            // descendant) belonging to the tree being canceled sends cancellation request of 
            // its own around the same time.
        if ( !ctx.my_cancellation_requested )
            ctx.propagate_cancellation_from_ancestors();
        node = node->my_next;
        __TBB_ASSERT( is_alive(ctx.my_version_and_traits), "Walked into a destroyed context while propagating cancellation" );
    }
    // Sync up local cancelation epoch with the global one. Release fence prevents 
    // reordering of possible store to my_cancellation_requested after the sync point.
    __TBB_store_with_release(local_cancel_count, global_cancel_count);
}

#if !__TBB_ARENA_PER_MASTER
/** Propagates cancellation down the tree of dependent contexts by walking each 
    thread's local list of contexts **/
void generic_scheduler::propagate_cancellation ( task_group_context& ctx ) {
    __TBB_ASSERT ( ctx.my_cancellation_requested, "No cancellation request in the context" );
    // The whole propagation algorithm is under the lock in order to ensure correctness 
    // in case of parallel cancellations at the different levels of the context tree.
    // See the note 2 at the bottom of the file.
    mutex::scoped_lock lock(the_scheduler_list_mutex);
    // Advance global cancellation state
    __TBB_FetchAndAddWrelease(&global_cancel_count, 1);
    // First propagate to workers using arena to access their context lists
    size_t num_workers = my_arena->prefix().number_of_workers;
    for ( size_t i = 0; i < num_workers; ++i ) {
        // No fence is necessary here since the context list of worker's scheduler 
        // can contain anything of interest only after the first stealing was done
        // by that worker. And doing it applies the necessary fence
        generic_scheduler *s = my_arena->prefix().worker_list[i].scheduler;
        // If the worker is in the middle of its startup sequence, skip it.
        if ( s )
            s->propagate_cancellation();
    }
    // Then propagate to masters using the global list of master's schedulers
    scheduler_list_node_t *node = the_scheduler_list_head.my_next;
    while ( node != &the_scheduler_list_head ) {
        __TBB_get_object_ref(generic_scheduler, my_node, node).propagate_cancellation();
        node = node->my_next;
    }
}
#endif /* !__TBB_ARENA_PER_MASTER */
#endif /* __TBB_TASK_GROUP_CONTEXT */


void generic_scheduler::init_stack_info () {
    // Stacks are growing top-down. Highest address is called "stack base", 
    // and the lowest is "stack limit".
#if __TBB_ARENA_PER_MASTER
    __TBB_ASSERT( !my_stealing_threshold, "Stealing threshold has already been calculated" );
    size_t  stack_size = my_market->worker_stack_size();
#else /* !__TBB_ARENA_PER_MASTER */
    size_t  stack_size = my_arena->prefix().stack_size;
#endif /* !__TBB_ARENA_PER_MASTER */
#if USE_WINTHREAD
#if defined(_MSC_VER)&&_MSC_VER<1400 && !_WIN64
    NT_TIB  *pteb = (NT_TIB*)__TBB_machine_get_current_teb();
#else
    NT_TIB  *pteb = (NT_TIB*)NtCurrentTeb();
#endif
    __TBB_ASSERT( &pteb < pteb->StackBase && &pteb > pteb->StackLimit, "invalid stack info in TEB" );
    __TBB_ASSERT( stack_size >0, "stack_size not initialized?" );
    // When a thread is created with the attribute STACK_SIZE_PARAM_IS_A_RESERVATION, stack limit 
    // in the TIB points to the committed part of the stack only. This renders the expression
    // "(uintptr_t)pteb->StackBase / 2 + (uintptr_t)pteb->StackLimit / 2" virtually useless.
    // Thus for worker threads we use the explicit stack size we used while creating them.
    // And for master threads we rely on the following fact and assumption:
    // - the default stack size of a master thread on Windows is 1M;
    // - if it was explicitly set by the application it is at least as large as the size of a worker stack.
    if ( is_worker() || stack_size < MByte )
        my_stealing_threshold = (uintptr_t)pteb->StackBase - stack_size / 2;
    else
        my_stealing_threshold = (uintptr_t)pteb->StackBase - MByte / 2;
#else /* USE_PTHREAD */
    // There is no portable way to get stack base address in Posix, so we use 
    // non-portable method (on all modern Linux) or the simplified approach 
    // based on the common sense assumptions. The most important assumption 
    // is that the main thread's stack size is not less than that of other threads.
    void    *stack_base = &stack_size;
#if __TBB_ipf
    void    *rsb_base = __TBB_get_bsp();
#endif
#if __linux__
    size_t  np_stack_size = 0;
    void    *stack_limit = NULL;
    pthread_attr_t  attr_stack, np_attr_stack;
    if( 0 == pthread_getattr_np(pthread_self(), &np_attr_stack) ) {
        if ( 0 == pthread_attr_getstack(&np_attr_stack, &stack_limit, &np_stack_size) ) {
            if ( 0 == pthread_attr_init(&attr_stack) ) {
                if ( 0 == pthread_attr_getstacksize(&attr_stack, &stack_size) )
                {
                    stack_base = (char*)stack_limit + np_stack_size;
                    if ( np_stack_size < stack_size ) {
                        // We are in a secondary thread. Use reliable data.
#if __TBB_ipf
                        // IA64 stack is split into RSE backup and memory parts
                        rsb_base = stack_limit;
                        stack_size = np_stack_size/2;
#else
                        stack_size = np_stack_size;
#endif /* !__TBB_ipf */
                    }
                    // We are either in the main thread or this thread stack 
                    // is bigger that that of the main one. As we cannot discern
                    // these cases we fall back to the default (heuristic) values.
                }
                pthread_attr_destroy(&attr_stack);
            }
        }
        pthread_attr_destroy(&np_attr_stack);
    }
#endif /* __linux__ */
    __TBB_ASSERT( stack_size>0, "stack size must be positive" );
    my_stealing_threshold = (uintptr_t)((char*)stack_base - stack_size/2);
#if __TBB_ipf
    my_rsb_stealing_threshold = (uintptr_t)((char*)rsb_base + stack_size/2);
#endif
#endif /* USE_PTHREAD */
}

/** The function uses synchronization scheme similar to the one in the destructor
    of task_group_context augmented with interlocked state change of each context
    object. The purpose of this algo is to prevent threads doing nonlocal context
    destruction from accessing destroyed owner-scheduler instance still pointed to 
    by the context object. **/
void generic_scheduler::cleanup_local_context_list () {
    // Detach contexts remaining in the local list
    bool wait_for_concurrent_destroyers_to_leave = false;
    uintptr_t local_count_snapshot = local_cancel_count;
    local_ctx_list_update = 1;
    {
        // This is just a definition. Actual lock is acquired only in case of conflict.
        spin_mutex::scoped_lock lock;
        // Full fence prevents reordering of store to local_ctx_list_update with 
        // load from nonlocal_ctx_list_update.
        __TBB_full_memory_fence();
        // Check for the conflict with concurrent destroyer or cancelation propagator
        if ( nonlocal_ctx_list_update || local_count_snapshot != global_cancel_count )
            lock.acquire(context_list_mutex);
        // No acquire fence is necessary for loading context_list_head.my_next,
        // as the list can be updated by this thread only.
        context_list_node_t *node = context_list_head.my_next;
        while ( node != &context_list_head ) {
            task_group_context &ctx = __TBB_get_object_ref(task_group_context, my_node, node);
            __TBB_ASSERT( ctx.my_kind != task_group_context::binding_required, "Only a context bound to a root task can be detached" );
            node = node->my_next;
            __TBB_ASSERT( is_alive(ctx.my_version_and_traits), "Walked into a destroyed context while detaching contexts from the local list" );
            // On 64-bit systems my_kind can be a 32-bit value padded with 32 uninitialized bits.
            // So the cast below is necessary to throw off the higher bytes containing garbage
            if ( (task_group_context::kind_type)(uintptr_t)__TBB_FetchAndStoreW(&ctx.my_kind, task_group_context::detached) == task_group_context::dying )
                wait_for_concurrent_destroyers_to_leave = true;
        }
    }
    __TBB_store_with_release( local_ctx_list_update, 0 );
    // Wait until other threads referencing this scheduler object finish with it
    if ( wait_for_concurrent_destroyers_to_leave )
        spin_wait_until_eq( nonlocal_ctx_list_update, 0u );
}

void generic_scheduler::free_scheduler() {
    if( in_arena() ) {
        acquire_task_pool();
        leave_arena();
    }
#if __TBB_TASK_GROUP_CONTEXT
    cleanup_local_context_list();
#if !__TBB_ARENA_PER_MASTER
    task_group_context* default_context = dummy_task->prefix().context;
    if ( default_context != &dummy_context) {
        // Only master thread's dummy task has a dynamically allocated context
        default_context->task_group_context::~task_group_context();
        NFS_Free(default_context);
        {
            mutex::scoped_lock lock(the_scheduler_list_mutex);
            my_node.my_next->my_prev = my_node.my_prev;
            my_node.my_prev->my_next = my_node.my_next;
        }
    }
#endif /* !__TBB_ARENA_PER_MASTER */
#endif /* __TBB_TASK_GROUP_CONTEXT */
    free_task<small_local_task>( *dummy_task );

    // k accounts for a guard reference and each task that we deallocate.
    intptr_t k = 1;
    for(;;) {
        while( task* t = free_list ) {
            free_list = t->prefix().next;
            deallocate_task(*t);
            ++k;
        }
        if( return_list==plugged_return_list() ) 
            break;
        free_list = (task*)__TBB_FetchAndStoreW( &return_list, (intptr_t)plugged_return_list() );
    }
#if __TBB_COUNT_TASK_NODES
#if __TBB_ARENA_PER_MASTER
    my_market->update_task_node_count( task_node_count );
#else /* !__TBB_ARENA_PER_MASTER */
    my_arena->prefix().task_node_count += task_node_count;
#endif /* !__TBB_ARENA_PER_MASTER */
#endif /* __TBB_COUNT_TASK_NODES */
#if !__TBB_ARENA_PER_MASTER && __TBB_STATISTICS
    dump_statistics(my_counters, arena_index < my_arena->prefix().number_of_workers ? arena_index + 1 : 0 );
#endif /* !__TBB_ARENA_PER_MASTER && __TBB_STATISTICS */
    free_task_pool( dummy_slot.task_pool );
    dummy_slot.task_pool = NULL;
    // Update small_task_count last.  Doing so sooner might cause another thread to free *this.
    __TBB_ASSERT( small_task_count>=k, "small_task_count corrupted" );
    governor::sign_off(this);
    if( __TBB_FetchAndAddW( &small_task_count, -k )==k )
        NFS_Free( this );
}

task& generic_scheduler::allocate_task( size_t number_of_bytes, 
                                            __TBB_CONTEXT_ARG(task* parent, task_group_context* context) ) {
    GATHER_STATISTIC(++my_counters.active_tasks);
    task* t = free_list;
    if( number_of_bytes<=quick_task_size ) {
        if( t ) {
            GATHER_STATISTIC(--my_counters.free_list_length);
            __TBB_ASSERT( t->state()==task::freed, "free list of tasks is corrupted" );
            free_list = t->prefix().next;
        } else if( return_list ) {
            // No fence required for read of return_list above, because __TBB_FetchAndStoreW has a fence.
            t = (task*)__TBB_FetchAndStoreW( &return_list, 0 );
            __TBB_ASSERT( t, "another thread emptied the return_list" );
            __TBB_ASSERT( t->prefix().origin==this, "task returned to wrong return_list" );
            ITT_NOTIFY( sync_acquired, &return_list );
            free_list = t->prefix().next;
        } else {
            t = (task*)((char*)NFS_Allocate( task_prefix_reservation_size+quick_task_size, 1, NULL ) + task_prefix_reservation_size );
#if __TBB_COUNT_TASK_NODES
            ++task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */
            t->prefix().origin = this;
            ++small_task_count;
        }
    } else {
        GATHER_STATISTIC(++my_counters.big_tasks);
        t = (task*)((char*)NFS_Allocate( task_prefix_reservation_size+number_of_bytes, 1, NULL ) + task_prefix_reservation_size );
#if __TBB_COUNT_TASK_NODES
        ++task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */
        t->prefix().origin = NULL;
    }
    task_prefix& p = t->prefix();
#if __TBB_TASK_GROUP_CONTEXT
    p.context = context;
#endif /* __TBB_TASK_GROUP_CONTEXT */
    p.owner = this;
    p.ref_count = 0;
    // Assign some not outrageously out-of-place value for a while
    p.depth = 0;
    p.parent = parent;
    // In TBB 2.1 and later, the constructor for task sets extra_state to indicate the version of the tbb/task.h header.
    // In TBB 2.0 and earlier, the constructor leaves extra_state as zero.
    p.extra_state = 0;
    p.affinity = 0;
    p.state = task::allocated;
    return *t;
}

void generic_scheduler::free_nonlocal_small_task( task& t ) {
    __TBB_ASSERT( t.state()==task::freed, NULL );
    generic_scheduler& s = *static_cast<generic_scheduler*>(t.prefix().origin);
    __TBB_ASSERT( &s!=this, NULL );
    for(;;) {
        task* old = s.return_list;
        if( old==plugged_return_list() ) 
            break;
        // Atomically insert t at head of s.return_list
        t.prefix().next = old; 
        ITT_NOTIFY( sync_releasing, &s.return_list );
        if( __TBB_CompareAndSwapW( &s.return_list, (intptr_t)&t, (intptr_t)old )==(intptr_t)old ) {
            GATHER_STATISTIC(++my_counters.free_list_length);
            return;
        }
    }
    deallocate_task(t);
    if( __TBB_FetchAndDecrementWrelease( &s.small_task_count )==1 ) {
        // We freed the last task allocated by scheduler s, so it's our responsibility
        // to free the scheduler.
        NFS_Free( &s );
    }
}

task** generic_scheduler::allocate_task_pool( size_t n ) {
    __TBB_ASSERT( n > task_pool_size, "Cannot shrink the task pool" );
    size_t byte_size = ((n * sizeof(task*) + NFS_MaxLineSize - 1) / NFS_MaxLineSize) * NFS_MaxLineSize;
    task_pool_size = byte_size / sizeof(task*);
    task** new_pool = (task**)NFS_Allocate( byte_size, 1, NULL );
    // No need to clear the fresh deque since valid items are designated by the head and tail members.
#if TBB_USE_ASSERT>=2
    // But clear it in the high vigilance debug mode
    memset( new_pool, reinterpret_cast<int>(poisoned_ptr), n );
#endif /* TBB_USE_ASSERT>=2 */
    return new_pool;
}

void generic_scheduler::grow_task_pool( size_t new_size ) {
    assert_task_pool_valid();
    if ( new_size < 2 * task_pool_size )
        new_size = 2 * task_pool_size;
    task** new_pool = allocate_task_pool( new_size ); // updates task_pool_size
    task** old_pool = dummy_slot.task_pool;
    acquire_task_pool();    // requires the old dummy_slot.task_pool value
    my_arena_slot->tail -= my_arena_slot->head;
    __TBB_ASSERT( my_arena_slot->tail <= task_pool_size, "new task pool is too short" );
    memcpy( new_pool, old_pool + my_arena_slot->head, my_arena_slot->tail * sizeof(task*) );
    my_arena_slot->head = 0;
    dummy_slot.task_pool = new_pool;
    release_task_pool();    // updates the task pool pointer in our arena slot
    free_task_pool( old_pool );
    assert_task_pool_valid();
}

/** ATTENTION: 
    This method is mostly the same as generic_scheduler::lock_task_pool(), with 
    a little different logic of slot state checks (slot is either locked or points 
    to our task pool).
    Thus if either of them is changed, consider changing the counterpart as well. **/
inline void generic_scheduler::acquire_task_pool() const {
    if ( !in_arena() )
        return; // we are not in arena - nothing to lock
    atomic_backoff backoff;
    bool sync_prepare_done = false;
    for(;;) {
#if TBB_USE_ASSERT
        __TBB_ASSERT( my_arena_slot == my_arena->slot + arena_index, "invalid arena slot index" );
        // Local copy of the arena slot task pool pointer is necessary for the next 
        // assertion to work correctly to exclude asynchronous state transition effect.
        task** tp = my_arena_slot->task_pool;
        __TBB_ASSERT( tp == LockedTaskPool || tp == dummy_slot.task_pool, "slot ownership corrupt?" );
#endif
        if( my_arena_slot->task_pool != LockedTaskPool && 
            __TBB_CompareAndSwapW( &my_arena_slot->task_pool, (intptr_t)LockedTaskPool, 
                                   (intptr_t)dummy_slot.task_pool ) == (intptr_t)dummy_slot.task_pool )
        {
            // We acquired our own slot
            ITT_NOTIFY(sync_acquired, my_arena_slot);
            break;
        } 
        else if( !sync_prepare_done ) {
            // Start waiting
            ITT_NOTIFY(sync_prepare, my_arena_slot);
            sync_prepare_done = true;
        }
        // Someone else acquired a lock, so pause and do exponential backoff.
        backoff.pause();
    }
    __TBB_ASSERT( my_arena_slot->task_pool == LockedTaskPool, "not really acquired task pool" );
} // generic_scheduler::acquire_task_pool

inline void generic_scheduler::release_task_pool() const {
    if ( !in_arena() )
        return; // we are not in arena - nothing to unlock
    __TBB_ASSERT( my_arena_slot, "we are not in arena" );
    __TBB_ASSERT( my_arena_slot->task_pool == LockedTaskPool, "arena slot is not locked" );
    ITT_NOTIFY(sync_releasing, my_arena_slot);
    __TBB_store_with_release( my_arena_slot->task_pool, dummy_slot.task_pool );
}

/** ATTENTION: 
    This method is mostly the same as generic_scheduler::acquire_task_pool(), 
    with a little different logic of slot state checks (slot can be empty, locked 
    or point to any task pool other than ours, and asynchronous transitions between 
    all these states are possible).
    Thus if any of them is changed, consider changing the counterpart as well **/
inline task** generic_scheduler::lock_task_pool( arena_slot* victim_arena_slot ) const {
    task** victim_task_pool;
    atomic_backoff backoff;
    bool sync_prepare_done = false;
    for(;;) {
        victim_task_pool = victim_arena_slot->task_pool;
        // NOTE: Do not use comparison of head and tail indices to check for
        // the presence of work in the victim's task pool, as they may give
        // incorrect indication because of task pool relocations and resizes.
        if ( victim_task_pool == EmptyTaskPool ) {
            // The victim thread emptied its task pool - nothing to lock
            if( sync_prepare_done )
                ITT_NOTIFY(sync_cancel, victim_arena_slot);
            break;
        }
        if( victim_task_pool != LockedTaskPool && 
            __TBB_CompareAndSwapW( &victim_arena_slot->task_pool, 
                (intptr_t)LockedTaskPool, (intptr_t)victim_task_pool ) == (intptr_t)victim_task_pool )
        {
            // We've locked victim's task pool
            ITT_NOTIFY(sync_acquired, victim_arena_slot);
            break;
        }
        else if( !sync_prepare_done ) {
            // Start waiting
            ITT_NOTIFY(sync_prepare, victim_arena_slot);
            sync_prepare_done = true;
        }
        GATHER_STATISTIC( ++my_counters.thieves_conflicts );
        // Someone else acquired a lock, so pause and do exponential backoff.
        backoff.pause();
    }
    __TBB_ASSERT( victim_task_pool == EmptyTaskPool || 
                  (victim_arena_slot->task_pool == LockedTaskPool && victim_task_pool != LockedTaskPool), 
                  "not really locked victim's task pool?" );
    return victim_task_pool;
} // generic_scheduler::lock_task_pool

inline void generic_scheduler::unlock_task_pool( arena_slot* victim_arena_slot, 
                                                task** victim_task_pool ) const {
    __TBB_ASSERT( victim_arena_slot, "empty victim arena slot pointer" );
    __TBB_ASSERT( victim_arena_slot->task_pool == LockedTaskPool, "victim arena slot is not locked" );
    ITT_NOTIFY(sync_releasing, victim_arena_slot);
    __TBB_store_with_release( victim_arena_slot->task_pool, victim_task_pool );
}


inline task* generic_scheduler::prepare_for_spawning( task* t ) {
    __TBB_ASSERT( t->state()==task::allocated, "attempt to spawn task that is not in 'allocated' state" );
    t->prefix().owner = this;
    t->prefix().state = task::ready;
#if TBB_USE_ASSERT
    if( task* parent = t->parent() ) {
        internal::reference_count ref_count = parent->prefix().ref_count;
        __TBB_ASSERT( ref_count>=0, "attempt to spawn task whose parent has a ref_count<0" );
        __TBB_ASSERT( ref_count!=0, "attempt to spawn task whose parent has a ref_count==0 (forgot to set_ref_count?)" );
        parent->prefix().extra_state |= es_ref_count_active;
    }
#endif /* TBB_USE_ASSERT */
    affinity_id dst_thread = t->prefix().affinity;
    __TBB_ASSERT( dst_thread == 0 || is_version_3_task(*t), "backwards compatibility to TBB 2.0 tasks is broken" );
    if( dst_thread != 0 && dst_thread != my_affinity_id ) {
        task_proxy& proxy = (task_proxy&)allocate_task( sizeof(task_proxy), 
                                                      __TBB_CONTEXT_ARG(NULL, NULL) );
        // Mark as a proxy
        proxy.prefix().extra_state = es_task_proxy;
        proxy.outbox = &my_arena->mailbox(dst_thread);
        proxy.task_and_tag = intptr_t(t)|3;
        ITT_NOTIFY( sync_releasing, proxy.outbox );
        // Mail the proxy - after this point t may be destroyed by another thread at any moment.
        proxy.outbox->push(proxy);
        return &proxy;
    }
    return t;
}

/** Conceptually, this method should be a member of class scheduler.
    But doing so would force us to publish class scheduler in the headers. */
void generic_scheduler::local_spawn( task& first, task*& next ) {
    __TBB_ASSERT( governor::is_set(this), NULL );
    assert_task_pool_valid();
    if ( &first.prefix().next == &next ) {
        // Single task is being spawned
        if ( my_arena_slot->tail == task_pool_size ) {
            // If the free space at the beginning of the task pool is too short
            // we are likely facing a pathological single-producer-multiple-consumers
            // scenario, and thus it's better to expand the task pool
            if ( my_arena_slot->head > min_task_pool_size/4 ) {
                // Move the busy part of the deque to the beginning of the allocated space
                acquire_task_pool();
                my_arena_slot->tail -= my_arena_slot->head;
                memmove( dummy_slot.task_pool, dummy_slot.task_pool + my_arena_slot->head, my_arena_slot->tail * sizeof(task*) );
                my_arena_slot->head = 0;
                release_task_pool();
            }
            else {
                grow_task_pool( task_pool_size + 1 );
            }
        }
        dummy_slot.task_pool[my_arena_slot->tail] = prepare_for_spawning( &first );
        ITT_NOTIFY(sync_releasing, my_arena_slot);
        // The following store with release is required on ia64 only
        size_t new_tail = my_arena_slot->tail + 1;
        __TBB_store_with_release( my_arena_slot->tail, new_tail );
        __TBB_ASSERT ( my_arena_slot->tail <= task_pool_size, "task deque end was overwritten" );
    }
    else {
        // Task list is being spawned
        const size_t initial_capacity = 64;
        task *arr[initial_capacity];
        fast_reverse_vector<task*> tasks(arr, initial_capacity);
        task *t_next = NULL;
        for( task* t = &first; ; t = t_next ) {
            // After prepare_for_spawning returns t may already have been destroyed. 
            // So milk it while it is alive.
            bool end = &t->prefix().next == &next;
            t_next = t->prefix().next;
            tasks.push_back( prepare_for_spawning(t) );
            if( end )
                break;
        }
        size_t num_tasks = tasks.size();
        __TBB_ASSERT ( arena_index != null_arena_index, "invalid arena slot index" );
        if ( my_arena_slot->tail + num_tasks > task_pool_size ) {
            // 1 compensates for head possibly temporarily incremented by a thief
            size_t new_size = my_arena_slot->tail - my_arena_slot->head + num_tasks + 1;
            if ( new_size <= task_pool_size ) {
                // Move the busy part of the deque to the beginning of the allocated space
                acquire_task_pool();
                my_arena_slot->tail -= my_arena_slot->head;
                memmove( dummy_slot.task_pool, dummy_slot.task_pool + my_arena_slot->head, my_arena_slot->tail * sizeof(task*) );
                my_arena_slot->head = 0;
                release_task_pool();
            }
            else {
                grow_task_pool( new_size );
            }
        }
#if DO_ITT_NOTIFY
        else {
            // The preceding if-branch issues the same ittnotify inside release_task_pool() or grow_task_pool() methods
            ITT_NOTIFY(sync_releasing, my_arena_slot);
        }
#endif /* DO_ITT_NOTIFY */
        tasks.copy_memory( dummy_slot.task_pool + my_arena_slot->tail );
        // The following store with release is required on ia64 only
        size_t new_tail = my_arena_slot->tail + num_tasks;
        __TBB_store_with_release( my_arena_slot->tail, new_tail );
        __TBB_ASSERT ( my_arena_slot->tail <= task_pool_size, "task deque end was overwritten" );
    }
#if __TBB_ARENA_PER_MASTER
    if ( !in_arena() )
        enter_arena();
    my_arena->advertise_new_work</*Spawned=*/true>();
#else /* !__TBB_ARENA_PER_MASTER */
    if ( !in_arena() ) {
        if ( is_worker() )
            enter_arena();
        else
            try_enter_arena();
    }
    my_arena->mark_pool_full();
#endif /* !__TBB_ARENA_PER_MASTER */
    assert_task_pool_valid();
}

void generic_scheduler::local_spawn_root_and_wait( task& first, task*& next ) {
    __TBB_ASSERT( governor::is_set(this), NULL );
    __TBB_ASSERT( &first, NULL );
    auto_empty_task dummy( __TBB_CONTEXT_ARG(this, first.prefix().context) );
    internal::reference_count n = 0;
    for( task* t=&first; ; t=t->prefix().next ) {
        ++n;
        __TBB_ASSERT( !t->prefix().parent, "not a root task, or already running" );
        t->prefix().parent = &dummy;
        if( &t->prefix().next==&next ) break;
#if __TBB_TASK_GROUP_CONTEXT
        __TBB_ASSERT( t->prefix().context == t->prefix().next->prefix().context, 
                    "all the root tasks in list must share the same context");
#endif /* __TBB_TASK_GROUP_CONTEXT */
    }
    dummy.prefix().ref_count = n+1;
    if( n>1 )
        local_spawn( *first.prefix().next, next );
    local_wait_for_all( dummy, &first );
}

inline task* generic_scheduler::get_mailbox_task() {
    __TBB_ASSERT( my_affinity_id>0, "not in arena" );
    task* result = NULL;
    while( task_proxy* t = inbox.pop() ) {
        intptr_t tat = __TBB_load_with_acquire(t->task_and_tag);
        __TBB_ASSERT( tat==task_proxy::mailbox_bit || (tat==(tat|3)&&tat!=3), NULL );
        if( tat!=task_proxy::mailbox_bit && __TBB_CompareAndSwapW( &t->task_and_tag, task_proxy::pool_bit, tat )==tat ) {
            // Successfully grabbed the task, and left pool seeker with job of freeing the proxy.
            ITT_NOTIFY( sync_acquired, inbox.outbox() );
            result = (task*)(tat & ~3);
            result->prefix().owner = this;
            break;
        }
        free_task_proxy( *t );
    }
    return result;
}

inline task* generic_scheduler::strip_proxy( task_proxy* tp ) {
    __TBB_ASSERT( tp->prefix().extra_state==es_task_proxy, NULL );
    intptr_t tat = __TBB_load_with_acquire(tp->task_and_tag);
    if( (tat&3)==3 ) {
        // proxy is shared by a pool and a mailbox.
        // Attempt to transition it to "empty proxy in mailbox" state.
        if( __TBB_CompareAndSwapW( &tp->task_and_tag, task_proxy::mailbox_bit, tat )==tat ) {
            // Successfully grabbed the task, and left the mailbox with the job of freeing the proxy.
            return (task*)(tat&~3);
        }
        __TBB_ASSERT( tp->task_and_tag==task_proxy::pool_bit, NULL );
    } else {
        // We have exclusive access to the proxy
        __TBB_ASSERT( (tat&3)==task_proxy::pool_bit, "task did not come from pool?" );
        __TBB_ASSERT ( !(tat&~3), "Empty proxy in the pool contains non-zero task pointer" );
    }
#if TBB_USE_ASSERT
    tp->prefix().state = task::allocated;
#endif
    free_task_proxy( *tp );
    // Another thread grabbed the underlying task via their mailbox
    return NULL;
}

#if __TBB_ARENA_PER_MASTER
void generic_scheduler::local_enqueue( task& t ) {
    __TBB_ASSERT( governor::is_set(this), NULL );
    __TBB_ASSERT( t.state()==task::allocated, "attempt to enqueue task that is not in 'allocated' state" );
    t.prefix().owner = this;
    t.prefix().state = task::ready;

#if TBB_USE_ASSERT
    if( task* parent = t.parent() ) {
        internal::reference_count ref_count = parent->prefix().ref_count;
        __TBB_ASSERT( ref_count>=0, "attempt to enqueue task whose parent has a ref_count<0" );
        __TBB_ASSERT( ref_count!=0, "attempt to enqueue task whose parent has a ref_count==0 (forgot to set_ref_count?)" );
        parent->prefix().extra_state |= es_ref_count_active;
    }
    __TBB_ASSERT(t.prefix().affinity==affinity_id(0), "affinity is ignored for enqueued tasks");
#endif /* TBB_USE_ASSERT */

    __TBB_ASSERT( my_arena, "thread is not in any arena" );
    ITT_NOTIFY(sync_releasing, &my_arena->my_task_stream);
    my_arena->my_task_stream.push( &t, my_arena_slot->hint_for_push );
    my_arena->advertise_new_work< /*Spawned=*/ false >();
    assert_task_pool_valid();
}

inline task* generic_scheduler::dequeue_task() {
    task* result = NULL;
    my_arena->my_task_stream.pop(result, my_arena_slot->hint_for_pop);
    if (result) ITT_NOTIFY(sync_acquired, &my_arena->my_task_stream);
    return result;
}
#endif /* __TBB_ARENA_PER_MASTER */

inline task* generic_scheduler::get_task() {
    task* result = NULL;
retry:
    --my_arena_slot->tail;
    __TBB_full_memory_fence();
    if ( (intptr_t)my_arena_slot->head > (intptr_t)my_arena_slot->tail ) {
        acquire_task_pool();
        if ( (intptr_t)my_arena_slot->head <= (intptr_t)my_arena_slot->tail ) {
            // The thief backed off - grab the task
            result = dummy_slot.task_pool[my_arena_slot->tail];
            __TBB_ASSERT( !is_poisoned(result), NULL );
            poison_pointer( dummy_slot.task_pool[my_arena_slot->tail] );
        }
        else {
            __TBB_ASSERT ( my_arena_slot->head == my_arena_slot->tail + 1, "victim/thief arbitration algorithm failure" );
        }
        if ( (intptr_t)my_arena_slot->head < (intptr_t)my_arena_slot->tail ) {
            release_task_pool();
        }
        else {
            // In any case the deque is empty now, so compact it
            my_arena_slot->head = my_arena_slot->tail = 0;
            if ( in_arena() )
                leave_arena();
        }
    }
    else {
        result = dummy_slot.task_pool[my_arena_slot->tail];
        __TBB_ASSERT( !is_poisoned(result), NULL );
        poison_pointer( dummy_slot.task_pool[my_arena_slot->tail] );
    }
    if( result && is_proxy(*result) ) {
        result = strip_proxy((task_proxy*)result);
        if( !result ) {
            goto retry;
        }
        GATHER_STATISTIC( ++my_counters.proxies_executed );
        // Following assertion should be true because TBB 2.0 tasks never specify affinity, and hence are not proxied.
        __TBB_ASSERT( is_version_3_task(*result), "backwards compatibility with TBB 2.0 broken" );
        // Task affinity has changed.
        innermost_running_task = result;
        result->note_affinity(my_affinity_id);
    }
    return result;
} // generic_scheduler::get_task

task* generic_scheduler::steal_task( arena_slot& victim_slot ) {
    task** victim_pool = lock_task_pool( &victim_slot );
    if ( !victim_pool )
        return NULL;
    const size_t none = ~size_t(0);
    size_t first_skipped_proxy = none;
    task* result = NULL;
retry:
    ++victim_slot.head;
    __TBB_full_memory_fence();
    if ( (intptr_t)victim_slot.head > (intptr_t)victim_slot.tail ) {
        --victim_slot.head;
    }
    else {
        result = victim_pool[victim_slot.head - 1];
        __TBB_ASSERT( !is_poisoned(result), NULL );
        if( is_proxy(*result) ) {
            task_proxy& tp = *static_cast<task_proxy*>(result);
            // If task will likely be grabbed by whom it was mailed to, skip it.
            if( (tp.task_and_tag & 3) == 3 && tp.outbox->recipient_is_idle() ) {
                GATHER_STATISTIC( ++my_counters.proxies_bypassed );
                if ( first_skipped_proxy == none )
                    first_skipped_proxy = victim_slot.head - 1;
                result = NULL;
                goto retry;
            }
        }
        poison_pointer(victim_pool[victim_slot.head - 1]);
    }
    if ( first_skipped_proxy != none ) {
        if ( result ) {
            victim_pool[victim_slot.head - 1] = victim_pool[first_skipped_proxy];
            poison_pointer( victim_pool[first_skipped_proxy] );
            __TBB_store_with_release( victim_slot.head, first_skipped_proxy + 1 );
        }
        else
            __TBB_store_with_release( victim_slot.head, first_skipped_proxy );
    }
    unlock_task_pool( &victim_slot, victim_pool );
    return result;
}

inline void generic_scheduler::do_enter_arena() {
    my_arena_slot = &my_arena->slot[arena_index];
    __TBB_ASSERT ( my_arena_slot->head == my_arena_slot->tail, "task deque of a free slot must be empty" );
    __TBB_ASSERT ( dummy_slot.head < dummy_slot.tail, "entering arena without tasks to share" );
    my_arena_slot->head = dummy_slot.head;
    my_arena_slot->tail = dummy_slot.tail;
    // Release signal on behalf of previously spawned tasks (when this thread was not in arena yet)
    ITT_NOTIFY(sync_releasing, my_arena_slot);
    __TBB_store_with_release( my_arena_slot->task_pool, dummy_slot.task_pool );
    // We'll leave arena only when it's empty, so clean up local instances of indices.
    dummy_slot.head = dummy_slot.tail = 0;
}

void generic_scheduler::enter_arena() {
    __TBB_ASSERT ( my_arena, "no arena: initialization not completed?" );
#if __TBB_ARENA_PER_MASTER
    __TBB_ASSERT ( !in_arena(), "thread is already in arena?" );
    __TBB_ASSERT ( arena_index < my_arena->my_num_slots, "arena slot index is out-of-bound" );
#else /* !__TBB_ARENA_PER_MASTER */
    __TBB_ASSERT ( is_worker(), "only workers should use enter_arena()" );
    __TBB_ASSERT ( !in_arena(), "worker already in arena?" );
    __TBB_ASSERT ( arena_index < my_arena->prefix().number_of_workers, "invalid worker arena slot index" );
#endif /* !__TBB_ARENA_PER_MASTER */
    __TBB_ASSERT ( my_arena->slot[arena_index].task_pool == EmptyTaskPool, "someone else grabbed my arena slot?" );
    do_enter_arena();
}

#if !__TBB_ARENA_PER_MASTER
void generic_scheduler::try_enter_arena() {
    __TBB_ASSERT ( !is_worker(), "only masters should use try_enter_arena()" );
    __TBB_ASSERT ( my_arena, "no arena: initialization not completed?" );
    __TBB_ASSERT ( !in_arena(), "master already in arena?" );
    __TBB_ASSERT ( arena_index >= my_arena->prefix().number_of_workers && 
                   arena_index < my_arena->prefix().number_of_slots, "invalid arena slot hint value" );

    size_t h = arena_index;
    // We do not lock task pool upon successful entering arena
    if( my_arena->slot[h].task_pool != EmptyTaskPool || 
        __TBB_CompareAndSwapW( &my_arena->slot[h].task_pool, (intptr_t)LockedTaskPool, 
                                                          (intptr_t)EmptyTaskPool ) != (intptr_t)EmptyTaskPool )
    {
        // Hinted arena slot is already busy, try some of the others at random
        unsigned first = my_arena->prefix().number_of_workers,
                 last = my_arena->prefix().number_of_slots;
        unsigned n = last - first - 1;
        /// \todo Is this limit reasonable?
        size_t max_attempts = last - first;
        for (;;) {
            size_t k = first + random.get() % n;
            if( k >= h )
                ++k;    // Adjusts random distribution to exclude previously tried slot
            h = k;
            if( my_arena->slot[h].task_pool == EmptyTaskPool && 
                __TBB_CompareAndSwapW( &my_arena->slot[h].task_pool, (intptr_t)LockedTaskPool, 
                                                                  (intptr_t)EmptyTaskPool ) == (intptr_t)EmptyTaskPool )
            {
                break;
            }
            if ( --max_attempts == 0 ) {
                // After so many attempts we are still unable to find a vacant arena slot.
                // Cease the vain effort and work outside of arena for a while.
                return;
            }
        }
    }
    // Successfully claimed a slot in the arena.
    ITT_NOTIFY(sync_acquired, &my_arena->slot[h]);
    __TBB_ASSERT ( my_arena->slot[h].task_pool == LockedTaskPool, "arena slot is not actually acquired" );
    arena_index = h;
    do_enter_arena();
    attach_mailbox( affinity_id(h+1) );
}
#endif /* !__TBB_ARENA_PER_MASTER */

void generic_scheduler::leave_arena() {
    __TBB_ASSERT( in_arena(), "Not in arena" );
    // Do not reset arena_index. It will be used to (attempt to) re-acquire the slot next time
    __TBB_ASSERT( &my_arena->slot[arena_index] == my_arena_slot, "arena slot and slot index mismatch" );
    __TBB_ASSERT ( my_arena_slot->task_pool == LockedTaskPool, "Task pool must be locked when leaving arena" );
    __TBB_ASSERT ( my_arena_slot->head == my_arena_slot->tail, "Cannot leave arena when the task pool is not empty" );
#if !__TBB_ARENA_PER_MASTER
    if ( !is_worker() ) {
        my_affinity_id = 0;
        inbox.detach();
    }
#endif /* !__TBB_ARENA_PER_MASTER */
    ITT_NOTIFY(sync_releasing, &my_arena->slot[arena_index]);
    __TBB_store_with_release( my_arena_slot->task_pool, EmptyTaskPool );
    my_arena_slot = &dummy_slot;
}

#if __TBB_ARENA_PER_MASTER
generic_scheduler* generic_scheduler::create_worker( market& m, size_t index ) {
    generic_scheduler* s = allocate_scheduler( NULL, index );
#if __TBB_TASK_GROUP_CONTEXT
    s->dummy_task->prefix().context = &dummy_context;
    // Sync up the local cancellation state with the global one. No need for fence here.
    s->local_cancel_count = global_cancel_count;
#endif /* __TBB_TASK_GROUP_CONTEXT */
    s->my_market = &m;
    s->init_stack_info();
    return s;
}

#else /* !__TBB_ARENA_PER_MASTER */

generic_scheduler* generic_scheduler::create_worker( arena& a, size_t index ) {
    generic_scheduler* s = allocate_scheduler( &a, index );

    // Put myself into the arena
#if __TBB_TASK_GROUP_CONTEXT
    s->dummy_task->prefix().context = &dummy_context;
    // Sync up the local cancellation state with the global one. No need for fence here.
    s->local_cancel_count = global_cancel_count;
#endif /* __TBB_TASK_GROUP_CONTEXT */
    s->attach_mailbox( index+1 );
    s->init_stack_info();

    __TBB_store_with_release( a.prefix().worker_list[index].scheduler, s );
    return s;
}
#endif /* !__TBB_ARENA_PER_MASTER */

generic_scheduler* generic_scheduler::create_master( arena& a ) {
    generic_scheduler* s = allocate_scheduler( &a,
#if __TBB_ARENA_PER_MASTER
        0                   // Master thread always occupies the first slot
#else /* !__TBB_ARENA_PER_MASTER */
        null_arena_index    // Master thread will have to search for a vacant slot
#endif /* !__TBB_ARENA_PER_MASTER */
        );
    task& t = *s->dummy_task;
    s->innermost_running_task = &t;
    t.prefix().ref_count = 1;
    governor::sign_on(s);
    __TBB_ASSERT( &task::self()==&t, "governor::sign_on failed?" );
#if __TBB_ARENA_PER_MASTER
#if __TBB_TASK_GROUP_CONTEXT
    // Context to be used by root tasks by default (if the user has not specified one).
    // Allocation is done by NFS allocator because we cannot reuse memory allocated 
    // for task objects since the free list is empty at the moment.
    t.prefix().context = a.my_master_default_ctx = 
        new ( NFS_Allocate(sizeof(task_group_context), 1, NULL) ) task_group_context(task_group_context::isolated);
#endif
    s->my_market = a.my_market;
    __TBB_ASSERT( s->arena_index == 0, "Master thread must occupy the first slot in its arena" );
    s->attach_mailbox(1);
    a.slot[0].my_scheduler = s;
#if _WIN32|_WIN64
    __TBB_ASSERT( s->my_market, NULL );
    s->my_market->register_master( s->master_exec_resource );
#endif /* _WIN32|_WIN64 */
#else /* !__TBB_ARENA_PER_MASTER */
#if _WIN32|_WIN64
    s->register_master();
#endif 
#if __TBB_TASK_GROUP_CONTEXT
    // Context to be used by root tasks by default (if the user has not specified one).
    // Allocation is done by NFS allocator because we cannot reuse memory allocated 
    // for task objects since the free list is empty at the moment.
    t.prefix().context = new ( NFS_Allocate(sizeof(task_group_context), 1, NULL) ) task_group_context(task_group_context::isolated);
    scheduler_list_node_t &node = s->my_node;
    {
        mutex::scoped_lock lock(the_scheduler_list_mutex);
        node.my_next = the_scheduler_list_head.my_next;
        node.my_prev = &the_scheduler_list_head;
        the_scheduler_list_head.my_next->my_prev = &node;
        the_scheduler_list_head.my_next = &node;
#endif /* __TBB_TASK_GROUP_CONTEXT */
        unsigned last = a.prefix().number_of_slots,
                 cur_limit = a.prefix().limit;
        // This slot index assignment is just a hint to ...
        if ( cur_limit < last ) {
            // ... to prevent competition between the first few masters.
            s->arena_index = cur_limit++;
            // In the absence of exception handling this code is a subject to data 
            // race in case of multiple masters concurrently entering empty arena.
            // But it does not affect correctness, and can only result in a few 
            // masters competing for the same arena slot during the first acquisition.
            // The cost of competition is low in comparison to that of oversubscription.
            a.prefix().limit = cur_limit;
        }
        else {
            // ... to minimize the probability of competition between multiple masters.
            unsigned first = a.prefix().number_of_workers;
            s->arena_index = first + s->random.get() % (last - first);
        }
#if __TBB_TASK_GROUP_CONTEXT
    }
#endif
#endif /* !__TBB_ARENA_PER_MASTER */
    s->init_stack_info();
#if __TBB_TASK_GROUP_CONTEXT
    // Sync up the local cancellation state with the global one. No need for fence here.
    s->local_cancel_count = global_cancel_count;
#endif
#if __TBB_SCHEDULER_OBSERVER
    // Process any existing observers.
    s->notify_entry_observers();
#endif /* __TBB_SCHEDULER_OBSERVER */
    return s;
}

void generic_scheduler::cleanup_worker( void* arg, bool is_worker ) {
    generic_scheduler& s = *(generic_scheduler*)arg;
    __TBB_ASSERT( s.dummy_slot.task_pool, "cleaning up worker with missing task pool" );
// APM TODO: Decide how observers should react to each entry/leave to/from arena
#if __TBB_SCHEDULER_OBSERVER
    s.notify_exit_observers( is_worker );
#endif /* __TBB_SCHEDULER_OBSERVER */
    // When comparing "head" and "tail" indices ">=" is used because this worker's
    // task pool may still be published in the arena, and thieves can optimistically
    // bump "head" (and then roll back).
    __TBB_ASSERT( s.my_arena_slot->task_pool == EmptyTaskPool || s.my_arena_slot->head >= s.my_arena_slot->tail, 
                  "worker has unfinished work at run down" );
    s.free_scheduler();
}

void generic_scheduler::cleanup_master() {
    generic_scheduler& s = *this; // for similarity with cleanup_worker
    __TBB_ASSERT( s.dummy_slot.task_pool, "cleaning up master with missing task pool" );
#if __TBB_SCHEDULER_OBSERVER
    s.notify_exit_observers(/*is_worker=*/false);
#endif /* __TBB_SCHEDULER_OBSERVER */
    if ( !local_task_pool_empty() ) {
        __TBB_ASSERT ( governor::is_set(this), "TLS slot is cleared before the task pool cleanup" );
        s.local_wait_for_all( *s.dummy_task, NULL );
        __TBB_ASSERT ( governor::is_set(this), "Other thread reused our TLS key during the task pool cleanup" );
    }
#if __TBB_ARENA_PER_MASTER
#if _WIN32|_WIN64
    __TBB_ASSERT( s.my_market, NULL );
    s.my_market->unregister_master( s.master_exec_resource );
#endif /* _WIN32|_WIN64 */
    arena* a = s.my_arena;
#if __TBB_STATISTICS
    *a->slot[0].my_counters += s.my_counters;
#endif /* __TBB_STATISTICS */
#else /* !__TBB_ARENA_PER_MASTER */
#if _WIN32|_WIN64
    s.unregister_master();
#endif /* _WIN32|_WIN64 */
#endif /* __TBB_ARENA_PER_MASTER */
    s.free_scheduler();
#if __TBB_ARENA_PER_MASTER
    a->slot[0].my_scheduler = NULL;
    // Do not close arena if some fire-and-forget tasks remain; workers should care of it.
    if( a->my_task_stream.empty() && a->pool_state.fetch_and_store(arena::SNAPSHOT_EMPTY)!=arena::SNAPSHOT_EMPTY )
        a->my_market->adjust_demand( *a, -(int)a->my_max_num_workers );
#if __TBB_STATISTICS_EARLY_DUMP
    GATHER_STATISTIC( a->dump_arena_statistics() );
#endif
    if ( --a->my_num_threads_active==0 && a->pool_state==arena::SNAPSHOT_EMPTY )
        a->close_arena();
#else /* !__TBB_ARENA_PER_MASTER */
    governor::finish_with_arena();
#endif /* !__TBB_ARENA_PER_MASTER */
}

#if __TBB_SCHEDULER_OBSERVER
    void generic_scheduler::notify_entry_observers() {
        local_last_observer_proxy = observer_proxy::process_list(local_last_observer_proxy,is_worker(),/*is_entry=*/true);
    }

    void generic_scheduler::notify_exit_observers( bool is_worker ) {
        observer_proxy::process_list(local_last_observer_proxy,is_worker,/*is_entry=*/false);
    }
#endif /* __TBB_SCHEDULER_OBSERVER */

} // namespace internal
} // namespace tbb

