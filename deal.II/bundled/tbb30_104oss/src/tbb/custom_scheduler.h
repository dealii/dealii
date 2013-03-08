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

#ifndef _TBB_custom_scheduler_H
#define _TBB_custom_scheduler_H

#include "scheduler.h"
#include "observer_proxy.h"
#include "itt_notify.h"

namespace tbb {
namespace internal {

//! Amount of time to pause between steals.
/** The default values below were found to be best empirically for K-Means
    on the 32-way Altix and 4-way (*2 for HT) fxqlin04. */
#if __TBB_ipf
static const long PauseTime = 1500;
#else 
static const long PauseTime = 80;
#endif

//------------------------------------------------------------------------
//! Traits classes for scheduler
//------------------------------------------------------------------------

struct DefaultSchedulerTraits {
    static const bool itt_possible = true;
    static const bool has_slow_atomic = false;
};

struct IntelSchedulerTraits {
    static const bool itt_possible = false;
#if __TBB_x86_32||__TBB_x86_64
    static const bool has_slow_atomic = true;
#else
    static const bool has_slow_atomic = false;
#endif /* __TBB_x86_32||__TBB_x86_64 */
};

//------------------------------------------------------------------------
// custom_scheduler
//------------------------------------------------------------------------

//! A scheduler with a customized evaluation loop.
/** The customization can use SchedulerTraits to make decisions without needing a run-time check. */
template<typename SchedulerTraits>
class custom_scheduler: private generic_scheduler {
    typedef custom_scheduler<SchedulerTraits> scheduler_type;

    //! Scheduler loop that dispatches tasks.
    /** If child is non-NULL, it is dispatched first.
        Then, until "parent" has a reference count of 1, other task are dispatched or stolen. */
    /*override*/
    void local_wait_for_all( task& parent, task* child );

    //! Entry point from client code to the scheduler loop that dispatches tasks. 
    /** The method is virtual, but the *this object is used only for sake of dispatching on the correct vtable,
        not necessarily the correct *this object.  The correct *this object is looked up in TLS. */
    /*override*/
    void wait_for_all( task& parent, task* child ) {
        static_cast<custom_scheduler*>(governor::local_scheduler())->scheduler_type::local_wait_for_all( parent, child );
    }

    //! Construct a custom_scheduler
    custom_scheduler( arena* a, size_t index ) : generic_scheduler(a, index) {}

    //! Decrements ref_count of a predecessor.
    /** If it achieves 0, the predecessor is scheduled for execution.
        When changing, remember that this is a hot path function. */
    void tally_completion_of_predecessor( task& s, task*& bypass_slot ) {
        task_prefix& p = s.prefix();
        if( SchedulerTraits::itt_possible )
            ITT_NOTIFY(sync_releasing, &p.ref_count);
        if( SchedulerTraits::has_slow_atomic && p.ref_count==1 ) {
            p.ref_count=0;
        } else {
            if( __TBB_FetchAndDecrementWrelease(&p.ref_count) > 1 ) // more references exist
                return;
        }
        __TBB_ASSERT(p.ref_count==0, "completion of task caused predecessor's reference count to underflow");
        if( SchedulerTraits::itt_possible )
            ITT_NOTIFY(sync_acquired, &p.ref_count);
#if TBB_USE_ASSERT
        p.extra_state &= ~es_ref_count_active;
#endif /* TBB_USE_ASSERT */

        if( bypass_slot==NULL )
            bypass_slot = &s;
        else
            local_spawn( s, s.prefix().next );
    }

public:
    static generic_scheduler* allocate_scheduler( arena* a, size_t index ) {
#if !__TBB_ARENA_PER_MASTER
        __TBB_ASSERT( a, "missing arena" );
#endif /* !__TBB_ARENA_PER_MASTER */
        scheduler_type* s = (scheduler_type*)NFS_Allocate(sizeof(scheduler_type),1,NULL);
        new( s ) scheduler_type( a, index );
        s->assert_task_pool_valid();
        ITT_SYNC_CREATE(s, SyncType_Scheduler, SyncObj_TaskPoolSpinning);
        return s;
    }

    //! Try getting a task from the mailbox or stealing from another scheduler.
    /** Returns the stolen task or NULL if all attempts fail. */
    /* override */ task* receive_or_steal_task( reference_count&, bool );

}; // class custom_scheduler<>

//------------------------------------------------------------------------
// custom_scheduler methods
//------------------------------------------------------------------------

template<typename SchedulerTraits>
task* custom_scheduler<SchedulerTraits>::receive_or_steal_task( reference_count& completion_ref_count,
                                                                bool return_if_no_work ) {
    task* t = NULL;
    inbox.set_is_idle( true );
    // The state "failure_count==-1" is used only when itt_possible is true,
    // and denotes that a sync_prepare has not yet been issued.
    for( int failure_count = -static_cast<int>(SchedulerTraits::itt_possible);; ++failure_count) {
        if( completion_ref_count==1 ) {
            if( SchedulerTraits::itt_possible ) {
                if( failure_count!=-1 ) {
                    ITT_NOTIFY(sync_prepare, &completion_ref_count);
                    // Notify Intel(R) Thread Profiler that thread has stopped spinning.
                    ITT_NOTIFY(sync_acquired, this);
                }
                ITT_NOTIFY(sync_acquired, &completion_ref_count);
            }
            inbox.set_is_idle( false );
            return NULL;
        }
#if __TBB_ARENA_PER_MASTER
        size_t n = my_arena->my_limit;
        __TBB_ASSERT( arena_index < n, NULL );
#else /* !__TBB_ARENA_PER_MASTER */
        size_t n = my_arena->prefix().limit;
#endif /* !__TBB_ARENA_PER_MASTER */
        if( n>1 ) {
            if( my_affinity_id && (t=get_mailbox_task()) ) {
                GATHER_STATISTIC( ++my_counters.mails_received );
            }
#if __TBB_ARENA_PER_MASTER
            // Check if there are tasks in starvation-resistant stream.
            // Only allowed for workers with empty stack, which is identified by return_if_no_work.
            else if ( return_if_no_work && (t=dequeue_task()) ) {
                // just proceed with the obtained task
            }
            // Check if the resource manager requires our arena to relinquish some threads 
            else if ( return_if_no_work && (my_arena->my_num_workers_allotted < my_arena->num_workers_active()) ) {
                if( SchedulerTraits::itt_possible ) {
                    if( failure_count!=-1 )
                        ITT_NOTIFY(sync_cancel, this);
                }
                return NULL;
            }
#endif /* __TBB_ARENA_PER_MASTER */
            else {
                // Try to steal a task from a random victim.
                if ( !can_steal() )
                    goto fail;
                size_t k = random.get() % (n-1);
                arena_slot* victim = &my_arena->slot[k];
                // The following condition excludes the master that might have 
                // already taken our previous place in the arena from the list .
                // of potential victims. But since such a situation can take 
                // place only in case of significant oversubscription, keeping
                // the checks simple seems to be preferable to complicating the code.
                if( k >= arena_index )
                    ++victim;               // Adjusts random distribution to exclude self
                t = steal_task( *victim );
                if( !t ) goto fail;
                if( is_proxy(*t) ) {
                    t = strip_proxy((task_proxy*)t);
                    if( !t ) goto fail;
                    GATHER_STATISTIC( ++my_counters.proxies_stolen );
                }
                t->prefix().extra_state |= es_task_is_stolen;
                if( is_version_3_task(*t) ) {
                    innermost_running_task = t;
                    t->note_affinity( my_affinity_id );
                }
                GATHER_STATISTIC( ++my_counters.steals_committed );
            }
            __TBB_ASSERT(t,NULL);
#if __TBB_SCHEDULER_OBSERVER
            // No memory fence required for read of global_last_observer_proxy, because prior fence on steal/mailbox suffices.
            if( local_last_observer_proxy!=global_last_observer_proxy ) {
                notify_entry_observers();
            }
#endif /* __TBB_SCHEDULER_OBSERVER */
            if( SchedulerTraits::itt_possible ) {
                if( failure_count!=-1 ) {
                    // FIXME - might be victim, or might be selected from a mailbox
                    // Notify Intel(R) Thread Profiler that thread has stopped spinning.
                    ITT_NOTIFY(sync_acquired, this);
                }
            }
            inbox.set_is_idle( false );
            break; // jumps to: return t;
        }
fail:
        GATHER_STATISTIC( ++my_counters.steals_failed );
        if( SchedulerTraits::itt_possible && failure_count==-1 ) {
            // The first attempt to steal work failed, so notify Intel(R) Thread Profiler that
            // the thread has started spinning.  Ideally, we would do this notification
            // *before* the first failed attempt to steal, but at that point we do not
            // know that the steal will fail.
            ITT_NOTIFY(sync_prepare, this);
            failure_count = 0;
        }
        // Pause, even if we are going to yield, because the yield might return immediately.
        __TBB_Pause(PauseTime);
        int yield_threshold = 2*int(n);
        if( failure_count>=yield_threshold ) {
            __TBB_Yield();
            if( failure_count>=yield_threshold+100 ) {
                // When a worker thread has nothing to do, return it to RML.
                // For purposes of affinity support, the thread is considered idle while in RML.
                if( return_if_no_work && my_arena->is_out_of_work() ) {
                    if( SchedulerTraits::itt_possible ) {
                        if( failure_count!=-1 )
                            ITT_NOTIFY(sync_cancel, this);
                    }
                    return NULL;
                }
                failure_count = yield_threshold;
            }
        }
    }
    return t;
}

template<typename SchedulerTraits>
void custom_scheduler<SchedulerTraits>::local_wait_for_all( task& parent, task* child ) {
    __TBB_ASSERT( governor::is_set(this), NULL );
    if( child ) {
        child->prefix().owner = this;
    }
    __TBB_ASSERT( parent.ref_count() >= (child && child->parent() == &parent ? 2 : 1), "ref_count is too small" );
    assert_task_pool_valid();
    // Using parent's refcount in sync_prepare (in the stealing loop below) is 
    // a workaround for TP. We need to name it here to display correctly in Ampl.
    if( SchedulerTraits::itt_possible )
        ITT_SYNC_CREATE(&parent.prefix().ref_count, SyncType_Scheduler, SyncObj_TaskStealingLoop);
#if __TBB_TASK_GROUP_CONTEXT
    __TBB_ASSERT( parent.prefix().context || (is_worker() && &parent == dummy_task), "parent task does not have context" );
#endif /* __TBB_TASK_GROUP_CONTEXT */
    task* t = child;
    // Constant all_local_work_done is an unreacheable refcount value that prevents
    // early quitting the dispatch loop. It is defined to be in the middle of the range 
    // of negative values representable by the reference_count type.
    static const reference_count 
        // For normal dispatch loops
        parents_work_done = 1,
        // For termination dispatch loops in masters
        all_local_work_done = (reference_count)3 << (sizeof(reference_count) * 8 - 2);
    reference_count quit_point;
    if( innermost_running_task == dummy_task ) {
        // We are in the outermost task dispatch loop of a master thread,
        __TBB_ASSERT( !is_worker(), NULL );
        quit_point = &parent == dummy_task ? all_local_work_done : parents_work_done;
    } else {
        quit_point = parents_work_done;
    }
    task* old_innermost_running_task = innermost_running_task;
#if __TBB_TASK_GROUP_CONTEXT && TBB_USE_EXCEPTIONS
exception_was_caught:
    try {
#endif /* __TBB_TASK_GROUP_CONTEXT && TBB_USE_EXCEPTIONS */
    // Outer loop steals tasks when necessary.
    for(;;) {
        // Middle loop evaluates tasks that are pulled off "array".
        do {
            // Inner loop evaluates tasks that are handed directly to us by other tasks.
            while(t) {
                __TBB_ASSERT( inbox.is_idle_state(false), NULL );
#if TBB_USE_ASSERT
                __TBB_ASSERT(!is_proxy(*t),"unexpected proxy");
                __TBB_ASSERT( t->prefix().owner==this, NULL );
#if __TBB_TASK_GROUP_CONTEXT
                if ( !t->prefix().context->my_cancellation_requested ) 
#endif
                    __TBB_ASSERT( 1L<<t->state() & (1L<<task::allocated|1L<<task::ready|1L<<task::reexecute), NULL );
                assert_task_pool_valid();
#endif /* TBB_USE_ASSERT */
                task* t_next = NULL;
                innermost_running_task = t;
                t->prefix().state = task::executing;
#if __TBB_TASK_GROUP_CONTEXT
                if ( !t->prefix().context->my_cancellation_requested )
#endif
                {
                    GATHER_STATISTIC( ++my_counters.tasks_executed );
#if __TBB_TASK_GROUP_CONTEXT
                    if( SchedulerTraits::itt_possible )
                        ITT_STACK(callee_enter, t->prefix().context->itt_caller);
#endif
                    t_next = t->execute();
#if __TBB_TASK_GROUP_CONTEXT
                    if( SchedulerTraits::itt_possible )
                        ITT_STACK(callee_leave, t->prefix().context->itt_caller);
#endif
                    if (t_next) {
                        __TBB_ASSERT( t_next->state()==task::allocated,
                                "if task::execute() returns task, it must be marked as allocated" );
                        t_next->prefix().extra_state &= ~es_task_is_stolen;
#if TBB_USE_ASSERT
                        affinity_id next_affinity=t_next->prefix().affinity;
                        if (next_affinity != 0 && next_affinity != my_affinity_id)
                            GATHER_STATISTIC( ++my_counters.affinity_ignored );
#endif
                    }
                }
                assert_task_pool_valid();
                switch( task::state_type(t->prefix().state) ) {
                    case task::executing: {
                        task* s = t->parent();
                        __TBB_ASSERT( innermost_running_task==t, NULL );
                        __TBB_ASSERT( t->prefix().ref_count==0, "Task still has children after it has been executed" );
                        t->~task();
                        if( s )
                            tally_completion_of_predecessor(*s, t_next);
                        free_task<no_hint>( *t );
                        assert_task_pool_valid();
                        break;
                    }

                    case task::recycle: // set by recycle_as_safe_continuation()
                        t->prefix().state = task::allocated;
                        __TBB_ASSERT( t_next != t, "a task returned from method execute() can not be recycled in another way" );
                        t->prefix().extra_state &= ~es_task_is_stolen;
                        // for safe continuation, need atomically decrement ref_count;
                        tally_completion_of_predecessor(*t, t_next);
                        assert_task_pool_valid();
                        break;

                    case task::reexecute: // set by recycle_to_reexecute()
                        __TBB_ASSERT( t_next, "reexecution requires that method execute() return another task" );
                        __TBB_ASSERT( t_next != t, "a task returned from method execute() can not be recycled in another way" );
                        t->prefix().state = task::allocated;
                        t->prefix().extra_state &= ~es_task_is_stolen;
                        local_spawn( *t, t->prefix().next );
                        assert_task_pool_valid();
                        break;
                    case task::allocated:
                        t->prefix().extra_state &= ~es_task_is_stolen;
                        break;
#if TBB_USE_ASSERT
                    case task::ready:
                        __TBB_ASSERT( false, "task is in READY state upon return from method execute()" );
                        break;
                    default:
                        __TBB_ASSERT( false, "illegal state" );
#else
                    default: // just to shut up some compilation warnings
                        break;
#endif /* TBB_USE_ASSERT */
                }

                if( t_next ) {
                    // The store here has a subtle secondary effect - it fetches *t_next into cache.
                    t_next->prefix().owner = this;
                    GATHER_STATISTIC( ++my_counters.spawns_bypassed );
                }
                t = t_next;
            } // end of scheduler bypass loop
            assert_task_pool_valid();

            if ( parent.prefix().ref_count == quit_point )
                break;
            t = get_task();
            __TBB_ASSERT(!t || !is_proxy(*t),"unexpected proxy");
#if TBB_USE_ASSERT
            assert_task_pool_valid();
            if(t) {
                assert_task_valid(*t);
                __TBB_ASSERT( t->prefix().owner==this, "thread got task that it does not own" );
            }
#endif /* TBB_USE_ASSERT */
        } while( t ); // end of local task array processing loop

        if ( quit_point == all_local_work_done ) {
            __TBB_ASSERT( my_arena_slot == &dummy_slot && my_arena_slot->head == 0 && my_arena_slot->tail == 0, NULL );
            innermost_running_task = old_innermost_running_task;
            return;
        }
#if __TBB_ARENA_PER_MASTER
        __TBB_ASSERT( my_arena->my_max_num_workers > 0 || parent.prefix().ref_count == 1, "deadlock detected" );
#else /* !__TBB_ARENA_PER_MASTER */
        __TBB_ASSERT( my_arena->prefix().number_of_workers>0||parent.prefix().ref_count==1, "deadlock detected" );
#endif /* !__TBB_ARENA_PER_MASTER */
        // old_innermost_running_task is NULL *iff* a worker thread is in its "inborn" dispath loop
        // (i.e. its execution stack is empty), and it should return from there if no work is available.
        t = receive_or_steal_task( parent.prefix().ref_count, !old_innermost_running_task );
        if (!t) {
            if( parent.prefix().ref_count==1 ) goto done;
            __TBB_ASSERT( is_worker() && !old_innermost_running_task, "a thread exits dispatch loop prematurely" );
            innermost_running_task = NULL;
            return;
        }
        __TBB_ASSERT(t,NULL);
        __TBB_ASSERT(!is_proxy(*t),"unexpected proxy");
        t->prefix().owner = this;
    } // end of stealing loop
#if __TBB_TASK_GROUP_CONTEXT && TBB_USE_EXCEPTIONS
    } TbbCatchAll( t->prefix().context );

    if( task::state_type(t->prefix().state) == task::recycle ) { // state set by recycle_as_safe_continuation()
        t->prefix().state = task::allocated;
        // for safe continuation, need to atomically decrement ref_count;
        if( SchedulerTraits::itt_possible )
            ITT_NOTIFY(sync_releasing, &t->prefix().ref_count);
        if( __TBB_FetchAndDecrementWrelease(&t->prefix().ref_count)==1 ) {
            if( SchedulerTraits::itt_possible )
                ITT_NOTIFY(sync_acquired, &t->prefix().ref_count);
        }else{
            t = NULL;
        }
    }
    goto exception_was_caught;
#endif /* __TBB_TASK_GROUP_CONTEXT && TBB_USE_EXCEPTIONS */
done:
    if ( !ConcurrentWaitsEnabled(parent) )
        parent.prefix().ref_count = 0;
#if TBB_USE_ASSERT
    parent.prefix().extra_state &= ~es_ref_count_active;
#endif /* TBB_USE_ASSERT */
    innermost_running_task = old_innermost_running_task;
#if __TBB_TASK_GROUP_CONTEXT
    __TBB_ASSERT(parent.prefix().context && dummy_task->prefix().context, NULL);
    task_group_context* parent_ctx = parent.prefix().context;
    if ( parent_ctx->my_cancellation_requested ) {
        task_group_context::exception_container_type *pe = parent_ctx->my_exception;
        if ( innermost_running_task == dummy_task && parent_ctx == dummy_task->prefix().context ) {
            // We are in the outermost task dispatch loop of a master thread, and 
            // the whole task tree has been collapsed. So we may clear cancellation data.
            parent_ctx->my_cancellation_requested = 0;
            __TBB_ASSERT(dummy_task->prefix().context == parent_ctx || !CancellationInfoPresent(*dummy_task), 
                         "Unexpected exception or cancellation data in the dummy task");
            // If possible, add assertion that master's dummy task context does not have children
        }
        if ( pe )
            pe->throw_self();
    }
    __TBB_ASSERT(!is_worker() || !CancellationInfoPresent(*dummy_task), 
                 "Worker's dummy task context modified");
    __TBB_ASSERT(innermost_running_task != dummy_task || !CancellationInfoPresent(*dummy_task), 
                 "Unexpected exception or cancellation data in the master's dummy task");
#endif /* __TBB_TASK_GROUP_CONTEXT */
    assert_task_pool_valid();
}

} // namespace internal
} // namespace tbb

#endif /* _TBB_custom_scheduler_H */
