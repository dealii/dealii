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

#include "arena.h"
#include "governor.h"
#include "scheduler.h"
#include "itt_notify.h"

#if __TBB_STATISTICS_STDOUT
#include <cstdio>
#endif

namespace tbb {
namespace internal {

#if !__TBB_ARENA_PER_MASTER
//------------------------------------------------------------------------
// UnpaddedArenaPrefix
//------------------------------------------------------------------------
inline arena& UnpaddedArenaPrefix::Arena() {
    return *static_cast<tbb::internal::arena*>(static_cast<void*>( static_cast<ArenaPrefix*>(this)+1 ));
}

void UnpaddedArenaPrefix::process( job& j ) {
    generic_scheduler& s = static_cast<generic_scheduler&>(j);
    __TBB_ASSERT( governor::is_set(&s), NULL );
    __TBB_ASSERT( !s.innermost_running_task, NULL );
    // Try to steal a task.
    // Passing reference count is technically unnecessary in this context,
    // but omitting it here would add checks inside the function.
    task* t = s.receive_or_steal_task( s.dummy_task->prefix().ref_count, /*return_if_no_work=*/true );
    if (t) {
        // A side effect of receive_or_steal_task is that innermost_running_task can be set.
        // But for the outermost dispatch loop of a worker it has to be NULL.
        s.innermost_running_task = NULL;
        s.local_wait_for_all(*s.dummy_task,t);
    }
    __TBB_ASSERT( s.inbox.is_idle_state(true), NULL );
    __TBB_ASSERT( !s.innermost_running_task, NULL );
}

void UnpaddedArenaPrefix::cleanup( job& j ) {
    generic_scheduler& s = static_cast<generic_scheduler&>(j);
    if( !governor::is_set( &s ) ) {
        bool is_master = governor::is_set( NULL );
        governor::assume_scheduler( &s );
        generic_scheduler::cleanup_worker( &s, !is_master );
        governor::assume_scheduler( NULL );
    } else {
        generic_scheduler::cleanup_worker( &s, true );
    }
}

void UnpaddedArenaPrefix::acknowledge_close_connection() {
    Arena().free_arena();
}

::rml::job* UnpaddedArenaPrefix::create_one_job() {
    generic_scheduler* s = generic_scheduler::create_worker( Arena(), next_job_index++ );
    governor::sign_on(s);
    return s;
}
#endif /* !__TBB_ARENA_PER_MASTER */

//------------------------------------------------------------------------
// arena
//------------------------------------------------------------------------

#if __TBB_ARENA_PER_MASTER

void arena::process( generic_scheduler& s ) {
    __TBB_ASSERT( is_alive(my_guard), NULL );
    __TBB_ASSERT( governor::is_set(&s), NULL );
    __TBB_ASSERT( !s.innermost_running_task, NULL );

    __TBB_ASSERT( my_num_slots != 1, NULL );
    // Start search for an empty slot from the one we occupied the last time
    unsigned index = s.arena_index < my_num_slots ? s.arena_index : s.random.get() % (my_num_slots - 1) + 1,
             end = index;
    __TBB_ASSERT( index != 0, "A worker cannot occupy slot 0" );
    __TBB_ASSERT( index < my_num_slots, NULL );

    // Find a vacant slot
    for ( ;; ) {
        if ( !slot[index].my_scheduler && __TBB_CompareAndSwapW( &slot[index].my_scheduler, (intptr_t)&s, 0 ) == 0 )
            break;
        if ( ++index == my_num_slots )
            index = 1;
        if ( index == end ) {
            // Likely this arena is already saturated
            if ( --my_num_threads_active == 0 )
                close_arena();
            return;
        }
    }
    ITT_NOTIFY(sync_acquired, &slot[index]);
    s.my_arena = this;
    s.arena_index = index;
    s.attach_mailbox( affinity_id(index+1) );

    slot[index].hint_for_push = index ^ unsigned(&s-(generic_scheduler*)NULL)>>16; // randomizer seed
    slot[index].hint_for_pop  = index; // initial value for round-robin

    unsigned new_limit = index + 1;
    unsigned old_limit = my_limit;
    while ( new_limit > old_limit ) {
        if ( my_limit.compare_and_swap(new_limit, old_limit) == old_limit )
            break;
        old_limit = my_limit;
    }

    for ( ;; ) {
        // Try to steal a task.
        // Passing reference count is technically unnecessary in this context,
        // but omitting it here would add checks inside the function.
        __TBB_ASSERT( is_alive(my_guard), NULL );
        task* t = s.receive_or_steal_task( s.dummy_task->prefix().ref_count, /*return_if_no_work=*/true );
        if (t) {
            // A side effect of receive_or_steal_task is that innermost_running_task can be set.
            // But for the outermost dispatch loop of a worker it has to be NULL.
            s.innermost_running_task = NULL;
            s.local_wait_for_all(*s.dummy_task,t);
        }
        ++my_num_threads_leaving;
        __TBB_ASSERT ( slot[index].head == slot[index].tail, "Worker cannot leave arena while its task pool is not empty" );
        __TBB_ASSERT( slot[index].task_pool == EmptyTaskPool, "Empty task pool is not marked appropriately" );
        // Revalidate quitting condition
        // This check prevents relinquishing more than necessary workers because 
        // of the non-atomicity of the decision making procedure
        if ( num_workers_active() >= my_num_workers_allotted || !my_num_workers_requested )
            break;
        --my_num_threads_leaving;
        __TBB_ASSERT( !slot[0].my_scheduler || my_num_threads_active > 0, "Who requested more workers after the last one left the dispatch loop and the master's gone?" );
    }
#if __TBB_STATISTICS
    ++s.my_counters.arena_roundtrips;
    *slot[index].my_counters += s.my_counters;
    s.my_counters.reset();
#endif /* __TBB_STATISTICS */
    __TBB_store_with_release( slot[index].my_scheduler, (generic_scheduler*)NULL );
    s.inbox.detach();
    __TBB_ASSERT( s.inbox.is_idle_state(true), NULL );
    __TBB_ASSERT( !s.innermost_running_task, NULL );
    __TBB_ASSERT( is_alive(my_guard), NULL );
    // Decrementing my_num_threads_active first prevents extra workers from leaving
    // this arena prematurely, but can result in some workers returning back just
    // to repeat the escape attempt. If instead my_num_threads_leaving is decremented
    // first, the result is the opposite - premature leaving is allowed and gratuitous
    // return is prevented. Since such a race has any likelihood only when multiple
    // workers are in the stealing loop, and consequently there is a lack of parallel
    // work in this arena, we'd rather let them go out and try get employment in 
    // other arenas (before returning into this one again).
    --my_num_threads_leaving;
    if ( !--my_num_threads_active )
        close_arena();
}

arena::arena ( market& m, unsigned max_num_workers ) {
    __TBB_ASSERT( !my_guard, "improperly allocated arena?" );
    __TBB_ASSERT( sizeof(slot[0]) % NFS_GetLineSize()==0, "arena::slot size not multiple of cache line size" );
    __TBB_ASSERT( (uintptr_t)this % NFS_GetLineSize()==0, "arena misaligned" );
    my_market = &m;
    my_limit = 1;
    // Two slots are mandatory: for the master, and for 1 worker (required to support starvation resistant tasks).
    my_num_slots = max(2u, max_num_workers + 1);
    my_max_num_workers = max_num_workers;
    my_num_threads_active = 1; // accounts for the master
    __TBB_ASSERT ( my_max_num_workers < my_num_slots, NULL );
    // Construct mailboxes. Mark internal synchronization elements for the tools.
    for( unsigned i = 0; i < my_num_slots; ++i ) {
        __TBB_ASSERT( !slot[i].my_scheduler && !slot[i].task_pool, NULL );
        ITT_SYNC_CREATE(slot + i, SyncType_Scheduler, SyncObj_WorkerTaskPool);
        mailbox(i+1).construct();
        ITT_SYNC_CREATE(&mailbox(i+1), SyncType_Scheduler, SyncObj_Mailbox);
#if __TBB_STATISTICS
        slot[i].my_counters = new ( NFS_Allocate(sizeof(statistics_counters), 1, NULL) ) statistics_counters;
#endif /* __TBB_STATISTICS */
    }
    my_task_stream.initialize(my_num_slots);
    ITT_SYNC_CREATE(&my_task_stream, SyncType_Scheduler, SyncObj_TaskStream);
    my_mandatory_concurrency = false;
#if __TBB_TASK_GROUP_CONTEXT
    my_master_default_ctx = NULL;
#endif
}

arena& arena::allocate_arena( market& m, unsigned max_num_workers ) {
    __TBB_ASSERT( sizeof(base_type) + sizeof(arena_slot) == sizeof(arena), "All arena data fields must go to arena_base" );
    __TBB_ASSERT( sizeof(base_type) % NFS_GetLineSize() == 0, "arena slots area misaligned: wrong padding" );
    __TBB_ASSERT( sizeof(mail_outbox) == NFS_MaxLineSize, "Mailbox padding is wrong" );

    unsigned num_slots = max(2u, max_num_workers + 1);
    size_t n = sizeof(base_type) + num_slots * (sizeof(mail_outbox) + sizeof(arena_slot));

    unsigned char* storage = (unsigned char*)NFS_Allocate( n, 1, NULL );
    // Zero all slots to indicate that they are empty
    memset( storage, 0, n );
    return *new( storage + num_slots * sizeof(mail_outbox) ) arena(m, max_num_workers);
}

void arena::free_arena () {
    __TBB_ASSERT( !my_num_threads_active, "There are threads in the dying arena" );
    poison_value( my_guard );
    intptr_t drained = 0;
    for ( unsigned i = 1; i <= my_num_slots; ++i )
        drained += mailbox(i).drain();
    __TBB_ASSERT(my_task_stream.empty() && my_task_stream.drain()==0, "Not all enqueued tasks were executed");
#if __TBB_COUNT_TASK_NODES
    my_market->update_task_node_count( -drained );
#endif /* __TBB_COUNT_TASK_NODES */
    my_market->release();
#if __TBB_TASK_GROUP_CONTEXT
    __TBB_ASSERT( my_master_default_ctx, "Master thread never entered the arena?" );
    my_master_default_ctx->~task_group_context();
    NFS_Free(my_master_default_ctx);
#endif /* __TBB_TASK_GROUP_CONTEXT */
#if __TBB_STATISTICS
    for( unsigned i = 0; i < my_num_slots; ++i )
        NFS_Free( slot[i].my_counters );
#endif /* __TBB_STATISTICS */
    void* storage  = &mailbox(my_num_slots);
    this->~arena();
    NFS_Free( storage );
}

void arena::close_arena () {
#if !__TBB_STATISTICS_EARLY_DUMP
    GATHER_STATISTIC( dump_arena_statistics() );
#endif
    my_market->detach_arena( *this );
    free_arena();
}

#if __TBB_STATISTICS
void arena::dump_arena_statistics () {
    statistics_counters total;
    for( unsigned i = 0; i < my_num_slots; ++i ) {
#if __TBB_STATISTICS_EARLY_DUMP
        generic_scheduler* s = slot[i].my_scheduler;
        if ( s )
            *slot[i].my_counters += s->my_counters;
#else
        __TBB_ASSERT( !slot[i].my_scheduler, NULL );
#endif
        if ( i != 0 ) {
            total += *slot[i].my_counters;
            dump_statistics( *slot[i].my_counters, i );
        }
    }
    dump_statistics( *slot[0].my_counters, 0 );
#if __TBB_STATISTICS_STDOUT
    printf( "----------------------------------------------\n" );
    dump_statistics( total, workers_counters_total );
    total += *slot[0].my_counters;
    dump_statistics( total, arena_counters_total );
    printf( "==============================================\n" );
#endif /* __TBB_STATISTICS_STDOUT */
}
#endif /* __TBB_STATISTICS */

#else /* !__TBB_ARENA_PER_MASTER */

arena* arena::allocate_arena( unsigned number_of_slots, unsigned number_of_workers, stack_size_type stack_size ) {
    __TBB_ASSERT( sizeof(ArenaPrefix) % NFS_GetLineSize()==0, "ArenaPrefix not multiple of cache line size" );
    __TBB_ASSERT( sizeof(mail_outbox)==NFS_MaxLineSize, NULL );
    __TBB_ASSERT( stack_size>0, NULL );

    size_t n = sizeof(ArenaPrefix) + number_of_slots*(sizeof(mail_outbox)+sizeof(arena_slot));

    unsigned char* storage = (unsigned char*)NFS_Allocate( n, 1, NULL );
    // Zero all slots to indicate that they are empty
    memset( storage, 0, n );
    arena* a = (arena*)(storage + sizeof(ArenaPrefix)+ number_of_slots*(sizeof(mail_outbox)));
    __TBB_ASSERT( sizeof(a->slot[0]) % NFS_GetLineSize()==0, "arena::slot size not multiple of cache line size" );
    __TBB_ASSERT( (uintptr_t)a % NFS_GetLineSize()==0, NULL );
    new( &a->prefix() ) ArenaPrefix( number_of_slots, number_of_workers );

    // Allocate the worker_list
    WorkerDescriptor * w = new WorkerDescriptor[number_of_workers];
    memset( w, 0, sizeof(WorkerDescriptor)*(number_of_workers));
    a->prefix().worker_list = w;

    // Construct mailboxes.
    for( unsigned j=1; j<=number_of_slots; ++j ) 
        a->mailbox(j).construct();

    a->prefix().stack_size = stack_size;
    size_t k;
    // Mark each internal sync element for the tools
    for( k=0; k<number_of_workers; ++k ) {
        ITT_SYNC_CREATE(a->slot + k, SyncType_Scheduler, SyncObj_WorkerTaskPool);
        ITT_SYNC_CREATE(&w[k].scheduler, SyncType_Scheduler, SyncObj_WorkerLifeCycleMgmt);
        ITT_SYNC_CREATE(&a->mailbox(k+1), SyncType_Scheduler, SyncObj_Mailbox);
    }
    for( ; k<number_of_slots; ++k ) {
        ITT_SYNC_CREATE(a->slot + k, SyncType_Scheduler, SyncObj_MasterTaskPool);
        ITT_SYNC_CREATE(&a->mailbox(k+1), SyncType_Scheduler, SyncObj_Mailbox);
    }

    return a;
}

void arena::free_arena () {
    // Drain mailboxes
    // TODO: each scheduler should plug-and-drain its own mailbox when it terminates.
    intptr_t drain_count = 0;
    for( unsigned i=1; i<=prefix().number_of_slots; ++i )
        drain_count += mailbox(i).drain();
#if __TBB_COUNT_TASK_NODES
    prefix().task_node_count -= drain_count;
    if( prefix().task_node_count ) {
        runtime_warning( "Leaked %ld task objects\n", long(prefix().task_node_count) );
    }
#endif /* __TBB_COUNT_TASK_NODES */
    void* storage  = &mailbox(prefix().number_of_slots);
    delete[] prefix().worker_list;
    prefix().~ArenaPrefix();
    NFS_Free( storage );
}

void arena::close_arena () {
    for(;;) {
        pool_state_t snapshot = prefix().pool_state;
        if( snapshot==SNAPSHOT_SERVER_GOING_AWAY ) 
            break;
        if( prefix().pool_state.compare_and_swap( SNAPSHOT_SERVER_GOING_AWAY, snapshot )==snapshot ) {
            if( snapshot!=SNAPSHOT_EMPTY )
                prefix().server->adjust_job_count_estimate( -int(prefix().number_of_workers) );
            break;
        }
    }
    prefix().server->request_close_connection();
}

#endif /* !__TBB_ARENA_PER_MASTER */

bool arena::is_out_of_work() {
    // TODO: rework it to return at least a hint about where a task was found; better if the task itself.
    for(;;) {
        pool_state_t snapshot = prefix().pool_state;
        switch( snapshot ) {
            case SNAPSHOT_EMPTY:
#if !__TBB_ARENA_PER_MASTER
            case SNAPSHOT_SERVER_GOING_AWAY:
#endif /* !__TBB_ARENA_PER_MASTER */
                return true;
            case SNAPSHOT_FULL: {
                // Use unique id for "busy" in order to avoid ABA problems.
                const pool_state_t busy = pool_state_t(this);
                // Request permission to take snapshot
                if( prefix().pool_state.compare_and_swap( busy, SNAPSHOT_FULL )==SNAPSHOT_FULL ) {
                    // Got permission.  Take the snapshot.
#if __TBB_ARENA_PER_MASTER
                    size_t n = my_limit;
#else /* !__TBB_ARENA_PER_MASTER */
                    size_t n = prefix().limit;
#endif /* !__TBB_ARENA_PER_MASTER */
                    size_t k; 
                    for( k=0; k<n; ++k ) 
                        if( slot[k].task_pool != EmptyTaskPool && slot[k].head < slot[k].tail )
                            break;
                    bool work_absent = k>=n;
#if __TBB_ARENA_PER_MASTER
                    work_absent = work_absent && my_task_stream.empty();
#endif /* __TBB_ARENA_PER_MASTER */
                    // Test and test-and-set.
                    if( prefix().pool_state==busy ) {
                        if( work_absent ) {
#if __TBB_ARENA_PER_MASTER
                            // save current demand value before setting SNAPSHOT_EMPTY,
                            // to avoid race with advertise_new_work.
                            int current_demand = (int)my_max_num_workers;
#endif
                            if( prefix().pool_state.compare_and_swap( SNAPSHOT_EMPTY, busy )==busy ) {
                                // This thread transitioned pool to empty state, and thus is responsible for
                                // telling RML that there is no other work to do.
#if __TBB_ARENA_PER_MASTER
                                my_market->adjust_demand( *this, -current_demand );
#else /* !__TBB_ARENA_PER_MASTER */
                                prefix().server->adjust_job_count_estimate( -int(prefix().number_of_workers) );
#endif /* !__TBB_ARENA_PER_MASTER */
                                return true;
                            }
                        } else {
                            // Undo previous transition SNAPSHOT_FULL-->busy, unless another thread undid it.
                            prefix().pool_state.compare_and_swap( SNAPSHOT_FULL, busy );
                        }
                    }
                } 
                return false;
            }
            default:
                // Another thread is taking a snapshot.
                return false;
        }
    }
}

#if __TBB_COUNT_TASK_NODES 
intptr_t arena::workers_task_node_count() {
    intptr_t result = 0;
#if __TBB_ARENA_PER_MASTER
    for( unsigned i = 1; i < my_num_slots; ++i ) {
        generic_scheduler* s = slot[i].my_scheduler;
#else /* !__TBB_ARENA_PER_MASTER */
    for( unsigned i=0; i<prefix().number_of_workers; ++i ) {
        generic_scheduler* s = prefix().worker_list[i].scheduler;
#endif /* !__TBB_ARENA_PER_MASTER */
        if( s )
            result += s->task_node_count;
    }
    return result;
}
#endif /* __TBB_COUNT_TASK_NODES */

} // namespace internal
} // namespace tbb
