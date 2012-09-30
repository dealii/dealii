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

#include "tbb/tbb_stddef.h"

#if __TBB_ARENA_PER_MASTER

#include "market.h"
#include "tbb_main.h"
#include "governor.h"
#include "scheduler.h"
#include "itt_notify.h"

namespace tbb {
namespace internal {

//------------------------------------------------------------------------
// market
//------------------------------------------------------------------------

market::market ( unsigned max_num_workers, size_t stack_size )
    : my_ref_count(1)
    , my_stack_size(stack_size)
    , my_max_num_workers(max_num_workers)
{
    my_next_arena = my_arenas.begin();

    // Once created RML server will start initializing workers that will need 
    // global market instance to get worker stack size
    my_server = governor::create_rml_server( *this );
    __TBB_ASSERT( my_server, "Failed to create RML server" );
}


market& market::global_market ( unsigned max_num_workers, size_t stack_size ) {
    global_market_mutex_type::scoped_lock lock( theMarketMutex );
    market *m = theMarket;
    if ( m ) {
        ++m->my_ref_count;
        if ( m->my_stack_size < stack_size )
            runtime_warning( "Newer master request for larger stack cannot be satisfied\n" );
    }
    else {
        max_num_workers = max( governor::default_num_threads() - 1, max_num_workers );
        // at least 1 worker is required to support starvation resistant tasks
        if( max_num_workers==0 ) max_num_workers = 1;
        // Create the global market instance
        size_t size = sizeof(market);
#if __TBB_TASK_GROUP_CONTEXT
        __TBB_ASSERT( __TBB_offsetof(market, my_workers) + sizeof(generic_scheduler*) == sizeof(market),
                      "my_workers must be the last data field of the market class");
        size += sizeof(generic_scheduler*) * (max_num_workers - 1);
#endif /* __TBB_TASK_GROUP_CONTEXT */
        __TBB_InitOnce::add_ref();
        void* storage = NFS_Allocate(size, 1, NULL);
        memset( storage, 0, size );
        // Initialize and publish global market
        m = new (storage) market( max_num_workers, stack_size );
        theMarket = m;
    }
    return *m;
}

void market::destroy () {
#if __TBB_COUNT_TASK_NODES
    if ( my_task_node_count )
        runtime_warning( "Leaked %ld task objects\n", (intptr_t)my_task_node_count );
#endif /* __TBB_COUNT_TASK_NODES */
    this->~market();
    NFS_Free( this );
    __TBB_InitOnce::remove_ref();
}

void market::release () {
    __TBB_ASSERT( theMarket == this, "Global market instance was destroyed prematurely?" );
    bool do_release = false;
    {
        global_market_mutex_type::scoped_lock lock(theMarketMutex);
        if ( --my_ref_count == 0 ) {
            do_release = true;
            theMarket = NULL;
        }
    }
    if( do_release )
        my_server->request_close_connection();
}

arena& market::create_arena ( unsigned max_num_workers, size_t stack_size ) {
    market &m = global_market( max_num_workers, stack_size ); // increases market's ref count
    arena& a = arena::allocate_arena( m, min(max_num_workers, m.my_max_num_workers) );
    // Add newly created arena into the existing market's list.
    spin_mutex::scoped_lock lock(m.my_arenas_list_mutex);
    m.my_arenas.push_front( a );
    if ( m.my_arenas.size() == 1 )
        m.my_next_arena = m.my_arenas.begin();
    return a;
}

void market::detach_arena ( arena& a ) {
    __TBB_ASSERT( theMarket == this, "Global market instance was destroyed prematurely?" );
    spin_mutex::scoped_lock lock(my_arenas_list_mutex);
    __TBB_ASSERT( my_next_arena != my_arenas.end(), NULL );
    if ( &*my_next_arena == &a )
        if ( ++my_next_arena == my_arenas.end() && my_arenas.size() > 1 )
            my_next_arena = my_arenas.begin();
    my_arenas.remove( a );
}

arena* market::arena_in_need () {
    spin_mutex::scoped_lock lock(my_arenas_list_mutex);
    if ( my_arenas.empty() )
        return NULL;
    __TBB_ASSERT( my_next_arena != my_arenas.end(), NULL );
    arena_list_type::iterator it = my_next_arena;
    do {
        arena& a = *it;
        if ( ++it == my_arenas.end() )
            it = my_arenas.begin();
        if ( a.num_workers_active() < a.my_num_workers_allotted ) {
            ++a.my_num_threads_active;
            my_next_arena = it;
            return &a;
        }
    } while ( it != my_next_arena );
    return NULL;
}

void market::update_allotment ( int max_workers ) {
    unsigned carry = 0;
    spin_mutex::scoped_lock lock(my_arenas_list_mutex);
    arena_list_type::iterator it = my_arenas.begin();
    int total_demand = my_total_demand;
    max_workers = min(max_workers, total_demand);
    if ( total_demand > 0 ) {
        for ( ; it != my_arenas.end(); ++it ) {
            arena& a = *it;
            int tmp = a.my_num_workers_requested * max_workers + carry;
            int allotted = tmp / total_demand;
            carry = tmp % total_demand;
            a.my_num_workers_allotted = min( allotted, (int)a.my_max_num_workers );
        }
    }
    else {
        for ( ; it != my_arenas.end(); ++it ) {
            it->my_num_workers_allotted = 0;
        }
    }
}

/** The balancing algorithm may be liable to data races. However the aberrations 
    caused by the races are not fatal and generally only temporarily affect fairness 
    of the workers distribution among arenas. **/
void market::adjust_demand ( arena& a, int delta ) {
    __TBB_ASSERT( theMarket, "market instance was destroyed prematurely?" );
    a.my_num_workers_requested += delta;
    my_total_demand += delta;
    update_allotment( my_max_num_workers );
    // Must be called outside of any locks
    my_server->adjust_job_count_estimate( delta );
    GATHER_STATISTIC( governor::local_scheduler_if_initialized() ? ++governor::local_scheduler_if_initialized()->my_counters.gate_switches : 0 );
}

void market::process( job& j ) {
    generic_scheduler& s = static_cast<generic_scheduler&>(j);
    while ( arena *a = arena_in_need() )
        a->process(s);
    GATHER_STATISTIC( ++s.my_counters.market_roundtrips );
}

void market::cleanup( job& j ) {
    __TBB_ASSERT( theMarket != this, NULL );
    generic_scheduler& s = static_cast<generic_scheduler&>(j);
    generic_scheduler* mine = governor::local_scheduler_if_initialized();
    __TBB_ASSERT( !mine || mine->arena_index!=0, NULL );
    if( mine!=&s ) {
        governor::assume_scheduler( &s );
        generic_scheduler::cleanup_worker( &s, mine!=NULL );
        governor::assume_scheduler( mine );
    } else {
        generic_scheduler::cleanup_worker( &s, true );
    }
}

void market::acknowledge_close_connection() {
    destroy();
}

::rml::job* market::create_one_job() {
    unsigned index = ++my_num_workers;
    __TBB_ASSERT( index > 0, NULL );
    ITT_THREAD_SET_NAME(_T("TBB Worker Thread"));
    // index serves as a hint decreasing conflicts between workers when they migrate between arenas
    generic_scheduler* s = generic_scheduler::create_worker( *this, index );
#if __TBB_TASK_GROUP_CONTEXT
    __TBB_ASSERT( !my_workers[index - 1], NULL );
    my_workers[index - 1] = s;
#endif /* __TBB_TASK_GROUP_CONTEXT */
    governor::sign_on(s);
    return s;
}

#if __TBB_TASK_GROUP_CONTEXT
/** Propagates cancellation down the tree of dependent contexts by walking each 
    thread's local list of contexts **/
void market::propagate_cancellation ( task_group_context& ctx ) {
    __TBB_ASSERT ( ctx.my_cancellation_requested, "No cancellation request in the context" );
    // The whole propagation algorithm is under the lock in order to ensure correctness 
    // in case of parallel cancellations at the different levels of the context tree.
    // See the note 1 at the bottom of this file.
    global_market_mutex_type::scoped_lock lock(theMarketMutex);
    // Advance global cancellation epoch
    __TBB_FetchAndAddWrelease(&global_cancel_count, 1);
    // Propagate to all workers and masters and sync up their local epochs with the global one
    unsigned num_workers = my_num_workers;
    for ( unsigned i = 0; i < num_workers; ++i ) {
        generic_scheduler *s = my_workers[i];
        // If the worker is only about to be registered, skip it.
        if ( s )
            s->propagate_cancellation();
    }
    arena_list_type::iterator it = my_arenas.begin();
    for ( ; it != my_arenas.end(); ++it ) {
        generic_scheduler *s = it->slot[0].my_scheduler;
        // If the master is under construction, skip it.
        if ( s )
            s->propagate_cancellation();
    }
}
#endif /* __TBB_TASK_GROUP_CONTEXT */

#if __TBB_COUNT_TASK_NODES 
intptr_t market::workers_task_node_count() {
    intptr_t result = 0;
    spin_mutex::scoped_lock lock(my_arenas_list_mutex);
    for ( arena_list_type::iterator it = my_arenas.begin(); it != my_arenas.end(); ++it )
        result += it->workers_task_node_count();
    return result;
}
#endif /* __TBB_COUNT_TASK_NODES */

} // namespace internal
} // namespace tbb

#endif /* __TBB_ARENA_PER_MASTER */

/*
    Notes:

1.  Consider parallel cancellations at the different levels of the context tree:

        Ctx1 <- Cancelled by Thread1            |- Thread2 started processing
         |                                      |
        Ctx2                                    |- Thread1 started processing
         |                                   T1 |- Thread2 finishes and syncs up local counters
        Ctx3 <- Cancelled by Thread2            |
         |                                      |- Ctx5 is bound to Ctx2
        Ctx4                                    |
                                             T2 |- Thread1 reaches Ctx2
                                             
    Thread-propagator of each cancellation increments global counter. However the thread 
    propagating the cancellation from the outermost context (Thread1) may be the last 
    to finish. Which means that the local counters may be synchronized earlier (by Thread2, 
    at Time1) than it propagated cancellation into Ctx2 (at time Time2). If a new context 
    (Ctx5) is created and bound to Ctx2 between Time1 and Time2, checking its parent only 
    (Ctx2) may result in cancellation request being lost.

    This issue is solved by doing the whole propagation under the lock (the_scheduler_list_mutex).

    If we need more concurrency while processing parallel cancellations, we could try 
    the following modification of the propagation algorithm:

    advance global counter and remember it
    for each thread:
        scan thread's list of contexts
    for each thread:
        sync up its local counter only if the global counter has not been changed

    However this version of the algorithm requires more analysis and verification.
*/
