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

#ifndef _TBB_market_H
#define _TBB_market_H

#include "tbb/tbb_stddef.h"

#if __TBB_ARENA_PER_MASTER

#include "tbb/atomic.h"
#include "tbb/spin_mutex.h"
#include "../rml/include/rml_tbb.h"

#include "intrusive_list.h"

#if defined(_MSC_VER) && defined(_Wp64)
    // Workaround for overzealous compiler warnings in /Wp64 mode
    #pragma warning (push)
    #pragma warning (disable: 4244)
#endif

namespace tbb {

class task_group_context;

namespace internal {

class arena;
class generic_scheduler;

//------------------------------------------------------------------------
// Class market
//------------------------------------------------------------------------

class market : no_copy, rml::tbb_client {
    friend void ITT_DoUnsafeOneTimeInitialization ();

    typedef intrusive_list<arena> arena_list_type;

    //! Currently active global market
    static market* theMarket;

    typedef spin_mutex global_market_mutex_type;

    //! Mutex guarding creation/destruction of theMarket, insertions/deletions in my_arenas, and cancellation propagation
    static global_market_mutex_type  theMarketMutex;

    //! Reference count controlling market object lifetime
    intptr_t my_ref_count;

    //! List of active arenas
    arena_list_type my_arenas;

    //! The first arena to be checked when idle worker seeks for an arena to enter
    /** The check happens in round-robin fashion. **/
    arena_list_type::iterator my_next_arena;

    //! Lightweight mutex guarding accounting operations with arenas list
    spin_mutex  my_arenas_list_mutex;

    //! Number of workers that were requested by all arenas
    atomic<int> my_total_demand;

    //! Pointer to the RML server object that services this TBB instance.
    rml::tbb_server* my_server;

    //! Stack size of worker threads
    size_t my_stack_size;

    //! Number of workers requested from the underlying resource manager
    unsigned my_max_num_workers;

#if __TBB_COUNT_TASK_NODES
    //! Net number of nodes that have been allocated from heap.
    /** Updated each time a scheduler or arena is destroyed. */
    atomic<intptr_t> my_task_node_count;
#endif /* __TBB_COUNT_TASK_NODES */

    //! Number of workers that have been delivered by RML
    atomic<unsigned> my_num_workers;

    //! Constructor
    market ( unsigned max_num_workers, size_t stack_size );

    //! Factory method creating new market object
    static market& global_market ( unsigned max_num_workers, size_t stack_size );

    //! Destroys and deallocates market object created by market::create()
    void destroy ();

    //! Returns next arena that needs more workers, or NULL.
    arena* arena_in_need ();

    //! Recalculates the number of workers assigned to each arena.
    /** The actual number of workers servicing a particular arena may temporarily 
        deviate from the calculated value. **/
    void update_allotment ( int max_workers );

    //! Returns number of masters doing computational (CPU-intensive) work
    int num_active_masters () { return 1; }  // APM TODO: replace with a real mechanism

    // // //
    // Implementation of rml::tbb_client interface methods

    /*override*/ version_type version () const { return 0; }

    /*override*/ unsigned max_job_count () const { return my_max_num_workers; }

    /*override*/ size_t min_stack_size () const { return worker_stack_size(); }

    /*override*/ policy_type policy () const { return throughput; }

    /*override*/ job* create_one_job ();

    /*override*/ void cleanup( job& j );

    /*override*/ void acknowledge_close_connection ();

    /*override*/ void process( job& j );

public:
    //! Creates an arena object
    /** If necessary, also creates global market instance, and boosts its ref count.
        Each call to create_arena() must be matched by the call to arena::free_arena(). **/
    static arena& create_arena ( unsigned max_num_workers, size_t stack_size );

    //! Removes the arena from the market's list
    void detach_arena ( arena& );

    //! Decrements market's refcount and destroys it in the end
    void release ();

    //! Request that arena's need in workers should be adjusted.
    /** Concurrent invocations are possible only on behalf of different arenas. **/
    void adjust_demand ( arena&, int delta );

    //! Returns the requested stack size of worker threads.
    size_t worker_stack_size () const { return my_stack_size; }

#if __TBB_COUNT_TASK_NODES
    //! Returns the number of task objects "living" in worker threads
    intptr_t workers_task_node_count();

    //! Net number of nodes that have been allocated from heap.
    /** Updated each time a scheduler or arena is destroyed. */
    void update_task_node_count( intptr_t delta ) { my_task_node_count += delta; }
#endif /* __TBB_COUNT_TASK_NODES */

#if __TBB_TASK_GROUP_CONTEXT
    //! Propagates cancellation request to all descendants of the context.
    void propagate_cancellation ( task_group_context& ctx );

    //! Array of pointers to the registered workers
    /** Used by cancellation propagation mechanism.
        Must be the last data member of the class market. **/
    generic_scheduler* my_workers[1];
#endif /* __TBB_TASK_GROUP_CONTEXT */
#if __TBB_ARENA_PER_MASTER && ( _WIN32||_WIN64 )
    //! register master with the resource manager
    void register_master( ::rml::server::execution_resource_t& rsc_handle ) {
        __TBB_ASSERT( my_server, "RML server not defined?" );
        // the server may ignore registration and set master_exec_resource to NULL.
        my_server->register_master( rsc_handle );
    }

    //! unregister master with the resource manager
    void unregister_master( ::rml::server::execution_resource_t& rsc_handle ) const {
        my_server->unregister_master( rsc_handle );
    }
#endif /* !__TBB_ARENA_PER_MASTER && ( _WIN32||_WIN64 ) */

}; // class market

} // namespace internal
} // namespace tbb

#if defined(_MSC_VER) && defined(_Wp64)
    // Workaround for overzealous compiler warnings in /Wp64 mode
    #pragma warning (pop)
#endif // warning 4244 is back

#endif /* __TBB_ARENA_PER_MASTER */

#endif /* _TBB_market_H */
