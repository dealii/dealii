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

#include "governor.h"
#include "tbb_main.h"
#include "scheduler.h"
#if __TBB_ARENA_PER_MASTER
#include "market.h"
#endif /* __TBB_ARENA_PER_MASTER */
#include "arena.h"

#include "tbb/task_scheduler_init.h"

#if __TBB_SURVIVE_THREAD_SWITCH
#include "dynamic_link.h"
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

namespace tbb {
namespace internal {

//------------------------------------------------------------------------
// governor
//------------------------------------------------------------------------

#if __TBB_SURVIVE_THREAD_SWITCH

#if _WIN32
#define CILKLIB_NAME "cilkrts20.dll"
#else
#define CILKLIB_NAME "libcilkrts.so"
#endif

//! Handler for memory allocation
static __cilk_tbb_retcode (*watch_stack_handler)(struct __cilk_tbb_unwatch_thunk* u,
                                                 struct __cilk_tbb_stack_op_thunk o);

//! Table describing the how to link the handlers.
static const dynamic_link_descriptor CilkLinkTable[] = {
    DLD(__cilkrts_watch_stack, watch_stack_handler)
};

void initialize_survive_thread_switch() {
    dynamic_link( CILKLIB_NAME, CilkLinkTable, 1 );
}
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

namespace rml {
    tbb_server* make_private_server( tbb_client& client );
}

void governor::acquire_resources () {
#if USE_PTHREAD
    int status = theTLS.create(auto_terminate);
#else
    int status = theTLS.create();
#endif
    if( status )
        handle_perror(status, "TBB failed to initialize TLS storage\n");

    ::rml::factory::status_type res = theRMLServerFactory.open(); 
    UsePrivateRML = res != ::rml::factory::st_success;
}

void governor::release_resources () {
    theRMLServerFactory.close();
#if TBB_USE_ASSERT
    if( __TBB_InitOnce::initialization_done() && theTLS.get() ) 
        runtime_warning( "TBB is unloaded while tbb::task_scheduler_init object is alive?" );
#endif
    int status = theTLS.destroy();
    if( status )
        handle_perror(status, "TBB failed to destroy TLS storage");
}

rml::tbb_server* governor::create_rml_server ( rml::tbb_client& client ) {
    rml::tbb_server* server = NULL;
    if( !UsePrivateRML ) {
        ::rml::factory::status_type status = theRMLServerFactory.make_server( server, client );
        if( status != ::rml::factory::st_success ) {
            UsePrivateRML = true;
            runtime_warning( "rml::tbb_factorymake_server failed with status %x, falling back on private rml", status );
        }
    }
    if ( !server ) {
        __TBB_ASSERT( UsePrivateRML, NULL );
        server = rml::make_private_server( client );
    }
    __TBB_ASSERT( server, "Failed to create RML server" );
    return server;
}

#if !__TBB_ARENA_PER_MASTER

arena* governor::obtain_arena( int number_of_threads, stack_size_type thread_stack_size ) {
    mutex::scoped_lock lock( theArenaMutex );
    arena* arena = theArena;
    if( arena ) {
        arena->prefix().number_of_masters += 1;
    } else {
        __TBB_ASSERT( number_of_threads > 0, NULL );
        arena = arena::allocate_arena( 2*number_of_threads, number_of_threads-1,
                                   thread_stack_size ? thread_stack_size : ThreadStackSize );
        __TBB_ASSERT( arena->prefix().number_of_masters==1, NULL );
        NumWorkers = arena->prefix().number_of_workers;

        arena->prefix().server = create_rml_server( arena->prefix() );

        // Publish the arena.  
        // A memory release fence is not required here, because workers have not started yet,
        // and concurrent masters inspect theArena while holding theArenaMutex.
        __TBB_ASSERT( !theArena, NULL );
        theArena = arena;
    }
    return arena;
}

void governor::finish_with_arena() {
    mutex::scoped_lock lock( theArenaMutex );
    arena* a = theArena;
    __TBB_ASSERT( a, "theArena is missing" );
    if( --(a->prefix().number_of_masters) )
        a = NULL;
    else {
        theArena = NULL;
        // Must do this while holding lock, otherwise terminate message might reach
        // RML thread *after* initialize message reaches it for the next arena,
        // which causes TLS to be set to new value before old one is erased!
        a->close_arena();
    }
}
#endif /* !__TBB_ARENA_PER_MASTER */

void governor::sign_on(generic_scheduler* s) {
    __TBB_ASSERT( !s->is_registered, NULL );  
    s->is_registered = true;
#if !__TBB_ARENA_PER_MASTER
    __TBB_InitOnce::add_ref();
#endif /* !__TBB_ARENA_PER_MASTER */
    theTLS.set(s);
#if __TBB_SURVIVE_THREAD_SWITCH
    __cilk_tbb_stack_op_thunk o;
    o.routine = &stack_op_handler;
    o.data = s;
    if( watch_stack_handler ) {
        if( (*watch_stack_handler)(&s->my_cilk_unwatch_thunk, o) ) {
            // Failed to register with Cilk, make sure we are clean
            s->my_cilk_unwatch_thunk.routine = NULL;
        }
#if TBB_USE_ASSERT
        else
            s->my_cilk_state = generic_scheduler::cs_running;
#endif /* TBB_USE_ASSERT */
    }
#endif /* __TBB_SURVIVE_THREAD_SWITCH */
}

void governor::sign_off(generic_scheduler* s) {
    if( s->is_registered ) {
        __TBB_ASSERT( theTLS.get()==s || (!s->is_worker() && !theTLS.get()), "attempt to unregister a wrong scheduler instance" );
        theTLS.set(NULL);
        s->is_registered = false;
#if __TBB_SURVIVE_THREAD_SWITCH
        __cilk_tbb_unwatch_thunk &ut = s->my_cilk_unwatch_thunk;
        if ( ut.routine )
           (*ut.routine)(ut.data);
#endif /* __TBB_SURVIVE_THREAD_SWITCH */
#if !__TBB_ARENA_PER_MASTER
        __TBB_InitOnce::remove_ref();
#endif /* !__TBB_ARENA_PER_MASTER */
    }
}

generic_scheduler* governor::init_scheduler( unsigned num_threads, stack_size_type stack_size, bool auto_init ) {
    if( !__TBB_InitOnce::initialization_done() )
        DoOneTimeInitializations();
    generic_scheduler* s = theTLS.get();
    if( s ) {
        s->ref_count += 1;
        return s;
    }
    if( (int)num_threads == task_scheduler_init::automatic )
        num_threads = default_num_threads();
#if __TBB_ARENA_PER_MASTER
    s = generic_scheduler::create_master( 
            market::create_arena( num_threads - 1, stack_size ? stack_size : ThreadStackSize ) );
#else /* !__TBB_ARENA_PER_MASTER */
    s = generic_scheduler::create_master( *obtain_arena(num_threads, stack_size) );
#endif /* !__TBB_ARENA_PER_MASTER */
    __TBB_ASSERT(s, "Somehow a local scheduler creation for a master thread failed");
    s->is_auto_initialized = auto_init;
    return s;
}

void governor::terminate_scheduler( generic_scheduler* s ) {
    __TBB_ASSERT( s == theTLS.get(), "Attempt to terminate non-local scheduler instance" );
    if( !--(s->ref_count) )
        s->cleanup_master();
}

void governor::auto_terminate(void* arg){
    generic_scheduler* s = static_cast<generic_scheduler*>(arg);
    if( s && s->is_auto_initialized ) {
        if( !--(s->ref_count) ) {
            if ( !theTLS.get() && !s->local_task_pool_empty() ) {
                // This thread's TLS slot is already cleared. But in order to execute
                // remaining tasks cleanup_master() will need TLS correctly set.
                // So we temporarily restore its value.
                theTLS.set(s);
                s->cleanup_master();
                theTLS.set(NULL);
            }
            else
                s->cleanup_master();
        }
    }
}

void governor::print_version_info () {
    if ( UsePrivateRML )
        PrintExtraVersionInfo( "RML", "private" );
    else {
        PrintExtraVersionInfo( "RML", "shared" );
        theRMLServerFactory.call_with_server_info( PrintRMLVersionInfo, (void*)"" );
    }
#if __TBB_SURVIVE_THREAD_SWITCH
    if( watch_stack_handler )
        PrintExtraVersionInfo( "CILK", CILKLIB_NAME );
#endif /* __TBB_SURVIVE_THREAD_SWITCH */
}

#if __TBB_SURVIVE_THREAD_SWITCH
__cilk_tbb_retcode governor::stack_op_handler( __cilk_tbb_stack_op op, void* data ) {
    __TBB_ASSERT(data,NULL);
    generic_scheduler* s = static_cast<generic_scheduler*>(data);
#if TBB_USE_ASSERT
    void* current = theTLS.get();
#if _WIN32||_WIN64
    unsigned thread_id = GetCurrentThreadId();
#else
    unsigned thread_id = unsigned(pthread_self());
#endif

#endif /* TBB_USE_ASSERT */
    switch( op ) {
        default:
            __TBB_ASSERT( 0, "invalid op" );
        case CILK_TBB_STACK_ADOPT: {
            __TBB_ASSERT( !current && s->my_cilk_state==generic_scheduler::cs_limbo || 
                          current==s && s->my_cilk_state==generic_scheduler::cs_running, "invalid adoption" );
#if TBB_USE_ASSERT
            if( current==s ) 
                runtime_warning( "redundant adoption of %p by thread %x\n", s, thread_id );
            s->my_cilk_state = generic_scheduler::cs_running;
#endif /* TBB_USE_ASSERT */
            theTLS.set(s);
            break;
        }
        case CILK_TBB_STACK_ORPHAN: {
            __TBB_ASSERT( current==s && s->my_cilk_state==generic_scheduler::cs_running, "invalid orphaning" ); 
#if TBB_USE_ASSERT
            s->my_cilk_state = generic_scheduler::cs_limbo;
#endif /* TBB_USE_ASSERT */
            theTLS.set(NULL);
            break;
        }
        case CILK_TBB_STACK_RELEASE: {
            __TBB_ASSERT( !current && s->my_cilk_state==generic_scheduler::cs_limbo || 
                          current==s && s->my_cilk_state==generic_scheduler::cs_running, "invalid release" );
#if TBB_USE_ASSERT
            s->my_cilk_state = generic_scheduler::cs_freed;
#endif /* TBB_USE_ASSERT */
            s->my_cilk_unwatch_thunk.routine = NULL;
            auto_terminate( s );
        } 
    }
    return 0;
}
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

} // namespace internal

//------------------------------------------------------------------------
// task_scheduler_init
//------------------------------------------------------------------------

using namespace internal;

/** Left out-of-line for the sake of the backward binary compatibility **/
void task_scheduler_init::initialize( int number_of_threads ) {
    initialize( number_of_threads, 0 );
}

void task_scheduler_init::initialize( int number_of_threads, stack_size_type thread_stack_size ) {
    if( number_of_threads!=deferred ) {
        __TBB_ASSERT( !my_scheduler, "task_scheduler_init already initialized" );
        __TBB_ASSERT( number_of_threads==-1 || number_of_threads>=1,
                    "number_of_threads for task_scheduler_init must be -1 or positive" );
        my_scheduler = governor::init_scheduler( number_of_threads, thread_stack_size );
    } else {
        __TBB_ASSERT( !thread_stack_size, "deferred initialization ignores stack size setting" );
    }
}

void task_scheduler_init::terminate() {
    generic_scheduler* s = static_cast<generic_scheduler*>(my_scheduler);
    my_scheduler = NULL;
    __TBB_ASSERT( s, "task_scheduler_init::terminate without corresponding task_scheduler_init::initialize()");
    governor::terminate_scheduler(s);
}

int task_scheduler_init::default_num_threads() {
    return governor::default_num_threads();
}

} // namespace tbb
