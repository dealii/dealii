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

#ifndef _TBB_governor_H
#define _TBB_governor_H

#include "tbb/task_scheduler_init.h"
#if !__TBB_ARENA_PER_MASTER
#include "tbb/mutex.h"
#endif /* !__TBB_ARENA_PER_MASTER */
#include "../rml/include/rml_tbb.h"

#include "tbb_misc.h" // for DetectNumberOfWorkers and ThreadStackSize
#include "tls.h"

#if __TBB_SURVIVE_THREAD_SWITCH
#include "cilk-tbb-interop.h"
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

namespace tbb {
namespace internal {

#if __TBB_ARENA_PER_MASTER
class market;
#else /* !__TBB_ARENA_PER_MASTER */
class arena;
#endif /* !__TBB_ARENA_PER_MASTER */
class generic_scheduler;
class __TBB_InitOnce;

//------------------------------------------------------------------------
// Class governor
//------------------------------------------------------------------------

#if __TBB_ARENA_PER_MASTER
//! The class handles access to the single instance of market, and to TLS to keep scheduler instances.
#else /* !__TBB_ARENA_PER_MASTER */
//! The class handles access to the single instance of arena, and to TLS to keep scheduler instances.
#endif /* !__TBB_ARENA_PER_MASTER */
/** It also supports automatic on-demand initialization of the TBB scheduler.
    The class contains only static data members and methods.*/
class governor {
    friend class __TBB_InitOnce;
#if __TBB_ARENA_PER_MASTER
    friend class market;
#else /* !__TBB_ARENA_PER_MASTER */
    friend void ITT_DoUnsafeOneTimeInitialization ();
#endif /* __TBB_ARENA_PER_MASTER */

    //! TLS for scheduler instances associated with individual threads
    static basic_tls<generic_scheduler*> theTLS;

#if !__TBB_ARENA_PER_MASTER
    //! Currently active arena
    static arena* theArena;

    //! Mutex guarding creation/destruction of theArena
    static mutex  theArenaMutex;

    //! Caches the number of workers in the currently active arena
    static unsigned NumWorkers;
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Caches the maximal level of paralellism supported by the hardware 
    static unsigned DefaultNumberOfThreads;
    
    static rml::tbb_factory theRMLServerFactory;

    static bool UsePrivateRML;

    //! Create key for thread-local storage and initialize RML.
    static void acquire_resources ();

    //! Destroy the thread-local storage key and deinitialize RML.
    static void release_resources ();

    static rml::tbb_server* create_rml_server ( rml::tbb_client& );

#if !__TBB_ARENA_PER_MASTER
    //! Obtain the instance of arena to register a new master thread
    /** If there is no active arena, create one. */
    static arena* obtain_arena( int number_of_threads, stack_size_type thread_stack_size );
#endif /* !__TBB_ARENA_PER_MASTER */

    //! The internal routine to undo automatic initialization.
    /** The signature is written with void* so that the routine
        can be the destructor argument to pthread_key_create. */
    static void auto_terminate(void* scheduler);

public:
    static unsigned default_num_threads () {
        // No memory fence required, because at worst each invoking thread calls DetectNumberOfWorkers once.
        return DefaultNumberOfThreads ? DefaultNumberOfThreads : 
                                        DefaultNumberOfThreads = DetectNumberOfWorkers();
    }
    //! Processes scheduler initialization request (possibly nested) in a master thread
    /** If necessary creates new instance of arena and/or local scheduler.
        The auto_init argument specifies if the call is due to automatic initialization. **/
    static generic_scheduler* init_scheduler( unsigned num_threads, stack_size_type stack_size, bool auto_init = false );

    //! Processes scheduler termination request (possibly nested) in a master thread
    static void terminate_scheduler( generic_scheduler* s );

#if __TBB_ARENA_PER_MASTER
    //! Returns number of worker threads in the currently active arena.
    inline static unsigned max_number_of_workers ();

#else /* !__TBB_ARENA_PER_MASTER */
    //! Dereference arena when a master thread stops using TBB.
    /** If no more masters in the arena, terminate workers and destroy it. */
    static void finish_with_arena();

    static unsigned max_number_of_workers() {
        __TBB_ASSERT( theArena, "thread did not activate a task_scheduler_init object?" );
        return NumWorkers;
    }
#endif /* !__TBB_ARENA_PER_MASTER */

    //! Register TBB scheduler instance in thread local storage.
    static void sign_on(generic_scheduler* s);

    //! Unregister TBB scheduler instance from thread local storage.
    static void sign_off(generic_scheduler* s);

    //! Used to check validity of the local scheduler TLS contents.
    static bool is_set ( generic_scheduler* s ) { return theTLS.get() == s; }

    //! Temporarily set TLS slot to the given scheduler
    static void assume_scheduler( generic_scheduler* s ) { 
#if !__TBB_ARENA_PER_MASTER
        // should be called by a Master
        __TBB_ASSERT( !s || !theTLS.get(), "should be called by master" );
#endif
        theTLS.set( s ); 
    }

    //! Obtain the thread local instance of the TBB scheduler.
    /** If the scheduler has not been initialized yet, initialization is done automatically.
        Note that auto-initialized scheduler instance is destroyed only when its thread terminates. **/
    static generic_scheduler* local_scheduler () {
        generic_scheduler* s = theTLS.get();
        return s ? s : init_scheduler( (unsigned)task_scheduler_init::automatic, 0, true );
    }

    static generic_scheduler* local_scheduler_if_initialized () {
        return theTLS.get();
    }

    //! Undo automatic initialization if necessary; call when a thread exits.
    static void terminate_auto_initialized_scheduler() {
        auto_terminate( theTLS.get() );
    }

    static void print_version_info ();

#if __TBB_SURVIVE_THREAD_SWITCH
    static __cilk_tbb_retcode stack_op_handler( __cilk_tbb_stack_op op, void* );
#endif /* __TBB_SURVIVE_THREAD_SWITCH */
}; // class governor

} // namespace internal
} // namespace tbb

#if __TBB_ARENA_PER_MASTER
#include "scheduler.h"

inline unsigned tbb::internal::governor::max_number_of_workers () {
    return local_scheduler()->number_of_workers_in_my_arena();
}
#endif /* __TBB_ARENA_PER_MASTER */

#endif /* _TBB_governor_H */
