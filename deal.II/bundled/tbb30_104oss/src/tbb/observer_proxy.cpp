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

#include "tbb/tbb_config.h"

#if __TBB_SCHEDULER_OBSERVER

#include "tbb/spin_rw_mutex.h"
#include "tbb/aligned_space.h"

#include "observer_proxy.h"
#include "tbb_main.h"
#include "governor.h"
#include "scheduler.h"

namespace tbb {
namespace internal {

typedef spin_rw_mutex::scoped_lock task_scheduler_observer_mutex_scoped_lock;

/** aligned_space used here to shut up warnings when mutex destructor is called while threads are still using it. */
static aligned_space<spin_rw_mutex,1> the_task_scheduler_observer_mutex;
static observer_proxy* global_first_observer_proxy;
observer_proxy* global_last_observer_proxy;


#if TBB_USE_ASSERT
static atomic<int> observer_proxy_count;

struct check_observer_proxy_count {
    ~check_observer_proxy_count() {
        if( observer_proxy_count!=0 ) {
            runtime_warning( "Leaked %ld observer_proxy objects\n", long(observer_proxy_count) );
        }
    }
};

static check_observer_proxy_count the_check_observer_proxy_count;
#endif /* TBB_USE_ASSERT */

observer_proxy::observer_proxy( task_scheduler_observer_v3& tso ) : next(NULL), observer(&tso) {
#if TBB_USE_ASSERT
    ++observer_proxy_count;
#endif /* TBB_USE_ASSERT */
    // 1 for observer
    gc_ref_count = 1;
    {
        // Append to the global list
        task_scheduler_observer_mutex_scoped_lock lock(the_task_scheduler_observer_mutex.begin()[0],/*is_writer=*/true);
        observer_proxy* p = global_last_observer_proxy;
        prev = p;
        if( p ) 
            p->next=this;
        else 
            global_first_observer_proxy = this;
        global_last_observer_proxy = this;
    }
}

void observer_proxy::remove_from_list() {
    // Take myself off the global list.  
    if( next ) 
        next->prev = prev;
    else 
        global_last_observer_proxy = prev;
    if( prev )
        prev->next = next;
    else 
        global_first_observer_proxy = next;
#if TBB_USE_ASSERT
    poison_pointer(prev);
    poison_pointer(next);
    gc_ref_count = -666;
#endif /* TBB_USE_ASSERT */
}

void observer_proxy::remove_ref_slow() {
    int r = gc_ref_count;
    while(r>1) {
        __TBB_ASSERT( r!=0, NULL );
        int r_old = gc_ref_count.compare_and_swap(r-1,r);
        if( r_old==r ) {
            // Successfully decremented count.
            return;
        } 
        r = r_old;
    } 
    __TBB_ASSERT( r==1, NULL );
    // Reference count might go to zero
    {
        task_scheduler_observer_mutex_scoped_lock lock(the_task_scheduler_observer_mutex.begin()[0],/*is_writer=*/true);
        r = --gc_ref_count;
        if( !r ) {
            remove_from_list();
        } 
    }
    if( !r ) {
        __TBB_ASSERT( gc_ref_count == -666, NULL );
#if TBB_USE_ASSERT
        --observer_proxy_count;
#endif /* TBB_USE_ASSERT */
        delete this;
    }
}

observer_proxy* observer_proxy::process_list( observer_proxy* local_last, bool is_worker, bool is_entry ) {
    // Pointer p marches though the list.
    // If is_entry, start with our previous list position, otherwise start at beginning of list.
    observer_proxy* p = is_entry ? local_last : NULL;
    for(;;) { 
        task_scheduler_observer* tso=NULL;
        // Hold lock on list only long enough to advance to next proxy in list.
        { 
            task_scheduler_observer_mutex_scoped_lock lock(the_task_scheduler_observer_mutex.begin()[0],/*is_writer=*/false);
            do {
                if( local_last && local_last->observer ) {
                    // 2 = 1 for observer and 1 for local_last
                    __TBB_ASSERT( local_last->gc_ref_count>=2, NULL );  
                    // Can decrement count quickly, because it cannot become zero here.
                    --local_last->gc_ref_count;
                    local_last = NULL;
                } else {
                    // Use slow form of decrementing the reference count, after lock is released.
                }  
                if( p ) {
                    // We were already processing the list.
                    if( observer_proxy* q = p->next ) {
                        // Step to next item in list.
                        p=q;
                    } else {
                        // At end of list.
                        if( is_entry ) {  
                            // Remember current position in the list, so we can start at on the next call.
                            ++p->gc_ref_count;
                        } else {
                            // Finishin running off the end of the list
                            p=NULL;
                        }
                        goto done;
                    }
                } else {
                    // Starting pass through the list
                    p = global_first_observer_proxy;
                    if( !p ) 
                        goto done;
                } 
                tso = p->observer;
            } while( !tso );
            ++p->gc_ref_count;
            ++tso->my_busy_count;
        }
        __TBB_ASSERT( !local_last || p!=local_last, NULL );
        if( local_last )
            local_last->remove_ref_slow();
        // Do not hold any locks on the list while calling user's code.
        __TBB_TRY {    
            if( is_entry )
                tso->on_scheduler_entry( is_worker );
            else
                tso->on_scheduler_exit( is_worker );
        } __TBB_CATCH(...) {
            // Suppress exception, because user routines are supposed to be observing, not changing
            // behavior of a master or worker thread.
#if TBB_USE_ASSERT
            runtime_warning( "%s threw exception\n", is_entry ? "on_scheduler_entry" : "on_scheduler_exit"); 
#endif /* __TBB_USE_ASSERT */        
        }
        intptr_t bc = --tso->my_busy_count;
        __TBB_ASSERT_EX( bc>=0, "my_busy_count underflowed" );
        local_last = p;
    }
done:
    // Return new value to be used as local_last next time.
    if( local_last )
        local_last->remove_ref_slow();
    __TBB_ASSERT( !p || is_entry, NULL );
    return p;
}

void task_scheduler_observer_v3::observe( bool state ) {
    if( state ) {
        if( !my_proxy ) {
            if( !__TBB_InitOnce::initialization_done() )
                DoOneTimeInitializations();
            my_busy_count = 0;
            my_proxy = new observer_proxy(*this);
            if( generic_scheduler* s = governor::local_scheduler_if_initialized() ) {
                // Notify newly created observer of its own thread.
                // Any other pending observers are notified too.
                s->notify_entry_observers();
            }
        } 
    } else {
        if( observer_proxy* proxy = my_proxy ) {
            my_proxy = NULL;
            __TBB_ASSERT( proxy->gc_ref_count>=1, "reference for observer missing" );
            {
                task_scheduler_observer_mutex_scoped_lock lock(the_task_scheduler_observer_mutex.begin()[0],/*is_writer=*/true);
                proxy->observer = NULL;
            }
            proxy->remove_ref_slow();
            while( my_busy_count ) {
                __TBB_Yield();
            }
        }
    }
}

} // namespace internal
} // namespace tbb

#endif /* __TBB_SCHEDULER_OBSERVER */
