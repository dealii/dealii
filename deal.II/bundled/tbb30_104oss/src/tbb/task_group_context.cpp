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

#include "scheduler.h"

#include "tbb/task.h"
#include "tbb/tbb_exception.h"
#include "tbb/cache_aligned_allocator.h"
#include "itt_notify.h"

namespace tbb {

#if __TBB_TASK_GROUP_CONTEXT

using namespace internal;

//------------------------------------------------------------------------
// captured_exception
//------------------------------------------------------------------------

inline char* duplicate_string ( const char* src ) {
    char* dst = NULL;
    if ( src ) {
        size_t len = strlen(src) + 1;
        dst = (char*)allocate_via_handler_v3(len);
        strncpy (dst, src, len);
    }
    return dst;
}

void captured_exception::set ( const char* name, const char* info ) throw() {
    my_exception_name = duplicate_string( name );
    my_exception_info = duplicate_string( info );
}

void captured_exception::clear () throw() {
    deallocate_via_handler_v3 (const_cast<char*>(my_exception_name));
    deallocate_via_handler_v3 (const_cast<char*>(my_exception_info));
}

captured_exception* captured_exception::move () throw() {
    captured_exception *e = (captured_exception*)allocate_via_handler_v3(sizeof(captured_exception));
    if ( e ) {
        ::new (e) captured_exception();
        e->my_exception_name = my_exception_name;
        e->my_exception_info = my_exception_info;
        e->my_dynamic = true;
        my_exception_name = my_exception_info = NULL;
    }
    return e;
}

void captured_exception::destroy () throw() {
    __TBB_ASSERT ( my_dynamic, "Method destroy can be used only on objects created by clone or allocate" );
    if ( my_dynamic ) {
        this->captured_exception::~captured_exception();
        deallocate_via_handler_v3 (this);
    }
}

captured_exception* captured_exception::allocate ( const char* name, const char* info ) {
    captured_exception *e = (captured_exception*)allocate_via_handler_v3( sizeof(captured_exception) );
    if ( e ) {
        ::new (e) captured_exception(name, info);
        e->my_dynamic = true;
    }
    return e;
}

const char* captured_exception::name() const throw() {
    return my_exception_name;
}

const char* captured_exception::what() const throw() {
    return my_exception_info;
}


//------------------------------------------------------------------------
// tbb_exception_ptr
//------------------------------------------------------------------------

#if !TBB_USE_CAPTURED_EXCEPTION

namespace internal {

template<typename T>
tbb_exception_ptr* AllocateExceptionContainer( const T& src ) {
    tbb_exception_ptr *eptr = (tbb_exception_ptr*)allocate_via_handler_v3( sizeof(tbb_exception_ptr) );
    if ( eptr )
        new (eptr) tbb_exception_ptr(src);
    return eptr;
}

tbb_exception_ptr* tbb_exception_ptr::allocate () {
    return AllocateExceptionContainer( std::current_exception() );
}

tbb_exception_ptr* tbb_exception_ptr::allocate ( const tbb_exception& ) {
    return AllocateExceptionContainer( std::current_exception() );
}

tbb_exception_ptr* tbb_exception_ptr::allocate ( captured_exception& src ) {
    tbb_exception_ptr *res = AllocateExceptionContainer( src );
    src.destroy();
    return res;
}

void tbb_exception_ptr::destroy () throw() {
    this->tbb_exception_ptr::~tbb_exception_ptr();
    deallocate_via_handler_v3 (this);
}

} // namespace internal
#endif /* !TBB_USE_CAPTURED_EXCEPTION */


//------------------------------------------------------------------------
// task_group_context
//------------------------------------------------------------------------

task_group_context::~task_group_context () {
    if ( my_kind != isolated ) {
        generic_scheduler *s = (generic_scheduler*)my_owner;
        if ( governor::is_set(s) ) {
            // Local update of the context list 
            uintptr_t local_count_snapshot = s->local_cancel_count;
            s->local_ctx_list_update = 1;
            __TBB_full_memory_fence();
            if ( s->nonlocal_ctx_list_update ) {
                spin_mutex::scoped_lock lock(s->context_list_mutex);
                my_node.my_prev->my_next = my_node.my_next;
                my_node.my_next->my_prev = my_node.my_prev;
                s->local_ctx_list_update = 0;
            }
            else {
                my_node.my_prev->my_next = my_node.my_next;
                my_node.my_next->my_prev = my_node.my_prev;
                __TBB_store_with_release( s->local_ctx_list_update, 0 );
                if ( local_count_snapshot != global_cancel_count ) {
                    // Another thread was propagating cancellation request when we removed
                    // ourselves from the list. We must ensure that it is not accessing us 
                    // when this destructor finishes. We'll be able to acquire the lock 
                    // below only after the other thread finishes with us.
                    spin_mutex::scoped_lock lock(s->context_list_mutex);
                }
            }
        }
        else {
            // Nonlocal update of the context list 
            if ( __TBB_FetchAndStoreW(&my_kind, dying) == detached ) {
                my_node.my_prev->my_next = my_node.my_next;
                my_node.my_next->my_prev = my_node.my_prev;
            }
            else {
                __TBB_FetchAndAddW(&s->nonlocal_ctx_list_update, 1);
                spin_wait_until_eq( s->local_ctx_list_update, 0u );
                s->context_list_mutex.lock();
                my_node.my_prev->my_next = my_node.my_next;
                my_node.my_next->my_prev = my_node.my_prev;
                s->context_list_mutex.unlock();
                __TBB_FetchAndAddW(&s->nonlocal_ctx_list_update, -1);
            }
        }
    }
#if TBB_USE_DEBUG
    my_version_and_traits = 0xDeadBeef;
#endif /* TBB_USE_DEBUG */
    if ( my_exception )
        my_exception->destroy();
    if (itt_caller != ITT_CALLER_NULL) ITT_STACK(caller_destroy, itt_caller);
}

void task_group_context::init () {
    __TBB_ASSERT ( sizeof(uintptr_t) < 32, "Layout of my_version_and_traits must be reconsidered on this platform" );
    __TBB_ASSERT ( sizeof(task_group_context) == 2 * NFS_MaxLineSize, "Context class has wrong size - check padding and members alignment" );
    __TBB_ASSERT ( (uintptr_t(this) & (sizeof(my_cancellation_requested) - 1)) == 0, "Context is improperly aligned" );
    __TBB_ASSERT ( my_kind == isolated || my_kind == bound, "Context can be created only as isolated or bound" );
    my_parent = NULL;
    my_cancellation_requested = 0;
    my_exception = NULL;
    itt_caller = ITT_CALLER_NULL;
    if ( my_kind == bound ) {
        generic_scheduler *s = governor::local_scheduler();
        my_owner = s;
        __TBB_ASSERT ( my_owner, "Thread has not activated a task_scheduler_init object?" );
        // Backward links are used by this thread only, thus no fences are necessary
        my_node.my_prev = &s->context_list_head;
        s->context_list_head.my_next->my_prev = &my_node;
        my_node.my_next = s->context_list_head.my_next;
        // Thread local list of contexts allows concurrent traversal by another 
        // thread while propagating cancellation request. Release fence ensures 
        // visibility of my_node's members in the traversing thread.
        __TBB_store_with_release(s->context_list_head.my_next, &my_node);
    }
}

bool task_group_context::cancel_group_execution () {
    __TBB_ASSERT ( my_cancellation_requested == 0 || my_cancellation_requested == 1, "Invalid cancellation state");
    if ( my_cancellation_requested || __TBB_CompareAndSwapW(&my_cancellation_requested, 1, 0) ) {
        // This task group has already been canceled
        return false;
    }
#if __TBB_ARENA_PER_MASTER
    governor::local_scheduler()->my_arena->propagate_cancellation( *this );
#else /* !__TBB_ARENA_PER_MASTER */
    governor::local_scheduler()->propagate_cancellation( *this );
#endif /* !__TBB_ARENA_PER_MASTER */
    return true;
}

bool task_group_context::is_group_execution_cancelled () const {
    return my_cancellation_requested != 0;
}

// IMPORTANT: It is assumed that this method is not used concurrently!
void task_group_context::reset () {
    //! \todo Add assertion that this context does not have children
    // No fences are necessary since this context can be accessed from another thread
    // only after stealing happened (which means necessary fences were used).
    if ( my_exception )  {
        my_exception->destroy();
        my_exception = NULL;
    }
    my_cancellation_requested = 0;
}

void task_group_context::propagate_cancellation_from_ancestors () {
    task_group_context *ancestor = my_parent;
    while ( ancestor && !ancestor->my_cancellation_requested )
        ancestor = ancestor->my_parent;
    if ( ancestor ) {
        // One of my ancestor groups was canceled. Cancel all its descendants in my heritage line.
        task_group_context *ctx = this;
        do {
            ctx->my_cancellation_requested = 1;
            ctx = ctx->my_parent;
        } while ( ctx != ancestor );
    }
}

void task_group_context::register_pending_exception () {
    if ( my_cancellation_requested )
        return;
#if TBB_USE_EXCEPTIONS
    try {
        throw;
    } TbbCatchAll( this );
#endif /* TBB_USE_EXCEPTIONS */
}

#endif /* __TBB_TASK_GROUP_CONTEXT */

} // namespace tbb
