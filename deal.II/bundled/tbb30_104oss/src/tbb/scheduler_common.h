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

#ifndef _TBB_scheduler_common_H
#define _TBB_scheduler_common_H

#include "tbb/tbb_stddef.h"

#include <string.h>  // for memset, memcpy, memmove

#include "tbb_statistics.h"

/* Temporarily change "private" to "public" while including "tbb/task.h".
   This hack allows us to avoid publishing internal types and methods
   in the public header files just for sake of friend declarations. */
#define private public
#include "tbb/task.h"
#include "tbb/tbb_exception.h"
#undef private

// This macro is an attempt to get rid of ugly ifdefs in the shared parts of the code. 
// It drops the second argument depending on whether the controlling macro is defined. 
// The first argument is just a convenience allowing to keep comma before the macro usage.
#if __TBB_TASK_GROUP_CONTEXT
    #define __TBB_CONTEXT_ARG(arg1, context) arg1, context
#else /* !__TBB_TASK_GROUP_CONTEXT */
    #define __TBB_CONTEXT_ARG(arg1, context) arg1
#endif /* !__TBB_TASK_GROUP_CONTEXT */

#if DO_TBB_TRACE
#include <cstdio>
#define TBB_TRACE(x) ((void)std::printf x)
#else
#define TBB_TRACE(x) ((void)(0))
#endif /* DO_TBB_TRACE */

#if _MSC_VER && !defined(__INTEL_COMPILER)
    // Workaround for overzealous compiler warnings
    // These particular warnings are so ubiquitous that no attempt is made to narrow 
    // the scope of the warnings.
    #pragma warning (disable: 4100 4127 4312 4244 4267 4706)
#endif

namespace tbb {
namespace internal {

/** Defined in scheduler.cpp **/
extern uintptr_t global_cancel_count;

//! Alignment for a task object
const size_t task_alignment = 16;

//! Number of bytes reserved for a task prefix
/** If not exactly sizeof(task_prefix), the extra bytes *precede* the task_prefix. */
const size_t task_prefix_reservation_size = ((sizeof(internal::task_prefix)-1)/task_alignment+1)*task_alignment;

//! Definitions for bits in task_prefix::extra_state
enum task_extra_state {
    //! Tag for v1 tasks (i.e. tasks in TBB 1.0 and 2.0)
    es_version_1_task = 0,
    //! Tag for v3 tasks (i.e. tasks in TBB 2.1-2.2)
    es_version_3_task = 1,
    //! Tag for v3 task_proxy.
    es_task_proxy = 0x20,
    //! Set if ref_count might be changed by another thread.  Used for debugging.
    es_ref_count_active = 0x40,
    //! Set if the task has been stolen
    es_task_is_stolen = 0x80
};

//! Optimization hint to free_task that enables it omit unnecessary tests and code.
enum free_task_hint {
    //! No hint 
    no_hint=0,
    //! Task is known to have been allocated by this scheduler
    local_task=1,
    //! Task is known to be a small task.
    /** Task should be returned to the free list of *some* scheduler, possibly not this scheduler. */
    small_task=2,
    //! Bitwise-OR of local_task and small_task.  
    /** Task should be returned to free list of this scheduler. */
    small_local_task=3
};

//------------------------------------------------------------------------
// Debugging support
//------------------------------------------------------------------------

#if TBB_USE_ASSERT

static const uintptr_t venom = 
#if __TBB_WORDSIZE == 8
        0xDDEEAADDDEADBEEF;
#else
        0xDEADBEEF;
#endif


/** In contrast to poison_pointer() and assert_task_valid() poison_value() is a macro 
    because the variable used as its argument may be undefined in release builds. **/
#define poison_value(g) (g = venom)

/** Expected to be used in assertions only, thus no empty form is defined. **/
inline bool is_alive( uintptr_t v ) { return v != venom; }

/** Logically, this method should be a member of class task.
    But we do not want to publish it, so it is here instead. */
inline void assert_task_valid( const task& task ) {
    __TBB_ASSERT( &task!=NULL, NULL );
    __TBB_ASSERT( !is_poisoned(&task), NULL );
    __TBB_ASSERT( (uintptr_t)&task % task_alignment == 0, "misaligned task" );
    __TBB_ASSERT( (unsigned)task.state()<=(unsigned)task::recycle, "corrupt task (invalid state)" );
}

#else /* !TBB_USE_ASSERT */

#define poison_value(g) ((void)0)

inline void assert_task_valid( const task& ) {}

#endif /* !TBB_USE_ASSERT */

//------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------

inline bool ConcurrentWaitsEnabled ( task& t ) {
    return (t.prefix().context->my_version_and_traits & task_group_context::concurrent_wait) != 0;
}

inline bool CancellationInfoPresent ( task& t ) {
    return t.prefix().context->my_cancellation_requested != 0;
}

#if __TBB_TASK_GROUP_CONTEXT
#if TBB_USE_CAPTURED_EXCEPTION
    inline tbb_exception* TbbCurrentException( task_group_context*, tbb_exception* src) { return src->move(); }
    inline tbb_exception* TbbCurrentException( task_group_context*, captured_exception* src) { return src; }
#else
    // Using macro instead of an inline function here allows to avoid evaluation of the 
    // TbbCapturedException expression when exact propagation is enabled for the context.
    #define TbbCurrentException(context, TbbCapturedException) \
        context->my_version_and_traits & task_group_context::exact_exception    \
            ? tbb_exception_ptr::allocate()    \
            : tbb_exception_ptr::allocate( *(TbbCapturedException) );
#endif /* !TBB_USE_CAPTURED_EXCEPTION */

#define TbbRegisterCurrentException(context, TbbCapturedException) \
    if ( context->cancel_group_execution() ) {  \
        /* We are the first to signal cancellation, so store the exception that caused it. */  \
        context->my_exception = TbbCurrentException( context, TbbCapturedException ); \
    }

#define TbbCatchAll(context)  \
    catch ( tbb_exception& exc ) {  \
        TbbRegisterCurrentException( context, &exc );   \
    } catch ( std::exception& exc ) {   \
        TbbRegisterCurrentException( context, captured_exception::allocate(typeid(exc).name(), exc.what()) ); \
    } catch ( ... ) {   \
        TbbRegisterCurrentException( context, captured_exception::allocate("...", "Unidentified exception") );\
    }
#endif /* __TBB_TASK_GROUP_CONTEXT */

} // namespace internal
} // namespace tbb

#endif /* _TBB_scheduler_common_H */
