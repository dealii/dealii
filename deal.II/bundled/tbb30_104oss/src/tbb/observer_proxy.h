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

#ifndef _TBB_observer_proxy_H
#define _TBB_observer_proxy_H

#if __TBB_SCHEDULER_OBSERVER

#include "tbb/task_scheduler_observer.h"

namespace tbb {
namespace internal {

class observer_proxy {
    friend class task_scheduler_observer_v3;
    //! Reference count used for garbage collection.
    /** 1 for reference from my task_scheduler_observer.
        1 for each local_last_observer_proxy that points to me. 
        No accounting for predecessor in the global list. 
        No accounting for global_last_observer_proxy that points to me. */
    atomic<int> gc_ref_count;
    //! Pointer to next task_scheduler_observer 
    /** Valid even when *this has been removed from the global list. */
    observer_proxy* next; 
    //! Pointer to previous task_scheduler_observer in global list.
    observer_proxy* prev; 
    //! Associated observer
    task_scheduler_observer* observer;
    //! Account for removing reference from p.  No effect if p is NULL.
    void remove_ref_slow();
    void remove_from_list(); 
    observer_proxy( task_scheduler_observer_v3& wo ); 
public:
    static observer_proxy* process_list( observer_proxy* local_last, bool is_worker, bool is_entry );
};

extern observer_proxy* global_last_observer_proxy;

} // namespace internal
} // namespace tbb

#endif /* __TBB_SCHEDULER_OBSERVER */

#endif /* _TBB_observer_proxy_H */
